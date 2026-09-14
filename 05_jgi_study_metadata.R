#!/usr/bin/env Rscript
# 05_jgi_study_metadata.R
#
# Builds Tables S11 (study metadata) and S12 (sequencing statistics) for the six
# matched JGI studies used in the viral-versus-total comparison. Joins the
# targeted gene counts to GSVA sample metadata (location, soil type, ecosystem)
# and to a keyed table of sequencing depth, platform and management history.
#
# In:  data/jgi_targeted_counts.csv
#      data/GSVA_sample_metadata_5.csv
#      data/jgi_study_metadata.csv
# Out: results/Tables_S11_S12_JGI_matched_comparison_metadata.xlsx
#      results/05_sessionInfo.txt

library(openxlsx)

set.seed(42)

counts <- read.csv("data/jgi_targeted_counts.csv", comment.char = "#",
                   stringsAsFactors = FALSE)

table2 <- data.frame(
  JGI_Study_ID = as.character(counts$img_taxon_oid),
  Study_Num    = counts$study_label,
  Function     = counts$function_label,
  Target_IDs   = paste0(counts$ko_term, "/", counts$pfam_id),
  Total_KEGG   = counts$total_kegg_genes,
  Viral_KEGG   = counts$viral_kegg_genes,
  stringsAsFactors = FALSE
)
table2$Viral_pct <- round(100 * table2$Viral_KEGG / table2$Total_KEGG, 2)

gsva <- read.csv("data/GSVA_sample_metadata_5.csv")
gsva$IMG.taxon_oid <- as.character(gsva$IMG.taxon_oid)

jgi_ids <- table2$JGI_Study_ID
gsva_sub <- gsva[gsva$IMG.taxon_oid %in% jgi_ids,
                 c("IMG.taxon_oid", "Genome.Name...Sample.Name",
                   "Geographic.Location", "Latitude", "Longitude",
                   "GOLD.Ecosystem.Type", "GOLD.Ecosystem.Subtype",
                   "GOLD.Specific.Ecosystem", "Habitat",
                   "Genome.Size.assembled", "Gene.Count.assembled")]

names(gsva_sub)[names(gsva_sub) == "IMG.taxon_oid"] <- "JGI_Study_ID"
names(gsva_sub)[names(gsva_sub) == "Genome.Name...Sample.Name"] <- "Sample_Name"

# Raw_Gbp for HiSeq studies is estimated as raw reads x 151 bp; NovaSeq READMEs
# report raw base counts directly (Raw_Gbp_source). Same_Assembly = viral contigs
# were identified within the same assembled metagenome as the total annotations.

study_meta <- read.csv("data/jgi_study_metadata.csv", comment.char = "#",
                       stringsAsFactors = FALSE)
study_meta$JGI_Study_ID <- as.character(study_meta$JGI_Study_ID)

expected_ids <- as.character(counts$img_taxon_oid)
stopifnot(
  length(unique(expected_ids)) == 6,
  setequal(study_meta$JGI_Study_ID, expected_ids),
  !anyDuplicated(study_meta$JGI_Study_ID),
  setequal(gsva_sub$JGI_Study_ID, expected_ids),
  !anyDuplicated(gsva_sub$JGI_Study_ID),
  all(table2$Total_KEGG > 0),
  all(table2$Viral_KEGG >= 0 & table2$Viral_KEGG <= table2$Total_KEGG)
)

merged <- merge(table2, gsva_sub, by = "JGI_Study_ID", all.x = TRUE)
merged <- merge(merged, study_meta, by = "JGI_Study_ID", all.x = TRUE)
stopifnot(
  nrow(merged) == 6,
  !anyNA(merged$Sample_Name),
  !anyNA(merged$GOLD.Ecosystem.Type),
  !anyNA(merged$Platform)
)

merged$Ecosystem <- ifelse(
  merged$GOLD.Ecosystem.Type %in% c("Rhizosphere", "Roots", "Rhizoplane"),
    "Plant-associated",
  ifelse(merged$GOLD.Ecosystem.Type == "Peat", "Peat",
  ifelse(merged$GOLD.Specific.Ecosystem %in%
           c("Agricultural land", "Agricultural", "Agricultural soil"),
           "Agricultural",
  ifelse(merged$GOLD.Specific.Ecosystem %in% c("Forest Soil", "Forest soil"),
           "Forest",
           "Unclassified"))))

tab_s11 <- data.frame(
  `JGI Study ID`         = merged$JGI_Study_ID,
  `Study #`              = merged$Study_Num,
  `Function`             = merged$Function,
  `Target IDs (KO/Pfam)` = merged$Target_IDs,
  `Ecosystem`            = merged$Ecosystem,
  `Soil type`            = merged$Soil_Type,
  `Geographic location`  = merged$Geographic.Location,
  `Latitude`             = merged$Latitude,
  `Longitude`            = merged$Longitude,
  `Management history`   = merged$Management_History,
  check.names = FALSE, stringsAsFactors = FALSE
)

tab_s12 <- data.frame(
  `JGI Study ID`         = merged$JGI_Study_ID,
  `Sequencing platform`  = merged$Platform,
  `Raw reads (M)`        = round(merged$Raw_Reads / 1e6, 1),
  `Raw Gbp`              = merged$Raw_Gbp,
  `Gbp source`           = merged$Raw_Gbp_source,
  `Filtered reads (M)`   = round(merged$Filtered_Reads / 1e6, 1),
  `Assembler`            = merged$Assembler,
  `Assembly date`        = merged$Assembly_Date,
  `Assembled size (Gbp)` = round(merged$Genome.Size.assembled / 1e9, 2),
  `Gene count`           = merged$Gene.Count.assembled,
  `Total KEGG (target)`  = merged$Total_KEGG,
  `Viral KEGG (target)`  = merged$Viral_KEGG,
  `Viral %`              = merged$Viral_pct,
  `Same assembly`        = merged$Same_Assembly,
  `Sample name`          = merged$Sample_Name,
  check.names = FALSE, stringsAsFactors = FALSE
)

tab_s11 <- tab_s11[order(tab_s11$`JGI Study ID`), ]
tab_s12 <- tab_s12[order(tab_s12$`JGI Study ID`), ]

dir.create("results", showWarnings = FALSE)
hdr <- createStyle(textDecoration = "bold", border = "Bottom", borderStyle = "thin")

wb <- createWorkbook()

addWorksheet(wb, "Table S11")
writeData(wb, "Table S11",
  x = "Table S11. Metadata for the six matched JGI studies used in the viral-versus-total metagenomic comparison.",
  startRow = 1, startCol = 1)
addStyle(wb, "Table S11", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, "Table S11", tab_s11, startRow = 2, headerStyle = hdr)
setColWidths(wb, "Table S11", cols = 1:ncol(tab_s11), widths = "auto")

addWorksheet(wb, "Table S12")
writeData(wb, "Table S12",
  x = "Table S12. Sequencing statistics for the six matched JGI studies used in the viral-versus-total metagenomic comparison.",
  startRow = 1, startCol = 1)
addStyle(wb, "Table S12", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, "Table S12", tab_s12, startRow = 2, headerStyle = hdr)
setColWidths(wb, "Table S12", cols = 1:ncol(tab_s12), widths = "auto")

footnote_row <- nrow(tab_s12) + 4
writeData(wb, "Table S12", startRow = footnote_row, startCol = 1,
  x = paste0(
    "Notes: Raw Gbp for HiSeq studies estimated as raw reads x 151 bp ",
    "(interleaved paired-end fastq). NovaSeq studies report raw base counts ",
    "directly from JGI README files. 'Same assembly' indicates that viral ",
    "contigs were identified within the same assembled metagenome as the total ",
    "metagenomic annotations; viral and total annotations derive from the same ",
    "sequencing run, not from separate virome libraries. LTER = Long-Term ",
    "Ecological Research site (Kellogg Biological Station). SPRUCE = Spruce ",
    "and Peatland Responses Under Changing Environments."
  ))

saveWorkbook(wb, "results/Tables_S11_S12_JGI_matched_comparison_metadata.xlsx",
             overwrite = TRUE)

cat("\n=== Tables S11 / S12 summary ===\n")
cat("Studies:", nrow(tab_s12), "\n\n")

cat("Sequencing depth range:\n")
cat("  Raw Gbp:", range(tab_s12$`Raw Gbp`), "(min/max)\n")
cat("  Fold difference:", round(max(tab_s12$`Raw Gbp`) / min(tab_s12$`Raw Gbp`), 1), "\n\n")

cat("Platforms:", paste(unique(tab_s12$`Sequencing platform`), collapse = ", "), "\n")
cat("Assemblers:", paste(unique(tab_s12$`Assembler`), collapse = ", "), "\n")
cat("All same assembly (viral from total):", all(tab_s12$`Same assembly` == "Yes"), "\n\n")

print(tab_s12[, c("JGI Study ID", "Raw Gbp",
                  "Total KEGG (target)", "Viral KEGG (target)", "Viral %")])

sink("results/05_sessionInfo.txt")
cat("05_jgi_study_metadata.R\n")
cat("Run:", format(Sys.time()), "\n\n")
sessionInfo()
sink()

cat("\nDone. Output: results/Tables_S11_S12_JGI_matched_comparison_metadata.xlsx\n")
