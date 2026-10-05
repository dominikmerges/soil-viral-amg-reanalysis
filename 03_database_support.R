#!/usr/bin/env Rscript
# 03_database_support.R
#
# Counts supporting databases (1-3) per functional category and builds Table 1.
#
# CAZy indexes only carbohydrate-active enzymes, so the maximum is 3 databases
# for carbon cycling but 2 for nitrogen cycling and AMR. Support counts are
# therefore reported within category, not compared across categories.
#
# Uses both filtering sets: `annot` (enzyme-resolved) for all reported values,
# `annot_any` (post-lysis) for the annotated-at-all rate and the tool
# comparison, where filtering one arm only would confound tool with filter.
#
# In:  data/metagenome_confidence_annotations.csv
# Out: results/03_database_support_by_category.tsv
#      results/03_Table_S_database_support.xlsx
#      results/03_sessionInfo.txt

source("_common.R")

library(tidyverse)
library(data.table)

set.seed(42)

data_dir    <- "data_10kb"
dirs        <- mfd_init_dirs()
results_dir <- dirs$results_dir

annot_any <- mfd_load_annotations(data_dir, set = "post_lysis")
annot     <- annot_any[evidence_type == "specific"]

counts    <- mfd_counts(annot_any)

cat("post_lysis:", nrow(annot_any), " specific:", nrow(annot), "\n")

summarise_db_support <- function(dt) {
  dt[, .(
    n_genes = .N,
    n_1db   = sum(n_db == 1),
    n_2db   = sum(n_db == 2),
    n_3db   = sum(n_db == 3),
    pct_1db = round(100 * sum(n_db == 1) / .N, 1),
    pct_2db = round(100 * sum(n_db == 2) / .N, 1),
    pct_3db = round(100 * sum(n_db == 3) / .N, 1),
    max_possible_db = fifelse(
      first(functional_category) %in% CAZY_ELIGIBLE_CATEGORIES, 3L, 2L
    )
  ), by = functional_category]
}

db_by_cat     <- summarise_db_support(annot)
db_by_cat_any <- summarise_db_support(annot_any)

cat("\n--- Database support by functional category (enzyme-resolved) ---\n")
print(db_by_cat)
fwrite(db_by_cat,
       file.path(results_dir, "03_database_support_by_category.tsv"),
       sep = "\t")

library(openxlsx)

# Helper: prettify category labels (same as script 02)
prettify <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)
  x <- gsub("Amr", "AMR", x)
  x <- gsub("Kegg", "KEGG", x)
  x <- gsub("Pfam", "Pfam", x)
  x <- gsub("Cazy", "CAZy", x)
  x
}

# Styles
title_style    <- createStyle(fontSize = 11, textDecoration = "bold")
header_style   <- createStyle(
  fontSize = 10, textDecoration = "bold",
  border = "TopBottom", borderStyle = "thin",
  halign = "center", wrapText = TRUE
)
body_style     <- createStyle(fontSize = 10, halign = "center")
body_left      <- createStyle(fontSize = 10, halign = "left")
bottom_border  <- createStyle(border = "bottom", borderStyle = "thin")
category_style <- createStyle(fontSize = 10, halign = "left",
                              textDecoration = "italic")
note_style     <- createStyle(fontSize = 9, halign = "left",
                              wrapText = TRUE)

# Compute values for all sheets

all_total <- nrow(annot)      # enzyme-resolved, canonical
any_total <- nrow(annot_any)  # post-lysis, for the rarity rows only

# Every enzyme-resolved gene has 1-3 supporting databases, and the support-level
# counts must partition the analysis set exactly.
stopifnot(
  all(annot$n_db %in% 1:3),
  sum(annot$n_db == 1) + sum(annot$n_db == 2) + sum(annot$n_db == 3) == all_total
)

# Database-specific gene counts (canonical set, for revised Table 1)
n_kegg_genes <- annot[has_kegg == TRUE, .N]
n_pfam_genes <- annot[has_pfam == TRUE, .N]
n_cazy_genes <- annot[has_cazy == TRUE, .N]

# Same counts on the post-lysis set, needed for the original-vs-MFD comparison
n_kegg_any <- annot_any[has_kegg == TRUE, .N]
n_pfam_any <- annot_any[has_pfam == TRUE, .N]
n_cazy_any <- annot_any[has_cazy == TRUE, .N]

pct_cat <- function(cat) round(100 * db_by_cat[functional_category == cat,
                                               n_genes] / all_total, 1)

# Sheet 1: Table 1
# Both denominators are carried explicitly so they cannot be conflated.

t1 <- data.table(
  Metric = c(
    "Total viral genes analysed",
    "Total soil samples with ≥1 viral contig",
    "Ecosystem types represented",
    "Functionally annotated viral genes (any annotation)",
    "  of which broad (no enzyme resolved)",
    "Enzyme-resolved annotations (analysis set)",
    "Carbon cycling genes",
    "Nitrogen cycling genes",
    "Antibiotic resistance-associated genes",
    "Genes supported by 1 database",
    "Genes supported by 2 databases",
    "Genes supported by 3 databases",
    "KEGG annotated genes",
    "Pfam annotated genes",
    "CAZy annotated genes"
  ),
  Value = c(
    formatC(GSVA_TOTAL_VIRAL_GENES, big.mark = ","),
    formatC(GSVA_N_SAMPLES, big.mark = ","),
    as.character(GSVA_N_ECOSYSTEMS),
    formatC(any_total, big.mark = ","),
    formatC(counts$broad, big.mark = ","),
    formatC(all_total, big.mark = ","),
    formatC(db_by_cat[functional_category == "carbon_cycling", n_genes],
            big.mark = ","),
    formatC(db_by_cat[functional_category == "nitrogen_cycling", n_genes],
            big.mark = ","),
    formatC(db_by_cat[functional_category == "antibiotics", n_genes],
            big.mark = ","),
    formatC(sum(annot$n_db == 1), big.mark = ","),
    formatC(sum(annot$n_db == 2), big.mark = ","),
    formatC(sum(annot$n_db == 3), big.mark = ","),
    formatC(n_kegg_genes, big.mark = ","),
    formatC(n_pfam_genes, big.mark = ","),
    formatC(n_cazy_genes, big.mark = ",")
  ),
  Notes = c(
    "GSV Atlas viral gene catalogue",
    "Samples represented in the catalogue",
    "Ecosystem categories",
    paste0(round(100 * any_total / GSVA_TOTAL_VIRAL_GENES, 2),
           "% of all viral genes; peptidoglycanase excluded"),
    "Umbrella GO terms or non-specific text patterns; no enzyme subcategory",
    paste0(round(100 * all_total / GSVA_TOTAL_VIRAL_GENES, 2),
           "% of all viral genes; basis for all results below"),
    paste0(pct_cat("carbon_cycling"),   "% of enzyme-resolved annotations"),
    paste0(pct_cat("nitrogen_cycling"), "% of enzyme-resolved annotations"),
    paste0(pct_cat("antibiotics"),      "% of enzyme-resolved annotations"),
    paste0(round(100 * sum(annot$n_db == 1) / all_total, 1),
           "% of enzyme-resolved annotations"),
    paste0(round(100 * sum(annot$n_db == 2) / all_total, 1),
           "% of enzyme-resolved annotations"),
    paste0(round(100 * sum(annot$n_db == 3) / all_total, 1),
           "% of enzyme-resolved annotations; carbon cycling only (see note)"),
    paste0(round(100 * n_kegg_genes / GSVA_TOTAL_VIRAL_GENES, 3),
           "% of all viral genes"),
    paste0(round(100 * n_pfam_genes / GSVA_TOTAL_VIRAL_GENES, 3),
           "% of all viral genes"),
    paste0(round(100 * n_cazy_genes / GSVA_TOTAL_VIRAL_GENES, 3),
           "% of all viral genes")
  )
)

# Sheet 2: Database support within each functional category

s2 <- copy(db_by_cat)
s2[, `Functional category` := prettify(functional_category)]
s2_out <- s2[, .(
  `Functional category`,
  `No. genes`       = n_genes,
  `1 database`      = n_1db,
  `2 databases`     = n_2db,
  `3 databases`     = n_3db,
  `% 1 database`    = pct_1db,
  `% 2 databases`   = pct_2db,
  `% 3 databases`   = pct_3db,
  `Max. achievable`  = max_possible_db
)]

# Sheet 3: manual keyword matching vs MetaFuncDecoder
kw <- read.csv(file.path(data_dir, "original_keyword_counts.csv"),
               comment.char = "#", stringsAsFactors = FALSE)
kw_val <- setNames(kw$count, kw$metric)
kw_keys <- c("functionally_annotated", "carbon_cycling", "nitrogen_cycling",
             "antibiotics", "db3", "db2", "db1", "kegg", "pfam", "cazy")
stopifnot(all(kw_keys %in% names(kw_val)))

comp <- data.table(
  Metric = c(
    "Functionally annotated viral genes",
    "Carbon cycling genes",
    "Nitrogen cycling genes",
    "AMR-associated genes",
    "3-database support (\"high\")",
    "2-database support (\"medium\")",
    "1-database support (\"low\")",
    "KEGG annotated genes",
    "Pfam annotated genes",
    "CAZy annotated genes"
  ),
  `Manual keyword matching` = unname(kw_val[kw_keys]),
  `MetaFuncDecoder v1.0.0` = c(
    any_total,
    db_by_cat_any[functional_category == "carbon_cycling", n_genes],
    db_by_cat_any[functional_category == "nitrogen_cycling", n_genes],
    db_by_cat_any[functional_category == "antibiotics", n_genes],
    sum(annot_any$n_db == 3), sum(annot_any$n_db == 2), sum(annot_any$n_db == 1),
    n_kegg_any, n_pfam_any, n_cazy_any
  )
)

comp[, Difference := `MetaFuncDecoder v1.0.0` - `Manual keyword matching`]
comp[, `% change` := round(100 * Difference / `Manual keyword matching`, 1)]
comp[, Direction := fifelse(
  Difference == 0, "Identical",
  fifelse(Difference < 0, "Fewer (stricter)", "More")
)]

# Build workbook
wb <- createWorkbook()

addWorksheet(wb, "Revised Table 1")

title1 <- paste0(
  "Table 1 (revised) | Global soil viral functional landscape overview. ",
  "All values derived from MetaFuncDecoder v1.0.0. Peptidoglycanase (lysis) ",
  "genes are excluded throughout; functional and database-support rows refer ",
  "to the enzyme-resolved analysis set (broad, non-enzyme-resolved matches ",
  "excluded)."
)
writeData(wb, 1, title1, startRow = 1, startCol = 1)
addStyle(wb, 1, title_style, rows = 1, cols = 1)
mergeCells(wb, 1, cols = 1:3, rows = 1)

writeData(wb, 1, t1, startRow = 3, headerStyle = header_style)
nt1 <- nrow(t1)
addStyle(wb, 1, body_left, rows = 4:(3 + nt1), cols = 1, gridExpand = TRUE)
addStyle(wb, 1, body_style, rows = 4:(3 + nt1), cols = 2, gridExpand = TRUE)
addStyle(wb, 1, body_left, rows = 4:(3 + nt1), cols = 3, gridExpand = TRUE)
addStyle(wb, 1, bottom_border, rows = 3 + nt1, cols = 1:3, stack = TRUE)

# Table note about 3-database ceiling
note_row <- 3 + nt1 + 2
note_text <- paste0(
  "Note: The maximum number of supporting databases is 3 for carbon-cycling ",
  "annotations (CAZy, KEGG, Pfam) but 2 for nitrogen-cycling and AMR ",
  "annotations (KEGG and Pfam only), because CAZy exclusively indexes ",
  "carbohydrate-active enzyme families."
)
writeData(wb, 1, note_text, startRow = note_row, startCol = 1)
addStyle(wb, 1, note_style, rows = note_row, cols = 1)
mergeCells(wb, 1, cols = 1:3, rows = note_row)

setColWidths(wb, 1, cols = 1, widths = 42)
setColWidths(wb, 1, cols = 2, widths = 14)
setColWidths(wb, 1, cols = 3, widths = 52)

addWorksheet(wb, "Database support by category")

title2 <- paste0(
  "Table S14 | Number of supporting databases per functional annotation, ",
  "stratified by functional category."
)
writeData(wb, 2, title2, startRow = 1, startCol = 1)
addStyle(wb, 2, title_style, rows = 1, cols = 1)
mergeCells(wb, 2, cols = 1:ncol(s2_out), rows = 1)

writeData(wb, 2, s2_out, startRow = 3, headerStyle = header_style)
ns2 <- nrow(s2_out)
addStyle(wb, 2, category_style, rows = 4:(3 + ns2), cols = 1, gridExpand = TRUE)
addStyle(wb, 2, body_style, rows = 4:(3 + ns2), cols = 2:ncol(s2_out),
         gridExpand = TRUE)
addStyle(wb, 2, bottom_border, rows = 3 + ns2, cols = 1:ncol(s2_out),
         stack = TRUE)

setColWidths(wb, 2, cols = 1, widths = 22)
setColWidths(wb, 2, cols = 2:ncol(s2_out), widths = 15)

addWorksheet(wb, "Keyword matching vs MFD")

title3 <- paste0(
  "Comparison | Manual keyword matching vs. MetaFuncDecoder v1.0.0 ",
  "(peptidoglycanase excluded in both; broad matches retained in both, so ",
  "the comparison isolates the change of annotation tool)."
)
writeData(wb, 3, title3, startRow = 1, startCol = 1)
addStyle(wb, 3, title_style, rows = 1, cols = 1)
mergeCells(wb, 3, cols = 1:ncol(comp), rows = 1)

writeData(wb, 3, comp, startRow = 3, headerStyle = header_style)
nc <- nrow(comp)
addStyle(wb, 3, body_left, rows = 4:(3 + nc), cols = 1, gridExpand = TRUE)
addStyle(wb, 3, body_style, rows = 4:(3 + nc), cols = 2:ncol(comp),
         gridExpand = TRUE)
addStyle(wb, 3, bottom_border, rows = 3 + nc, cols = 1:ncol(comp),
         stack = TRUE)

setColWidths(wb, 3, cols = 1, widths = 38)
setColWidths(wb, 3, cols = 2:3, widths = 22)
setColWidths(wb, 3, cols = 4, widths = 12)
setColWidths(wb, 3, cols = 5, widths = 10)
setColWidths(wb, 3, cols = 6, widths = 18)

# Save
saveWorkbook(wb, file.path(results_dir, "03_Table_S_database_support.xlsx"),
             overwrite = TRUE)
cat("Publication-ready Excel saved to",
    file.path(results_dir, "03_Table_S_database_support.xlsx"), "\n")

mfd_write_session_info("03", results_dir)

cat("\nDone. All outputs written to", results_dir, "\n")
