#!/usr/bin/env Rscript
# 02_database_decomposition.R
#
# Decomposes annotations by supporting database (CAZy, KEGG, Pfam) and by
# evidence specificity (specific = enzyme resolved, broad = no subcategory),
# then re-tests category proportions under three filters:
#   all annotations / specific only / KEGG+Pfam only (drops CAZy-only hits).
#
# CAZy indexes only carbohydrate-active enzymes, so it can inflate the carbon
# share. These tests establish whether that drives the observed composition.
#
# In:  data/metagenome_confidence_annotations.csv
# Out: results/02_database_decomposition.tsv
#      results/02_specific_vs_broad.tsv
#      results/02_enrichment_specific_only.tsv
#      results/02_enrichment_kegg_pfam_only.tsv
#      results/02_supplementary_table_S_database_bias.tsv
#      results/02_Table_S_database_bias.xlsx
#      results/figures/02_database_coverage_bias.{pdf,png,tiff}
#      results/02_sessionInfo.txt

source("_common.R")

library(tidyverse)
library(data.table)

set.seed(42)

data_dir    <- "data_10kb"
dirs        <- mfd_init_dirs()
results_dir <- dirs$results_dir
fig_dir     <- dirs$fig_dir

# Post-lysis set: this script reports the specific/broad split, so it needs the
# broad matches that the analysis set drops. evidence_type and the has_* /
# n_db flags are attached by mfd_load_annotations().

annot <- mfd_load_annotations(data_dir, set = "post_lysis")

print(annot[, .N, by = functional_category])

spec_broad <- annot[, .N, by = .(functional_category, evidence_type)]
spec_broad <- dcast(spec_broad, functional_category ~ evidence_type,
                    value.var = "N", fill = 0)
cat("\n--- Specific vs broad matches ---\n")
print(spec_broad)
fwrite(spec_broad, file.path(results_dir, "02_specific_vs_broad.tsv"), sep = "\t")

# Database support summary per category
db_decomp <- annot[, .(
  n_total      = .N,
  n_cazy       = sum(has_cazy),
  n_kegg       = sum(has_kegg),
  n_pfam       = sum(has_pfam),
  n_cazy_only  = sum(has_cazy & !has_kegg & !has_pfam),
  n_kegg_only  = sum(!has_cazy & has_kegg & !has_pfam),
  n_pfam_only  = sum(!has_cazy & !has_kegg & has_pfam),
  n_kegg_or_pfam = sum(has_kegg | has_pfam),
  n_multi_db   = sum((has_cazy + has_kegg + has_pfam) >= 2)
), by = functional_category]

cat("\n--- Database decomposition ---\n")
print(db_decomp)
fwrite(db_decomp, file.path(results_dir, "02_database_decomposition.tsv"), sep = "\t")

# Compare category proportions across three filters (all, specific, KEGG/Pfam)
# using chi-squared goodness-of-fit: do proportions change when we remove
# broad/CAZy-only hits? If not, the bias does not drive the result.

annot_specific <- annot[evidence_type == "specific"]

cat("\n--- Specific-only counts ---\n")
counts_specific <- annot_specific[, .N, by = functional_category][order(-N)]
print(counts_specific)

# Chi-squared on mutually exclusive sets: does the category composition of the
# retained (specific) set differ from the excluded (broad) set? A goodness-of-fit
# test against proportions estimated from the superset that contains the specific
# genes would not be independent, so the two arms are compared directly instead.
tab_spec <- dcast(annot[, .N, by = .(functional_category, evidence_type)],
                  functional_category ~ evidence_type, value.var = "N", fill = 0)
chi_specific <- chisq.test(as.matrix(tab_spec[, .(specific, broad)]))
cat("\n--- Chi-squared: specific vs broad category composition ---\n")
print(chi_specific)

enrich_specific <- data.table(
  filter = "specific_only",
  n_carbon   = counts_specific[functional_category == "carbon_cycling", N],
  n_nitrogen = counts_specific[functional_category == "nitrogen_cycling", N],
  n_amr      = counts_specific[functional_category == "antibiotics", N],
  n_total    = nrow(annot_specific),
  pct_carbon = round(100 * counts_specific[functional_category == "carbon_cycling", N] /
                       nrow(annot_specific), 1),
  chi_sq_stat = round(chi_specific$statistic, 2),
  chi_sq_p    = chi_specific$p.value
)
fwrite(enrich_specific,
       file.path(results_dir, "02_enrichment_specific_only.tsv"), sep = "\t")

# Removes genes where the ONLY evidence is CAZy — the database-neutral test.

annot_kegg_pfam <- annot[has_kegg | has_pfam]

cat("\n--- KEGG/Pfam-supported counts ---\n")
counts_kp <- annot_kegg_pfam[, .N, by = functional_category][order(-N)]
print(counts_kp)

# Mutually exclusive arms again: KEGG/Pfam-supported vs CAZy-only genes partition
# the post-lysis set (every gene has at least one supporting database).
annot[, db_arm := fifelse(has_kegg | has_pfam, "kegg_pfam", "cazy_only")]
tab_kp <- dcast(annot[, .N, by = .(functional_category, db_arm)],
                functional_category ~ db_arm, value.var = "N", fill = 0)
chi_kp <- chisq.test(as.matrix(tab_kp[, .(kegg_pfam, cazy_only)]))
cat("\n--- Chi-squared: KEGG/Pfam vs CAZy-only category composition ---\n")
print(chi_kp)

enrich_kp <- data.table(
  filter = "kegg_pfam_only",
  n_carbon   = counts_kp[functional_category == "carbon_cycling", N],
  n_nitrogen = counts_kp[functional_category == "nitrogen_cycling", N],
  n_amr      = counts_kp[functional_category == "antibiotics", N],
  n_total    = nrow(annot_kegg_pfam),
  pct_carbon = round(100 * counts_kp[functional_category == "carbon_cycling", N] /
                       nrow(annot_kegg_pfam), 1),
  chi_sq_stat = round(chi_kp$statistic, 2),
  chi_sq_p    = chi_kp$p.value
)
fwrite(enrich_kp,
       file.path(results_dir, "02_enrichment_kegg_pfam_only.tsv"), sep = "\t")

# Rows: functional_category x evidence_type
# Columns: database counts, totals

supp_table <- annot[, .(
  n_genes        = .N,
  n_cazy         = sum(has_cazy),
  n_kegg         = sum(has_kegg),
  n_pfam         = sum(has_pfam),
  n_cazy_only    = sum(has_cazy & !has_kegg & !has_pfam),
  n_kegg_or_pfam = sum(has_kegg | has_pfam),
  n_multi_db     = sum((has_cazy + has_kegg + has_pfam) >= 2),
  pct_1db        = round(100 * sum(n_db == 1) / .N, 1),
  pct_2db        = round(100 * sum(n_db == 2) / .N, 1),
  pct_3db        = round(100 * sum(n_db == 3) / .N, 1)
), by = .(functional_category, evidence_type)]

supp_table <- supp_table[order(functional_category, -n_genes)]

cat("\n--- Supplementary table ---\n")
print(supp_table)
fwrite(supp_table,
       file.path(results_dir, "02_supplementary_table_S_database_bias.tsv"),
       sep = "\t")

# One row per category x subcategory. Specific matches only.

subcat_table <- annot_specific[, .(
  n_genes      = .N,
  n_cazy       = sum(has_cazy),
  n_kegg       = sum(has_kegg),
  n_pfam       = sum(has_pfam),
  n_multi_db   = sum((has_cazy + has_kegg + has_pfam) >= 2),
  n_1db        = sum(n_db == 1),
  n_2db        = sum(n_db == 2),
  n_3db        = sum(n_db == 3)
), by = .(functional_category, subcategories)]

subcat_table <- subcat_table[order(functional_category, -n_genes)]

cat("\n--- Subcategory breakdown (specific matches) ---\n")
print(subcat_table)

library(openxlsx)

# Prettify category and subcategory labels
prettify <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)  # title case
  # Fix specific terms
  x <- gsub("Amr", "AMR", x)
  x <- gsub("Dnra", "DNRA", x)
  x <- gsub("Beta ", "\u03B2-", x)       # β-
  x <- gsub("Alpha ", "\u03B1-", x)      # α-
  x <- gsub("Kegg", "KEGG", x)
  x <- gsub("Pfam", "Pfam", x)           # already correct
  x <- gsub("Cazy", "CAZy", x)
  x <- gsub("Db", "database", x)
  x
}

# Sheet 1: Subcategory breakdown
s1 <- copy(subcat_table)
s1[, `Functional category` := prettify(functional_category)]
s1[, Subcategory := prettify(subcategories)]
# Multi-subcategory entries: semicolon -> " / "
s1[, Subcategory := gsub(";", " / ", Subcategory)]

# GH16 recategorization: in soil, GH16 is β-1,3-glucanase (not a marine-algal
# substrate) and the algal_polysaccharide-only rows are alginate lyases (PL6/7/14).
# Rename the three affected subcategory labels to match Fig 1b / Table S7 and the
# manuscript text.
beta <- "β"
s1[Subcategory == paste0("Algal Polysaccharide / ", beta, "-Glucanase"),
   Subcategory := paste0(beta, "-1,3-Glucanase (GH16)")]
s1[Subcategory == "Algal Polysaccharide",
   Subcategory := "Alginate/polysaccharide lyase"]
s1[Subcategory == paste0(beta, "-Glucanase"),
   Subcategory := paste0(beta, "-1,3-Glucanase (other GH)")]

s1_out <- s1[, .(
  `Functional category` = `Functional category`,
  Subcategory,
  `No. genes` = n_genes,
  CAZy   = n_cazy,
  KEGG   = n_kegg,
  Pfam   = n_pfam,
  `Multi-database` = n_multi_db,
  `1 database`  = n_1db,
  `2 databases` = n_2db,
  `3 databases` = n_3db
)]

# Sheet 2: Database decomposition
s2 <- copy(supp_table)
s2[, `Functional category` := prettify(functional_category)]
s2[, `Evidence type` := fifelse(evidence_type == "specific", "Specific", "Broad")]

s2_out <- s2[, .(
  `Functional category` = `Functional category`,
  `Evidence type`,
  `No. genes` = n_genes,
  CAZy       = n_cazy,
  KEGG       = n_kegg,
  Pfam       = n_pfam,
  `CAZy only` = n_cazy_only,
  `KEGG or Pfam` = n_kegg_or_pfam,
  `Multi-database` = n_multi_db,
  `1 database (%)`  = pct_1db,
  `2 databases (%)` = pct_2db,
  `3 databases (%)` = pct_3db
)]

# Build workbook with formatting
wb <- createWorkbook()

# Styles
title_style <- createStyle(fontSize = 11, textDecoration = "bold")
header_style <- createStyle(
  fontSize = 10, textDecoration = "bold",
  border = "TopBottom", borderStyle = "thin",
  halign = "center", wrapText = TRUE
)
body_style <- createStyle(fontSize = 10, halign = "center")
body_left  <- createStyle(fontSize = 10, halign = "left")
bottom_border <- createStyle(border = "bottom", borderStyle = "thin")
category_style <- createStyle(fontSize = 10, halign = "left", textDecoration = "italic")

# Write Sheet 1
addWorksheet(wb, "Subcategory breakdown")

title1 <- "Table S3 | Subcategory breakdown of specific functional annotations by database support."
writeData(wb, 1, title1, startRow = 1, startCol = 1)
addStyle(wb, 1, title_style, rows = 1, cols = 1)
mergeCells(wb, 1, cols = 1:ncol(s1_out), rows = 1)

writeData(wb, 1, s1_out, startRow = 3, headerStyle = header_style)

# Body formatting
n1 <- nrow(s1_out)
addStyle(wb, 1, body_left, rows = 4:(3 + n1), cols = 1, gridExpand = TRUE)
addStyle(wb, 1, body_left, rows = 4:(3 + n1), cols = 2, gridExpand = TRUE)
addStyle(wb, 1, body_style, rows = 4:(3 + n1), cols = 3:ncol(s1_out), gridExpand = TRUE)

# Italic for functional category column
addStyle(wb, 1, category_style, rows = 4:(3 + n1), cols = 1, gridExpand = TRUE)

# Bottom border on last row
addStyle(wb, 1, bottom_border, rows = 3 + n1, cols = 1:ncol(s1_out), stack = TRUE)

# Column widths
setColWidths(wb, 1, cols = 1, widths = 20)
setColWidths(wb, 1, cols = 2, widths = 42)
setColWidths(wb, 1, cols = 3:ncol(s1_out), widths = 14)

# Merge repeated category labels
cats1 <- s1_out$`Functional category`
run_start <- 1
for (i in 2:length(cats1)) {
  if (cats1[i] != cats1[run_start]) {
    if (i - run_start > 1) {
      mergeCells(wb, 1, cols = 1, rows = (run_start + 3):(i + 2))
    }
    run_start <- i
  }
}
if (length(cats1) - run_start >= 1) {
  mergeCells(wb, 1, cols = 1, rows = (run_start + 3):(length(cats1) + 3))
}

# Write Sheet 2
addWorksheet(wb, "Database decomposition")

title2 <- "Table S1 | Database decomposition of functional annotations by evidence specificity."
writeData(wb, 2, title2, startRow = 1, startCol = 1)
addStyle(wb, 2, title_style, rows = 1, cols = 1)
mergeCells(wb, 2, cols = 1:ncol(s2_out), rows = 1)

writeData(wb, 2, s2_out, startRow = 3, headerStyle = header_style)

n2 <- nrow(s2_out)
addStyle(wb, 2, body_left, rows = 4:(3 + n2), cols = 1:2, gridExpand = TRUE)
addStyle(wb, 2, body_style, rows = 4:(3 + n2), cols = 3:ncol(s2_out), gridExpand = TRUE)
addStyle(wb, 2, category_style, rows = 4:(3 + n2), cols = 1, gridExpand = TRUE)
addStyle(wb, 2, bottom_border, rows = 3 + n2, cols = 1:ncol(s2_out), stack = TRUE)

setColWidths(wb, 2, cols = 1, widths = 20)
setColWidths(wb, 2, cols = 2, widths = 16)
setColWidths(wb, 2, cols = 3:ncol(s2_out), widths = 14)

# Merge repeated category labels
cats2 <- s2_out$`Functional category`
run_start <- 1
for (i in 2:length(cats2)) {
  if (cats2[i] != cats2[run_start]) {
    if (i - run_start > 1) {
      mergeCells(wb, 2, cols = 1, rows = (run_start + 3):(i + 2))
    }
    run_start <- i
  }
}
if (length(cats2) - run_start >= 1) {
  mergeCells(wb, 2, cols = 1, rows = (run_start + 3):(length(cats2) + 3))
}

# Save
saveWorkbook(wb, file.path(results_dir, "02_Table_S_database_bias.xlsx"),
             overwrite = TRUE)
cat("Publication-ready Excel saved to",
    file.path(results_dir, "02_Table_S_database_bias.xlsx"), "\n")

# Panel A: stacked bar — database support per category (specific matches only).
# The three groups partition each category: single-database CAZy, single-database
# KEGG/Pfam (n_db == 1, else a KEGG+Pfam gene would also count in Multi_database),
# and multi-database (n_db >= 2).
panel_a_data <- annot_specific[, .(
  CAZy_only      = sum(n_db == 1 & has_cazy),
  KEGG_or_Pfam   = sum(n_db == 1 & !has_cazy & (has_kegg | has_pfam)),
  Multi_database = sum(n_db >= 2)
), by = functional_category]

# The three plotted groups must sum to the functional-category total.
panel_a_check <- merge(panel_a_data,
                       annot_specific[, .(n_cat = .N), by = functional_category],
                       by = "functional_category")
stopifnot(all(panel_a_check[, CAZy_only + KEGG_or_Pfam + Multi_database] ==
                panel_a_check$n_cat))

panel_a_long <- melt(panel_a_data, id.vars = "functional_category",
                     variable.name = "db_support", value.name = "n")

cat_labels <- c(carbon_cycling = "Carbon cycling",
                nitrogen_cycling = "Nitrogen cycling",
                antibiotics = "AMR")
panel_a_long[, category_label := factor(cat_labels[functional_category],
                                         levels = cat_labels)]
panel_a_long[, db_support := factor(db_support,
                                     levels = c("Multi_database",
                                                "KEGG_or_Pfam",
                                                "CAZy_only"))]

# Data-label positions matching position_stack(): the first factor level
# (Multi_database) is drawn at the top, so cumulate from the bottom-most level.
lab_dat <- copy(panel_a_long)[n > 0]
lab_dat[, stack_rank := as.integer(db_support)]
setorder(lab_dat, category_label, -stack_rank)
lab_dat[, y_center := cumsum(n) - n / 2, by = category_label]
lab_in  <- lab_dat[functional_category == "carbon_cycling"]
lab_out <- lab_dat[functional_category != "carbon_cycling"]
lab_out[, bar_total := sum(n), by = category_label]
lab_out[, rank_desc := frank(-n, ties.method = "first"), by = category_label]
# pmax() floors the in-slice label so a very short bottom slice (AMR: 11) is
# not clipped by the x-axis.
lab_out[, y_label := fifelse(rank_desc == 1L, pmax(y_center, 16), bar_total + 45)]

p_a <- ggplot(panel_a_long, aes(x = category_label, y = n, fill = db_support)) +
  geom_col(width = 0.6) +
  geom_text(data = lab_in,
            aes(x = category_label, y = y_center, label = n),
            inherit.aes = FALSE, size = 3) +
  geom_text(data = lab_out,
            aes(x = category_label, y = y_label, label = n),
            inherit.aes = FALSE, size = 3) +
  scale_fill_manual(
    values = c("Multi_database" = "#2166AC",
               "KEGG_or_Pfam"   = "#67A9CF",
               "CAZy_only"      = "#D1E5F0"),
    labels = c("Multi_database" = "≥2 databases",
               "KEGG_or_Pfam"   = "Single database: KEGG or Pfam",
               "CAZy_only"      = "CAZy only")
  ) +
  # Stack the legend to one column so every label prints in full: the long
  # "Single database: KEGG or Pfam" entry overflows a single horizontal row
  # under the narrower left panel and clips the "CAZy only" swatch.
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL,
       y = "Number of genes (specific matches)",
       fill = "Database support",
       tag = "A") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 14))

# Panel B: comparison of all vs specific vs KEGG/Pfam-only proportions
filter_labels <- c("All annotations", "Specific only", "KEGG/Pfam only")

comparison_data <- rbind(
  annot[,           .(filter = "All annotations", n = .N),
        by = .(category = functional_category)],
  annot_specific[,  .(filter = "Specific only", n = .N),
        by = .(category = functional_category)],
  annot_kegg_pfam[, .(filter = "KEGG/Pfam only", n = .N),
        by = .(category = functional_category)]
)[, .(filter, category, n)]

comparison_data[, total := sum(n), by = filter]
comparison_data[, pct := 100 * n / total]
comparison_data[, category_label := factor(cat_labels[category], levels = cat_labels)]
comparison_data[, filter := factor(filter, levels = filter_labels)]

p_b <- ggplot(comparison_data, aes(x = filter, y = pct, fill = category_label)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = paste0("n=", n)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Carbon cycling"    = "#2166AC",
                                "Nitrogen cycling" = "#B2182B",
                                "AMR"              = "#4DAF4A")) +
  labs(x = "Annotation filter",
       y = "Proportion of annotated genes (%)",
       fill = NULL,
       tag = "B") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 15, hjust = 1))

pdf(file.path(fig_dir, "02_database_coverage_bias.pdf"),
    width = 10, height = 5)
gridExtra::grid.arrange(p_a, p_b, ncol = 2, widths = c(1, 1.3))
dev.off()

png(file.path(fig_dir, "02_database_coverage_bias.png"),
    width = 10, height = 5, units = "in", res = 300)
gridExtra::grid.arrange(p_a, p_b, ncol = 2, widths = c(1, 1.3))
dev.off()

# TIFF at 300 dpi (AEM preferred format for figure upload)
tiff(file.path(fig_dir, "02_database_coverage_bias.tiff"),
     width = 10, height = 5, units = "in", res = 300,
     compression = "lzw", type = "cairo")
gridExtra::grid.arrange(p_a, p_b, ncol = 2, widths = c(1, 1.3))
dev.off()

cat("\nFigures saved to", fig_dir, "\n")

mfd_write_session_info("02", results_dir)

cat("\nDone. All outputs written to", results_dir, "\n")
