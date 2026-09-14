#!/usr/bin/env Rscript
# 01_contig_length_stratification.R
#
# Stratifies annotation rate and functional composition by contig-length bin
# (1-5 kb, 5-10 kb, >=10 kb). 
#
# Headline results use the enzyme-resolved set; the post-lysis set is carried
# through as a sensitivity analysis in the same output files.
#
# In:  data/GSVA_soil_viruses_genome_metadata_2.tsv.gz  (contig_id, length)
#      data/GSVA_soil_viruses_gene_metadata_4.tsv.gz    (gene_id -> contig_id)
#      data/metagenome_confidence_annotations.csv
# Out: results/01_contig_bin_summary.tsv
#      results/01_annotation_rate_by_bin.tsv
#      results/01_functional_composition_by_bin.tsv
#      results/01_confidence_by_bin.tsv
#      results/01_enrichment_tests.tsv
#      results/figures/01_contig_length_stratification.pdf
#      results/01_sessionInfo.txt

source("_common.R")

library(tidyverse)
library(data.table)

set.seed(42)

data_dir    <- "data"
dirs        <- mfd_init_dirs()
results_dir <- dirs$results_dir
fig_dir     <- dirs$fig_dir

# Genome metadata: contig_id + contig_length
genome_meta <- fread(
  file.path(data_dir, "GSVA_soil_viruses_genome_metadata_2.tsv.gz"),
  select = c("contig_id", "contig_length")
)

# Gene metadata: gene_id -> contig_id (1.4M rows, keep lean)
gene_meta <- fread(
  file.path(data_dir, "GSVA_soil_viruses_gene_metadata_4.tsv.gz"),
  select = c("gene_id", "contig_id")
)

# Peptidoglycanase already removed; enzyme-resolved subset taken after the join.
annotations <- mfd_load_annotations(data_dir, set = "post_lysis")

cat("Loaded:", nrow(genome_meta), "contigs,",
    nrow(gene_meta), "genes,",
    nrow(annotations), "annotated genes (post-lysis)\n")

genome_meta[, length_bin := fcase(
  contig_length >= 1000  & contig_length < 5000,  "1-5 kb",
  contig_length >= 5000  & contig_length < 10000, "5-10 kb",
  contig_length >= 10000,                         ">=10 kb"
)]

# Order factor for consistent plotting
bin_levels <- c("1-5 kb", "5-10 kb", ">=10 kb")
genome_meta[, length_bin := factor(length_bin, levels = bin_levels)]

stopifnot(!anyDuplicated(genome_meta$contig_id))
# Inner join on contig_id: every gene's contig must be in genome_meta, so guard
# against silently dropping genes from the length-bin denominators.
n_genes_before <- nrow(gene_meta)
gene_meta <- merge(gene_meta, genome_meta, by = "contig_id", all.x = FALSE)
stopifnot(nrow(gene_meta) == n_genes_before)

# gene_id is the unique join key, so the join must not multiply annotation rows.
stopifnot(!anyDuplicated(gene_meta$gene_id))
n_annot_in <- nrow(annotations)
annot_any <- merge(annotations, gene_meta[, .(gene_id, contig_length, length_bin)],
                   by = "gene_id", all.x = TRUE)
stopifnot(nrow(annot_any) == n_annot_in)

# Every annotated gene is a catalogue gene, so complete matching is expected.
stopifnot(sum(is.na(annot_any$length_bin)) == 0)

cat("Annotated genes with contig info (post-lysis):", nrow(annot_any), "\n")

annot <- annot_any[evidence_type == "specific"]
cat("Analysis set (enzyme-resolved):", nrow(annot), "of", nrow(annot_any),
    "post-lysis annotations (",
    nrow(annot_any) - nrow(annot), "broad matches removed )\n")

# Both rates reported per bin: any annotation, and enzyme-resolved.

bin_summary <- gene_meta[, .(
  n_contigs            = uniqueN(contig_id),
  n_total_genes        = .N,
  n_annotated_any      = sum(gene_id %in% annot_any$gene_id),
  n_annotated_specific = sum(gene_id %in% annot$gene_id)
), by = length_bin][order(length_bin)]

stopifnot(all(bin_summary$n_total_genes > 0))
bin_summary[, annotation_rate_any      := n_annotated_any / n_total_genes]
bin_summary[, annotation_rate_specific := n_annotated_specific / n_total_genes]

cat("\n--- Contig-length bin summary ---\n")
print(bin_summary)

fwrite(bin_summary, file.path(results_dir, "01_contig_bin_summary.tsv"), sep = "\t")

# Chi-squared: does annotation rate differ across bins? Run on both denominators.

test_rate_by_bin <- function(n_annotated, n_total) {
  m <- cbind(annotated = n_annotated, unannotated = n_total - n_annotated)
  # Bins with no genes at all carry no information and break the test.
  m <- m[n_total > 0, , drop = FALSE]
  chisq.test(m)
}

chi_rate_any  <- test_rate_by_bin(bin_summary$n_annotated_any,
                                  bin_summary$n_total_genes)
chi_rate_spec <- test_rate_by_bin(bin_summary$n_annotated_specific,
                                  bin_summary$n_total_genes)

cat("\n--- Annotation rate chi-squared (any annotation) ---\n")
print(chi_rate_any)
cat("\n--- Annotation rate chi-squared (enzyme-resolved) ---\n")
print(chi_rate_spec)

annotation_rate_out <- copy(bin_summary)
annotation_rate_out[, chi_sq_p_any      := chi_rate_any$p.value]
annotation_rate_out[, chi_sq_p_specific := chi_rate_spec$p.value]

fwrite(annotation_rate_out,
       file.path(results_dir, "01_annotation_rate_by_bin.tsv"), sep = "\t")

composition_by_bin <- function(dt) {
  x <- dt[, .N, by = .(length_bin, functional_category)]
  dcast(x, length_bin ~ functional_category, value.var = "N", fill = 0)
}

func_comp     <- composition_by_bin(annot)
func_comp_any <- composition_by_bin(annot_any)

func_comp[,     analysis_set := "specific"]
func_comp_any[, analysis_set := "post_lysis"]

cat("\n--- Functional composition by bin (enzyme-resolved) ---\n")
print(func_comp)
cat("\n--- Functional composition by bin (post-lysis, sensitivity) ---\n")
print(func_comp_any)

fwrite(rbind(func_comp, func_comp_any, fill = TRUE),
       file.path(results_dir, "01_functional_composition_by_bin.tsv"), sep = "\t")

# For each bin: Fisher's exact test — is carbon_cycling enriched relative to
# the other categories combined?

enrichment_by_bin <- function(dt) {
  overall_carbon <- sum(dt$functional_category == "carbon_cycling")
  overall_other  <- sum(dt$functional_category != "carbon_cycling")

  dt[, {
    n_carbon <- sum(functional_category == "carbon_cycling")
    n_other  <- sum(functional_category != "carbon_cycling")

    # Carbon vs other in this bin against the OTHER bins only, so the background
    # excludes the bin under test (overall_carbon - n_carbon).
    fisher_mat <- matrix(c(n_carbon, n_other,
                           overall_carbon - n_carbon,
                           overall_other  - n_other),
                         nrow = 2)
    fisher_res <- fisher.test(fisher_mat)

    .(n_carbon      = n_carbon,
      n_other       = n_other,
      pct_carbon    = round(100 * n_carbon / (n_carbon + n_other), 1),
      fisher_OR     = round(fisher_res$estimate, 3),
      fisher_CI_lo  = round(fisher_res$conf.int[1], 3),
      fisher_CI_hi  = round(fisher_res$conf.int[2], 3),
      fisher_p      = fisher_res$p.value)
  }, by = length_bin][order(length_bin)]
}

enrichment_tests     <- enrichment_by_bin(annot)
enrichment_tests_any <- enrichment_by_bin(annot_any)

cat("\n--- Carbon-cycling enrichment per bin, enzyme-resolved (Fisher) ---\n")
print(enrichment_tests)
cat("\n--- Carbon-cycling enrichment per bin, post-lysis sensitivity ---\n")
print(enrichment_tests_any)

fwrite(rbind(
  cbind(analysis_set = "specific",   enrichment_tests),
  cbind(analysis_set = "post_lysis", enrichment_tests_any)
), file.path(results_dir, "01_enrichment_tests.tsv"), sep = "\t")

# Number of supporting databases, not the high/medium/low labels: the ceiling
# differs by category (see 03_database_support.R).

conf_by_bin <- annot[, .N, by = .(length_bin, functional_category, n_db)]
conf_by_bin <- conf_by_bin[order(length_bin, functional_category, n_db)]

cat("\n--- Database support distribution by bin and category ---\n")
print(conf_by_bin)

fwrite(conf_by_bin,
       file.path(results_dir, "01_confidence_by_bin.tsv"), sep = "\t")

# Prepare annotation data for plotting
plot_data <- annot[, .N, by = .(length_bin, functional_category)]
plot_data[, pct := 100 * N / sum(N), by = length_bin]

plot_data[, category_label := factor(CATEGORY_LABELS[functional_category],
                                      levels = CATEGORY_LABELS)]

# Panel A: functional category proportions per bin
p_a <- ggplot(plot_data, aes(x = length_bin, y = pct, fill = category_label)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = paste0("n=", N)),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = CATEGORY_COLOURS) +
  labs(x = "Contig-length bin",
       y = "Proportion of enzyme-resolved annotations (%)",
       fill = NULL,
       tag = "A") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 14))

# Panel B: forest plot of ORs for carbon enrichment
or_data <- enrichment_tests[, .(length_bin, fisher_OR, fisher_CI_lo,
                                 fisher_CI_hi, fisher_p)]
or_data[, sig_label := ifelse(fisher_p < 0.05, "*", "n.s.")]

p_b <- ggplot(or_data, aes(x = length_bin, y = fisher_OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = fisher_CI_lo, ymax = fisher_CI_hi),
                  size = 0.7, colour = "#2166AC", fatten = 3) +
  geom_text(aes(label = sig_label, y = fisher_CI_hi + 0.2),
            size = 3.5, colour = "grey30") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Contig-length bin",
       y = "Odds ratio\n(carbon cycling vs. other categories)",
       tag = "B") +
  theme_bw(base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 14))

# Combine panels
pdf(file.path(fig_dir, "01_contig_length_stratification.pdf"),
    width = 9, height = 4.5)
gridExtra::grid.arrange(p_a, p_b, ncol = 2, widths = c(1.2, 1))
dev.off()

png(file.path(fig_dir, "01_contig_length_stratification.png"),
    width = 9, height = 4.5, units = "in", res = 300)
gridExtra::grid.arrange(p_a, p_b, ncol = 2, widths = c(1.2, 1))
dev.off()

cat("\nFigures saved to", fig_dir, "\n")

mfd_write_session_info("01", results_dir)

cat("\nDone. All outputs written to", results_dir, "\n")
