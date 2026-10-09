#!/usr/bin/env Rscript
# =============================================================================
# run_all.R - runs every analysis script in order.
#
# Working directory must be the repository root, with data/ populated as
# described in README.md. Outputs are written to results/.
#
# 00_filter_10kb.R runs first: it builds data_10kb/ (contigs >=10 kb) from
# data/. The primary analysis scripts (02, 03, 04, 07) read data_10kb/; the
# contig-length sensitivity analysis (01) reads the full data/ (see README).
# =============================================================================

scripts <- c(
  "00_filter_10kb.R",
  "01_contig_length_stratification.R",
  "02_database_decomposition.R",
  "03_database_support.R",
  "04_amg_distribution.R",
  "05_jgi_study_metadata.R",
  "06_jgi_viral_contribution.R",
  "07_lysis_gene_decomposition.R"
)

for (s in scripts) {
  cat("\n==== ", s, " ====\n", sep = "")
  t0 <- Sys.time()
  source(s, echo = FALSE)
  cat("---- done in ", round(difftime(Sys.time(), t0, units = "secs")), "s\n",
      sep = "")
}

cat("\nAll scripts completed.\n")
