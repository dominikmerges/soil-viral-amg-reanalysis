#!/usr/bin/env Rscript
# _common.R - shared loading and filtering for all analysis scripts.
#
# Filtering chain applied to the MetaFuncDecoder output:
#   raw          5,760  all calls
#     - peptidoglycanase (lysis genes, not metabolic)
#   post_lysis   1,858  any annotation          0.13% of catalogue
#     - broad matches (empty subcategories: no enzyme resolved)
#   specific     1,169  enzyme-resolved         0.08% of catalogue
#
# `specific` is the analysis set for all reported results. `post_lysis` is
# used only for the annotated-at-all rate and the tool-comparison table in
# script 03, where filtering one arm would confound tool with filter.


suppressPackageStartupMessages({
  library(data.table)
})

# Catalogue-level constants (GSVA; Graham et al. 2024)
# Denominators for all "% of all viral genes" values.
GSVA_TOTAL_VIRAL_GENES <- 1238728L
GSVA_TOTAL_CONTIGS     <-   31344L
# Samples that carry >=1 viral gene: the source of the gene catalogue and the
# denominator distinct-sample count, NOT the full sample set. Equals the number
# of distinct sample IDs among the gene metadata.
GSVA_N_SAMPLES         <-    1093L
# Full soil-sample metadata (GSV Atlas File 5). Asserted against
# nrow(GSVA_sample_metadata_5.csv) in scripts 04/05 where that file is loaded.
GSVA_N_TOTAL_SAMPLES   <-    2953L
# Harmonized ecosystem bins used in the analysis (mapping in script 04 / Table
# S16). Cross-checked against the harmonized sample_meta in script 04.
GSVA_N_ECOSYSTEMS      <-      12L

# Presentation constants
CATEGORY_LABELS <- c(
  carbon_cycling   = "Carbon cycling",
  nitrogen_cycling = "Nitrogen cycling",
  antibiotics      = "AMR"
)

CATEGORY_COLOURS <- c(
  "Carbon cycling"   = "#2166AC",
  "Nitrogen cycling" = "#B2182B",
  "AMR"              = "#4DAF4A"
)

# CAZy indexes only carbohydrate-active enzymes, so other categories have a
# ceiling of 2 supporting databases (KEGG + Pfam).
CAZY_ELIGIBLE_CATEGORIES <- "carbon_cycling"

# Loading and filtering

#' Read the raw MetaFuncDecoder output.
#'
#' @param data_dir Directory holding metagenome_confidence_annotations.csv.
#' @return data.table of all MFD calls, unfiltered.
mfd_load_raw <- function(data_dir = "data") {
  f <- file.path(data_dir, "metagenome_confidence_annotations.csv")
  if (!file.exists(f)) {
    stop("MetaFuncDecoder output not found at '", f, "'.\n",
         "See README.md for how to obtain the GSVA inputs and regenerate ",
         "this file with MetaFuncDecoder v1.0.0 ",
         "(DOI:10.5281/zenodo.19009291).")
  }
  fread(f)
}

#' Drop peptidoglycanase (phage lysis) genes.
mfd_drop_lysis <- function(dt) {
  n_before <- nrow(dt)
  # Keep NA subcategories: they are broad (no enzyme resolved), not lysis genes.
  # A bare `!=` would drop them because NA != x is NA.
  out <- dt[is.na(subcategories) | subcategories != "peptidoglycanase"]
  attr(out, "n_lysis_dropped") <- n_before - nrow(out)
  out
}

#' Tag each annotation as specific (enzyme-resolved) or broad.
#'
#' MetaFuncDecoder assigns a subcategory only when a specific enzyme is
#' resolved. Broad matches (umbrella GO terms, non-specific text patterns)
#' carry an empty `subcategories` field.
mfd_add_evidence_type <- function(dt) {
  dt[, evidence_type := fifelse(subcategories == "" | is.na(subcategories),
                                "broad", "specific")]
  dt[]
}

#' Add per-database support flags and the supporting-database count.
#'
#' `supporting_databases` occasionally suffers comma-in-quoted-field parsing
#' artefacts, so presence is confirmed against the per-database ID columns.
mfd_add_db_flags <- function(dt) {
  dt[, has_cazy := grepl("CAZy", supporting_databases, fixed = TRUE) |
                   (cazy_families != "" & !is.na(cazy_families))]
  dt[, has_kegg := grepl("KEGG", supporting_databases, fixed = TRUE) |
                   (ko_terms != "" & !is.na(ko_terms))]
  dt[, has_pfam := grepl("Pfam", supporting_databases, fixed = TRUE) |
                   (pfam_ids != "" & !is.na(pfam_ids))]
  dt[, n_db := as.integer(has_cazy) + as.integer(has_kegg) + as.integer(has_pfam)]
  dt[]
}

#' Load annotations with the full canonical chain applied.
#'
#' @param data_dir Directory holding the MFD output.
#' @param set "specific" (default, n = 1,169), "post_lysis" (n = 1,858) or
#'   "raw" (n = 5,760).
#' @param verbose Print the filtering chain to stdout (default TRUE).
mfd_load_annotations <- function(data_dir = "data",
                                 set = c("specific", "post_lysis", "raw"),
                                 verbose = TRUE) {
  set <- match.arg(set)

  raw <- mfd_load_raw(data_dir)
  n_raw <- nrow(raw)

  post_lysis <- mfd_drop_lysis(raw)
  n_lysis <- attr(post_lysis, "n_lysis_dropped")
  post_lysis <- mfd_add_evidence_type(post_lysis)
  post_lysis <- mfd_add_db_flags(post_lysis)

  n_post <- nrow(post_lysis)
  n_broad <- post_lysis[evidence_type == "broad", .N]
  n_spec <- n_post - n_broad

  # Fail-fast guards on the documented GSV Atlas totals. Two input states are
  # sanctioned: the full QA/QC catalogue (read by 01, the contig-length
  # sensitivity analysis) and the >=10 kb primary catalogue (read by 02/03).
  # The loaded set must match one of them exactly; any other count means the
  # input or filtering logic moved silently, and the pipeline stops.
  n_carbon <- post_lysis[evidence_type == "specific" &
                           functional_category == "carbon_cycling",   .N]
  n_nitro  <- post_lysis[evidence_type == "specific" &
                           functional_category == "nitrogen_cycling", .N]
  n_amr    <- post_lysis[evidence_type == "specific" &
                           functional_category == "antibiotics",      .N]
  got      <- c(n_raw, n_post, n_spec, n_carbon, n_nitro, n_amr)
  full_set <- c(5760L, 1858L, 1169L, 1130L, 25L, 14L)   # all QA/QC contigs
  ge10_set <- c(5020L, 1698L, 1067L, 1032L, 22L, 13L)   # >=10 kb contigs
  stopifnot(identical(got, full_set) || identical(got, ge10_set))

  if (verbose) {
    cat("--- MetaFuncDecoder filtering chain -------------------------------\n")
    cat(sprintf("  raw MFD calls                        %7s\n",
                format(n_raw, big.mark = ",")))
    cat(sprintf("  - peptidoglycanase (lysis)           %7s\n",
                format(-n_lysis, big.mark = ",")))
    cat(sprintf("  = post-lysis ('annotated at all')    %7s   (%.2f%% of catalogue)\n",
                format(n_post, big.mark = ","),
                100 * n_post / GSVA_TOTAL_VIRAL_GENES))
    cat(sprintf("  - broad (no enzyme resolved)         %7s\n",
                format(-n_broad, big.mark = ",")))
    cat(sprintf("  = ANALYSIS SET ('enzyme-resolved')   %7s   (%.2f%% of catalogue)\n",
                format(n_spec, big.mark = ","),
                100 * n_spec / GSVA_TOTAL_VIRAL_GENES))
    cat(sprintf("  returning: %s\n", set))
    cat("-------------------------------------------------------------------\n")
  }

  out <- switch(set,
    raw        = mfd_add_db_flags(mfd_add_evidence_type(raw)),
    post_lysis = post_lysis,
    specific   = post_lysis[evidence_type == "specific"]
  )

  # Chain counts travel with the table so callers need not re-derive them.
  attr(out, "mfd_counts") <- list(
    raw = n_raw, lysis_dropped = n_lysis, post_lysis = n_post,
    broad = n_broad, specific = n_spec, set = set
  )
  out
}

#' Retrieve the filtering-chain counts attached by mfd_load_annotations().
mfd_counts <- function(dt) {
  x <- attr(dt, "mfd_counts")
  if (is.null(x)) stop("No mfd_counts attribute; use mfd_load_annotations().")
  x
}

#' Ensure the results/ and results/figures/ directories exist.
mfd_init_dirs <- function(results_dir = "results") {
  fig_dir <- file.path(results_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(list(results_dir = results_dir, fig_dir = fig_dir))
}

#' Write sessionInfo() for a script.
mfd_write_session_info <- function(script_tag, results_dir = "results") {
  writeLines(capture.output(sessionInfo()),
             file.path(results_dir, paste0(script_tag, "_sessionInfo.txt")))
}
