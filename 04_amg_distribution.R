#!/usr/bin/env Rscript
# 04_amg_distribution.R
#
# Distribution of annotated AMGs across soil ecosystems and functional groups.
# Builds all of Figure 1:
#   A  AMG rate per ecosystem, normalized per 10,000 viral genes
#   B  carbon-cycling AMGs by enzyme group
#   C  nitrogen-cycling AMGs by mechanism
#   D  AMR AMGs by mechanism
# plus the composite and the ecosystem supplement workbook.
#
# Panel A is produced in two versions: all ecosystem bins (supplement) and
# bins with n >= 30 samples plus "Other" (main figure).
#
# In:  data/GSVA_sample_metadata_5.csv                (sample -> ecosystem)
#      data/GSVA_soil_viruses_gene_metadata_4.tsv.gz  (gene_id -> contig_id)
#      data/metagenome_confidence_annotations.csv
# Out: results/04_ecosystem_summary.tsv
#      results/04_ecosystem_normalized.tsv
#      results/04_ecosystem_sample_counts.tsv
#      results/04_ecosystem_supplement.xlsx
#      results/figures/04_Figure1_composite.pdf
#      results/04_sessionInfo.txt

library(tidyverse)
library(data.table)

set.seed(42)
source("_common.R")

data_dir    <- "data_10kb"
results_dir <- "results"
fig_dir     <- file.path(results_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Sample metadata: IMG.taxon_oid -> GOLD ecosystem columns
sample_meta <- fread(file.path(data_dir, "GSVA_sample_metadata_5.csv"),
                     select = c("IMG.taxon_oid",
                                "GOLD.Ecosystem.Type",
                                "GOLD.Specific.Ecosystem"))
# File 5 is the full soil-sample metadata; guard the documented total.
stopifnot(nrow(sample_meta) == GSVA_N_TOTAL_SAMPLES)

# Gene metadata: gene_id -> contig_id (need gene_id to extract sample prefix)
# pfam/cazyme/kegg_ortholog carried through for the broad (any-hit) annotation
# rate reported alongside the AMG-specific rate.
gene_meta <- fread(file.path(data_dir, "GSVA_soil_viruses_gene_metadata_4.tsv.gz"),
                   select = c("gene_id", "contig_id", "pfam", "cazyme", "kegg_ortholog"))

annotations <- fread(file.path(data_dir, "metagenome_confidence_annotations.csv"))

cat("Loaded:", nrow(sample_meta), "samples,",
    nrow(gene_meta), "genes,",
    nrow(annotations), "annotated genes\n")

# gene_id format: "2124908025.a:MBSR2a_Contig_1912_60"
# prefix before "." is the IMG.taxon_oid

gene_meta[, taxon_oid := bit64::as.integer64(sub("\\..*", "", gene_id))]
# The taxon prefix is the join key to sample metadata; a missing parse would
# silently unmatch the gene downstream, so require every gene_id to yield one.
stopifnot(!any(is.na(gene_meta$taxon_oid)))

# Map Ecosystem.Type + Specific.Ecosystem -> clean ecosystem bin

sample_meta[, ecosystem := fcase(
  # Plant-associated (Rhizosphere, Roots, Rhizoplane)
  GOLD.Ecosystem.Type %in% c("Rhizosphere", "Roots", "Rhizoplane"),
    "Plant-associated",

  # Peat
  GOLD.Ecosystem.Type == "Peat",
    "Peat",

  # Volcanic
  GOLD.Ecosystem.Type == "Volcanic",
    "Volcanic",

  # Marine/Salt marsh
  GOLD.Ecosystem.Type == "Marine",
    "Salt marsh",

  # Soil types: bin by Specific.Ecosystem
  GOLD.Specific.Ecosystem %in% c("Agricultural land", "Agricultural",
                                   "Agricultural soil"),
    "Agricultural",

  GOLD.Specific.Ecosystem %in% c("Forest Soil", "Forest soil"),
    "Forest",

  GOLD.Specific.Ecosystem == "Grasslands",
    "Grassland",

  GOLD.Specific.Ecosystem == "Permafrost",
    "Permafrost",

  GOLD.Specific.Ecosystem == "Desert",
    "Desert",

  GOLD.Specific.Ecosystem == "Tropical rainforest",
    "Tropical rainforest",

  GOLD.Specific.Ecosystem == "Shrubland",
    "Shrubland",

  default = "Unclassified"
)]

cat("\n--- Ecosystem bin counts (samples) ---\n")
print(sample_meta[, .N, by = ecosystem][order(-N)])

stopifnot(uniqueN(sample_meta$ecosystem) == GSVA_N_ECOSYSTEMS)

stopifnot(!anyDuplicated(sample_meta$IMG.taxon_oid))
gene_eco <- merge(gene_meta,
                  sample_meta[, .(taxon_oid = IMG.taxon_oid, ecosystem)],
                  by = "taxon_oid", all.x = FALSE)
stopifnot(!anyDuplicated(gene_eco$gene_id))

stopifnot(nrow(gene_eco) == GSVA_TOTAL_VIRAL_GENES)

cat("\nGenes matched to an ecosystem:", nrow(gene_eco), "of", nrow(gene_meta), "\n")

# Enzyme-resolved set: excludes peptidoglycanase and broad (empty subcategory).

annot <- annotations[subcategories != "peptidoglycanase" & subcategories != ""]
cat("Annotations after peptidoglycanase + broad-match exclusion:", nrow(annot),
    "(", sum(annotations$subcategories != "peptidoglycanase") - nrow(annot),
    "broad matches removed)\n")

# Add ecosystem to annotations via gene_id (unique key, so no row multiplication)
n_annot_in <- nrow(annot)
annot <- merge(annot,
               gene_eco[, .(gene_id, ecosystem)],
               by = "gene_id", all.x = TRUE)
stopifnot(nrow(annot) == n_annot_in)

n_unmatched <- sum(is.na(annot$ecosystem))
if (n_unmatched > 0) {
  stop(n_unmatched, " annotated genes could not be matched to an ecosystem")
}

# Recompute the supporting-database count the same way as _common.R (presence
# confirmed against the per-database ID columns, since the raw supporting_databases
# string can suffer comma-in-quoted-field parsing artefacts). Export n_db, not the
# raw n_databases, and assert the two agree so the exported column is verified.
annot[, has_cazy := grepl("CAZy", supporting_databases, fixed = TRUE) |
                    (cazy_families != "" & !is.na(cazy_families))]
annot[, has_kegg := grepl("KEGG", supporting_databases, fixed = TRUE) |
                    (ko_terms != "" & !is.na(ko_terms))]
annot[, has_pfam := grepl("Pfam", supporting_databases, fixed = TRUE) |
                    (pfam_ids != "" & !is.na(pfam_ids))]
annot[, n_db := as.integer(has_cazy) + as.integer(has_kegg) + as.integer(has_pfam)]
stopifnot(all(annot$n_db == annot$n_databases))

# Total viral genes per ecosystem (denominator for normalization)
total_genes_eco <- gene_eco[, .(n_total_genes = .N), by = ecosystem]

samples_eco <- sample_meta[, .(n_samples = .N), by = ecosystem]

# Normalization denominators must be strictly positive.
stopifnot(all(total_genes_eco$n_total_genes > 0), all(samples_eco$n_samples > 0))

amg_counts <- annot[, .(n_annotated = .N),
                    by = .(ecosystem, functional_category)]

# Build the full ecosystem x functional-category grid before summarising, so an
# ecosystem with zero annotated AMGs (e.g. Salt marsh) is retained in the written
# summaries rather than dropped. The grid uses the ecosystem list
# (samples_eco) crossed with the functional categories present in the catalogue.
func_categories <- sort(unique(amg_counts$functional_category))
eco_grid <- CJ(ecosystem = samples_eco$ecosystem,
               functional_category = func_categories)
eco_summary <- merge(eco_grid, amg_counts,
                     by = c("ecosystem", "functional_category"), all.x = TRUE)
eco_summary[is.na(n_annotated), n_annotated := 0L]
eco_summary <- merge(eco_summary, total_genes_eco, by = "ecosystem", all.x = TRUE)
eco_summary <- merge(eco_summary, samples_eco, by = "ecosystem", all.x = TRUE)
eco_summary[is.na(n_total_genes), n_total_genes := 0L]

# Every ecosystem with samples must be present in the summary.
stopifnot(uniqueN(eco_summary$ecosystem) == nrow(samples_eco))

eco_summary[, `:=`(
  per_10k_genes = fifelse(n_total_genes > 0,
                          (n_annotated / n_total_genes) * 10000, 0),
  per_sample    = n_annotated / n_samples
)]

cat("\n--- Ecosystem summary (all bins) ---\n")
print(eco_summary[order(ecosystem, functional_category)])

fwrite(eco_summary[order(ecosystem, functional_category)],
       file.path(results_dir, "04_ecosystem_normalized.tsv"), sep = "\t")

fwrite(samples_eco[order(-n_samples)],
       file.path(results_dir, "04_ecosystem_sample_counts.tsv"), sep = "\t")

raw_wide <- dcast(eco_summary, ecosystem + n_samples + n_total_genes ~
                    functional_category, value.var = "n_annotated", fill = 0)
fwrite(raw_wide[order(-n_samples)],
       file.path(results_dir, "04_ecosystem_summary.tsv"), sep = "\t")

brewer_set2 <- RColorBrewer::brewer.pal(3, "Set2")
cat_colours <- c(
  "carbon_cycling"   = brewer_set2[1],
  "nitrogen_cycling"  = brewer_set2[2],
  "antibiotics"       = brewer_set2[3]
)

cat_labels <- c(
  "carbon_cycling"   = "Carbon Cycling",
  "nitrogen_cycling"  = "Nitrogen Cycling",
  "antibiotics"       = "Antibiotic Resistance"
)

# Build the grid from the ecosystem list (samples_eco).
all_combos <- CJ(
  ecosystem = samples_eco$ecosystem,
  functional_category = names(cat_colours)
)
plot_data_full <- merge(all_combos, eco_summary,
                        by = c("ecosystem", "functional_category"), all.x = TRUE)
plot_data_full[is.na(per_10k_genes), per_10k_genes := 0]
plot_data_full[is.na(n_annotated), n_annotated := 0]

# Fill sample counts for combos that had no annotations
plot_data_full <- merge(plot_data_full[, !"n_samples"],
                        samples_eco, by = "ecosystem", all.x = TRUE)

eco_order_full <- plot_data_full[, .(total_rate = sum(per_10k_genes)),
                                 by = ecosystem][order(-total_rate)]$ecosystem

eco_labels_full <- samples_eco[match(eco_order_full, ecosystem),
                                paste0(ecosystem, "\n(n=", n_samples, ")")]

plot_data_full[, ecosystem := factor(ecosystem, levels = eco_order_full,
                                      labels = eco_labels_full)]
plot_data_full[, functional_category := factor(functional_category,
                                                levels = names(cat_colours))]

p_full <- ggplot(plot_data_full,
                 aes(x = ecosystem, y = per_10k_genes,
                     fill = functional_category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = cat_colours, labels = cat_labels, name = NULL) +
  labs(
    title = "A. Viral AMGs Across Ecosystems (Normalized)",
    x = "Ecosystem",
    y = "Annotated AMGs per 10,000 viral genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "04_fig1a_full.pdf"), p_full,
       width = 10, height = 5.5)
cat("\nSaved: figures/04_fig1a_full.pdf\n")

log_baseline <- 0.01
n_cats <- 3
dodge_total <- 0.8
bar_w <- dodge_total / n_cats * 0.85

plot_data_full_log <- copy(plot_data_full)
plot_data_full_log <- plot_data_full_log[per_10k_genes > 0]
plot_data_full_log[, eco_idx := as.numeric(ecosystem)]
plot_data_full_log[, cat_idx := as.numeric(functional_category)]
plot_data_full_log[, xpos := eco_idx +
                     (cat_idx - (n_cats + 1) / 2) * dodge_total / n_cats]

p_full_log <- ggplot(plot_data_full_log) +
  geom_rect(aes(xmin = xpos - bar_w / 2, xmax = xpos + bar_w / 2,
                ymin = log_baseline, ymax = per_10k_genes,
                fill = functional_category)) +
  scale_x_continuous(breaks = seq_along(eco_labels_full),
                     labels = eco_labels_full) +
  scale_y_log10(limits = c(log_baseline, 50),
                expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = cat_colours, labels = cat_labels, name = NULL) +
  labs(
    title = "A. Viral AMGs Across Ecosystems (Normalized, log scale)",
    x = "Ecosystem",
    y = "Annotated AMGs per 10,000 viral genes (log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "04_fig1a_full_log.pdf"), p_full_log,
       width = 10, height = 5.5)
cat("Saved: figures/04_fig1a_full_log.pdf\n")

min_samples <- 30

large_ecosystems <- samples_eco[n_samples >= min_samples]$ecosystem
small_ecosystems <- samples_eco[n_samples < min_samples]$ecosystem

# Collapse small ecosystems into "Other".

amg_pub <- copy(amg_counts)
amg_pub[ecosystem %in% small_ecosystems, ecosystem := "Other"]
amg_pub <- amg_pub[, .(n_annotated = sum(n_annotated)),
                   by = .(ecosystem, functional_category)]

tot_pub <- copy(total_genes_eco)
tot_pub[ecosystem %in% small_ecosystems, ecosystem := "Other"]
tot_pub <- tot_pub[, .(n_total_genes = sum(n_total_genes)), by = ecosystem]

samp_pub_denom <- copy(samples_eco)
samp_pub_denom[ecosystem %in% small_ecosystems, ecosystem := "Other"]
samp_pub_denom <- samp_pub_denom[, .(n_samples = sum(n_samples)), by = ecosystem]

eco_summary_pub <- merge(amg_pub, tot_pub, by = "ecosystem", all.x = TRUE)
eco_summary_pub <- merge(eco_summary_pub, samp_pub_denom, by = "ecosystem", all.x = TRUE)

eco_summary_pub[, `:=`(
  per_10k_genes = (n_annotated / n_total_genes) * 10000,
  per_sample    = n_annotated / n_samples
)]

# Build from the pooled ecosystem list (tot_pub), so a bar with
# zero AMGs is still drawn rather than dropped (same guard as the full figure).
all_combos_pub <- CJ(
  ecosystem = tot_pub$ecosystem,
  functional_category = names(cat_colours)
)
plot_data_pub <- merge(all_combos_pub, eco_summary_pub,
                       by = c("ecosystem", "functional_category"), all.x = TRUE)
plot_data_pub[is.na(per_10k_genes), per_10k_genes := 0]
plot_data_pub[is.na(n_annotated), n_annotated := 0]

samples_pub <- copy(samples_eco)
samples_pub[ecosystem %in% small_ecosystems, ecosystem := "Other"]
samples_pub <- samples_pub[, .(n_samples = sum(n_samples)), by = ecosystem]
plot_data_pub <- merge(plot_data_pub[, !"n_samples"],
                       samples_pub, by = "ecosystem", all.x = TRUE)

# Order: named ecosystems by rate, then "Other" and "Unclassified" at end
eco_order_pub <- plot_data_pub[!ecosystem %in% c("Other", "Unclassified"),
                               .(total_rate = sum(per_10k_genes)),
                               by = ecosystem][order(-total_rate)]$ecosystem
eco_order_pub <- c(eco_order_pub,
                   intersect(c("Unclassified", "Other"),
                             unique(plot_data_pub$ecosystem)))

eco_labels_pub <- samples_pub[match(eco_order_pub, ecosystem),
                               paste0(ecosystem, "\n(n=", n_samples, ")")]

plot_data_pub[, ecosystem := factor(ecosystem, levels = eco_order_pub,
                                     labels = eco_labels_pub)]
plot_data_pub[, functional_category := factor(functional_category,
                                               levels = names(cat_colours))]

p_pub <- ggplot(plot_data_pub,
                aes(x = ecosystem, y = per_10k_genes,
                    fill = functional_category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = cat_colours, labels = cat_labels, name = NULL) +
  labs(
    tag = "A",
    x = "Ecosystem",
    y = "Annotated AMGs per 10,000 viral genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top",
    plot.tag = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "04_fig1a_publication.pdf"), p_pub,
       width = 8, height = 5)
cat("Saved: figures/04_fig1a_publication.pdf\n")

# Log-scale version (geom_rect for proper baseline)
plot_data_pub_log <- copy(plot_data_pub)
plot_data_pub_log <- plot_data_pub_log[per_10k_genes > 0]
plot_data_pub_log[, eco_idx := as.numeric(ecosystem)]
plot_data_pub_log[, cat_idx := as.numeric(functional_category)]
plot_data_pub_log[, xpos := eco_idx +
                    (cat_idx - (n_cats + 1) / 2) * dodge_total / n_cats]

p_pub_log <- ggplot(plot_data_pub_log) +
  geom_rect(aes(xmin = xpos - bar_w / 2, xmax = xpos + bar_w / 2,
                ymin = log_baseline, ymax = per_10k_genes,
                fill = functional_category)) +
  scale_x_continuous(breaks = seq_along(eco_labels_pub),
                     labels = eco_labels_pub) +
  scale_y_log10(limits = c(log_baseline, 50),
                expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = cat_colours, labels = cat_labels, name = NULL) +
  labs(
    title = "A. Viral AMGs Across Ecosystems (Normalized, log scale)",
    x = "Ecosystem",
    y = "Annotated AMGs per 10,000 viral genes (log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "04_fig1a_publication_log.pdf"), p_pub_log,
       width = 8, height = 5)
cat("Saved: figures/04_fig1a_publication_log.pdf\n")

brewer_greens <- RColorBrewer::brewer.pal(7, "Set2")[1]  # carbon
brewer_orange <- RColorBrewer::brewer.pal(7, "Set2")[2]  # nitrogen
brewer_purple <- RColorBrewer::brewer.pal(7, "Set2")[3]  # AMR

carbon <- annot[functional_category == "carbon_cycling"]

# Map fine-grained MFD subcategories to display groups
carbon[, enzyme_group := fcase(
  subcategories == "chitinase", "Chitinase",

  subcategories %like% "beta_glucanase", "β-1,3-glucanase",

  subcategories %like% "algal_polysaccharide", "Polysaccharide lyase",

  subcategories == "pectinase", "Pectinase",

  subcategories %in% c("mannanase", "other_hemicellulase", "xylanase") |
    subcategories %like% "arabinofuranosidase", "Hemicellulase",

  subcategories %in% c("cellulase", "beta_glucosidase") |
    subcategories %like% "cellulase.*mannanase|arabinogalactanase.*cellulase",
    "Cellulase",

  default = "Other"
)]

carbon_summary <- carbon[, .N, by = enzyme_group][order(-N)]
cat("\n--- Carbon cycling enzyme groups ---\n")
print(carbon_summary)

carbon_summary[, enzyme_group := factor(enzyme_group,
                                         levels = enzyme_group)]

p_panelB <- ggplot(carbon_summary, aes(x = enzyme_group, y = N)) +
  geom_col(fill = brewer_greens, width = 0.7) +
  geom_text(aes(label = N), vjust = -0.3, size = 3, colour = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    tag = "B",
    x = "Functional group",
    y = "Gene count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.tag = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "04_fig1b_carbon_enzymes.pdf"), p_panelB,
       width = 6, height = 4.5, device = cairo_pdf)
cat("Saved: figures/04_fig1b_carbon_enzymes.pdf\n")

nitrogen <- annot[functional_category == "nitrogen_cycling"]

nitrogen[, mechanism := fcase(
  subcategories == "nitrogen_fixation", "Nitrogen fixation",
  subcategories == "denitrification", "Denitrification",
  subcategories == "ammonia_assimilation", "Ammonia assimilation",
  default = "Other"
)]

nitrogen_summary <- nitrogen[, .N, by = mechanism][order(-N)]
cat("\n--- Nitrogen cycling mechanisms ---\n")
print(nitrogen_summary)

nitrogen_summary[, mechanism := factor(mechanism, levels = mechanism)]

p_panelC <- ggplot(nitrogen_summary, aes(x = mechanism, y = N)) +
  geom_col(fill = brewer_orange, width = 0.7) +
  geom_text(aes(label = N), vjust = -0.3, size = 3, colour = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    tag = "C",
    x = "Functional group",
    y = "Gene count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.tag = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "04_fig1c_nitrogen_mechanisms.pdf"), p_panelC,
       width = 5, height = 4.5)
cat("Saved: figures/04_fig1c_nitrogen_mechanisms.pdf\n")

amr <- annot[functional_category == "antibiotics"]

amr[, mechanism := fcase(
  subcategories == "antibiotic_efflux", "Efflux",
  subcategories == "antibiotic_target_alteration", "Target alteration",
  subcategories == "antibiotic_inactivation", "Inactivation",
  subcategories %like% "inactivation.*alteration|alteration.*inactivation",
    "Inactivation + Target alteration",
  default = "Other"
)]

amr_summary <- amr[, .N, by = mechanism][order(-N)]
cat("\n--- AMR mechanisms ---\n")
print(amr_summary)

amr_summary[, mechanism := factor(mechanism, levels = mechanism)]

p_panelD <- ggplot(amr_summary, aes(x = mechanism, y = N)) +
  geom_col(fill = brewer_purple, width = 0.7) +
  geom_text(aes(label = N), vjust = -0.3, size = 3, colour = "grey30") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    tag = "D",
    x = "Functional group",
    y = "Gene count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.tag = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "04_fig1d_amr_mechanisms.pdf"), p_panelD,
       width = 5, height = 4.5)
cat("Saved: figures/04_fig1d_amr_mechanisms.pdf\n")

library(openxlsx)

prettify <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)
  x <- gsub("Amr", "AMR", x)
  x
}

title_style    <- createStyle(fontSize = 11, textDecoration = "bold")
header_style   <- createStyle(
  fontSize = 10, textDecoration = "bold",
  border = "TopBottom", borderStyle = "thin",
  halign = "center", wrapText = TRUE
)
body_style     <- createStyle(fontSize = 10, halign = "center")
body_left      <- createStyle(fontSize = 10, halign = "left")
body_right     <- createStyle(fontSize = 10, halign = "right", numFmt = "0.00")
bottom_border  <- createStyle(border = "bottom", borderStyle = "thin")
note_style     <- createStyle(fontSize = 9, halign = "left", wrapText = TRUE)

# Compute per-ecosystem statistics including zero-gene sample breakdown

gene_eco[, taxon_oid_chr := as.character(taxon_oid)]
genes_per_sample <- gene_eco[, .(n_viral_genes = .N), by = taxon_oid]

# All samples with ecosystem assignment (incl. those with 0 viral genes)
all_samples <- sample_meta[, .(taxon_oid = IMG.taxon_oid, ecosystem)]
all_samples <- merge(all_samples,
                     genes_per_sample[, .(taxon_oid = taxon_oid, n_viral_genes)],
                     by = "taxon_oid", all.x = TRUE)
all_samples[is.na(n_viral_genes), n_viral_genes := 0]

eco_detail <- all_samples[, .(
  n_samples_total     = .N,
  n_samples_with_genes = sum(n_viral_genes > 0),
  n_samples_zero_genes = sum(n_viral_genes == 0),
  pct_with_genes       = round(100 * sum(n_viral_genes > 0) / .N, 1),
  total_viral_genes    = sum(n_viral_genes),
  mean_genes_per_sample = round(mean(n_viral_genes), 1),
  median_genes_per_sample = as.double(median(n_viral_genes))
), by = ecosystem]

# Broad annotation rate: genes with >= 1 hit in Pfam, KEGG or CAZy, independent
# of the three AMG categories in "Normalized AMG rates" below. Gives the
# unannotated-gene fraction (total functional diversity) per ecosystem.
gene_eco[, any_hit := (!is.na(pfam) & nzchar(pfam)) |
                       (!is.na(cazyme) & nzchar(cazyme)) |
                       (!is.na(kegg_ortholog) & nzchar(kegg_ortholog))]
annot_any_eco <- gene_eco[, .(n_annotated_any = sum(any_hit)), by = ecosystem]

eco_detail <- merge(eco_detail, annot_any_eco, by = "ecosystem", all.x = TRUE)
eco_detail[is.na(n_annotated_any), n_annotated_any := 0]
eco_detail[, pct_annotated_any := round(100 * n_annotated_any / total_viral_genes, 1)]

amg_wide <- dcast(eco_summary, ecosystem ~ functional_category,
                  value.var = "n_annotated", fill = 0)
setnames(amg_wide,
         c("carbon_cycling", "nitrogen_cycling", "antibiotics"),
         c("n_carbon_amgs", "n_nitrogen_amgs", "n_antibiotic_amgs"),
         skip_absent = TRUE)
amg_wide[, n_total_amgs := rowSums(.SD),
         .SDcols = intersect(c("n_carbon_amgs", "n_nitrogen_amgs", "n_antibiotic_amgs"),
                             names(amg_wide))]

eco_table <- merge(eco_detail, amg_wide, by = "ecosystem", all.x = TRUE)
# Fill NAs for ecosystems with no AMGs
for (col in c("n_carbon_amgs", "n_nitrogen_amgs", "n_antibiotic_amgs", "n_total_amgs"))
  if (col %in% names(eco_table)) eco_table[is.na(get(col)), (col) := 0]

eco_table[, `:=`(
  carbon_per_10k   = round(ifelse(total_viral_genes > 0,
                                   n_carbon_amgs / total_viral_genes * 10000, 0), 2),
  nitrogen_per_10k = round(ifelse(total_viral_genes > 0,
                                   n_nitrogen_amgs / total_viral_genes * 10000, 0), 2),
  abr_per_10k      = round(ifelse(total_viral_genes > 0,
                                   n_antibiotic_amgs / total_viral_genes * 10000, 0), 2),
  total_per_10k    = round(ifelse(total_viral_genes > 0,
                                   n_total_amgs / total_viral_genes * 10000, 0), 2)
)]

eco_table <- eco_table[order(-total_per_10k)]

wb <- createWorkbook()

addWorksheet(wb, "Ecosystem overview")

s1 <- eco_table[, .(
  Ecosystem                       = ecosystem,
  `Total samples`                 = n_samples_total,
  `Samples with viral genes`      = n_samples_with_genes,
  `Samples without viral genes`   = n_samples_zero_genes,
  `% with viral genes`            = pct_with_genes,
  `Total viral genes`             = total_viral_genes,
  `Mean viral genes per sample`   = mean_genes_per_sample,
  `Median viral genes per sample` = median_genes_per_sample,
  `Genes with any functional annotation (Pfam/KEGG/CAZy)` = n_annotated_any,
  `% annotated (any hit)`         = pct_annotated_any
)]

title1 <- paste0(
  "Table S13 | Ecosystem sampling overview. Number of metagenome samples, ",
  "viral gene recovery, zero-gene sample fraction, and broad functional ",
  "annotation rate (any Pfam/KEGG/CAZy hit) per ecosystem category."
)
writeData(wb, 1, title1, startRow = 1, startCol = 1)
addStyle(wb, 1, title_style, rows = 1, cols = 1)
mergeCells(wb, 1, cols = 1:ncol(s1), rows = 1)

writeData(wb, 1, s1, startRow = 3, headerStyle = header_style)
ns1 <- nrow(s1)
addStyle(wb, 1, body_left, rows = 4:(3 + ns1), cols = 1, gridExpand = TRUE)
for (j in 2:ncol(s1))
  addStyle(wb, 1, body_style, rows = 4:(3 + ns1), cols = j, gridExpand = TRUE)
addStyle(wb, 1, bottom_border, rows = 3 + ns1, cols = 1:ncol(s1), stack = TRUE)

note1 <- paste0(
  "Note: 'Unclassified' comprises samples for which GOLD.Specific.Ecosystem ",
  "was not curated in the IMG/M database. This is a metadata gap, not an ",
  "analytical artefact. Ecosystem labels were harmonized from GOLD ontology ",
  "(e.g., 'Forest Soil'/'Forest soil' merged as 'Forest'; 'Agricultural land'/",
  "'Agricultural'/'Agricultural soil' merged as 'Agricultural'). ",
  "Shrubland (n=7, all peptidoglycanase, zero AMGs) is a named GOLD ecosystem pooled into 'Other' by the <30-sample rule. ",
  "'Genes with any functional annotation' counts genes matching at least one ",
  "of Pfam, KEGG or CAZy, independent of the three AMG categories reported in ",
  "the 'Normalized AMG rates' sheet; this is the broad annotation rate behind ",
  "the catalogue-wide dark-matter estimate in the Discussion. Ecosystems with ",
  "fewer than 100 samples (Volcanic, Tropical rainforest, Salt marsh, Desert) ",
  "rest on 2,001-6,232 genes and their annotation rate should be read with caution."
)
writeData(wb, 1, note1, startRow = 3 + ns1 + 2, startCol = 1)
addStyle(wb, 1, note_style, rows = 3 + ns1 + 2, cols = 1)
mergeCells(wb, 1, cols = 1:ncol(s1), rows = 3 + ns1 + 2)

setColWidths(wb, 1, cols = 1, widths = 22)
setColWidths(wb, 1, cols = 2:ncol(s1), widths = 16)

addWorksheet(wb, "Normalized AMG rates")

s2 <- eco_table[, .(
  Ecosystem                            = ecosystem,
  `Total samples`                      = n_samples_total,
  `Samples with viral genes`           = n_samples_with_genes,
  `Total viral genes`                  = total_viral_genes,
  `Carbon cycling AMGs`                = n_carbon_amgs,
  `Nitrogen cycling AMGs`              = n_nitrogen_amgs,
  `Antibiotic resistance AMGs`         = n_antibiotic_amgs,
  `Total AMGs`                         = n_total_amgs,
  `Carbon per 10k viral genes`         = carbon_per_10k,
  `Nitrogen per 10k viral genes`       = nitrogen_per_10k,
  `Antibiotic res. per 10k viral genes` = abr_per_10k,
  `Total AMGs per 10k viral genes`     = total_per_10k
)]

title2 <- paste0(
  "Table S2 | Normalized AMG rates per ecosystem. Raw AMG counts and rates ",
  "per 10,000 viral genes, stratified by functional category. ",
  "Peptidoglycanase and broad (non-enzyme-resolved) matches excluded."
)
writeData(wb, 2, title2, startRow = 1, startCol = 1)
addStyle(wb, 2, title_style, rows = 1, cols = 1)
mergeCells(wb, 2, cols = 1:ncol(s2), rows = 1)

writeData(wb, 2, s2, startRow = 3, headerStyle = header_style)
ns2 <- nrow(s2)
addStyle(wb, 2, body_left, rows = 4:(3 + ns2), cols = 1, gridExpand = TRUE)
for (j in 2:8)
  addStyle(wb, 2, body_style, rows = 4:(3 + ns2), cols = j, gridExpand = TRUE)
for (j in 9:ncol(s2))
  addStyle(wb, 2, body_right, rows = 4:(3 + ns2), cols = j, gridExpand = TRUE)
addStyle(wb, 2, bottom_border, rows = 3 + ns2, cols = 1:ncol(s2), stack = TRUE)

note2 <- paste0(
  "Note: Normalization denominator is the total number of viral genes ",
  "(all genes on geNomad-identified viral contigs) per ecosystem, not the ",
  "number of samples. This controls for differences in sequencing depth ",
  "and viral gene recovery across ecosystems."
)
writeData(wb, 2, note2, startRow = 3 + ns2 + 2, startCol = 1)
addStyle(wb, 2, note_style, rows = 3 + ns2 + 2, cols = 1)
mergeCells(wb, 2, cols = 1:ncol(s2), rows = 3 + ns2 + 2)

setColWidths(wb, 2, cols = 1, widths = 22)
setColWidths(wb, 2, cols = 2:ncol(s2), widths = 16)

addWorksheet(wb, "Ecosystem label mapping")

s3 <- data.table(
  `Ecosystem bin`     = c("Agricultural", "Forest", "Grassland", "Permafrost",
                           "Desert", "Tropical rainforest", "Shrubland", "Peat",
                           "Plant-associated", "Volcanic", "Salt marsh",
                           "Unclassified"),
  `Source column`     = c(rep("GOLD.Specific.Ecosystem", 7),
                           rep("GOLD.Ecosystem.Type", 4),
                           "Various"),
  `Original labels`   = c("Agricultural land; Agricultural; Agricultural soil",
                           "Forest Soil; Forest soil",
                           "Grasslands",
                           "Permafrost",
                           "Desert",
                           "Tropical rainforest",
                           "Shrubland",
                           "Peat (all subtypes)",
                           "Rhizosphere; Roots; Rhizoplane",
                           "Volcanic",
                           "Marine (all = Salt marsh)",
                           "Soil/Unclassified; Unclassified/Unclassified"),
  `Rationale`         = c("Synonym harmonization",
                           "Case harmonization",
                           "Direct mapping",
                           "Direct mapping",
                           "Direct mapping",
                           "Direct mapping",
                           "Named ecosystem, <30 samples -> pooled as Other",
                           "Distinct ecosystem type in GOLD",
                           "All plant-associated soil types combined",
                           "Distinct ecosystem type in GOLD",
                           "All marine samples were salt marsh",
                           "No specific ecosystem metadata in IMG/GOLD")
)

title3 <- paste0(
  "Table S16 | Ecosystem label harmonization. Mapping of GOLD ontology ",
  "metadata to the ecosystem categories used in this study."
)
writeData(wb, 3, title3, startRow = 1, startCol = 1)
addStyle(wb, 3, title_style, rows = 1, cols = 1)
mergeCells(wb, 3, cols = 1:ncol(s3), rows = 1)

writeData(wb, 3, s3, startRow = 3, headerStyle = header_style)
ns3 <- nrow(s3)
addStyle(wb, 3, body_left, rows = 4:(3 + ns3), cols = 1:ncol(s3), gridExpand = TRUE)
addStyle(wb, 3, bottom_border, rows = 3 + ns3, cols = 1:ncol(s3), stack = TRUE)

setColWidths(wb, 3, cols = 1, widths = 22)
setColWidths(wb, 3, cols = 2, widths = 24)
setColWidths(wb, 3, cols = 3, widths = 52)
setColWidths(wb, 3, cols = 4, widths = 38)

addWorksheet(wb, "Subcategory summaries")

# Carbon cycling enzyme groups
carbon_out <- carbon_summary[, .(
  Category = "Carbon cycling",
  Subcategory = enzyme_group,
  `Gene count` = N,
  `%` = round(100 * N / sum(N), 1)
)]
# Nitrogen cycling
nitrogen_out <- nitrogen_summary[, .(
  Category = "Nitrogen cycling",
  Subcategory = mechanism,
  `Gene count` = N,
  `%` = round(100 * N / sum(N), 1)
)]
# AMR
amr_out <- amr_summary[, .(
  Category = "Antibiotic resistance",
  Subcategory = mechanism,
  `Gene count` = N,
  `%` = round(100 * N / sum(N), 1)
)]

s4 <- rbindlist(list(carbon_out, nitrogen_out, amr_out))

title4 <- paste0(
  "Table S4 | Subcategory breakdown for Figures 1B-D. Gene counts per ",
  "enzyme group / mechanism within each functional category. ",
  "Only specific (enzyme-resolved) annotations included."
)
writeData(wb, 4, title4, startRow = 1, startCol = 1)
addStyle(wb, 4, title_style, rows = 1, cols = 1)
mergeCells(wb, 4, cols = 1:ncol(s4), rows = 1)

writeData(wb, 4, s4, startRow = 3, headerStyle = header_style)
ns4 <- nrow(s4)
addStyle(wb, 4, body_left, rows = 4:(3 + ns4), cols = 1:2, gridExpand = TRUE)
addStyle(wb, 4, body_style, rows = 4:(3 + ns4), cols = 3:4, gridExpand = TRUE)
addStyle(wb, 4, bottom_border, rows = 3 + ns4, cols = 1:ncol(s4), stack = TRUE)

carbon_end <- 3 + nrow(carbon_out)
nitrogen_end <- carbon_end + nrow(nitrogen_out)
addStyle(wb, 4, bottom_border, rows = carbon_end, cols = 1:ncol(s4), stack = TRUE)
addStyle(wb, 4, bottom_border, rows = nitrogen_end, cols = 1:ncol(s4), stack = TRUE)

setColWidths(wb, 4, cols = 1, widths = 22)
setColWidths(wb, 4, cols = 2, widths = 30)
setColWidths(wb, 4, cols = 3:4, widths = 14)

addWorksheet(wb, "Gene annotations")

s5 <- annot[, .(
  `Gene ID`             = gene_id,
  `Functional category` = prettify(functional_category),
  `Subcategory`         = subcategories,
  `Ecosystem`           = ecosystem,
  `CAZy families`       = cazy_families,
  `KEGG ortholog`       = ko_terms,
  `Pfam IDs`            = pfam_ids,
  `Supporting databases` = supporting_databases,
  `No. databases`       = n_db
)]

title5 <- paste0(
  "Figure 1B-D source data | Complete gene annotation list. All specific ",
  "(enzyme-resolved) AMG annotations across Figures 1B-D. Peptidoglycanase and ",
  "broad matches excluded. (Provided as a dataset, not a numbered supplementary table.)"
)
writeData(wb, 5, title5, startRow = 1, startCol = 1)
addStyle(wb, 5, title_style, rows = 1, cols = 1)
mergeCells(wb, 5, cols = 1:ncol(s5), rows = 1)

writeData(wb, 5, s5, startRow = 3, headerStyle = header_style)
ns5 <- nrow(s5)
addStyle(wb, 5, body_left, rows = 4:(3 + ns5), cols = 1:ncol(s5), gridExpand = TRUE)

setColWidths(wb, 5, cols = 1, widths = 40)
setColWidths(wb, 5, cols = 2:3, widths = 22)
setColWidths(wb, 5, cols = 4, widths = 18)
setColWidths(wb, 5, cols = 5:9, widths = 16)

saveWorkbook(wb, file.path(results_dir, "04_ecosystem_supplement.xlsx"),
             overwrite = TRUE)
cat("Saved: results/04_ecosystem_supplement.xlsx\n")


library(patchwork)

# B (6 bars) gets more width than C (3 bars) and D (4 bars).
composite <- p_pub /
  (p_panelB + p_panelC + p_panelD +
     plot_layout(widths = c(3, 1.5, 2))) +
  plot_layout(heights = c(1, 1))

# TIFF at 300 dpi (AEM preferred format)
tiff(file.path(fig_dir, "04_Figure1_composite.tiff"),
     width = 12, height = 9, units = "in", res = 300,
     compression = "lzw", type = "cairo")
print(composite)
dev.off()
cat("Saved: figures/04_Figure1_composite.tiff\n")

cairo_pdf(file.path(fig_dir, "04_Figure1_composite.pdf"),
    width = 12, height = 9)
print(composite)
dev.off()
cat("Saved: figures/04_Figure1_composite.pdf\n")

sink(file.path(results_dir, "04_sessionInfo.txt"))
cat("Script: 04_amg_distribution.R\n")
cat("Date:", format(Sys.time()), "\n\n")
sessionInfo()
sink()

cat("\nDone. All outputs in", results_dir, "\n")
