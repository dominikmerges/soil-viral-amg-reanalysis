#!/usr/bin/env Rscript
# 07_lysis_gene_decomposition.R
#
# Figure 2: separating genes that serve core viral functions from candidate
# auxiliary metabolic genes in the soil viral CAZyme inventory. Glycoside
# hydrolase annotations in viral genomes may reflect host entry, virion
# decoration or lysis rather than nutrient acquisition (Martin et al. 2025).
#
#   A  what the peptidoglycanase (lysis) filter removes from the carbon
#      repertoire, and how the composition would look without it
#   B  CAZy-family decomposition of all six carbon-cycling groups of Fig. 1b,
#      classified by whether the family has a documented virion-structural or
#      lysis role, or no catalytic module at all
#   C  lysis-cassette membership test: is the gene adjacent to a holin/spanin,
#      as an endolysin would be?  GH24 phage lysozyme = positive control,
#      random gene on the same contigs = null.
#   D  length-matched carriage of a SECOND peptidoglycan hydrolase on
#      GH19-bearing contigs vs all viral contigs.
#
# In:  data/metagenome_confidence_annotations.csv
#      data/GSVA_soil_viruses_gene_metadata_4.tsv.gz
#      data/GSVA_soil_viruses_genome_metadata_2.tsv.gz
# Out: results/07_chitinase_family_decomposition.tsv
#      results/07_carbon_family_decomposition.tsv
#      results/07_lysis_filter_effect.tsv
#      results/07_cassette_adjacency.tsv
#      results/07_gh19_contig_context.tsv
#      results/07_length_matched_pgh.tsv
#      results/07_Table_S_lysis_decomposition.xlsx
#      results/figures/07_Figure2_composite.{pdf,tiff}
#      results/07_sessionInfo.txt

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(openxlsx)
})

set.seed(42)
source("_common.R")

dirs        <- mfd_init_dirs()
results_dir <- dirs$results_dir
fig_dir     <- dirs$fig_dir
data_dir    <- "data"

# Peptidoglycan hydrolase (endolysin) catalytic domains. Pfam accessions
# verified against InterPro; CAZy families are the canonical phage endolysin
# families (Cahill & Young 2019; Latka et al. 2017).
PGH_PFAM <- c(
  "PF00959",  # Phage_lysozyme
  "PF01183",  # Glyco_hydro_25
  "PF01464",  # SLT, transglycosylase
  "PF05838",  # Glyco_hydro_108
  "PF01832",  # Glucosaminidase (mannosyl-glycoprotein endo-beta-GlcNAc'ase)
  "PF01510",  # Amidase_2
  "PF01520",  # Amidase_3
  "PF05257",  # CHAP
  "PF00877",  # NlpC/P60
  "PF01551",  # Peptidase_M23
  "PF08291",  # Peptidase_M15_3
  "PF13539",  # Peptidase_M15_4
  "PF02557",  # VanY
  "PF00062"   # C-type lysozyme
)

# CAZy aligned to the MetaFuncDecoder peptidoglycanase subcategory
# (Lopez-Mondejar et al. 2022, Table S3). GH22, GH102, GH153 and CE9 contribute
# zero genes in this catalogue, so this reproduces the previous 7-family signal
# exactly while using a single, citable family list across the analysis.
PGH_CAZY <- c("GH22", "GH23", "GH24", "GH25", "GH73", "GH102",
              "GH103", "GH104", "GH108", "GH153", "CE9")

PGH_KO <- c("K01185", "K07273", "K02395", "K03642", "K12054",
            "K08309", "K08305", "K21508", "K01448")

# Holins: all Pfam families whose name contains "holin" (InterPro query).
HOLIN_PFAM <- c(
  "PF04020", "PF04531", "PF04550", "PF04688", "PF04971", "PF05102", "PF05105",
  "PF05106", "PF05449", "PF06946", "PF07332", "PF09682", "PF10746", "PF10960",
  "PF11031", "PF11351", "PF13272", "PF16079", "PF16080", "PF16081", "PF16082",
  "PF16083", "PF16084", "PF16085", "PF16931", "PF16935", "PF16936", "PF16938",
  "PF16945", "PF23778", "PF23809", "PF23987", "PF24205", "PF26844", "PF27028",
  "PF28015"
)

# Spanins (outer-membrane disruption, third component of the lysis cassette).
SPANIN_PFAM <- c("PF03245", "PF17531", "PF26831")

# CAZy family annotation for panel B
# Non-catalytic: carbohydrate-binding modules carry no hydrolytic activity.
CBM_NONCATALYTIC <- c("CBM9", "CBM12", "CBM13", "CBM32", "CBM35",
                      "CBM50", "CBM61", "CBM67")

# Families with a documented virion-structural, host-entry or lysis role in
# viruses. GH19: Orlando et al. 2021; Edvardsen et al. 2025; Meng et al. 2024.
# Polysaccharide lyases and GH28: virion-associated depolymerases / tailspikes
# (Latka et al. 2017; Martin et al. 2025).
DUAL_USE_CAZY <- c("GH19",
                   "PL1", "PL1_2", "PL6", "PL7", "PL9", "PL14",
                   "GH28")

cat("--- loading ---------------------------------------------------------\n")

annot_all <- mfd_load_raw(data_dir)                       # 5,760 MFD calls
stopifnot(nrow(annot_all) == 5760L)
annot_all[, evidence_type := fifelse(subcategories == "" | is.na(subcategories),
                                     "broad", "specific")]

genes <- fread(file.path(data_dir, "GSVA_soil_viruses_gene_metadata_4.tsv.gz"),
               select = c("gene_id", "contig_id", "start_coordinate",
                          "end_coordinate", "pfam", "cazyme", "kegg_ortholog"))
setnafill_chr <- function(x) fifelse(is.na(x), "", x)
for (j in c("pfam", "cazyme", "kegg_ortholog"))
  set(genes, j = j, value = setnafill_chr(genes[[j]]))

contigs <- fread(file.path(data_dir, "GSVA_soil_viruses_genome_metadata_2.tsv.gz"),
                 select = c("contig_id", "contig_length", "checkv_quality",
                            "votu", "genus_cluster", "family_cluster"))

cat(sprintf("  %s genes on %s contigs; %s MFD calls\n",
            format(nrow(genes), big.mark = ","),
            format(uniqueN(genes$contig_id), big.mark = ","),
            format(nrow(annot_all), big.mark = ",")))

# Pfam field is ";"-separated with version suffixes: "PF00959.22;PF01471.9".
pfam_hit <- function(col, acc) {
  pat <- paste0("(^|;)(", paste(acc, collapse = "|"), ")\\.")
  grepl(pat, col)
}
# CAZy/KEGG fields carry no version suffix but may be ";"-delimited; match whole
# tokens so a family inside a multi-family record is not missed.
token_hit <- function(col, ids) {
  pat <- paste0("(^|;)(", paste(ids, collapse = "|"), ")($|;)")
  grepl(pat, col)
}
# TRUE only when EVERY ";"-delimited token in col is one of ids (and there is at
# least one token). A record is non-catalytic only if it carries no catalytic
# family alongside its binding module(s), so "GH19;CBM50" is not non-catalytic.
all_tokens_in <- function(col, ids) {
  vapply(strsplit(col, ";", fixed = TRUE), function(toks) {
    toks <- toks[toks != ""]
    length(toks) > 0L && all(toks %in% ids)
  }, logical(1))
}

genes[, is_pgh    := pfam_hit(pfam, PGH_PFAM) |
                     token_hit(cazyme, PGH_CAZY) |
                     token_hit(kegg_ortholog, PGH_KO)]
genes[, is_holin  := pfam_hit(pfam, HOLIN_PFAM)]
genes[, is_spanin := pfam_hit(pfam, SPANIN_PFAM)]
genes[, is_cassette_marker := is_holin | is_spanin]
genes[, is_lysis_any := is_pgh | is_holin | is_spanin]

cat("\n--- lysis genes in the catalogue ------------------------------------\n")
for (fl in c("is_pgh", "is_holin", "is_spanin", "is_lysis_any"))
  cat(sprintf("  %-13s %6d genes on %5d contigs\n", fl,
              genes[get(fl) == TRUE, .N],
              genes[get(fl) == TRUE, uniqueN(contig_id)]))

# Gene order along the contig (for adjacency).
setorder(genes, contig_id, start_coordinate)
genes[, gidx := seq_len(.N) - 1L, by = contig_id]
genes[, aa_len := (end_coordinate - start_coordinate + 1) / 3 - 1]

carbon <- annot_all[functional_category == "carbon_cycling" &
                    evidence_type == "specific"]

# Fig. 1b grouping (identical rules to 04_amg_distribution.R), plus lysis.
carbon[, enzyme_group := fcase(
  subcategories == "peptidoglycanase", "Peptidoglycanase (lysis)",
  subcategories == "chitinase", "Chitinase",
  # Unresolved GH16 (combined "algal_polysaccharide;beta_glucanase") + dedicated
  # β-glucanase GH families; soil GH16 activity is β-1,3-glucanase, not algal.
  subcategories %like% "beta_glucanase", "β-1,3-glucanase",
  # Genuine alginate lyases (PL6/7/14) from the algal "alginate.*lyase" pattern.
  subcategories %like% "algal_polysaccharide", "Polysaccharide lyase",
  subcategories == "pectinase", "Pectinase",
  subcategories %in% c("mannanase", "other_hemicellulase", "xylanase") |
    subcategories %like% "arabinofuranosidase", "Hemicellulase",
  subcategories %in% c("cellulase", "beta_glucosidase") |
    subcategories %like% "cellulase.*mannanase|arabinogalactanase.*cellulase",
    "Cellulase",
  default = "Other"
)]

lysis_effect <- rbind(
  carbon[, .(n = .N), by = enzyme_group][
    , .(scenario = "Lysis genes retained", enzyme_group, n)],
  carbon[enzyme_group != "Peptidoglycanase (lysis)", .(n = .N), by = enzyme_group][
    , .(scenario = "Lysis genes excluded", enzyme_group, n)]
)
lysis_effect[, pct := 100 * n / sum(n), by = scenario]

GROUP_LEVELS <- c("Peptidoglycanase (lysis)", "Chitinase", "Hemicellulase",
                  "β-1,3-glucanase", "Cellulase", "Pectinase",
                  "Polysaccharide lyase", "Other")
lysis_effect[, enzyme_group := factor(enzyme_group, levels = GROUP_LEVELS)]
lysis_effect[, scenario := factor(scenario, levels = c(
  "Lysis genes retained", "Lysis genes excluded"))]

fwrite(dcast(lysis_effect, enzyme_group ~ scenario, value.var = c("n", "pct")),
       file.path(results_dir, "07_lysis_filter_effect.tsv"), sep = "\t")

cat("\n--- panel A: lysis filter effect ------------------------------------\n")
print(dcast(lysis_effect, enzyme_group ~ scenario, value.var = "pct"))

# Catalogue-wide share of GH-class CAZy annotations that are endolysin families.
gh_all <- genes[cazyme %like% "^GH"]
cat(sprintf(paste("\n  GH-class CAZy annotations in the catalogue: %d;",
                  "in endolysin families %d (%.1f%%); incl. GH19 %d (%.1f%%)\n"),
            nrow(gh_all),
            gh_all[token_hit(cazyme, PGH_CAZY), .N],
            100 * gh_all[token_hit(cazyme, PGH_CAZY), .N] / nrow(gh_all),
            gh_all[token_hit(cazyme, c(PGH_CAZY, "GH19")), .N],
            100 * gh_all[token_hit(cazyme, c(PGH_CAZY, "GH19")), .N] / nrow(gh_all)))

group_cols <- c("Peptidoglycanase (lysis)" = "#B2182B",   # red  — lysis
                "Chitinase"                = "#08519C",
                "Hemicellulase"            = "#2171B5",
                "β-1,3-glucanase"          = "#4292C6",
                "Cellulase"                = "#6BAED6",
                "Pectinase"                = "#9ECAE1",
                "Polysaccharide lyase"     = "#C6DBEF",
                "Other"                    = "#BABABA")   # grey


setorder(lysis_effect, scenario, enzyme_group)
lysis_effect[, ypos := 100 - (cumsum(pct) - pct / 2), by = scenario]
lysis_effect[, lab := fifelse(pct >= 3,
                sprintf("%s  (%.1f%%)", format(n, big.mark = ","), pct), "")]
lysis_effect[, txt_col := fifelse(enzyme_group == "Peptidoglycanase (lysis)",
                                  "white", "grey15")]

pA <- ggplot(lysis_effect, aes(x = scenario, y = pct, fill = enzyme_group)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(y = ypos, label = lab, colour = txt_col), size = 2.6) +
  scale_colour_identity() +
  scale_fill_manual(values = group_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "A", x = "Annotation set",
       y = "% of enzyme-resolved carbon-cycling calls") +
  theme_minimal(base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 14),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8),
        panel.grid.major.x = element_blank())

carb_noly <- carbon[enzyme_group != "Peptidoglycanase (lysis)"]
carb_noly[, cazy_fam := fifelse(cazy_families == "" | is.na(cazy_families),
                                "(no CAZy family)", cazy_families)]

# Classify by whole CAZy tokens. Dual-use catalytic families are evaluated before
# the CBM test so a record carrying both (e.g. "GH19;CBM50") is assigned by its
# catalytic family; "Non-catalytic binding module" requires ALL tokens to be CBMs.
carb_noly[, family_class := fcase(
  cazy_fam == "(no CAZy family)",             "No CAZy family (KEGG/Pfam only)",
  token_hit(cazy_fam, DUAL_USE_CAZY),         "Documented virion-structural or lysis role",
  all_tokens_in(cazy_fam, CBM_NONCATALYTIC),  "Non-catalytic binding module",
  default =                                   "Polysaccharide-degrading, no viral structural role reported"
)]

CLASS_LEVELS <- c("Documented virion-structural or lysis role",
                  "Non-catalytic binding module",
                  "No CAZy family (KEGG/Pfam only)",
                  "Polysaccharide-degrading, no viral structural role reported")
carb_noly[, family_class := factor(family_class, levels = CLASS_LEVELS)]

fam_decomp <- carb_noly[, .N, by = .(enzyme_group, cazy_fam, family_class)][order(-N)]

other_sub <- carb_noly[enzyme_group == "Other",
                       .(subcat = paste(sort(unique(
                          sub("^(.)", "\\U\\1", gsub("_", " ", subcategories), perl = TRUE)
                        )), collapse = "; ")), by = cazy_fam]
fam_decomp[other_sub, on = "cazy_fam", subcat := i.subcat]
fam_decomp[enzyme_group == "Other" & !is.na(subcat),
           enzyme_group := paste0("Other (", subcat, ")")]
fam_decomp[, subcat := NULL]
fwrite(fam_decomp, file.path(results_dir, "07_carbon_family_decomposition.tsv"),
       sep = "\t")

chit_decomp <- carb_noly[enzyme_group == "Chitinase",
                         .N, by = .(cazy_fam, family_class)][order(-N)]
chit_decomp[, pct := round(100 * N / sum(N), 1)]
fwrite(chit_decomp, file.path(results_dir, "07_chitinase_family_decomposition.tsv"),
       sep = "\t")

cat("\n--- panel B: the 487 chitinase-associated calls ----------------------\n")
print(chit_decomp)

class_by_group <- carb_noly[, .N, by = .(enzyme_group, family_class)]
class_by_group[, pct := 100 * N / sum(N), by = enzyme_group]
grp_tot <- carb_noly[, .(tot = .N), by = enzyme_group]
class_by_group <- merge(class_by_group, grp_tot, by = "enzyme_group")
class_by_group[, enzyme_group := factor(enzyme_group, levels = GROUP_LEVELS)]

cat("\n--- panel B: composition by carbon group ----------------------------\n")
print(dcast(class_by_group, enzyme_group ~ family_class, value.var = "pct"))

class_cols <- c("Documented virion-structural or lysis role"  = "#B2182B",
                "Non-catalytic binding module"                = "#762A83",
                "No CAZy family (KEGG/Pfam only)"             = "#BABABA",
                "Polysaccharide-degrading, no viral structural role reported"
                                                              = "#2166AC")

grp_lab <- setNames(sprintf("%s\n(n = %d)", grp_tot$enzyme_group, grp_tot$tot),
                    grp_tot$enzyme_group)


class_by_group[, lab := as.character(N)]   # label every present segment

pB <- ggplot(class_by_group, aes(x = enzyme_group, y = pct, fill = family_class)) +
  geom_col(width = 0.66, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = lab),
            position = position_stack(vjust = 0.5), size = 2.6,
            colour = "white") +
  scale_fill_manual(values = class_cols, name = NULL,
                    labels = function(x) {
                      x <- gsub(", no viral structural role reported",
                                " (no viral structural role reported)", x)
                      vapply(x, function(s)
                             paste(strwrap(s, width = 22), collapse = "\n"),
                             character(1))
                    }) +
  scale_x_discrete(labels = grp_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(tag = "B", x = "Functional group", y = "% of calls in group") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_minimal(base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 7.5, angle = 30, hjust = 1),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7.5),
        panel.grid.major.x = element_blank())

gh19_ids <- annot_all[token_hit(subcategories, "chitinase") &
                      token_hit(cazy_families, "GH19"), gene_id]
cat(sprintf("\nGH19 chitinase-assigned genes: %d\n", length(gh19_ids)))

# Index of every holin/spanin per contig, for a fast nearest-neighbour lookup.
cass_pos <- genes[is_cassette_marker == TRUE, .(pos = list(gidx)), by = contig_id]
setkey(cass_pos, contig_id)

nearest_cassette <- function(gene_ids) {
  q <- genes[gene_id %chin% gene_ids, .(gene_id, contig_id, gidx)]
  q <- merge(q, cass_pos, by = "contig_id", all.x = TRUE)
  q[, dist := mapply(function(i, p) {
        if (is.null(p)) return(NA_integer_)
        p <- p[p != i]
        if (!length(p)) return(NA_integer_)
        min(abs(p - i))
      }, gidx, pos)]
  q[, .(gene_id, dist)]
}

adjacency_row <- function(label, gene_ids) {
  d <- nearest_cassette(gene_ids)
  n_have <- sum(!is.na(d$dist))
  n_adj  <- sum(d$dist <= 3, na.rm = TRUE)
  ci <- if (n_have > 0) binom.test(n_adj, n_have)$conf.int else c(NA, NA)
  data.table(set = label, n_genes = nrow(d), n_with_cassette = n_have,
             n_adjacent = n_adj,
             pct_adjacent = 100 * n_adj / max(n_have, 1),
             ci_lo = 100 * ci[1], ci_hi = 100 * ci[2],
             median_dist = as.numeric(median(d$dist, na.rm = TRUE)))
}

# Null: one randomly chosen non-GH19 gene from each GH19-bearing contig.
gh19_contigs <- unique(genes[gene_id %chin% gh19_ids, contig_id])
null_ids <- genes[contig_id %chin% gh19_contigs & !(gene_id %chin% gh19_ids)][
  sample(.N)][, .SD[1], by = contig_id]$gene_id

adjacency <- rbindlist(list(
  adjacency_row("GH19\nchitinase-called", gh19_ids),
  adjacency_row("GH24\nlysozyme",         genes[token_hit(cazyme, "GH24"), gene_id]),
  adjacency_row("GH108\nendolysin",       genes[token_hit(cazyme, "GH108"), gene_id]),
  adjacency_row("GH25\nendolysin",        genes[token_hit(cazyme, "GH25"), gene_id]),
  adjacency_row("GH16\nβ-1,3-glucanase",  genes[token_hit(cazyme, "GH16"), gene_id]),
  adjacency_row("Random gene\nsame contigs", null_ids)
))
fwrite(adjacency, file.path(results_dir, "07_cassette_adjacency.tsv"), sep = "\t")

cat("\n--- panel C: adjacency to a holin/spanin (<=3 genes) ----------------\n")
print(adjacency)

# Fisher tests against GH19
gh19_row <- adjacency[set %like% "GH19"]
fisher_res <- rbindlist(lapply(adjacency[!(set %like% "GH19"), set], function(s) {
  r <- adjacency[set == s]
  ft <- fisher.test(matrix(c(gh19_row$n_adjacent,
                             gh19_row$n_with_cassette - gh19_row$n_adjacent,
                             r$n_adjacent, r$n_with_cassette - r$n_adjacent),
                           nrow = 2, byrow = TRUE))
  data.table(comparison = paste("GH19 vs", gsub("\n", " ", s)),
             odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}))
cat("\n--- Fisher tests ----------------------------------------------------\n")
print(fisher_res)

adjacency[, set := factor(set, levels = adjacency$set)]
adjacency[, role := fcase(
  set %like% "GH19",   "GH19",
  set %like% "Random", "Null",
  set %like% "GH16",   "Non-lysis CAZyme control",
  default =            "Endolysin control")]

role_cols <- c("GH19" = "#B2182B", "Endolysin control" = "#762A83",
               "Non-lysis CAZyme control" = "#2166AC", "Null" = "#BABABA")
adjacency[, role := factor(role, levels = names(role_cols))]

pC <- ggplot(adjacency, aes(x = set, y = pct_adjacent, fill = role)) +
  geom_col(width = 0.66, colour = "grey25", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.16,
                linewidth = 0.35, colour = "grey25") +
  geom_text(aes(y = ci_hi + 4,
                label = sprintf("%d/%d", n_adjacent, n_with_cassette)),
            size = 2.7, colour = "grey25") +
  scale_fill_manual(values = role_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0))) +
  labs(tag = "C", x = "Gene / control set",
       y = "% within 3 genes of a holin/spanin\n(of genes on contigs carrying one)") +
  theme_minimal(base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 7),
        legend.position = "top",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7.5),
        panel.grid.major.x = element_blank())

# "Separate PG hydrolase" means one OTHER than the focal GH19 gene, so the focal
# genes are excluded from has_pgh/has_lysis (four GH19 genes independently satisfy
# a PGH Pfam/KO and would otherwise self-count). Holin/spanin markers are never
# GH19, so is_cassette_marker needs no exclusion.
ctg_flags <- genes[, .(
  has_pgh      = any(is_pgh & !(gene_id %chin% gh19_ids)),
  has_lysis    = any(is_lysis_any & !(gene_id %chin% gh19_ids)),
  has_cassette = any(is_cassette_marker)), by = contig_id]
ctg <- merge(contigs, ctg_flags, by = "contig_id", all.x = TRUE)
ctg[is.na(has_pgh), `:=`(has_pgh = FALSE, has_lysis = FALSE, has_cassette = FALSE)]
ctg[, is_gh19 := contig_id %chin% gh19_contigs]

# right = FALSE makes the bins left-closed [lo, hi), so a contig of exactly 10 kb
# falls in the >=10 kb bin, matching the >= threshold in 01_contig_length_stratification.R.
# The last bin is [50 kb, Inf), so its label is >=50 (it includes exactly 50 kb).
LEN_BREAKS <- c(0, 10e3, 20e3, 30e3, 50e3, Inf)
LEN_LABS   <- c("<10", "10-20", "20-30", "30-50", ">=50")
ctg[, len_bin := cut(contig_length, breaks = LEN_BREAKS, labels = LEN_LABS,
                     right = FALSE)]

len_matched <- rbind(
  ctg[, .(set = "All viral contigs", n = .N, pgh = sum(has_pgh)), by = len_bin],
  ctg[is_gh19 == TRUE, .(set = "GH19-bearing contigs", n = .N,
                         pgh = sum(has_pgh)), by = len_bin]
)
len_matched <- len_matched[!is.na(len_bin) & n > 0]
len_matched[, pct := 100 * pgh / n]
len_matched[, c("ci_lo", "ci_hi") := {
  ci <- mapply(function(x, nn) binom.test(x, nn)$conf.int, pgh, n)
  list(100 * ci[1, ], 100 * ci[2, ])
}]
len_matched[, len_bin := factor(len_bin, levels = LEN_LABS)]
fwrite(len_matched, file.path(results_dir, "07_length_matched_pgh.tsv"), sep = "\t")

cat("\n--- panel D: separate PG hydrolase by contig length -----------------\n")
print(len_matched)

# Same comparison restricted to the highest-confidence genomes, where a missing
# lysis gene is least likely to be an assembly artefact.
hq <- ctg[checkv_quality %chin% c("High-quality", "Complete")]
cat(sprintf(paste("\n  CheckV high-quality/complete contigs: all n=%d, %.1f%% carry a",
                  "PG hydrolase; GH19-bearing n=%d, %.1f%%\n"),
            nrow(hq), 100 * mean(hq$has_pgh),
            nrow(hq[is_gh19 == TRUE]), 100 * mean(hq[is_gh19 == TRUE]$has_pgh)))
hq_row <- data.table(
  set = c("All viral contigs", "GH19-bearing contigs"),
  len_bin = "CheckV HQ/Complete",
  n = c(nrow(hq), nrow(hq[is_gh19 == TRUE])),
  pgh = c(sum(hq$has_pgh), sum(hq[is_gh19 == TRUE]$has_pgh)))
hq_row[, pct := 100 * pgh / n]
fwrite(hq_row, file.path(results_dir, "07_checkv_hq_pgh.tsv"), sep = "\t")

# GH19 contig context breakdown (reported in Results, tabulated here)
gh19_ctx <- ctg[is_gh19 == TRUE, .N, by = .(has_pgh, has_cassette)]
gh19_ctx[, context := fcase(
  !has_pgh & !has_cassette, "GH19 only: no other lysis gene detected",
  !has_pgh &  has_cassette, "Holin/spanin present, no other PG hydrolase",
   has_pgh & !has_cassette, "Separate PG hydrolase, no holin/spanin",
   has_pgh &  has_cassette, "Both a separate PG hydrolase and a holin/spanin"
)]
gh19_ctx[, pct := round(100 * N / sum(N), 1)]
setorder(gh19_ctx, -N)
fwrite(gh19_ctx[, .(context, N, pct)],
       file.path(results_dir, "07_gh19_contig_context.tsv"), sep = "\t")
cat("\n--- GH19 contig context (n = ", nrow(ctg[is_gh19 == TRUE]), ") ---\n", sep = "")
print(gh19_ctx[, .(context, N, pct)])

# Taxonomic spread, to show the signal is not one clade
gh19_tax <- ctg[is_gh19 == TRUE]
gh19_tax <- merge(gh19_tax[, .(contig_id)],
                  contigs[, .(contig_id, votu, genus_cluster, family_cluster)],
                  by = "contig_id")
cat(sprintf("\nGH19 contigs: %d contigs, %d vOTUs, %d genus clusters, %d family clusters\n",
            nrow(gh19_tax), uniqueN(gh19_tax$votu),
            uniqueN(gh19_tax$genus_cluster), uniqueN(gh19_tax$family_cluster)))

# Protein lengths, for the Results text
len_tab <- rbind(
  data.table(set = "GH19 (called chitinase)",
             aa = genes[gene_id %chin% gh19_ids, aa_len]),
  genes[token_hit(cazyme, c("GH18", "GH24", "GH108", "GH16")),
        .(set = cazyme, aa = aa_len)]
)[, .(n = .N, median_aa = round(median(aa)),
      q10 = round(quantile(aa, .1)), q90 = round(quantile(aa, .9))), by = set]
cat("\n--- protein length ---------------------------------------------------\n")
print(len_tab)

len_matched[, set := factor(set, levels = c("All viral contigs",
                                            "GH19-bearing contigs"))]

pD <- ggplot(len_matched, aes(x = len_bin, y = pct, fill = set)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.68,
           colour = "grey25", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                position = position_dodge(width = 0.74), width = 0.16,
                linewidth = 0.35, colour = "grey25") +
  geom_text(aes(y = ci_hi + 1.6, label = paste0("n=", n)),
            position = position_dodge(width = 0.74), size = 2.3,
            colour = "grey35") +
  scale_fill_manual(values = c("All viral contigs"    = "#BABABA",
                               "GH19-bearing contigs" = "#B2182B"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(tag = "D", x = "Contig length (kb)",
       y = "% of contigs carrying a separate\npeptidoglycan hydrolase") +
  theme_minimal(base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 14),
        legend.position = "top",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 8),
        panel.grid.major.x = element_blank())

fig2 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1))

ggsave(file.path(fig_dir, "07_Figure2_composite.pdf"), fig2,
       width = 13, height = 9.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "07_Figure2_composite.tiff"), fig2,
       width = 13, height = 9.5, dpi = 300, compression = "lzw", type = "cairo")
cat("\nSaved: figures/07_Figure2_composite.{pdf,tiff}\n")

wb <- createWorkbook()
add <- function(name, x) {
  addWorksheet(wb, name); writeData(wb, name, x); freezePane(wb, name, firstRow = TRUE)
}
add("A_lysis_filter_effect",  dcast(lysis_effect, enzyme_group ~ scenario,
                                    value.var = c("n", "pct")))
add("B_chitinase_families",   chit_decomp)
add("B_all_carbon_families",  fam_decomp)
add("C_cassette_adjacency",   adjacency[, .(set = gsub("\n", " ", set), n_genes,
                                            n_with_cassette, n_adjacent,
                                            pct_adjacent, ci_lo, ci_hi,
                                            median_dist)])
add("C_fisher_tests",         fisher_res)
add("D_length_matched_pgh",   len_matched)
add("D_gh19_contig_context",  gh19_ctx[, .(context, N, pct)])
add("protein_lengths",        len_tab)
add("lysis_gene_definitions",
    rbind(data.table(type = "PG hydrolase Pfam",  id = PGH_PFAM),
          data.table(type = "PG hydrolase CAZy",  id = PGH_CAZY),
          data.table(type = "PG hydrolase KEGG",  id = PGH_KO),
          data.table(type = "Holin Pfam",         id = HOLIN_PFAM),
          data.table(type = "Spanin Pfam",        id = SPANIN_PFAM),
          data.table(type = "Non-catalytic CBM",  id = CBM_NONCATALYTIC),
          data.table(type = "Dual-use CAZy",      id = DUAL_USE_CAZY)))
saveWorkbook(wb, file.path(results_dir, "07_Table_S_lysis_decomposition.xlsx"),
             overwrite = TRUE)
cat("Saved: results/07_Table_S_lysis_decomposition.xlsx\n")

mfd_write_session_info("07", results_dir)
cat("\nDone.\n")
