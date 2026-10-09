#!/usr/bin/env Rscript
# 00_filter_10kb.R
#
# Builds data_10kb/: the >=10 kb-filtered GSV Atlas catalogue used by the
# primary analysis (scripts 02, 03, 04, 07). A single contig-length threshold
# is applied once, here, to the input tables; the downstream scripts are
# otherwise unchanged and simply read from data_10kb/.
#
# Rationale: the GSV Atlas applied different viral-identification criteria above
# and below 10 kb, and MIUViG guidance cautions against interpreting <10 kb
# viral sequences; the primary inventory is therefore restricted to >=10 kb
# contigs. The full catalogue is retained in data/ and is read by
# 01_contig_length_stratification.R as the contig-length sensitivity analysis.
# Script 06 (JGI viral contribution, Figure 3) reads data_10kb/: its viral gene
# counts are restricted here to viral contigs >=10 kb (the six matched studies'
# viral contigs are in the GSV Atlas, keyed by IMG taxon-OID prefix, so the same
# contig_length used above also filters the Figure 3 numerator). The
# total_kegg_genes denominator is the whole assembled metagenome and is left
# unchanged, as the >=10 kb criterion concerns viral identification only.
#
# In:  data/GSVA_soil_viruses_genome_metadata_2.tsv.gz
#      data/GSVA_soil_viruses_gene_metadata_4.tsv.gz
#      data/metagenome_confidence_annotations.csv
#      (+ length-agnostic inputs copied through unchanged)
# Out: data_10kb/ with the same filenames.

suppressPackageStartupMessages(library(data.table))

in_dir  <- "data"
out_dir <- "data_10kb"
THRESH  <- 10000L
dir.create(out_dir, showWarnings = FALSE)

genome <- fread(file.path(in_dir, "GSVA_soil_viruses_genome_metadata_2.tsv.gz"))
genes  <- fread(file.path(in_dir, "GSVA_soil_viruses_gene_metadata_4.tsv.gz"))
ann    <- fread(file.path(in_dir, "metagenome_confidence_annotations.csv"))

keep_contigs <- genome[contig_length >= THRESH, contig_id]
genome_f <- genome[contig_length >= THRESH]
genes_f  <- genes[contig_id %in% keep_contigs]
ann_f    <- ann[gene_id %in% genes_f$gene_id]

# Guard: every retained annotation and gene must sit on a retained contig.
stopifnot(
  all(genome_f$contig_length >= THRESH),
  all(genes_f$contig_id %in% keep_contigs),
  all(ann_f$gene_id %in% genes_f$gene_id)
)

fwrite(genome_f, file.path(out_dir, "GSVA_soil_viruses_genome_metadata_2.tsv.gz"), sep = "\t")
fwrite(genes_f,  file.path(out_dir, "GSVA_soil_viruses_gene_metadata_4.tsv.gz"),  sep = "\t")
fwrite(ann_f,    file.path(out_dir, "metagenome_confidence_annotations.csv"))

# --- Figure 3 (JGI): restrict viral counts to viral contigs >=10 kb ---------
# The six matched studies' viral contigs are in the Atlas, keyed by the IMG
# taxon-OID prefix of contig_id; viral_kegg_genes is recomputed as the number of
# target-KO genes on a viral contig >=10 kb. total_kegg_genes (whole metagenome)
# is left unchanged. A guard verifies the unfiltered recount reproduces the
# transcribed viral_kegg_genes in data/jgi_targeted_counts.csv before filtering.
jgi <- read.csv(file.path(in_dir, "jgi_targeted_counts.csv"),
                comment.char = "#", stringsAsFactors = FALSE)
genes[, prefix := sub("[.:].*$", "", contig_id)]
vcount <- function(oid, ko, ge10) {
  sel <- genes[prefix == oid & kegg_ortholog == ko]
  if (ge10) sum(sel$contig_id %in% keep_contigs) else nrow(sel)
}
viral_all  <- mapply(function(o, k) vcount(as.character(o), k, FALSE),
                     jgi$img_taxon_oid, jgi$ko_term)
stopifnot(all(viral_all == jgi$viral_kegg_genes))   # reproduces the manuscript
jgi$viral_kegg_genes <- mapply(function(o, k) vcount(as.character(o), k, TRUE),
                               jgi$img_taxon_oid, jgi$ko_term)
jgi_hdr <- c(
  "# Gene counts for the six matched JGI studies, per targeted function.",
  "# total_kegg_genes = genes with the given KO in the whole assembled metagenome (not length-filtered).",
  "# viral_kegg_genes = subset on a viral contig >=10 kb (GSVA genome metadata contig_length; 2nd revision).",
  "# Viral percentage is derived, not stored.")
jgi_out <- file.path(out_dir, "jgi_targeted_counts.csv")
writeLines(jgi_hdr, jgi_out)
suppressWarnings(
  write.table(jgi, jgi_out, sep = ",", row.names = FALSE, col.names = TRUE,
              quote = FALSE, append = TRUE))
cat(sprintf("JGI viral (>=10 kb)     %s -> %s\n",
            paste(viral_all, collapse = "/"),
            paste(jgi$viral_kegg_genes, collapse = "/")))

# Length-agnostic inputs copied through unchanged.
for (f in c("GSVA_sample_metadata_5.csv",
            "jgi_study_metadata.csv", "original_keyword_counts.csv")) {
  src <- file.path(in_dir, f)
  if (file.exists(src)) file.copy(src, file.path(out_dir, f), overwrite = TRUE)
}

cat(sprintf("contigs      %7d -> %7d\n", nrow(genome), nrow(genome_f)))
cat(sprintf("genes        %7d -> %7d\n", nrow(genes),  nrow(genes_f)))
cat(sprintf("annotations  %7d -> %7d\n", nrow(ann),    nrow(ann_f)))
cat(sprintf("samples(>=1 viral gene) %d -> %d\n",
            uniqueN(sub("\\..*", "", genes$gene_id)),
            uniqueN(sub("\\..*", "", genes_f$gene_id))))
