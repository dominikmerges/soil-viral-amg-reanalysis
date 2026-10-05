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
# Script 06 (JGI viral contribution) is length-agnostic and also reads data/.
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

# Length-agnostic inputs copied through unchanged.
for (f in c("GSVA_sample_metadata_5.csv", "jgi_targeted_counts.csv",
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
