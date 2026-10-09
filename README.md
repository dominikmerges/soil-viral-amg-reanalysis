# Soil viral AMG analysis — Merges 2026

Analysis code for a re-analysis of the Global Soil Virus (GSV) Atlas
(Graham et al. 2024) quantifying the frequency and functional distribution of
virus-encoded metabolic genes in soils.

No sequencing data were generated. Functional category assignment uses
MetaFuncDecoder v1.0.0 (DOI:10.5281/zenodo.19009291).

## Input data

Two sources. Place everything under `data/`.

**1. GSV Atlas** — download from https://doi.org/10.25584/2229733 and keep the
original filenames:

| File | Contents |
|---|---|
| `GSVA_soil_viruses_genome_metadata_2.tsv.gz` | File 2, 49,649 QA/QC contigs |
| `GSVA_soil_viruses_gene_metadata_4.tsv.gz` | File 4, 1,432,147 genes |
| `GSVA_sample_metadata_5.csv` | File 5, 2,953 samples |

**2. Derived inputs** — download from [FILL: Figshare DOI]:

| File | Contents |
|---|---|
| `metagenome_confidence_annotations.csv` | MetaFuncDecoder v1.0.0 output |
| `jgi_targeted_counts.csv` | Targeted gene counts, six matched JGI studies |
| `jgi_study_metadata.csv` | Per-study sequencing stats and management history (S11/S12); values from the JGI assembly READMEs |
| `original_keyword_counts.csv` | Original-submission manual keyword-matching counts (tool-comparison table in script 03) |
| `jgi_readmes/` | JGI assembly READMEs (portal download 1977383); provenance for `jgi_study_metadata.csv`, not read by the scripts |

## Running

```
Rscript run_all.R
```

Scripts are also runnable individually and in any order, except that all
require `data/` to be populated, and 02/03/04/06/07 additionally require
`data_10kb/` (run `00_filter_10kb.R` first). Outputs go to `results/`.

## Contig-length threshold

The primary analysis is restricted to viral contigs ≥10 kb. `00_filter_10kb.R`
writes a filtered copy of the three GSV Atlas input tables to `data_10kb/`
(31,344 contigs; 1,238,728 gene calls; 1,093 samples), and scripts 02, 03, 04
and 07 read from it. `01_contig_length_stratification.R` reads the full `data/`
and provides the contig-length sensitivity analysis across the 1–5, 5–10 and
≥10 kb bins. `06_jgi_viral_contribution.R` (Figure 3) also reads `data_10kb/`:
`00_filter_10kb.R` restricts its viral gene counts to viral contigs ≥10 kb (the
six matched studies' viral contigs are in the Atlas, keyed by IMG taxon-OID
prefix), while the `total_kegg_genes` denominator — the whole assembled,
bacterial-dominated metagenome — is left unchanged, since the ≥10 kb criterion
concerns viral identification only.

## Scripts

| Script | Produces |
|---|---|
| `_common.R` | Shared loading and filtering; sourced by 01, 02, 03, 07 |
| `00_filter_10kb.R` | Builds `data_10kb/` (contigs ≥10 kb) read by the primary scripts 02/03/04/07 |
| `01_contig_length_stratification.R` | Fig. S1 |
| `02_database_decomposition.R` | Fig. S2, Tables S1, S3 |
| `03_database_support.R` | Table 1, Table S14 |
| `04_amg_distribution.R` | Figure 1 (all panels), Tables S2, S4, S13, S16 |
| `05_jgi_study_metadata.R` | Tables S11, S12 |
| `06_jgi_viral_contribution.R` | Figure 3 |
| `07_lysis_gene_decomposition.R` | Figure 2, Tables S5–S10, S15 |

## Annotation filtering

`_common.R` defines the filtering chain applied to the MetaFuncDecoder output:

```
raw          5,760   all calls
  - peptidoglycanase (lysis genes, not metabolic)
post_lysis   1,858   any annotation        0.13% of the catalogue
  - broad matches (no enzyme resolved)
specific     1,169   enzyme-resolved       0.08% of the catalogue
```

`specific` is the analysis set for all reported results. `post_lysis` is used
only for the annotated-at-all rate and for the annotation-method comparison in
`03_database_support.R`, where filtering one arm only would confound the
change of tool with the change of filter.

The counts above are the full QA/QC catalogue (read by `01`). Under the ≥10 kb
primary restriction (`data_10kb/`, read by 02/03/04/07) the same chain gives
raw 5,020 → post_lysis 1,698 → specific 1,067. The fail-fast guard in
`_common.R` accepts either sanctioned state and stops on any other.

## Environment

R 4.1.2 with `data.table` 1.14.8, `tidyverse` 2.0.0, `ggplot2` 3.4.2,
`gridExtra` 2.3, `openxlsx`, `patchwork`, `RColorBrewer`, `bit64`.
`sessionInfo()` for each script is written to `results/`.

## Citation

[FILL: manuscript citation once available]

## License

MIT — see `LICENSE`.
