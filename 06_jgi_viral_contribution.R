#!/usr/bin/env Rscript
# 06_jgi_viral_contribution.R
#
# Figure 3: viral share of targeted functional annotations in six matched JGI
# studies. Panel A gives total and viral gene counts per targeted function on a
# log scale; panel B gives the viral percentage.
#
# Styling follows 04_amg_distribution.R (theme_minimal, base 11, Brewer Set2).
#
# Viral gene counts are restricted to viral contigs >=10 kb (built by
# 00_filter_10kb.R into data_10kb/jgi_targeted_counts.csv); the total_kegg_genes
# denominator is the whole assembled metagenome and is unchanged.
#
# In:  data_10kb/jgi_targeted_counts.csv
# Out: results/figures/06_Figure3_composite.pdf
#      results/figures/06_Figure3_composite.tiff
#      results/06_sessionInfo.txt

library(ggplot2)
library(patchwork)
library(RColorBrewer)

set.seed(42)

results_dir <- "results"
fig_dir     <- file.path(results_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv("data_10kb/jgi_targeted_counts.csv", comment.char = "#",
                   stringsAsFactors = FALSE)

dat <- data.frame(
  Study      = counts$study_label,
  JGI_ID     = as.character(counts$img_taxon_oid),
  Function   = counts$function_label,
  Category   = counts$functional_category,
  Total_KEGG = counts$total_kegg_genes,
  Viral_KEGG = counts$viral_kegg_genes,
  stringsAsFactors = FALSE
)
dat$Viral_pct <- round(100 * dat$Viral_KEGG / dat$Total_KEGG, 2)

# Wrap long function labels onto two lines for the x axis.
wrap_label <- function(x) {
  x <- sub("^(Antibiotic) (inactivation)$", "\\1\n\\2", x)
  x <- sub(" (\\+ target alteration)$", "\n\\1", x)
  x <- sub("^(Ammonia) (assimilation)$", "\\1\n\\2", x)
  x <- sub("^(Nitrogen) (fixation)$", "\\1\n\\2", x)
  x
}

dat$x_label <- paste0(dat$Study, "\n(", wrap_label(dat$Function), ")")
dat$x_label <- factor(dat$x_label, levels = dat$x_label)

# Colour palette (same as Figure 1)

brewer_set2 <- brewer.pal(3, "Set2")
cat_colours <- c(
  "carbon_cycling"   = brewer_set2[1],
  "nitrogen_cycling"  = brewer_set2[2],
  "antibiotics"       = brewer_set2[3]
)

# Reshape to long format for grouped bars
dat_long <- rbind(
  data.frame(x_label = dat$x_label, Category = dat$Category,
             Type = "Total KEGG", Count = dat$Total_KEGG,
             stringsAsFactors = FALSE),
  data.frame(x_label = dat$x_label, Category = dat$Category,
             Type = "Viral KEGG", Count = dat$Viral_KEGG,
             stringsAsFactors = FALSE)
)
dat_long$Type <- factor(dat_long$Type, levels = c("Total KEGG", "Viral KEGG"))

log_baseline <- 0.8
n_type      <- 2
dodge_total <- 0.7
bar_w       <- dodge_total / n_type * 0.9
dat_long$x_idx    <- as.numeric(dat_long$x_label)
dat_long$type_idx <- as.numeric(dat_long$Type)
dat_long$xpos     <- dat_long$x_idx +
  (dat_long$type_idx - (n_type + 1) / 2) * dodge_total / n_type

# Zero-safe log axis: a viral count of 0 (two studies after the >=10 kb filter)
# is drawn as a zero-height bar at the baseline but still labelled "0".
dat_long$Count_plot <- ifelse(dat_long$Count < log_baseline, log_baseline, dat_long$Count)
dat_long$lab_y      <- ifelse(dat_long$Count < log_baseline, log_baseline, dat_long$Count)

p_a <- ggplot(dat_long) +
  geom_rect(aes(xmin = xpos - bar_w / 2, xmax = xpos + bar_w / 2,
                ymin = log_baseline, ymax = Count_plot, fill = Type)) +
  geom_text(aes(x = xpos, y = lab_y, label = Count),
            vjust = -0.3, size = 3, colour = "grey30") +
  scale_x_continuous(breaks = seq_along(levels(dat_long$x_label)),
                     labels = levels(dat_long$x_label)) +
  scale_y_log10(
    limits = c(log_baseline, 2000),
    breaks = c(1, 10, 100, 1000),
    labels = scales::comma
  ) +
  scale_fill_manual(
    values = c("Total KEGG" = "#6BAED6", "Viral KEGG" = "#956B9E"),
    name = NULL
  ) +
  labs(
    tag = "A",
    title = NULL,
    x = "Matched study (targeted function)",
    y = "Gene Count (log scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    plot.title = element_text(face = "bold", size = 12),
    plot.tag = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# Colour each bar by functional category
p_b <- ggplot(dat, aes(x = x_label, y = Viral_pct, fill = Category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(Viral_pct, "%")),
            vjust = -0.3, size = 3, colour = "grey30") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    limits = c(0, NA)
  ) +
  scale_fill_manual(
    values = cat_colours,
    labels = c("carbon_cycling" = "Carbon Cycling",
               "nitrogen_cycling" = "Nitrogen Cycling",
               "antibiotics" = "Antibiotic Resistance"),
    name = NULL
  ) +
  labs(
    tag = "B",
    title = NULL,
    x = "Matched study (targeted function)",
    y = "Viral Contribution (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    plot.title = element_text(face = "bold", size = 12),
    plot.tag = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

composite <- p_a + p_b + plot_layout(widths = c(1, 1))

# TIFF at 300 dpi (AEM preferred)
tiff(file.path(fig_dir, "06_Figure3_composite.tiff"),
     width = 12, height = 5.5, units = "in", res = 300,
     compression = "lzw")
print(composite)
dev.off()
cat("Saved:", file.path(fig_dir, "06_Figure3_composite.tiff"), "\n")

# PDF for convenience
pdf(file.path(fig_dir, "06_Figure3_composite.pdf"),
    width = 12, height = 5.5)
print(composite)
dev.off()
cat("Saved:", file.path(fig_dir, "06_Figure3_composite.pdf"), "\n")

sink(file.path(results_dir, "06_sessionInfo.txt"))
cat("06_jgi_viral_contribution.R\n")
cat("Run:", format(Sys.time()), "\n\n")
sessionInfo()
sink()

cat("Done.\n")
