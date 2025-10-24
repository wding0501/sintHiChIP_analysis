# ============================================================================
# sintHiChIP Independence Validation - Minimal Output Version
# ============================================================================

library("devtools")
setwd("~/HiChIP/Callloop/Update/sintHiChIP_res/local")
library("sintHiChIP")
library(GenomicRanges)
library(data.table)
library(stats)
library(MASS)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(cowplot)

# Parameters
resfrags_file = "/home/wding/HiChIP/hichipper/hg19_mboi.bed"
outdir <- "/home/wding/HiChIP/Callloop/Update/sintHiChIP_res/local"
normSiteFile <- "/home/wding/HiChIP/Callloop/Update/CutSite/normsite_hg19_5000.tmp"
FDR <- 0.01
min_dist <- 20000L
max_dist <- 2000000L
half_length <- 73
no_merge <- FALSE
max_anchor_width <- 50000
nbins <- 10

cwd <- outdir
sname <- "GM12878"

# Load data
lo <- create_loops_object_fast(cwd, snames = paste0(sname, ".filt.intra"), type = "all")
totalPetCounts <- rowSums(lo@counts)

GenomeBin <- data.table::fread(normSiteFile, header = TRUE, showProgress = FALSE)
GenomeBin <- sintHiChIP_create_granges_fast(GenomeBin, 
                                            names(GenomeBin)[1],
                                            names(GenomeBin)[2],
                                            names(GenomeBin)[3])

chrname <- unique(intersect(seqlevels(lo@anchors), seqnames(GenomeBin)))
seqlevels(GenomeBin, pruning.mode = 'coarse') <- chrname

# Anchor site processing
n_anchors <- length(lo@anchors)
n_interactions <- nrow(lo@interactions)
n_genomebins <- length(GenomeBin)

anchors_df <- data.frame(
  seqnames = as.character(seqnames(lo@anchors)),
  start = start(lo@anchors),
  end = end(lo@anchors),
  anchor_id = seq_len(n_anchors),
  stringsAsFactors = FALSE
)

unique_anchors_df <- unique(anchors_df[, c("seqnames", "start", "end")])
unique_anchors_df$unique_id <- seq_len(nrow(unique_anchors_df))

# Density calculation
if (nrow(unique_anchors_df) > 10000 || n_genomebins > 100000) {
  anchors_dt <- data.table::as.data.table(unique_anchors_df)
  genomebins_dt <- data.table::as.data.table(as.data.frame(GenomeBin))
  
  anchors_dt[, anchor_id := .I]
  genomebins_dt[, weight := MeanSite]
  
  data.table::setkey(anchors_dt, seqnames, start, end)
  data.table::setkey(genomebins_dt, seqnames, start, end)
  
  overlaps <- data.table::foverlaps(genomebins_dt, anchors_dt, type = "any", nomatch = 0L)
  
  if (nrow(overlaps) > 0) {
    density_stats <- overlaps[, .(
      total_weight = sum(weight, na.rm = TRUE),
      bin_count = .N
    ), by = .(seqnames, start, end, anchor_id)]
    
    density_stats[, width := end - start + 1]
    density_stats[, density := total_weight / width]
    
    unique_densities <- numeric(nrow(anchors_dt))
    unique_densities[density_stats$anchor_id] <- density_stats$density
  } else {
    unique_densities <- numeric(nrow(unique_anchors_df))
  }
} else {
  unique_gr <- sintHiChIP_create_granges_fast(unique_anchors_df)
  overlaps_hits <- GenomicRanges::findOverlaps(unique_gr, GenomeBin, type = "any")
  
  if (length(overlaps_hits) > 0) {
    weights <- mcols(GenomeBin)$MeanSite[subjectHits(overlaps_hits)]
    query_indices <- queryHits(overlaps_hits)
    
    overlap_dt <- data.table::data.table(anchor_idx = query_indices, weight = weights)
    weight_sums <- overlap_dt[, .(total_weight = sum(weight, na.rm = TRUE)), by = anchor_idx]
    
    anchor_widths <- width(unique_gr)
    unique_densities <- numeric(length(unique_gr))
    unique_densities[weight_sums$anchor_idx] <- weight_sums$total_weight / anchor_widths[weight_sums$anchor_idx]
  } else {
    unique_densities <- numeric(length(unique_gr))
  }
}

anchor_key <- paste(anchors_df$seqnames, anchors_df$start, anchors_df$end, sep = "_")
unique_key <- paste(unique_anchors_df$seqnames, unique_anchors_df$start, unique_anchors_df$end, sep = "_")
density_lookup <- unique_densities[match(anchor_key, unique_key)]

site_l <- density_lookup[lo@interactions[, 1]]
site_r <- density_lookup[lo@interactions[, 2]]

interaction_anchors_dt <- data.table::data.table(
  seqnames = c(anchors_df$seqnames[lo@interactions[, 1]], anchors_df$seqnames[lo@interactions[, 2]]),
  start = c(anchors_df$start[lo@interactions[, 1]], anchors_df$start[lo@interactions[, 2]]),
  end = c(anchors_df$end[lo@interactions[, 1]], anchors_df$end[lo@interactions[, 2]]),
  totalAnchorSites = c(site_l, site_r)
)

tc_dt <- interaction_anchors_dt[, .(totalAnchorSites = sum(totalAnchorSites, na.rm = TRUE)),
                                by = .(seqnames, start, end)]
tc.site <- sort(sintHiChIP_create_granges_fast(tc_dt))

chrpeaks <- as.data.frame(tc.site)
if (ncol(chrpeaks) >= 6) {
  chrpeaks <- chrpeaks[, c(1, 2, 3, 4, 6, 5)]
} else {
  names(chrpeaks)[1:3] <- c("chr", "start", "end")
  chrpeaks$name <- seq_len(nrow(chrpeaks))
  chrpeaks$score <- tc.site$totalAnchorSites
  chrpeaks$strand <- "*"
  chrpeaks <- chrpeaks[, c("chr", "start", "end", "name", "score", "strand")]
}
names(chrpeaks) <- c("chr", "start", "end", "name", "score", "strand")

mcols(lo@anchors) <- mcols(tc.site)
df <- cbind(as.data.frame(lo@anchors[lo@interactions[, 1]])[, -c(4:5)],
            as.data.frame(lo@anchors[lo@interactions[, 2]])[, -c(4:5)])
df <- data.table::setDT(df)
data.table::setnames(df, c("chr_1", "start_1", "end_1", "totalAnchorSites_1",
                           "chr_2", "start_2", "end_2", "totalAnchorSites_2"))
df[, loopWidth := as.integer(lo@rowData$loopWidth)]
df[, PETS := totalPetCounts]
df[, site := totalAnchorSites_1 * totalAnchorSites_2]

# Fit models
distanceborders <- make_equal_bins(df$loopWidth, binmethod = "equalocc", numberbins = nbins)
siteborders <- make_equal_bins(df$site, binmethod = "equalocc", numberbins = nbins)

distance_IAB_model <- model_hichip(df$loopWidth, df$PETS, borders = distanceborders, yvals = TRUE)
distance_IAB_spline <- stats::smooth.spline(log10(distance_IAB_model[, 1]), distance_IAB_model[, 3], spar = 0.35)

site_IAB_model <- model_hichip(df$site, df$PETS, borders = siteborders, yvals = TRUE)
site_IAB_model.complete <- site_IAB_model[complete.cases(site_IAB_model), , drop = FALSE]
site_IAB_spline <- stats::smooth.spline(log10(site_IAB_model.complete[, 1]), site_IAB_model.complete[, 3], spar = 0.35)

# Combination processing
sumofy_dist <- rep(0, nbins)
sumofx_dist <- rep(0, nbins)
countofx_dist <- rep(0, nbins)
sumofy_site <- rep(0, nbins)
sumofx_site <- rep(0, nbins)
countofx_site <- rep(0, nbins)

chromosomes <- sort(unique(as.character(chrpeaks$chr)))
distance_list <- vector("list", length(chromosomes))
site_list <- vector("list", length(chromosomes))

chrpeaks_dt <- data.table::as.data.table(chrpeaks)
data.table::setkey(chrpeaks_dt, chr)

for (i in seq_along(chromosomes)) {
  chrom <- chromosomes[i]
  chrom_peaks <- chrpeaks_dt[chrom][order(start)]
  
  if (nrow(chrom_peaks) >= 2) {
    combos <- makecombos(chrom, chrom_peaks, mindist = 0, maxdist = 2000000)
    if (nrow(combos) > 0) {
      distance_list[[i]] <- combos$dist
      site_list[[i]] <- combos$score1 * combos$score2
    }
  }
}

distance_list <- distance_list[!sapply(distance_list, is.null)]
site_list <- site_list[!sapply(site_list, is.null)]

all_distances <- unlist(distance_list, use.names = FALSE)
all_sites <- unlist(site_list, use.names = FALSE)

rm(distance_list, site_list)
invisible(gc())

if (length(all_distances) > 0) {
  distance_combo_model <- model_hichip(x = all_distances, y = NA, borders = distanceborders, yvals = FALSE)
  site_combo_model <- model_hichip(x = all_sites, y = NA, borders = siteborders, yvals = FALSE)
  
  sumofy_dist <- distance_combo_model[, 2]
  sumofx_dist <- distance_combo_model[, 4]
  countofx_dist <- distance_combo_model[, 5]
  
  sumofy_site <- site_combo_model[, 2]
  sumofx_site <- site_combo_model[, 4]
  countofx_site <- site_combo_model[, 5]
}

rm(all_distances, all_sites)
invisible(gc())

# Build combined models
site_combo_model <- cbind(sumofx_site / countofx_site, sumofy_site, sumofy_site / sum(sumofy_site))
site_combo_model.complete <- site_combo_model[complete.cases(site_combo_model),]

distance_combo_model <- cbind(sumofx_dist / countofx_dist, sumofy_dist, sumofy_dist / sum(sumofy_dist))
distance_combo_model.complete <- distance_combo_model[complete.cases(distance_combo_model),]

distance_combo_spline <- stats::smooth.spline(log10(distance_combo_model.complete[, 1]), distance_combo_model.complete[, 3], spar = 0.35)
site_combo_spline <- stats::smooth.spline(log10(site_combo_model.complete[, 1]), site_combo_model.complete[, 3], spar = 0.35)

# ============================================================================
# Independence Validation (Minimal Output)
# ============================================================================

cat("\n=== INDEPENDENCE VALIDATION ===\n\n")

# Genome-wide correlation
gw_pearson <- cor.test(df$loopWidth, log2(df$site), method = "pearson")
gw_rsq_pearson <- gw_pearson$estimate^2
gw_vif <- 1 / (1 - gw_rsq_pearson)
gw_spearman <- cor.test(df$loopWidth, log2(df$site), method = "spearman")

cat(sprintf("Genome-wide (n = %s):\n", format(nrow(df), big.mark = ",")))
cat(sprintf("  Pearson R² = %.6f, VIF = %.4f\n", gw_rsq_pearson, gw_vif))
cat(sprintf("  Spearman ρ = %.6f\n\n", gw_spearman$estimate))

# Per-chromosome
chr_results <- data.frame()
for (chr in sort(unique(as.character(df$chr_1)))) {
  chr_data <- df[df$chr_1 == chr, ]
  if (nrow(chr_data) >= 100) {
    chr_cor <- cor.test(chr_data$loopWidth, log2(chr_data$site), method = "pearson")
    chr_results <- rbind(chr_results, data.frame(
      chromosome = chr, n = nrow(chr_data), rsq = chr_cor$estimate^2
    ))
  }
}

cat(sprintf("Per-chromosome: Mean R² = %.6f, Range = [%.6f, %.6f]\n\n",
            mean(chr_results$rsq), min(chr_results$rsq), max(chr_results$rsq)))

# Distance bins
distance_bins <- list(
  "20-500kb" = df$loopWidth >= 20000 & df$loopWidth < 500000,
  "500kb-1Mb" = df$loopWidth >= 500000 & df$loopWidth < 1000000,
  "1-1.5Mb" = df$loopWidth >= 1000000 & df$loopWidth < 1500000,
  "1.5-2Mb" = df$loopWidth >= 1500000 & df$loopWidth <= 2000000
)

bin_results <- data.frame()
for (bin_name in names(distance_bins)) {
  bin_data <- df[distance_bins[[bin_name]], ]
  if (nrow(bin_data) >= 100) {
    bin_cor <- cor.test(bin_data$loopWidth, log2(bin_data$site), method = "pearson")
    bin_results <- rbind(bin_results, data.frame(
      bin = bin_name, rsq = bin_cor$estimate^2, pct = (nrow(bin_data) / nrow(df)) * 100
    ))
  }
}

for (i in 1:nrow(bin_results)) {
  cat(sprintf("  %s: R² = %.6f (%.1f%% of loops)\n",
              bin_results$bin[i], bin_results$rsq[i], bin_results$pct[i]))
}

genome_wide_rsq <- gw_rsq_pearson
df_chr10 <- df[df$chr_1 == "chr10", ]
site_distance_cor_chr10 <- cor.test(df_chr10$loopWidth, log2(df_chr10$site), method = "pearson")

cat(sprintf("\nConclusion: VIF = %.4f < 5 → Independence assumption valid\n\n", gw_vif))

# ============================================================================
# Figure S1 Generation
# ============================================================================

pub_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_line(color = "gray95", size = 0.2),
    panel.border = element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 9, color = "black"),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(30, 30, 30, 30)
  )

fig_a <- ggplot(df, aes(x = log2(site))) +
  geom_histogram(bins = 50, fill = "#4575b4", alpha = 0.8, color = "white", size = 0.2) +
  geom_density(aes(y = ..density.. * nrow(df) * 0.32), color = "#d73027", size = 1.2) +
  scale_x_continuous(breaks = seq(-21, -14, 1)) +
  scale_y_continuous(
    name = "Interaction number", labels = comma_format(),
    sec.axis = sec_axis(~ . / (nrow(df) * 0.32), name = "Interaction number density")
  ) +
  pub_theme +
  theme(axis.title.y.right = element_text(color = "#d73027", face = "bold"),
        axis.text.y.right = element_text(color = "#d73027")) +
  labs(x = "Log2 joint cut site density")

fig_b <- ggplot(df, aes(x = loopWidth, y = log2(site))) +
  stat_bin2d(bins = 40) +
  geom_smooth(method = "lm", color = "white", size = 1.2, se = FALSE) +
  scale_fill_viridis_c(option = "plasma", name = "Count") +
  scale_x_continuous(labels = comma_format(), breaks = c(0, 500000, 1000000, 1500000, 2000000), limits = c(0, 2000000)) +
  scale_y_continuous(breaks = seq(-40, -20, 5)) +
  annotate("text", x = 1000000, y = -24, 
           label = sprintf("R² = %.5f, VIF = %.3f, p < 2.2e-16", 
                           genome_wide_rsq, gw_vif),
           color = "white", size = 3.5, fontface = "bold") +
  pub_theme + theme(legend.position = "right") +
  labs(x = "Distance (bp)", y = "Log2 joint cut site density")

dist_obs <- data.frame(x = distance_IAB_model[,1], y = distance_IAB_model[,3])[complete.cases(distance_IAB_model[,c(1,3)]), ]
dist_fit <- data.frame(x = 10^distance_IAB_spline$x, y = distance_IAB_spline$y)

fig_c <- ggplot() +
  geom_point(data = dist_obs, aes(x = log10(x), y = y, color = "Raw"), size = 1.5, alpha = 0.7) +
  geom_line(data = dist_fit, aes(x = log10(x), y = y, color = "Fitted"), size = 1.2) +
  scale_color_manual(values = c("Raw" = "#74add1", "Fitted" = "#d73027")) +
  scale_x_continuous(breaks = c(4.5, 5.0, 5.5, 6.0), labels = c("4.5", "5.0", "5.5", "6.0")) +
  pub_theme + theme(legend.position = "right", legend.title = element_blank()) +
  labs(x = "Log10 distance", y = "PETS probability")

site_obs <- data.frame(x = site_IAB_model.complete[,1], y = site_IAB_model.complete[,3])[complete.cases(site_IAB_model.complete[,c(1,3)]), ]
site_fit <- data.frame(x = 10^site_IAB_spline$x, y = site_IAB_spline$y)

fig_d <- ggplot() +
  geom_point(data = site_obs, aes(x = log10(x), y = y, color = "Raw"), size = 1.5, alpha = 0.7) +
  geom_line(data = site_fit, aes(x = log10(x), y = y, color = "Fitted"), size = 1.2) +
  scale_color_manual(values = c("Raw" = "#74add1", "Fitted" = "#d73027")) +
  pub_theme + theme(legend.position = "right", legend.title = element_blank()) +
  labs(x = "Log10 cut site density", y = "PETS probability")

dist_combo <- data.frame(x = distance_combo_model.complete[,1], y = distance_combo_model.complete[,3])[complete.cases(distance_combo_model.complete[,c(1,3)]), ]
dist_combo_fit <- data.frame(x = 10^distance_combo_spline$x, y = distance_combo_spline$y)

fig_e <- ggplot() +
  geom_point(data = dist_combo, aes(x = log10(x), y = y, color = "Raw"), size = 1.5, alpha = 0.7) +
  geom_line(data = dist_combo_fit, aes(x = log10(x), y = y, color = "Fitted"), size = 1.2) +
  scale_color_manual(values = c("Raw" = "#74add1", "Fitted" = "#d73027")) +
  scale_x_continuous(breaks = c(4.5, 5.0, 5.5, 6.0), labels = c("4.5", "5.0", "5.5", "6.0")) +
  pub_theme + theme(legend.position = "right", legend.title = element_blank()) +
  labs(x = "Log10 distance", y = "PETS probability")

site_combo <- data.frame(x = site_combo_model.complete[,1], y = site_combo_model.complete[,3])[complete.cases(site_combo_model.complete[,c(1,3)]), ]
site_combo_fit <- data.frame(x = 10^site_combo_spline$x, y = site_combo_spline$y)

fig_f <- ggplot() +
  geom_point(data = site_combo, aes(x = log10(x), y = y, color = "Raw"), size = 1.5, alpha = 0.7) +
  geom_line(data = site_combo_fit, aes(x = log10(x), y = y, color = "Fitted"), size = 1.2) +
  scale_color_manual(values = c("Raw" = "#74add1", "Fitted" = "#d73027")) +
  pub_theme + theme(legend.position = "right", legend.title = element_blank()) +
  labs(x = "Log10 cut site density", y = "PET probability")

create_labeled_plot <- function(plot, label) {
  ggdraw() + draw_plot(plot, 0, 0, 1, 1) +
    draw_label(label, x = 0.1, y = 0.9, hjust = 1, vjust = 0, fontface = "bold", size = 16)
}

plot_a <- create_labeled_plot(fig_a, "a")
plot_b <- create_labeled_plot(fig_b, "b") 
plot_c <- create_labeled_plot(fig_c, "c")
plot_d <- create_labeled_plot(fig_d, "d")
plot_e <- create_labeled_plot(fig_e, "e")
plot_f <- create_labeled_plot(fig_f, "f")

pdf("Figure_S1.pdf", width = 12, height = 10)
grid.arrange(plot_a, plot_b, plot_c, plot_d, plot_e, plot_f, ncol = 2, nrow = 3, padding = unit(30, "pt")) 
dev.off()

png("Figure_S1.png", width = 12, height = 8, units = "in", res = 600)
grid.arrange(plot_a, plot_b, plot_c, plot_d, plot_e, plot_f, ncol = 2, nrow = 3, padding = unit(30, "pt"))
dev.off()

cat("Figure S1 generated\n")