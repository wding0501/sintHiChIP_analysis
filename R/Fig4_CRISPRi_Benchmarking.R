
#!/usr/bin/env Rscript

# ============================================================================
# Complete CRISPRi Analysis - Find sintHiChIP (global) Maximum Advantage Genes
# ============================================================================

library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)

if (!require("PRROC", character.only = TRUE, quietly = TRUE)) {
  install.packages("PRROC", dependencies = TRUE)
  library(PRROC, character.only = TRUE)
}

setwd("/home/wding/HiChIP/Callloop/Update/Result/CRISPR_K562")

# ============================================================================
# Configuration and Functions
# ============================================================================
tolerance <- 5000
q_thresholds <- c(0.01, 0.001, 1e-5, 1e-7)

# Color scheme (local=red, global=green) - CORRECTED
method_colors <- c(
  "sintHiChIP (global)" = "#DC143C",  # Crimson Red  
  "sintHiChIP (local)" = "#228B22",   # Forest Green
  "HiC-DC+" = "#FFD700",             # Gold
  "MAPS" = "#4169E1",                # Royal Blue
  "hichipper" = "#8B4513",           # Saddle Brown
  "FitHiChIP" = "#9932CC",           # Dark Orchid
  "cLoops2" = "#FF8C00",             # Dark Orange
  "cLoops" = "#708090"               # Slate Gray
)

# Nature journal theme
nature_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.line = element_line(color = "black", size = 0.8),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 1),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black")
  )

# auPR calculation using PRROC
calculate_aupr_consistent <- function(scores, labels) {
  scores <- as.numeric(scores)
  labels <- as.numeric(labels)
  
  if(length(unique(labels)) < 2 || all(scores == 0) || sum(labels == 1) == 0) {
    return(0)
  }
  
  tryCatch({
    pr_result <- pr.curve(
      scores.class0 = scores[labels == 1],
      scores.class1 = scores[labels == 0],
      curve = FALSE
    )
    return(pr_result$auc.integral)
  }, error = function(e) return(0))
}


# Corrected detection scoring function with proper scoring for all methods
check_detection_with_score <- function(chr, start, end, gene_tss, loops, method_name = "default") {
  enhancer_center <- (start + end) / 2
  same_chr_loops <- loops[chr1 == chr & chr2 == chr]
  if (nrow(same_chr_loops) == 0) return(0)
  
  # Find matching loops
  matches <- same_chr_loops[
    (start1 <= end & end1 >= start & start2 <= gene_tss+tolerance & end2 >= gene_tss-tolerance) |
      (start2 <= end & end2 >= start & start1 <= gene_tss+tolerance & end1 >= gene_tss-tolerance)
  ]
  
  if (nrow(matches) > 0) {
    if (method_name == "hichipper") {
      # For hichipper: use maximum PET count directly as score
      best_count <- max(matches$count, na.rm = TRUE)
      return(best_count)
    } else {
      # For all other methods: sort by qvalue first, then by PET count (descending)
      # This handles cases like cLoops/cLoops2 where q-values have low resolution
      matches_sorted <- matches[order(matches$qvalue, -matches$count)]
      best_match <- matches_sorted[1, ]
      
      best_qvalue <- best_match$qvalue
      best_count <- best_match$count
      
      if (best_qvalue == 0 || is.na(best_qvalue)) best_qvalue <- 1e-300
      
      # For methods with low q-value resolution (like cLoops), 
      # PET count becomes important for tie-breaking
      if (method_name %in% c("cLoops", "cLoops2")) {
        # For cLoops methods: q-value is primary, PET count is significant secondary
        primary_score <- -log10(best_qvalue)
        secondary_score <- log10(best_count + 1) * 0.1  # Higher weight for cLoops
        return(primary_score + secondary_score)
      } else {
        # For other methods with good q-value resolution: minimal PET influence
        primary_score <- -log10(best_qvalue)
        secondary_score <- log10(best_count + 1) * 0.001  # Small weight for tie-breaking
        return(primary_score + secondary_score)
      }
    }
  } else {
    return(0)
  }
}

# ============================================================================
# Load Data
# ============================================================================

# Ground truth
crispr_data <- read_excel("K562_CRISPR.xls")
crispr_dt <- as.data.table(crispr_data)

ground_truth <- crispr_dt[, .(
  chr = gsub("^chr", "", chr),
  start = as.numeric(start),
  end = as.numeric(end),
  gene = Gene,
  gene_tss = as.numeric(`Gene TSS`),
  element_class = if("class" %in% names(crispr_dt)) class else "enhancer",
  effect_size = as.numeric(`Fraction change in gene expr`),
  pvalue = as.numeric(`Adjusted p-value`),
  distance = abs((start + end)/2 - as.numeric(`Gene TSS`))
)]

ground_truth <- ground_truth[
  !is.na(gene_tss) & !is.na(start) & !is.na(end) & !is.na(pvalue) & !is.na(effect_size) &
    distance >= 20000 & distance <= 1500000 & element_class != "promoter"
]

ground_truth[, is_functional := (pvalue < 0.05 & abs(effect_size) > 0)]
ground_truth <- ground_truth[!is.na(is_functional)]

cat("Ground truth:", nrow(ground_truth), "enhancer-promoter pairs\n")

# Load methods (including hichipper with PET thresholds)
method_files <- c(
  "cLoops2" = "Sorted_cLoops2_qval_count_sorted.txt",
  "cLoops" = "Sorted_cLoops_qval_count_sorted.txt", 
  "HiC-DC+" = "Sorted_HiCDCPlus_qval_count_sorted.txt",
  "sintHiChIP (local)" = "Sorted_sintHiChIP.local_qval_count_sorted.txt",
  "sintHiChIP (global)" = "Sorted_sintHiChIP.global_qval_count_sorted.txt",
  "FitHiChIP" = "Sorted_FitHiChIP_qval_count_sorted.txt",
  "MAPS" = "Sorted_MAPS_qval_count_sorted.txt",
  "hichipper" = "Sorted_hichipper_qval_count_sorted.txt"
)

load_method_data <- function(file_path, method_name) {
  if (!file.exists(file_path)) return(NULL)
  
  loop_data <- fread(file_path)
  if (nrow(loop_data) == 0) return(NULL)
  
  # Standardize columns
  col_mapping <- list(
    "Chr1" = "chr1", "CHR1" = "chr1", "seqnames1" = "chr1", "V1" = "chr1",
    "Chr2" = "chr2", "CHR2" = "chr2", "seqnames2" = "chr2", "V4" = "chr2",
    "Start1" = "start1", "START1" = "start1", "start1" = "start1", "V2" = "start1",
    "Start2" = "start2", "START2" = "start2", "start2" = "start2", "V5" = "start2",
    "End1" = "end1", "END1" = "end1", "end1" = "end1", "V3" = "end1", 
    "End2" = "end2", "END2" = "end2", "end2" = "end2", "V6" = "end2",
    "ContactCount" = "count", "PETCount" = "count", "PETs" = "count", 
    "Contactcount" = "count", "counts" = "count", "V7" = "count",
    "Q-Value" = "qvalue", "FDR" = "qvalue", "qvalue" = "qvalue", 
    "Q-value" = "qvalue", "pvalue" = "qvalue", "V8" = "qvalue"
  )
  
  for (old_name in names(col_mapping)) {
    if (old_name %in% names(loop_data)) {
      setnames(loop_data, old_name, col_mapping[[old_name]])
    }
  }
  
  required_cols <- c("chr1", "chr2", "start1", "end1", "start2", "end2")
  if (!all(required_cols %in% names(loop_data))) return(NULL)
  
  if (!"count" %in% names(loop_data)) loop_data[, count := 1]
  if (!"qvalue" %in% names(loop_data)) loop_data[, qvalue := 0.01]
  
  loop_data[qvalue == 0, qvalue := 1e-300]
  loop_data[qvalue > 1, qvalue := 1]
  loop_data[qvalue < 0, qvalue := 1e-300]
  
  loop_data[, chr1 := gsub("^chr", "", chr1)]
  loop_data[, chr2 := gsub("^chr", "", chr2)]
  loop_data[, distance := abs((as.numeric(start1) + as.numeric(end1))/2 - (as.numeric(start2) + as.numeric(end2))/2)]
  
  loop_data <- loop_data[chr1 == chr2 & distance >= 20000 & distance <= 1500000]
  
  return(loop_data)
}

all_methods <- list()
for (method_name in names(method_files)) {
  result <- load_method_data(method_files[method_name], method_name)
  if (!is.null(result)) {
    all_methods[[method_name]] <- result
    cat("Loaded", method_name, ":", nrow(result), "loops\n")
  }
}

# Find valid genes
gene_summary <- ground_truth[, .(
  total_elements = .N,
  functional_elements = sum(is_functional),
  non_functional_elements = sum(!is_functional)
), by = gene]

valid_genes <- gene_summary[functional_elements > 0 & non_functional_elements > 0]$gene
cat("Valid genes for analysis:", length(valid_genes), "\n")

# ============================================================================
# Main Analysis - hichipper with percentile-based thresholds (based on PET≥2 data)
# ============================================================================

# Calculate percentile-based thresholds for hichipper based on PET≥2 filtered data
hichipper_loops <- all_methods[["hichipper"]]
hichipper_percentiles <- NULL

if(!is.null(hichipper_loops)) {
  # First filter by PET≥2 (baseline quality filter)
  hichipper_filtered <- hichipper_loops[count >= 2]
  
  # Calculate percentiles on the PET≥2 filtered data
  hichipper_percentiles <- quantile(hichipper_filtered$count, c(0.95, 0.99, 0.999), na.rm = TRUE)
  
  cat("hichipper thresholds (based on PET≥2 filtered data):\n")
  cat(sprintf("  Baseline: PET ≥ 2\n"))
  cat(sprintf("  Top 5%% (95th percentile): PET ≥ %.0f\n", hichipper_percentiles[1]))
  cat(sprintf("  Top 1%% (99th percentile): PET ≥ %.0f\n", hichipper_percentiles[2]))  
  cat(sprintf("  Top 0.1%% (99.9th percentile): PET ≥ %.0f\n", hichipper_percentiles[3]))
  
  # Define threshold mappings
  threshold_mapping <- data.frame(
    q_threshold = c(0.01, 0.001, 1e-5, 1e-7),
    pet_threshold = c(2, hichipper_percentiles[1], hichipper_percentiles[2], hichipper_percentiles[3]),
    percentile_label = c("PET≥2", "top 5%", "top 1%", "top 0.1%"),
    full_label = c("PET≥2 (baseline)", 
                   sprintf("top 5%% (PET≥%.0f)", hichipper_percentiles[1]),
                   sprintf("top 1%% (PET≥%.0f)", hichipper_percentiles[2]),
                   sprintf("top 0.1%% (PET≥%.0f)", hichipper_percentiles[3]))
  )
}
# Modified main analysis loop
# Replace the existing loop with this corrected version:

results_list <- list()
gene_aupr_list <- list()

for(method_name in names(all_methods)) {
  loops <- all_methods[[method_name]]
  if(is.null(loops)) next
  
  cat(sprintf("Evaluating %s\n", method_name))
  
  for(i in 1:length(q_thresholds)) {
    q_thresh <- q_thresholds[i]
    
    # Apply appropriate threshold based on method
    if(method_name == "hichipper" && !is.null(hichipper_percentiles)) {
      pet_thresh <- threshold_mapping$pet_threshold[i]
      thresh_loops <- loops[count >= pet_thresh]
      cat(sprintf("  %s (PET ≥ %.0f, corresponds to q < %g)\n", 
                  threshold_mapping$percentile_label[i], pet_thresh, q_thresh))
    } else if(method_name == "hichipper") {
      # Fallback if percentile calculation failed
      thresh_loops <- loops[count >= 2]
      cat(sprintf("  PET ≥ 2 (baseline fallback)\n"))
    } else {
      thresh_loops <- loops[qvalue < q_thresh]
      cat(sprintf("  q < %g\n", q_thresh))
    }
    
    gene_metrics <- list()
    
    for(target_gene in valid_genes) {
      gene_enhancers <- ground_truth[gene == target_gene]
      if(nrow(gene_enhancers) < 3) next
      
      gene_chr <- unique(gsub("^chr", "", gene_enhancers$chr))
      gene_tss <- unique(gene_enhancers$gene_tss)
      
      if(length(gene_chr) != 1 || length(gene_tss) != 1) next
      
      scores <- numeric(nrow(gene_enhancers))
      for(j in 1:nrow(gene_enhancers)) {
        row <- gene_enhancers[j]
        # IMPORTANT: Pass method_name to the function
        scores[j] <- check_detection_with_score(row$chr, row$start, row$end, 
                                                row$gene_tss, thresh_loops, method_name)
      }
      
      labels <- as.numeric(gene_enhancers$is_functional)
      
      if(sum(labels) > 0 && sum(scores > 0) > 0) {
        aupr <- calculate_aupr_consistent(scores, labels)
        
        binary_scores <- ifelse(scores > 0, 1, 0)
        TP <- sum(binary_scores == 1 & labels == 1)
        FP <- sum(binary_scores == 1 & labels == 0)
        FN <- sum(binary_scores == 0 & labels == 1)
        
        precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
        recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
        f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
        
        gene_metrics[[target_gene]] <- list(
          aupr = aupr,
          f1 = f1,
          precision = precision,
          recall = recall
        )
        
        gene_aupr_list[[paste(method_name, q_thresh, target_gene, sep = "_")]] <- list(
          method = method_name,
          q_threshold = q_thresh,
          gene = target_gene,
          aupr = aupr
        )
      }
    }
    
    if(length(gene_metrics) > 0) {
      results_list[[paste(method_name, q_thresh, sep = "_")]] <- data.frame(
        method = method_name,
        q_threshold = q_thresh,
        n_genes_aupr = length(gene_metrics),
        mean_aupr = mean(sapply(gene_metrics, function(x) x$aupr)),
        median_aupr = median(sapply(gene_metrics, function(x) x$aupr)),
        mean_f1 = mean(sapply(gene_metrics, function(x) x$f1)),
        median_f1 = median(sapply(gene_metrics, function(x) x$f1))
      )
    }
  }
}

# ============================================================================
# Follow Reference Pattern for Finding Advantage Genes
# ============================================================================
library(purrr)
# Combine results
final_results <- bind_rows(results_list)
gene_aupr_df <- bind_rows(gene_aupr_list)

# Save results
final_results <- final_results %>%
  mutate(q_threshold = formatC(q_threshold, format = "e", digits = 0))


# Calculate auPR differences for sintHiChIP (global) - EXACTLY as in reference
gene_aupr_diff <- gene_aupr_df %>%
  group_by(q_threshold, gene) %>%
  summarise(
    sintHiChIP_global_aupr = aupr[method == "sintHiChIP (global)"],
    max_other_aupr = max(aupr[method != "sintHiChIP (global)"], na.rm = TRUE),
    best_other_method = method[method != "sintHiChIP (global)"][which.max(aupr[method != "sintHiChIP (global)"])],
    aupr_diff = sintHiChIP_global_aupr - max_other_aupr,
    .groups = "drop"
  ) %>%
  filter(!is.na(sintHiChIP_global_aupr))
# Get top 5 genes based on maximum auPR difference (adjust n=5 as needed)
top_genes_max <- gene_aupr_diff %>%
  group_by(q_threshold) %>%
  arrange(desc(aupr_diff)) %>%
  slice_head(n = 17) %>%  # Changed from 17 to 5 for top 5; adjust if needed
  ungroup()

# Function to create detailed comparison table for a given threshold
create_detailed_comparison <- function(threshold) {
  # Get top 5 genes for this threshold
  top_genes_this_threshold <- top_genes_max[top_genes_max$q_threshold == threshold, ]$gene
  
  # Get auPR values for all methods on these genes
  detailed_results <- gene_aupr_df %>%
    filter(q_threshold == threshold & gene %in% top_genes_this_threshold) %>%
    select(gene, method, aupr) %>%
    pivot_wider(names_from = method, values_from = aupr) %>%
    # Reorder columns to put sintHiChIP (global) first
    select(gene, `sintHiChIP (global)`, everything()) %>%
    # Join with auPR difference for sorting
    left_join(
      top_genes_max %>%
        filter(q_threshold == threshold) %>%
        select(gene, aupr_diff),
      by = "gene"
    ) %>%
    arrange(desc(aupr_diff)) %>%
    select(-aupr_diff)  # Remove diff column from final table
  
  # Add threshold column for identification
  detailed_results$threshold <- threshold
  
  return(detailed_results)
}

# Use map_dfr to apply function over all thresholds and bind rows
all_detailed_tables <- map_dfr(q_thresholds, function(thresh) {
  create_detailed_comparison(thresh)
})

# Print summary
cat("\n=== Top 5 Genes - All Methods auPR Values Across Thresholds ===\n")
print(all_detailed_tables)

# Save combined table to CSV
write.csv(all_detailed_tables, "genes_all_methods_all_thresholds.csv", row.names = FALSE)

# Optional: Print individual thresholds for verification
for(thresh in q_thresholds) {
  cat(sprintf("\n=== Q-value < %g ===\n", thresh))
  detailed_table <- create_detailed_comparison(thresh)
  print(detailed_table)
  
  # Save individual CSV
  write.csv(detailed_table, 
            paste0("genes_all_methods_q", gsub("\\.", "", gsub("-", "", thresh)), ".csv"), 
            row.names = FALSE)
}

# Nature-format CRISPRi validation plots
# Required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)


data<-data.frame(read.csv("genes_all_methods_q001.csv",header = TRUE))
data<-data[,-ncol(data)]
names(data)<-c("gene", "sintHiChIP (global)", "cLoops2", "cLoops", "HiC-DC+", "sintHiChIP (local)", "FitHiChIP", "MAPS", "hichipper")
# # Create data matrix
# data <- tribble(
#   ~gene, ~`sintHiChIP (global)`, ~cLoops2, ~cLoops, ~`HiC-DC+`, ~`sintHiChIP (local)`, ~FitHiChIP, ~MAPS, ~hichipper,
#   "PQBP1", 0.340, 0.0256, 0.0274, 0.0245, 0.180, 0.0945, NA, 0.120,
#   "HNRNPA1", 0.213, 0.0305, NA, 0.169, 0.0379, 0.0607, 0.0276, 0.0232,
#   "JUNB", 0.0577, 0.0148, 0.0143, 0.0126, 0.0105, 0.0469, 0.0148, 0.00916,
#   "BAX", 0.333, NA, NA, 0.130, 0.117, 0.333, 0.333, 0.25,
#   "RAD23A", 0.5, NA, NA, 0.225, 0.5, 0.5, 0.5, 0.25,
#   "HDAC6", 0.0119, 0.0155, NA, 0.0107, 0.0143, 0.0109, NA, 0.0120,
#   "MYC", 0.697, 0.0530, 0.0393, 0.691, 0.688, 0.665, 0.703, 0.206,
#   "GATA1", 0.00717, NA, NA, 0.00722, 0.0151, 0.00707, NA, 0.00950,
#   "PRDX2", 0.545, NA, NA, NA, 0.242, 0.560, 0.249, 0.222,
#   "PLP2", 0.505, 0.0486, NA, 0.348, 0.351, 0.231, NA, 0.532,
#   "RPN1", 0.0138, 0.0115, 0.0131, 0.0123, 0.0242, 0.0393, 0.0122, 0.0405,
#   "H1FX", 0.0215, 0.0212, NA, 0.0226, 0.0198, 0.0636, 0.0218, 0.0211,
#   "KLF1", 0.131, NA, NA, 0.0469, 0.0642, 0.159, 0.209, 0.0508,
#   "FTL", 0.532, NA, 0.0315, 0.555, 0.645, 0.529, 0.396, 0.571,
#   "NUCB1", 0.353, NA, NA, 0.353, 0.389, 0.324, 0.0204, 0.75,
#   "FUT1", 0.107, 0.0129, NA, 0.0120, 0.0970, 0.166, 0.0128, 0.548,
#   "COPZ1", 0.0537, 0.00982, NA, 0.5, 0.00620, 0.0375, 0.00804, 0.0909
# )
# Nature journal style theme

# ============================================================================
# CRISPRi Visualization Code - Plot Generation Only
# ============================================================================

library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)

# Assumes you already have the following data objects from previous analysis:
# - bubble_plot_data: top genes data for bubble plot
# - box_data: all genes auPR data for box plot  
# - plot_data: line plot data with x_position, n_genes_aupr, median_aupr
# - method_colors: color scheme for methods

# ============================================================================
# Color Scheme and Theme Setup
# ============================================================================

# Consistent color scheme for all methods
method_colors <- c(
  "sintHiChIP (global)" = "#DC143C",  # Crimson Red  
  "sintHiChIP (local)" = "#228B22",   # Forest Green
  "HiC-DC+" = "#FFD700",             # Gold
  "MAPS" = "#4169E1",                # Royal Blue
  "hichipper" = "#8B4513",           # Saddle Brown
  "FitHiChIP" = "#9932CC",           # Dark Orchid
  "cLoops2" = "#FF8C00",             # Dark Orange
  "cLoops" = "#708090"               # Slate Gray
)

# Unified nature theme for all plots
nature_theme <- theme_minimal() +
  theme(
    # Text settings - consistent across all plots
    text = element_text(color = "black"),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    
    # Background and grid
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    
    # Borders and axes
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    
    # Consistent axis text margins for all plots
    axis.text.x = element_text(margin = margin(t = 10, unit = "pt")),
    axis.text.y = element_text(margin = margin(r = 5, unit = "pt")),
    axis.title.x = element_text(margin = margin(t = 15, unit = "pt")),
    axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
    
    # Legend defaults
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white"),
    
    # Consistent plot margins for all plots
    plot.margin = margin(50, 40, 50, 60)
  )

# ============================================================================
# Data Preparation for Plots
# ============================================================================
bubble_plot_data<-data
# Prepare bubble plot data
bubble_data <- bubble_plot_data %>%
  pivot_longer(-gene, names_to = "method", values_to = "aupr") %>%
  filter(!is.na(aupr)) %>%
  mutate(
    gene = factor(gene, levels = bubble_plot_data$gene),
    method = factor(method, levels = names(method_colors))
  )

# Calculate gene-wise scaled auPR values for bubble size
bubble_data <- bubble_data %>%
  group_by(gene) %>%
  mutate(
    max_aupr_in_gene = max(aupr, na.rm = TRUE),
    gene_scaled_aupr = aupr / max_aupr_in_gene
  ) %>%
  ungroup()

# Calculate gene ordering by sintHiChIP (global) advantage
gene_advantage <- bubble_data %>%
  filter(method != "sintHiChIP (global)") %>%
  group_by(gene) %>%
  summarise(max_other = max(aupr, na.rm = TRUE), .groups = 'drop') %>%
  left_join(
    bubble_data %>% filter(method == "sintHiChIP (global)") %>% 
      select(gene, sintHiChIP_aupr = aupr),
    by = "gene"
  ) %>%
  mutate(advantage = sintHiChIP_aupr - max_other) %>%
  arrange(desc(advantage))

# Reorder genes by advantage magnitude
bubble_data$gene <- factor(bubble_data$gene, levels = gene_advantage$gene)


# Prepare box_data for Figure B (auPR distribution boxplot, using q<0.01 threshold)
box_data <- gene_aupr_df %>%
  filter(q_threshold == 0.01) %>%  # Use the same threshold as bubble plot
  select(method, gene, aupr) %>%
  filter(!is.na(aupr) & aupr > 0) %>%  # Remove invalid values
  mutate(
    method = factor(method, levels = names(method_colors)),
    gene = factor(gene)  # For jitter positioning
  )

# Prepare plot_data for Figures C & D (line plots across thresholds)
plot_data <- final_results %>%
  mutate(
    x_position = match(q_threshold, q_thresholds),  # 1:4 mapping
    method = factor(method, levels = names(method_colors))
  ) %>%
  select(method, x_position, n_genes_aupr, median_aupr) %>%
  filter(!is.na(median_aupr))  # Remove any NA


# ============================================================================
# CRISPRi Visualization Code - Plot Generation Only
# ============================================================================



# ============================================================================
# Figure A: Bubble Plot with Dual Encoding
# ============================================================================
p_bubble <- ggplot(bubble_data, aes(x = gene, y = method)) +
  geom_point(aes(size = aupr, fill = gene_scaled_aupr), 
             shape = 21, stroke = 0.8, alpha = 0.8) +
  scale_size_continuous(
    name = "auPR",
    range = c(1, 20),
    breaks = c(0, 0.1, 0.2, 0.4, 0.6),
    labels = c("0", "0.1", "0.2", "0.4", "0.6")
  ) +
  scale_fill_gradient2(
    name = "Scaled auPR",
    low = "blue",
    mid = "yellow", 
    high = "red",
    midpoint = 0.5,
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    limits=c(0, 1)
  ) +
  labs(
    x = "Target Genes",
    y = "HiChIP Analysis Methods"
  ) +
  nature_theme +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.spacing.y = unit(5, "pt"),
    legend.box.spacing = unit(10, "pt"),
    legend.margin = margin(5, 5, 5, 5),
    plot.margin = margin(60, 60, 80, 60, unit = "pt"),
    legend.title = element_text(hjust = 0.5)  # Center legend titles
  ) +
  guides(
    size = guide_legend(
      override.aes = list(alpha = 1, fill = "white"),
      ncol = 1,
      order = 1,
      title.position = "top",
      keywidth = unit(0.6, "cm"),
      keyheight = unit(1.2, "cm")
    ),
    fill = guide_colorbar(
      barwidth = unit(1.2, "cm"),
      barheight = unit(6, "cm"),
      order = 2,
      title.position = "top"
    )
  )

# ============================================================================
# Figure B: auPR Distribution Boxplot
# ============================================================================

p_boxplot <- ggplot(box_data, aes(x = reorder(method, aupr, FUN = median, na.rm = TRUE), 
                                  y = aupr, fill = method)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5, 
               outlier.stroke = 0.8, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, stroke = 0.3, shape = 21) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = "HiChIP Analysis Methods",
    y = "auPR"
  ) +
  nature_theme +
  theme(
    legend.position = "none",  # Remove redundant legend since X-axis shows methods
    plot.margin = margin(80, 40, 50, 60)
  ) +
  stat_summary(fun = median, geom = "text", 
               aes(label = sprintf("%.3f", after_stat(y))),
               vjust = -0.5, size = 3.5, fontface = "bold")

# ============================================================================
# Figure C: Gene Detection Line Plot
# ============================================================================

p_genes <- ggplot(plot_data, aes(x = x_position, y = n_genes_aupr, 
                                 color = method, group = method)) +
  geom_line(size = 2, alpha = 0.8) +
  geom_point(size = 3) +
  geom_line(data = plot_data[plot_data$method == "hichipper",], 
            linetype = "dashed", size = 2) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("q<0.01 (PET>=2)", "q<0.001 (top 5%)", 
               "q<1e-5 (top 1%)", "q<1e-7 (top 0.1%)")
  ) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  nature_theme +
  theme(
    # aspect.ratio = 0.8,
    legend.position = "right",
    legend.key.size = unit(0.6, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 11)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) +
  labs(
    x = "PET Threshold", 
    y = "Genes with auPR",
    color = "Method"
  )

# ============================================================================
# Figure D: auPR Performance Line Plot
# ============================================================================

p_aupr <- ggplot(plot_data, aes(x = x_position, y = median_aupr, 
                                color = method, group = method)) +
  geom_line(size = 2, alpha = 0.8) +
  geom_point(size = 3) +
  geom_line(data = plot_data[plot_data$method == "hichipper",], 
            linetype = "dashed", size = 2) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("q<0.01 (PET≥>=2)", "q<0.001 (top 5%)", 
               "q<1e-5 (top 1%)", "q<1e-7 (top 0.1%)")
  ) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) +
  nature_theme +
  theme(
    # aspect.ratio = 0.8,
    legend.position = "right",
    legend.key.size = unit(0.6, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 11)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) +
  labs(
    x = "PET Threshold", 
    y = "Median auPR",
    color = "Method"
  )

# ============================================================================
# Combined Layout
# ============================================================================

# Create CD combined plot with equal heights
cd_combined <- grid.arrange(
  p_genes,
  p_aupr,
  nrow = 2,
  heights = c(1, 1),
  padding = unit(2, "cm")
)

# Create final combined plot with proper spacing
combined_plot <- grid.arrange(
  grid.rect(gp = gpar(col = "white", fill = "white")),  # Top spacing
  p_bubble,                    # A: Bubble plot
  arrangeGrob(                 # Bottom row
    p_boxplot,                 # B: Box plot
    cd_combined,               # CD: Line plots
    ncol = 2,
    widths = c(1.1, 0.9)       # Give CD section more space for legends
  ),
  grid.rect(gp = gpar(col = "white", fill = "white")),  # Bottom spacing
  nrow = 4,
  heights = c(0.2, 2.5, 2.5, 0.1)  # Balanced heights
)

# ============================================================================
# Save Plots
# ============================================================================

# Save the combined plot
ggsave("Figure4_Combined_Fixed.pdf", combined_plot, 
       width = 22, height = 20, dpi = 600, device = "pdf")

# Save individual plots
ggsave("Figure4A_Bubble.pdf", p_bubble, width = 16, height = 8, dpi = 600, device = "pdf")
ggsave("Figure4B_Boxplot.pdf", p_boxplot, width = 12, height = 6, dpi = 600, device = "pdf")
ggsave("Figure4C_Genes.pdf", p_genes, width = 10, height = 6, dpi = 600, device = "pdf")
ggsave("Figure4D_auPR.pdf", p_aupr, width = 10, height = 6, dpi = 600, device = "pdf")

# Display the combined plot
print(combined_plot)

# Print completion message
cat("Plot generation completed successfully!\n")
cat("Generated files:\n")
cat("- Figure4_Combined_Fixed.pdf (main combined figure)\n")
cat("- Figure4A_Bubble.pdf (bubble plot)\n")
cat("- Figure4B_Boxplot.pdf (boxplot)\n") 
cat("- Figure4C_Genes.pdf (gene detection)\n")
cat("- Figure4D_auPR.pdf (auPR performance)\n")

library(data.table)
library(dplyr)
library(readxl)

setwd("/home/wding/HiChIP/Callloop/Update/Result/CRISPR_K562")

# ============================================================================
# Configuration - Must match auPR analysis exactly
# ============================================================================
tolerance <- 5000  # Same as in auPR analysis
q_thresholds <- c(0.01, 0.001, 1e-5, 1e-7)  # Same thresholds

# ============================================================================
# Detection Functions - Unified with auPR analysis logic
# ============================================================================

# Main detection function - exactly matching auPR analysis


# Get all interactions for a given enhancer-promoter pair
# This returns the actual loop data for WashU track generation
get_enhancer_promoter_interactions <- function(chr, start, end, gene_tss, loops) {
  # Filter loops on the same chromosome
  same_chr_loops <- loops[chr1 == chr & chr2 == chr]
  if (nrow(same_chr_loops) == 0) return(data.table())
  
  # Use EXACT same matching logic as auPR analysis
  matches <- same_chr_loops[
    (start1 <= end & end1 >= start & start2 <= gene_tss+tolerance & end2 >= gene_tss-tolerance) |
      (start2 <= end & end2 >= start & start1 <= gene_tss+tolerance & end1 >= gene_tss-tolerance)
  ]
  
  return(matches)
}

# ============================================================================
# Load Data - Exactly matching auPR analysis
# ============================================================================

# Load CRISPR ground truth data
crispr_data <- read_excel("K562_CRISPR.xls")
crispr_dt <- as.data.table(crispr_data)

# Process ground truth exactly as in auPR analysis
ground_truth <- crispr_dt[, .(
  chr = gsub("^chr", "", chr),
  start = as.numeric(start),
  end = as.numeric(end),
  gene = Gene,
  gene_tss = as.numeric(`Gene TSS`),
  element_class = if("class" %in% names(crispr_dt)) class else "enhancer",
  effect_size = as.numeric(`Fraction change in gene expr`),
  pvalue = as.numeric(`Adjusted p-value`),
  distance = abs((start + end)/2 - as.numeric(`Gene TSS`))
)]

# Apply same filtering criteria as auPR analysis
ground_truth <- ground_truth[
  !is.na(gene_tss) & !is.na(start) & !is.na(end) & !is.na(pvalue) & !is.na(effect_size) &
    distance >= 20000 & distance <= 1500000 & element_class != "promoter"
]

ground_truth[, is_functional := (pvalue < 0.05 & abs(effect_size) > 0)]
ground_truth <- ground_truth[!is.na(is_functional)]

cat("Ground truth:", nrow(ground_truth), "enhancer-promoter pairs\n")

# Method files - same as auPR analysis
method_files <- c(
  "cLoops2" = "Sorted_cLoops2_qval_count_sorted.txt",
  "cLoops" = "Sorted_cLoops_qval_count_sorted.txt", 
  "HiC-DC+" = "Sorted_HiCDCPlus_qval_count_sorted.txt",
  "sintHiChIP (local)" = "Sorted_sintHiChIP.local_qval_count_sorted.txt",
  "sintHiChIP (global)" = "Sorted_sintHiChIP.global_qval_count_sorted.txt",
  "FitHiChIP" = "Sorted_FitHiChIP_qval_count_sorted.txt",
  "MAPS" = "Sorted_MAPS_qval_count_sorted.txt",
  "hichipper" = "Sorted_hichipper_qval_count_sorted.txt"
)

# Load method data function - exactly as in auPR analysis
load_method_data <- function(file_path, method_name) {
  if (!file.exists(file_path)) return(NULL)
  
  loop_data <- fread(file_path)
  if (nrow(loop_data) == 0) return(NULL)
  
  # Standardize column names - same mapping as auPR analysis
  col_mapping <- list(
    "Chr1" = "chr1", "CHR1" = "chr1", "seqnames1" = "chr1", "V1" = "chr1",
    "Chr2" = "chr2", "CHR2" = "chr2", "seqnames2" = "chr2", "V4" = "chr2",
    "Start1" = "start1", "START1" = "start1", "start1" = "start1", "V2" = "start1",
    "Start2" = "start2", "START2" = "start2", "start2" = "start2", "V5" = "start2",
    "End1" = "end1", "END1" = "end1", "end1" = "end1", "V3" = "end1", 
    "End2" = "end2", "END2" = "end2", "end2" = "end2", "V6" = "end2",
    "ContactCount" = "count", "PETCount" = "count", "PETs" = "count", 
    "Contactcount" = "count", "counts" = "count", "V7" = "count",
    "Q-Value" = "qvalue", "FDR" = "qvalue", "qvalue" = "qvalue", 
    "Q-value" = "qvalue", "pvalue" = "qvalue", "V8" = "qvalue"
  )
  
  for (old_name in names(col_mapping)) {
    if (old_name %in% names(loop_data)) {
      setnames(loop_data, old_name, col_mapping[[old_name]])
    }
  }
  
  required_cols <- c("chr1", "chr2", "start1", "end1", "start2", "end2")
  if (!all(required_cols %in% names(loop_data))) return(NULL)
  
  # Handle missing columns with same defaults as auPR analysis
  if (!"count" %in% names(loop_data)) loop_data[, count := 1]
  if (!"qvalue" %in% names(loop_data)) loop_data[, qvalue := 0.01]
  
  # Handle edge cases for q-values - same as auPR analysis
  loop_data[qvalue == 0, qvalue := 1e-300]
  loop_data[qvalue > 1, qvalue := 1]
  loop_data[qvalue < 0, qvalue := 1e-300]
  
  # Process chromosomes and calculate distances
  loop_data[, chr1 := gsub("^chr", "", chr1)]
  loop_data[, chr2 := gsub("^chr", "", chr2)]
  loop_data[, distance := abs((as.numeric(start1) + as.numeric(end1))/2 - 
                                (as.numeric(start2) + as.numeric(end2))/2)]
  
  # Apply same distance filtering as auPR analysis
  loop_data <- loop_data[chr1 == chr2 & distance >= 20000 & distance <= 1500000]
  
  return(loop_data)
}

# Load all methods
all_methods <- list()
for (method_name in names(method_files)) {
  result <- load_method_data(method_files[method_name], method_name)
  if (!is.null(result)) {
    all_methods[[method_name]] <- result
    cat("Loaded", method_name, ":", nrow(result), "loops\n")
  }
}

# Calculate hichipper percentiles based on PET≥2 data - exactly as in auPR analysis
hichipper_percentiles <- NULL
if(!is.null(all_methods[["hichipper"]])) {
  # First filter by PET≥2 (baseline quality filter)
  hichipper_filtered <- all_methods[["hichipper"]][count >= 2]
  
  # Calculate percentiles on the PET≥2 filtered data
  hichipper_percentiles <- quantile(hichipper_filtered$count, c(0.95, 0.99, 0.999), na.rm = TRUE)
  
  cat("\nhichipper thresholds (based on PET>=2 filtered data):\n")
  cat(sprintf("  Baseline: PET >= 2\n"))
  cat(sprintf("  Top 5%% (95th percentile): PET >= %.0f\n", hichipper_percentiles[1]))
  cat(sprintf("  Top 1%% (99th percentile): PET >= %.0f\n", hichipper_percentiles[2]))  
  cat(sprintf("  Top 0.1%% (99.9th percentile): PET >= %.0f\n", hichipper_percentiles[3]))
}

# ============================================================================
# Generate WashU Tracks for Target Genes
# ============================================================================

output_dir <- "washu_specific_genes"
dir.create(output_dir, showWarnings = FALSE)

# Target genes - confirmed as best performing genes
target_genes <- c("PQBP1", "FTL", "RAD23A")

cat("\n=== Generating WashU tracks for target genes ===\n")

for(target_gene in target_genes) {
  cat(sprintf("\nProcessing gene: %s\n", target_gene))
  
  # Get enhancers for this gene
  gene_enhancers <- ground_truth[gene == target_gene]
  if(nrow(gene_enhancers) < 3) {
    cat(sprintf("  Skipping %s: insufficient enhancers\n", target_gene))
    next
  }
  
  # Get gene information
  gene_chr <- unique(gsub("^chr", "", gene_enhancers$chr))
  gene_tss <- unique(gene_enhancers$gene_tss)
  
  if(length(gene_chr) != 1 || length(gene_tss) != 1) {
    cat(sprintf("  Skipping %s: inconsistent chromosome or TSS\n", target_gene))
    next
  }
  
  cat(sprintf("  Chr: %s, TSS: %d\n", gene_chr, gene_tss))
  cat(sprintf("  Total enhancers: %d (Functional: %d, Non-functional: %d)\n",
              nrow(gene_enhancers),
              sum(gene_enhancers$is_functional),
              sum(!gene_enhancers$is_functional)))
  
  # Process each method
  for(method_name in names(all_methods)) {
    loops <- all_methods[[method_name]]
    if(is.null(loops)) next
    
    # Collect all loops for this gene's enhancers
    all_gene_loops <- data.table()
    
    for(j in 1:nrow(gene_enhancers)) {
      enh <- gene_enhancers[j]
      
      # Get interactions using the unified detection function
      matches <- get_enhancer_promoter_interactions(
        enh$chr, enh$start, enh$end, enh$gene_tss, loops
      )
      
      if(nrow(matches) > 0) {
        # Add metadata about the enhancer
        matches[, gene := target_gene]
        matches[, enhancer_start := enh$start]
        matches[, enhancer_end := enh$end]
        matches[, enhancer_functional := enh$is_functional]
        matches[, enhancer_effect := enh$effect_size]
        matches[, enhancer_pvalue := enh$pvalue]
        all_gene_loops <- rbind(all_gene_loops, matches)
      }
    }
    
    if(nrow(all_gene_loops) == 0) {
      cat(sprintf("    %s: No interactions detected\n", method_name))
      next
    }
    
    # Get unique loops (remove duplicates)
    unique_loops <- unique(all_gene_loops, by = c("chr1", "start1", "end1", "chr2", "start2", "end2"))
    
    cat(sprintf("    %s: %d unique interactions\n", method_name, nrow(unique_loops)))
    
    # Add 'chr' prefix back for WashU format
    unique_loops$chr1 <- paste0("chr", unique_loops$chr1)
    unique_loops$chr2 <- paste0("chr", unique_loops$chr2)
    
    # Prepare interaction data with appropriate scoring
    if(method_name == "hichipper") {
      # For hichipper, use PET count as score
      interact_data <- unique_loops[, .(
        chr1, start1, end1,
        chr2, start2, end2,
        score = count,  # Use raw PET count for hichipper
        gene_info = gene
      )]
    } else {
      # For other methods, use -log10(q-value) as score
      interact_data <- unique_loops[, .(
        chr1, start1, end1,
        chr2, start2, end2,
        score = round(-log10(pmax(qvalue, 1e-300)), 2),
        gene_info = gene
      )]
    }
    
    # Clean method name for filename
    method_clean <- gsub("\\(|\\)", "", gsub("[+]", "plus", method_name))
    temp_file <- file.path(output_dir, paste0(target_gene, "_", method_clean, "_temp.txt"))
    
    # Format for WashU Browser
    interact_data_washu <- interact_data[, .(
      chr1, start1, end1,
      interaction = paste0(chr2, ":", start2, "-", end2, "|", gene_info),
      score
    )]
    
    # Write temporary file
    write.table(interact_data_washu[, .(chr1, start1, end1, interaction, score)], 
                temp_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Convert to WashU paired-end format using awk
    output_file <- file.path(output_dir, paste0(target_gene, "_", method_clean, ".washu.bed"))
    
    awk_cmd <- "awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\",\"$5\"\\t\"(NR*2-1)\"\\t.\\n\"substr($4,1,index($4,\":\")-1)\"\\t\"substr($4,index($4,\":\")+1,index($4,\"-\")-index($4,\":\")-1)\"\\t\"substr($4,index($4,\"-\")+1)\"\\t\"$1\":\"$2\"-\"$3\",\"$5\"\\t\"(NR*2)\"\\t.\"}'"
    system(paste(awk_cmd, shQuote(temp_file), "| sort -k1,1 -k2,2n >", shQuote(output_file)))
    
    # Compress and index
    if(file.exists(output_file) && file.size(output_file) > 0) {
      system(paste("bgzip -f", shQuote(output_file)))
      system(paste("tabix -p bed", shQuote(paste0(output_file, ".gz"))))
      cat(sprintf("      Created: %s.gz\n", basename(output_file)))
    }
    
    # Clean up temporary file
    unlink(temp_file)
  }
  
  # Generate enhancer tracks for this gene
  gene_enhancers[, chr := paste0("chr", chr)]
  
  # Functional enhancers track
  func_enh <- gene_enhancers[is_functional == TRUE]
  if(nrow(func_enh) > 0) {
    func_file <- file.path(output_dir, paste0(target_gene, "_functional_enhancers.bedGraph"))
    func_data <- func_enh[, .(chr, start, end, score = 1)]
    func_data <- func_data[order(chr, start)]
    write.table(func_data, func_file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created functional enhancers track: %d enhancers\n", nrow(func_enh)))
  }
  
  # Non-functional enhancers track
  nonfunc_enh <- gene_enhancers[is_functional == FALSE]
  if(nrow(nonfunc_enh) > 0) {
    nonfunc_file <- file.path(output_dir, paste0(target_gene, "_nonfunctional_enhancers.bedGraph"))
    nonfunc_data <- nonfunc_enh[, .(chr, start, end, score = 1)]
    nonfunc_data <- nonfunc_data[order(chr, start)]
    write.table(nonfunc_data, nonfunc_file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created non-functional enhancers track: %d enhancers\n", nrow(nonfunc_enh)))
  }
  
  # Gene TSS location track
  tss_file <- file.path(output_dir, paste0(target_gene, "_tss_location.bed"))
  tss_data <- data.table(
    chr = paste0("chr", gene_chr),
    start = gene_tss - 1000,
    end = gene_tss + 1000,
    name = target_gene,
    score = 1000
  )
  write.table(tss_data, tss_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Created TSS location track\n"))
}

# ============================================================================
# Generate comprehensive README with threshold information
# ============================================================================

readme_file <- file.path(output_dir, "README.txt")
readme_content <- c(
  "========================================================================",
  "WashU Browser Track Files for CRISPR-validated Enhancer-Promoter Interactions",
  "========================================================================",
  "",
  paste("Target genes:", paste(target_genes, collapse = ", ")),
  "These genes showed the highest auPR advantage for sintHiChIP (global)",
  "",
  "========================================================================",
  "FILE DESCRIPTIONS",
  "========================================================================",
  "",
  "For each gene, the following files are generated:",
  "",
  "1. Interaction Tracks:",
  "   {GENE}_{METHOD}.washu.bed.gz - All detected enhancer-promoter interactions",
  "   Score represents interaction confidence:",
  "   - For q-value methods: score = -log10(q-value)",
  "   - For hichipper: score = PET count",
  "",
  "2. Enhancer Annotation Tracks:",
  "   {GENE}_functional_enhancers.bedGraph - CRISPR-validated functional enhancers",
  "   {GENE}_nonfunctional_enhancers.bedGraph - Non-functional enhancers",
  "",
  "3. Gene Location:",
  "   {GENE}_tss_location.bed - Transcription start site +/-1kb",
  "",
  "========================================================================",
  "FILTERING THRESHOLDS IN WASHU BROWSER",
  "========================================================================",
  "",
  "To reproduce different stringency levels from the auPR analysis,",
  "apply the following score filters in WashU Browser:",
  "",
  "For Q-VALUE BASED METHODS (all except hichipper):",
  "  score >= 2  ->  q-value < 0.01",
  "  score >= 3  ->  q-value < 0.001",
  "  score >= 5  ->  q-value < 1e-5",
  "  score >= 7  ->  q-value < 1e-7",
  "",
  "For HICHIPPER (PET count based):",
  "  score >= 2  ->  Baseline (PET >= 2)"
)

# Add hichipper percentile thresholds if calculated
if(!is.null(hichipper_percentiles)) {
  readme_content <- c(
    readme_content,
    paste0("  score >= ", round(hichipper_percentiles[1]), "  ->  Top 5% (95th percentile)"),
    paste0("  score >= ", round(hichipper_percentiles[2]), "  ->  Top 1% (99th percentile)"),
    paste0("  score >= ", round(hichipper_percentiles[3]), "  ->  Top 0.1% (99.9th percentile)")
  )
}

readme_content <- c(
  readme_content,
  "",
  "========================================================================",
  "VISUALIZATION TIPS",
  "========================================================================",
  "",
  "1. Load functional and non-functional enhancer tracks with different colors:",
  "   - Functional: Green or Red (active enhancers)",
  "   - Non-functional: Gray (inactive/neutral enhancers)",
  "",
  "2. Compare methods by loading multiple interaction tracks:",
  "   - sintHiChIP (global) should show enrichment for functional enhancers",
  "   - Other methods may show more interactions with non-functional enhancers",
  "",
  "3. Use score filtering to see how methods perform at different thresholds:",
  "   - Higher thresholds should retain more functional interactions",
  "   - sintHiChIP (global) maintains better precision at stringent thresholds",
  "",
  "========================================================================",
  "ANALYSIS PARAMETERS",
  "========================================================================",
  "",
  "- Distance range: 20kb - 1.5Mb",
  "- Tolerance for overlap: ±5kb",
  "- Detection logic: Consistent with auPR analysis",
  "- Ground truth: CRISPR perturbation data from K562 cells",
  "",
  "Generated on:", format(Sys.Date(), "%B %d, %Y")
)

writeLines(readme_content, readme_file)

cat("\n=== WashU track generation completed ===\n")
cat(sprintf("Output directory: %s\n", output_dir))
