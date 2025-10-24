#!/usr/bin/env Rscript

# ==============================================================================
# Figure 5: Complete eQTL Analysis - Refactored Version
# Separated logic for complete figure vs individual panels for efficiency
# ==============================================================================

# Load required libraries
suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(data.table)
  library(tidyr)
  library(reshape2)
  library(RColorBrewer)
  library(pheatmap)
  library(tibble)
  library(grid)
  library(cowplot)
  library(scales)
  library(stringr)
})

# Global parameters
MIN_DISTANCE <- 20000
MAX_DISTANCE <- 1500000
EQTL_QVAL_THRESHOLD <- 0.05
MIN_TSS_DISTANCE <- 20000
LOOP_ANCHOR_WINDOW <- 10000
GENE_TOLERANCE <- 0
OUTPUT_DIR <- "eqtl_results"
QVALUE_THRES <- 0.01

# Configuration for 8 methods
method_files <- c(
  "HiC-DC+" = "Sorted_HiCDCPlus_qval_count_sorted_hg38.txt",
  "sintHiChIP (local)" = "Sorted_sintHiChIP.local_qval_count_sorted_hg38.txt",
  "sintHiChIP (global)" = "Sorted_sintHiChIP.global_qval_count_sorted_hg38.txt",
  "FitHiChIP" = "Sorted_FitHiChIP_qval_count_sorted_hg38.txt",
  "MAPS" = "Sorted_MAPS_qval_count_sorted_hg38.txt",
  "hichipper" = "Sorted_hichipper_qval_count_sorted_hg38.txt",
  "cLoops" = "Sorted_cLoops_qval_count_sorted_hg38.txt",
  "cLoops2" = "Sorted_cLoops2_qval_count_sorted_hg38.txt"
)

# Color scheme for 8 methods (Nature-compatible palette)
method_colors <- c(
  "HiC-DC+" = "#FFD700",           # Gold
  "sintHiChIP (local)" = "#228B22", # Forest Green
  "sintHiChIP (global)" = "#DC143C", # Crimson Red
  "FitHiChIP" = "#9932CC",         # Dark Orchid
  "MAPS" = "#4169E1",              # Royal Blue
  "hichipper" = "#8B4513",         # Saddle Brown
  "cLoops" = "#708090",            # Slate Gray
  "cLoops2" = "#FF8C00"            # Dark Orange
)

EQTL_FILE <- "Cells_EBV-transformed_lymphocytes.v10.eGenes.txt"
setwd("/home/wding/HiChIP/Callloop/Update/Result/Motif/GM12878_output/")

# Fix for dplyr select conflict
select <- dplyr::select 

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ==============================================================================
# Execution control parameters
# ==============================================================================

# Set what to generate
GENERATE_COMPLETE_FIGURE <- TRUE    # Generate complete combined figure
GENERATE_INDIVIDUAL_PANELS <- FALSE # Generate individual panel files
SAVE_DATA <- TRUE                   # Save data files

# Select specific panels to generate (only when GENERATE_INDIVIDUAL_PANELS=TRUE)
PANELS_TO_GENERATE <- c("A", "B", "C", "D", "E", "F")  # Or subset like c("A", "F")

# ==============================================================================
# Core data processing functions
# ==============================================================================

read_loop_data <- function(file_path, method_name, max_loops = Inf) {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  # Auto-detect header
  first_line <- readLines(file_path, n = 1)
  has_header <- grepl("chr|seqnames|start|end|count|qvalue|PETCount", first_line, ignore.case = TRUE)
  
  if(has_header) {
    loops <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Special handling for hichipper format
    if (tolower(method_name) == "hichipper") {
      if("PETCount" %in% colnames(loops)) {
        colnames(loops)[1:8] <- c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "count", "qvalue")
        loops$qvalue <- 0.001  # Set dummy qvalue for hichipper
      }
    } else {
      if(ncol(loops) >= 8) {
        colnames(loops)[1:8] <- c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "count", "qvalue")
      }
    }
  } else {
    loops <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(loops) <- c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "count", "qvalue")
  }
  
  # Clean chromosome naming
  loops$seqnames1 <- gsub("chr+", "", loops$seqnames1)
  loops$seqnames1 <- paste0("chr", loops$seqnames1)
  loops$seqnames2 <- gsub("chr+", "", loops$seqnames2)
  loops$seqnames2 <- paste0("chr", loops$seqnames2)
  
  # Convert to numeric
  numeric_cols <- c("start1", "end1", "start2", "end2", "count", "qvalue")
  for(col in numeric_cols) {
    loops[[col]] <- suppressWarnings(as.numeric(as.character(loops[[col]])))
  }
  
  # Clean and filter data
  loops <- loops[complete.cases(loops[numeric_cols]), ]
  loops$distance <- abs(loops$start2 - loops$start1)
  
  # Apply distance and chromosome filters
  if (tolower(method_name) == "hichipper") {
    loops_filtered <- loops[loops$distance >= MIN_DISTANCE & 
                              loops$distance <= MAX_DISTANCE, ]
    loops_filtered <- loops_filtered[loops_filtered$count >= 2, ]
    loops_sorted <- loops_filtered[order(-loops_filtered$count), ]
  } else {
    loops_filtered <- loops[loops$distance >= MIN_DISTANCE & 
                              loops$distance <= MAX_DISTANCE & 
                              loops$qvalue <= QVALUE_THRES, ]
    loops_filtered <- loops_filtered[loops_filtered$seqnames1 == loops_filtered$seqnames2, ]
    loops_sorted <- loops_filtered[order(loops_filtered$qvalue, -loops_filtered$count), ]
  }
  
  if(is.infinite(max_loops)) {
    top_loops <- loops_sorted
  } else {
    top_loops <- head(loops_sorted, min(max_loops, nrow(loops_sorted)))
  }
  
  return(top_loops)
}

load_eqtl_data <- function(eqtl_file) {
  eqtl_data <- fread(eqtl_file, sep = "\t", header = TRUE)
  
  # Parse variant information
  variant_parts <- strsplit(eqtl_data$variant_id, "_")
  eqtl_data$variant_chr <- sapply(variant_parts, function(x) x[1])
  eqtl_data$variant_pos <- as.numeric(sapply(variant_parts, function(x) x[2]))
  eqtl_data$variant_ref <- sapply(variant_parts, function(x) x[3])
  eqtl_data$variant_alt <- sapply(variant_parts, function(x) x[4])
  
  # Clean chromosome naming
  eqtl_data$variant_chr <- gsub("^chr", "", eqtl_data$variant_chr)
  eqtl_data$variant_chr <- paste0("chr", eqtl_data$variant_chr)
  
  if (!any(grepl("^chr", eqtl_data$gene_chr))) {
    eqtl_data$gene_chr <- paste0("chr", eqtl_data$gene_chr)
  }
  
  # Apply filters
  valid_records <- !is.na(eqtl_data$variant_pos) & !is.na(eqtl_data$slope) & 
    eqtl_data$variant_chr %in% paste0("chr", c(1:22, "X", "Y")) &
    eqtl_data$gene_chr %in% paste0("chr", c(1:22, "X", "Y"))
  
  eqtl_data_filtered <- eqtl_data[valid_records, ]
  eqtl_significant <- eqtl_data_filtered[eqtl_data_filtered$qval < EQTL_QVAL_THRESHOLD, ]
  
  # Filter by gene biotypes if available
  if("biotype" %in% colnames(eqtl_significant) && nrow(eqtl_significant) > 0) {
    coding_biotypes <- c("protein_coding", "lincRNA", "lncRNA", "antisense", 
                         "processed_transcript", "sense_intronic")
    eqtl_significant <- eqtl_significant[eqtl_significant$biotype %in% coding_biotypes, ]
  }
  
  return(as.data.frame(eqtl_significant))
}

create_eqtl_ranges <- function(eqtl_data) {
  variant_gr <- GRanges(
    seqnames = eqtl_data$variant_chr,
    ranges = IRanges(start = eqtl_data$variant_pos, end = eqtl_data$variant_pos),
    strand = "*",
    variant_id = eqtl_data$variant_id,
    gene_id = eqtl_data$gene_id,
    gene_name = eqtl_data$gene_name,
    gene_chr = eqtl_data$gene_chr,
    gene_start = eqtl_data$gene_start,
    gene_end = eqtl_data$gene_end,
    tss_distance = eqtl_data$tss_distance,
    pval_nominal = if("pval_nominal" %in% colnames(eqtl_data)) eqtl_data$pval_nominal else NA,
    qval = eqtl_data$qval,
    slope = eqtl_data$slope,
    slope_se = if("slope_se" %in% colnames(eqtl_data)) eqtl_data$slope_se else NA,
    af = if("af" %in% colnames(eqtl_data)) eqtl_data$af else NA
  )
  
  unique_genes <- unique(eqtl_data[, c("gene_id", "gene_name", "gene_chr", 
                                       "gene_start", "gene_end", "strand")])
  if("biotype" %in% colnames(eqtl_data)) {
    unique_genes$biotype <- eqtl_data$biotype[match(unique_genes$gene_id, eqtl_data$gene_id)]
  }
  
  gene_gr <- GRanges(
    seqnames = unique_genes$gene_chr,
    ranges = IRanges(start = unique_genes$gene_start, end = unique_genes$gene_end),
    strand = unique_genes$strand,
    gene_id = unique_genes$gene_id,
    gene_name = unique_genes$gene_name,
    biotype = if("biotype" %in% colnames(unique_genes)) unique_genes$biotype else NA
  )
  
  return(list(variants = variant_gr, genes = gene_gr))
}

calculate_loop_eqtl_overlap <- function(loops_data, eqtl_variants, method_name) {
  anchor1_regions <- GRanges(
    seqnames = loops_data$seqnames1,
    ranges = IRanges(
      start = pmax(1, (loops_data$start1 + loops_data$end1)/2 - LOOP_ANCHOR_WINDOW),
      end = (loops_data$start1 + loops_data$end1)/2 + LOOP_ANCHOR_WINDOW
    ),
    loop_id = 1:nrow(loops_data),
    distance = loops_data$distance,
    loop_count = loops_data$count,
    loop_qvalue = loops_data$qvalue
  )
  
  anchor2_regions <- GRanges(
    seqnames = loops_data$seqnames2,
    ranges = IRanges(
      start = pmax(1, (loops_data$start2 + loops_data$end2)/2 - LOOP_ANCHOR_WINDOW),
      end = (loops_data$start2 + loops_data$end2)/2 + LOOP_ANCHOR_WINDOW
    ),
    loop_id = 1:nrow(loops_data),
    distance = loops_data$distance,
    loop_count = loops_data$count,
    loop_qvalue = loops_data$qvalue
  )
  
  # Find overlaps
  nearest_anchor1 <- distanceToNearest(eqtl_variants, anchor1_regions)
  anchor1_overlaps_idx <- which(mcols(nearest_anchor1)$distance == 0)
  
  nearest_anchor2 <- distanceToNearest(eqtl_variants, anchor2_regions)
  anchor2_overlaps_idx <- which(mcols(nearest_anchor2)$distance == 0)
  
  # Extract overlapping records
  eqtl_associations <- data.frame()
  
  if(length(anchor1_overlaps_idx) > 0) {
    eqtl_idx <- queryHits(nearest_anchor1)[anchor1_overlaps_idx]
    loop_idx <- subjectHits(nearest_anchor1)[anchor1_overlaps_idx]
    
    anchor1_df <- data.frame(
      loop_id = anchor1_regions$loop_id[loop_idx],
      variant_anchor = "anchor1",
      variant_id = eqtl_variants$variant_id[eqtl_idx],
      gene_id = eqtl_variants$gene_id[eqtl_idx],
      gene_name = eqtl_variants$gene_name[eqtl_idx],
      gene_chr = eqtl_variants$gene_chr[eqtl_idx],
      gene_start = eqtl_variants$gene_start[eqtl_idx],
      gene_end = eqtl_variants$gene_end[eqtl_idx],
      tss_distance = eqtl_variants$tss_distance[eqtl_idx],
      loop_distance = anchor1_regions$distance[loop_idx],
      pval_nominal = eqtl_variants$pval_nominal[eqtl_idx],
      qval = eqtl_variants$qval[eqtl_idx],
      slope = eqtl_variants$slope[eqtl_idx],
      af = eqtl_variants$af[eqtl_idx],
      loop_count = anchor1_regions$loop_count[loop_idx],
      loop_qvalue = anchor1_regions$loop_qvalue[loop_idx],
      stringsAsFactors = FALSE
    )
    eqtl_associations <- rbind(eqtl_associations, anchor1_df)
  }
  
  if(length(anchor2_overlaps_idx) > 0) {
    eqtl_idx <- queryHits(nearest_anchor2)[anchor2_overlaps_idx]
    loop_idx <- subjectHits(nearest_anchor2)[anchor2_overlaps_idx]
    
    anchor2_df <- data.frame(
      loop_id = anchor2_regions$loop_id[loop_idx],
      variant_anchor = "anchor2",
      variant_id = eqtl_variants$variant_id[eqtl_idx],
      gene_id = eqtl_variants$gene_id[eqtl_idx],
      gene_name = eqtl_variants$gene_name[eqtl_idx],
      gene_chr = eqtl_variants$gene_chr[eqtl_idx],
      gene_start = eqtl_variants$gene_start[eqtl_idx],
      gene_end = eqtl_variants$gene_end[eqtl_idx],
      tss_distance = eqtl_variants$tss_distance[eqtl_idx],
      loop_distance = anchor2_regions$distance[loop_idx],
      pval_nominal = eqtl_variants$pval_nominal[eqtl_idx],
      qval = eqtl_variants$qval[eqtl_idx],
      slope = eqtl_variants$slope[eqtl_idx],
      af = eqtl_variants$af[eqtl_idx],
      loop_count = anchor2_regions$loop_count[loop_idx],
      loop_qvalue = anchor2_regions$loop_qvalue[loop_idx],
      stringsAsFactors = FALSE
    )
    eqtl_associations <- rbind(eqtl_associations, anchor2_df)
  }
  
  # Remove duplicates
  eqtl_associations <- eqtl_associations[
    !duplicated(paste(eqtl_associations$loop_id, eqtl_associations$variant_id, eqtl_associations$gene_id)), 
  ]
  
  return(eqtl_associations)
}

identify_regulatory_loops <- function(loops_data, eqtl_associations, eqtl_genes, method_name) {
  regulatory_loops <- data.frame()
  
  for(i in 1:nrow(eqtl_associations)) {
    assoc <- eqtl_associations[i, ]
    loop_data <- loops_data[assoc$loop_id, ]
    
    if(assoc$variant_anchor == "anchor1") {
      other_anchor_chr <- loop_data$seqnames2
      other_anchor_start <- loop_data$start2
      other_anchor_end <- loop_data$end2
    } else {
      other_anchor_chr <- loop_data$seqnames1
      other_anchor_start <- loop_data$start1
      other_anchor_end <- loop_data$end1
    }
    
    if(other_anchor_chr == assoc$gene_chr) {
      gene_overlap <- (other_anchor_start <= assoc$gene_end + GENE_TOLERANCE) & 
        (other_anchor_end >= assoc$gene_start - GENE_TOLERANCE)
      
      if(gene_overlap) {
        regulatory_loops <- rbind(regulatory_loops, data.frame(
          loop_id = assoc$loop_id,
          variant_id = assoc$variant_id,
          variant_anchor = assoc$variant_anchor,
          gene_id = assoc$gene_id,
          gene_name = assoc$gene_name,
          gene_anchor = ifelse(assoc$variant_anchor == "anchor1", "anchor2", "anchor1"),
          tss_distance = assoc$tss_distance,
          loop_distance = assoc$loop_distance,
          eqtl_pval = assoc$pval_nominal,
          eqtl_qval = assoc$qval,
          eqtl_slope = assoc$slope,
          eqtl_af = assoc$af,
          loop_count = assoc$loop_count,
          loop_qvalue = assoc$loop_qvalue,
          method = method_name,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(regulatory_loops)
}

compute_eqtl_metrics <- function(regulatory_loops) {
  if(missing(regulatory_loops) || is.null(regulatory_loops) || nrow(regulatory_loops) == 0) {
    stop("Error: regulatory_loops is missing, null, or empty")
  }
  
  required_cols <- c("method", "eqtl_slope", "eqtl_qval", "eqtl_af", "tss_distance", 
                     "gene_id", "variant_id")
  missing_cols <- setdiff(required_cols, colnames(regulatory_loops))
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  metrics_summary <- regulatory_loops %>%
    group_by(method) %>%
    summarise(
      # Discovery capacity
      total_regulatory_loops = n(),
      unique_genes = n_distinct(gene_id, na.rm = TRUE),
      unique_variants = n_distinct(variant_id, na.rm = TRUE),
      
      # Effect size metrics
      median_effect_size = median(abs(eqtl_slope), na.rm = TRUE),
      mean_effect_size = mean(abs(eqtl_slope), na.rm = TRUE),
      strong_effects = sum(abs(eqtl_slope) > 0.4, na.rm = TRUE),
      strong_effect_rate = strong_effects / n() * 100,
      very_strong_effects = sum(abs(eqtl_slope) > 0.6, na.rm = TRUE),
      very_strong_rate = very_strong_effects / n() * 100,
      
      # Statistical significance metrics
      median_qval = median(eqtl_qval, na.rm = TRUE),
      high_significance = sum(eqtl_qval < 0.001, na.rm = TRUE),
      high_significance_rate = high_significance / n() * 100,
      very_high_significance = sum(eqtl_qval < 1e-5, na.rm = TRUE),
      very_high_significance_rate = very_high_significance / n() * 100,
      
      # Distance metrics
      proximal_eqtls = sum(abs(tss_distance) < 100000, na.rm = TRUE),
      proximal_eqtl_rate = proximal_eqtls / n() * 100,
      median_tss_distance = median(abs(tss_distance), na.rm = TRUE),
      
      # Allele frequency metrics
      median_af = median(eqtl_af, na.rm = TRUE),
      rare_eqtl_rate = sum(eqtl_af < 0.05, na.rm = TRUE) / n() * 100,
      common_eqtl_rate = sum(eqtl_af > 0.2, na.rm = TRUE) / n() * 100,
      
      .groups = 'drop'
    )
  
  return(metrics_summary)
}

# ==============================================================================
# Publication theme
# ==============================================================================

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = 48, face = "bold", hjust = 0, color = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold", color = "black"),
      axis.text = element_text(size = base_size - 1, face = "bold", color = "black"),
      axis.line = element_line(color = "black", size = 1),
      axis.ticks = element_line(color = "black", size = 1),
      
      # Adjust spacing for x-axis text (distance from axis line)
      axis.text.x = element_text(
        size = base_size - 1, 
        face = "bold", 
        color = "black",
        margin = margin(t = 10, r = 0, b = 0, l = 0)  # t = top margin (distance from axis)
      ),
      
      # Adjust spacing for x-axis title (distance from axis text)
      axis.title.x = element_text(
        size = base_size + 2, 
        face = "bold", 
        color = "black",
        margin = margin(t = 15, r = 0, b = 0, l = 0)  # t = top margin (distance from axis.text)
      ),
      
      # Same for y-axis if needed
      axis.text.y = element_text(
        size = base_size - 1, 
        face = "bold", 
        color = "black",
        margin = margin(t = 0, r = 10, b = 0, l = 0)  # r = right margin
      ),
      
      axis.title.y = element_text(
        size = base_size + 2, 
        face = "bold", 
        color = "black",
        margin = margin(t = 0, r = 15, b = 0, l = 0)  # r = right margin
      ),
      
      legend.title = element_text(size = base_size + 1, face = "bold", color = "black"),
      legend.text = element_text(size = base_size - 1, face = "bold", color = "black"),
      legend.key.size = unit(0.5, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.margin = margin(t = 20, r = 20, b = 40, l = 20, unit = "pt")
    )
}
# ==============================================================================
# Panel generation functions with two-line method labels
# ==============================================================================

# Helper function to convert method names to two-line format
convert_method_labels <- function(method_names) {
  sapply(method_names, function(x) {
    if(grepl("sintHiChIP.*local", x, ignore.case = TRUE)) {
      return("sintHiChIP\n(local)")
    } else if(grepl("sintHiChIP.*global", x, ignore.case = TRUE)) {
      return("sintHiChIP\n(global)")
    }
    return(x)
  })
}

# ==============================================================================
# Panel generation functions - renamed to match final figure labels
# ==============================================================================

# Helper function to convert method names to two-line format
convert_method_labels <- function(method_names) {
  sapply(method_names, function(x) {
    if(grepl("sintHiChIP.*local", x, ignore.case = TRUE)) {
      return("sintHiChIP\n(local)")
    } else if(grepl("sintHiChIP.*global", x, ignore.case = TRUE)) {
      return("sintHiChIP\n(global)")
    }
    return(x)
  })
}

# Panel A: Progressive eQTL Loop Discovery
create_panel_a_progressive_discovery <- function(loop_data, eqtl_ranges, output_dir) {
  # Calculate minimum loops across all methods for fair comparison
  min_loops <- min(sapply(loop_data, nrow))
  step_size <- max(1000, floor(min_loops / 20))
  loop_counts <- seq(step_size, min_loops, by = step_size)
  if(max(loop_counts) < min_loops) {
    loop_counts <- c(loop_counts, min_loops)
  }
  
  progressive_results <- list()
  
  # Calculate progressive discovery metrics for each step
  for(i in seq_along(loop_counts)) {
    current_loop_count <- loop_counts[i]
    step_results <- data.frame()
    
    for(method in names(loop_data)) {
      loops <- head(loop_data[[method]], current_loop_count)
      eqtl_overlaps <- calculate_loop_eqtl_overlap(loops, eqtl_ranges$variants, method)
      regulatory_loops <- identify_regulatory_loops(loops, eqtl_overlaps, eqtl_ranges$genes, method)
      
      if(nrow(regulatory_loops) > 0) {
        step_metrics <- regulatory_loops %>%
          summarise(
            method = method,
            loop_count = current_loop_count,
            total_regulatory_loops = n(),
            unique_genes = n_distinct(gene_id),
            .groups = 'drop'
          )
      } else {
        step_metrics <- data.frame(
          method = method,
          loop_count = current_loop_count,
          total_regulatory_loops = 0,
          unique_genes = 0
        )
      }
      step_results <- rbind(step_results, step_metrics)
    }
    progressive_results[[i]] <- step_results
  }
  
  # Combine all progressive results
  all_progressive_results <- do.call(rbind, progressive_results)
  method_order <- names(method_colors)
  all_progressive_results$method <- factor(all_progressive_results$method, levels = method_order)
  
  # Create the plot
  p_a <- ggplot(all_progressive_results, aes(x = loop_count, y = total_regulatory_loops, 
                                             color = method)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1.5) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(labels = scales::comma_format()) +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(
      x = "Number of loops",
      y = "Total regulatory loops discovered",
      color = "Method"
    ) +
    theme_publication(base_size = 12) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm")
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  return(p_a)
}

# Panel B: High Significance Rate (originally Panel C)
create_panel_b_significance <- function(all_regulatory_loops, output_dir) {
  method_order <- names(method_colors)
  
  # Calculate high significance rates for each method
  summary_data <- all_regulatory_loops %>%
    mutate(
      method = factor(method, levels = method_order),
      high_significance = eqtl_qval < 1e-5
    ) %>%
    filter(method %in% method_order, !is.na(method)) %>%
    group_by(method) %>%
    summarise(
      high_significance_pct = mean(high_significance, na.rm = TRUE) * 100,
      .groups = 'drop'
    )
  
  # Convert method labels to two-line format
  summary_data$method_label <- convert_method_labels(as.character(summary_data$method))
  summary_data$method_label <- factor(summary_data$method_label, 
                                      levels = convert_method_labels(method_order))
  
  # Reorder by high significance percentage (descending)
  summary_data <- summary_data[order(-summary_data$high_significance_pct), ]
  summary_data$method_label <- reorder(summary_data$method_label, -summary_data$high_significance_pct)
  
  # Create the plot
  p_b <- ggplot(summary_data, aes(x = method_label, y = high_significance_pct, fill = method)) +
    geom_col(alpha = 0.8, width = 0.7) +
    scale_fill_manual(values = method_colors) +
    labs(
      y = "High significance rate (%)",
      x = "Method"
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 7.5, face = "bold", color = "black",
                                 margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 15, l = 5)
    ) +
    geom_text(aes(label = sprintf("%.1f", high_significance_pct)), 
              vjust = -0.5, size = 3, fontface = "bold", color = "black") +
    ylim(0, max(summary_data$high_significance_pct, na.rm = TRUE) * 1.15)
  
  return(p_b)
}

# Panel C: Strong Effect Rate (originally Panel B)
create_panel_c_strong_effects <- function(all_regulatory_loops, output_dir) {
  method_order <- names(method_colors)
  
  # Calculate strong effect rates for each method
  summary_data <- all_regulatory_loops %>%
    mutate(
      method = factor(method, levels = method_order),
      strong_effect = abs(eqtl_slope) > 0.4
    ) %>%
    filter(method %in% method_order, !is.na(method)) %>%
    group_by(method) %>%
    summarise(
      strong_effect_pct = mean(strong_effect, na.rm = TRUE) * 100,
      .groups = 'drop'
    )
  
  # Convert method labels to two-line format
  summary_data$method_label <- convert_method_labels(as.character(summary_data$method))
  summary_data$method_label <- factor(summary_data$method_label, 
                                      levels = convert_method_labels(method_order))
  
  # Reorder by strong effect percentage (descending)
  summary_data <- summary_data[order(-summary_data$strong_effect_pct), ]
  summary_data$method_label <- reorder(summary_data$method_label, -summary_data$strong_effect_pct)
  
  # Create the plot
  p_c <- ggplot(summary_data, aes(x = method_label, y = strong_effect_pct, fill = method)) +
    geom_col(alpha = 0.8, width = 0.7) +
    scale_fill_manual(values = method_colors) +
    labs(
      y = "Strong effect rate (%)",
      x = "Method"
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 7.5, face = "bold", color = "black",
                                 margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 15, l = 5)
    ) +
    geom_text(aes(label = sprintf("%.1f", strong_effect_pct)), 
              vjust = -0.5, size = 3, fontface = "bold", color = "black") +
    ylim(0, max(summary_data$strong_effect_pct, na.rm = TRUE) * 1.15)
  
  return(p_c)
}

# Panel D: Effect Size Distribution
create_panel_d_effect_size <- function(all_regulatory_loops, output_dir) {
  method_order <- names(method_colors)
  
  # Prepare data for effect size distribution analysis
  plot_data <- all_regulatory_loops %>%
    mutate(
      method = factor(method, levels = method_order),
      abs_effect_size = abs(eqtl_slope)
    ) %>%
    filter(method %in% method_order, !is.na(method))
  
  # Convert method labels to two-line format
  plot_data$method_label <- convert_method_labels(as.character(plot_data$method))
  plot_data$method_label <- factor(plot_data$method_label, 
                                   levels = convert_method_labels(method_order))
  
  # Reorder by median effect size (descending)
  median_order <- plot_data %>%
    group_by(method_label) %>%
    summarise(median_effect = median(abs_effect_size, na.rm = TRUE), .groups = 'drop') %>%
    arrange(desc(median_effect))
  
  plot_data$method_label <- factor(plot_data$method_label, levels = median_order$method_label)
  
  # Create the plot
  p_d <- ggplot(plot_data, aes(x = method_label, y = abs_effect_size, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = method_colors) +
    coord_cartesian(ylim = quantile(plot_data$abs_effect_size, c(0.05, 0.95), na.rm = TRUE)) +
    labs(
      y = "Absolute effect size",
      x = "Method"
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 7.5, face = "bold", color = "black",
                                 margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 15, l = 5)
    ) +
    stat_summary(fun = median, geom = "text", 
                 aes(label = sprintf("%.3f", after_stat(y))), 
                 vjust = -0.5, size = 3, fontface = "bold", color = "black")
  
  return(p_d)
}

# Panel E: TSS Distance Distribution
create_panel_e_tss_distance <- function(all_regulatory_loops, output_dir) {
  method_order <- names(method_colors)
  
  # Prepare data for TSS distance analysis
  plot_data <- all_regulatory_loops %>%
    mutate(
      method = factor(method, levels = method_order),
      abs_tss_distance = abs(tss_distance)
    ) %>%
    filter(method %in% method_order, !is.na(method), abs_tss_distance > 1000)
  
  # Convert method labels to two-line format
  plot_data$method_label <- convert_method_labels(as.character(plot_data$method))
  plot_data$method_label <- factor(plot_data$method_label, 
                                   levels = convert_method_labels(method_order))
  
  # Reorder by median TSS distance (ascending - closer to TSS is better)
  median_order <- plot_data %>%
    group_by(method_label) %>%
    summarise(median_distance = median(abs_tss_distance, na.rm = TRUE), .groups = 'drop') %>%
    arrange(median_distance)
  
  plot_data$method_label <- factor(plot_data$method_label, levels = median_order$method_label)
  
  # Create the plot
  p_e <- ggplot(plot_data, aes(x = method_label, y = abs_tss_distance/1000, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 0.5) +
    scale_fill_manual(values = method_colors) +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    labs(
      y = "Distance to TSS (kb, log scale)",
      x = "Method"
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 7.5, face = "bold", color = "black",
                                 margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
      legend.position = "none",
      plot.margin = margin(t = 5, r = 5, b = 15, l = 5)
    ) +
    stat_summary(fun = median, geom = "text", 
                 aes(label = sprintf("%.0f", after_stat(y))), 
                 vjust = -0.5, size = 3, fontface = "bold", color = "black")
  
  return(p_e)
}

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = 48, face = "bold", hjust = 0, color = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold", color = "black"),
      axis.text = element_text(size = base_size - 1, face = "bold", color = "black"),
      axis.line = element_line(color = "black", size = 1),
      axis.ticks = element_line(color = "black", size = 1),
      
      # Adjust spacing for x-axis text (distance from axis line)
      axis.text.x = element_text(
        size = base_size - 1, 
        face = "bold", 
        color = "black",
        margin = margin(t = 10, r = 0, b = 0, l = 0)  # t = top margin (distance from axis)
      ),
      
      # Adjust spacing for x-axis title (distance from axis text)
      axis.title.x = element_text(
        size = base_size + 2, 
        face = "bold", 
        color = "black",
        margin = margin(t = 15, r = 0, b = 0, l = 0)  # t = top margin (distance from axis.text)
      ),
      
      # Same for y-axis if needed
      axis.text.y = element_text(
        size = base_size - 1, 
        face = "bold", 
        color = "black",
        margin = margin(t = 0, r = 10, b = 0, l = 0)  # r = right margin
      ),
      
      axis.title.y = element_text(
        size = base_size + 2, 
        face = "bold", 
        color = "black",
        margin = margin(t = 0, r = 15, b = 0, l = 0)  # r = right margin
      ),
      
      legend.title = element_text(size = base_size + 1, face = "bold", color = "black"),
      legend.text = element_text(size = base_size - 1, face = "bold", color = "black"),
      legend.key.size = unit(0.5, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.margin = margin(t = 20, r = 20, b = 40, l = 20, unit = "pt")
    )
}
# Function to create panel f pathway analysis
create_panel_f_pathway_analysis <- function(output_dir) {
  # Define pathway order - group shared pathways at top
  pathway_order <- c(
    "Intracellular signal transduction",
    "Fc receptor signaling",
    "Circulating IL-6 level increased",
    "IgM level decreased",
    "Marginal zone B cell number decreased",
    "Humoral immune response decreased",
    "Transitional B cell number decreased",
    "IL-2 production increased",
    "Sialylation"
  )
  
  # Create pathway enrichment data (updated to match image and latest KEGG data)
  pathway_data <- data.frame(
    Method = c(rep("sintHiChIP\n(local)", 8),
               rep("sintHiChIP\n(global)", 3),
               rep("FitHiChIP", 2),
               "hichipper",
               rep("HiC-DC+", 0), rep("MAPS", 0), rep("cLoops", 0), rep("cLoops2", 0)), # Only methods with data
    Pathway = factor(c(
      "Intracellular signal transduction", "Marginal zone B cell number decreased",
      "Humoral immune response decreased", "Circulating IL-6 level increased",
      "Transitional B cell number decreased", "IgM level decreased",
      "IL-2 production increased", "Sialylation",
      "Intracellular signal transduction", "Fc receptor signaling", "IgM level decreased",
      "Intracellular signal transduction", "Circulating IL-6 level increased",
      "Intracellular signal transduction"
    ), levels = rev(pathway_order)),
    q_value = c(0.00003652, 0.02810, 0.02810, 0.02810, 0.02810, 0.02810, 0.08097, 0.08097,
                0.01141, 0.01390, 0.08181,
                0.0000968, 0.01417,
                0.03351),
    Category = c("BCR signaling", "B cell phenotype", "B cell phenotype", "Cytokine response",
                 "B cell phenotype", "B cell phenotype", "Cytokine response", "Glycosylation",
                 "BCR signaling", "Immune activation", "B cell phenotype",
                 "BCR signaling", "Cytokine response",
                 "BCR signaling"),
    stringsAsFactors = FALSE
  )
  
  # Define all methods (including empty ones for grid)
  all_methods <- c("sintHiChIP\n(local)", "sintHiChIP\n(global)", "FitHiChIP", "hichipper",
                   "HiC-DC+", "MAPS", "cLoops", "cLoops2")
  
  # Create full matrix
  full_matrix <- expand.grid(
    Method = factor(all_methods, levels = all_methods),
    Pathway = factor(rev(pathway_order), levels = rev(pathway_order)),
    stringsAsFactors = FALSE
  )
  
  # Join with pathway data and calculate significance levels
  full_data <- full_matrix %>%
    left_join(pathway_data, by = c("Method", "Pathway")) %>%
    mutate(
      neg_log_q = -log10(ifelse(is.na(q_value), 1, q_value)),  # Avoid NA in log
      sig_level = case_when(
        q_value < 0.001 ~ "***",
        q_value < 0.01 ~ "**",
        q_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # Filter data with actual q-values for plotting
  point_data <- full_data %>% filter(!is.na(q_value))
  
  # Create the plot
  p_f <- ggplot(full_data, aes(x = Method, y = Pathway)) +
    annotate("rect", xmin = -Inf, xmax = Inf, 
             ymin = 0.5, ymax = 5.5,
             fill = "#E8F4F8", alpha = 0.3) +
    geom_tile(color = "grey95", fill = NA, size = 0.3) +
    geom_point(data = point_data,
               aes(size = neg_log_q, color = Category)) +
    geom_text(data = point_data,
              aes(label = sig_level), 
              size = 10, vjust = 1.0, hjust = 0.5, fontface = "bold", color = "black") +
    scale_size_continuous(
      range = c(3, 9),
      breaks = -log10(c(0.1, 0.05, 0.01, 0.001)),
      labels = c("0.1", "0.05", "0.01", "0.001"),
      name = "q-value"
    ) +
    scale_color_manual(
      values = c(
        "BCR signaling" = "#D62728",
        "B cell phenotype" = "#1F77B4",
        "Cytokine response" = "#2CA02C",
        "Immune activation" = "#FF7F0E",
        "Glycosylation" = "#9467BD"
      ),
      name = "Pathway category",
      breaks = c("BCR signaling", "B cell phenotype", "Cytokine response", 
                 "Immune activation", "Glycosylation"),
      na.translate = FALSE
    ) +
    scale_x_discrete(name = "Method") +
    scale_y_discrete(
      name = "GM12878-specific pathways (q < 0.1)",
      labels = function(x) stringr::str_wrap(x, width = 22)
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 8, face = "bold", color = "black",
                                 margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(size = 9, face = "bold", color = "black"),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11, color = "black"),
      legend.text = element_text(size = 10, color = "black"),
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(t = 5, r = 5, b = 15, l = 5)
    ) 
  # Save as PDF
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggsave(filename = file.path(output_dir, "Figure5_f.pdf"), plot = p_f, 
         width = 10, height = 6, device = "pdf")
  
  return(p_f)
}

# ==============================================================================
# Main figure generation function
# ==============================================================================

generate_figure5 <- function(prepared_data, output_directory, 
                             mode = "both",
                             panels_to_generate = c("a", "b", "c", "d", "e", "f"),
                             label_size = 48,
                             individual_width = 8.5,
                             individual_height = 7,
                             combined_width = 18,
                             combined_height = 18) {
  
  # Extract data from prepared dataset
  loop_data <- prepared_data$loop_data
  eqtl_ranges <- prepared_data$eqtl_ranges
  all_regulatory_loops <- prepared_data$all_regulatory_loops
  
  # Create individual panel plots with corrected function names
  panel_plots <- list()
  panel_plots[["A"]] <- create_panel_a_progressive_discovery(loop_data, eqtl_ranges, output_directory)
  panel_plots[["B"]] <- create_panel_b_significance(all_regulatory_loops, output_directory)
  panel_plots[["C"]] <- create_panel_c_strong_effects(all_regulatory_loops, output_directory)
  panel_plots[["D"]] <- create_panel_d_effect_size(all_regulatory_loops, output_directory)
  panel_plots[["E"]] <- create_panel_e_tss_distance(all_regulatory_loops, output_directory)
  panel_plots[["F"]] <- create_panel_f_pathway_analysis(output_directory)
  
  results <- list()
  
  # Generate individual panels if requested
  if (mode %in% c("individual", "both")) {
    cat("Generating individual panel PDFs...\n")
    
    individual_panels <- list()
    
    for(panel in panels_to_generate) {
      if (!panel %in% names(panel_plots)) {
        warning(paste("Panel", panel, "not found. Skipping."))
        next
      }
      
      # Apply final theme without labels for individual panels
      panel_final <- panel_plots[[panel]] + 
        theme(plot.margin = margin(20, 20, 40, 20))
      
      # Save individual panel
      output_file <- file.path(output_directory, sprintf("Figure5%s_%s.pdf", panel, 
                                                         switch(panel,
                                                                "A" = "progressive_discovery",
                                                                "B" = "significance",
                                                                "C" = "strong_effects", 
                                                                "D" = "effect_size",
                                                                "E" = "tss_distance",
                                                                "F" = "pathway_analysis"
                                                         )))
      
      ggsave(output_file, panel_final, 
             width = individual_width, 
             height = individual_height, 
             dpi = 600, 
             device = "pdf")
      
      individual_panels[[panel]] <- panel_final
      cat(paste("Saved:", output_file, "\n"))
    }
    
    results$individual_panels <- individual_panels
  }
  
  # Generate combined figure if requested
  if (mode %in% c("combined", "both")) {
    cat("Generating combined figure PDF...\n")
    
    # Apply final themes with panel labels
    panel_a_final <- panel_plots[["A"]] + 
      theme(
        plot.margin = margin(5, 5, 5, 5),
        plot.title = element_text(size = label_size, face = "bold", hjust = 0, color = "black")
      ) +
      labs(title = "a")
    
    panel_b_final <- panel_plots[["B"]] +
      theme(
        plot.margin = margin(5, 5, 5, 5),
        plot.title = element_text(size = label_size, face = "bold", hjust = 0, color = "black"),
        axis.text.x = element_text(size = 7.5, face = "bold", color = "black",
                                   margin = margin(t = 5, r = 0, b = 0, l = 0))
      ) +
      labs(title = "b")
    
    panel_c_final <- panel_plots[["C"]] +
      theme(
        plot.margin = margin(5, 5, 5, 5),
        plot.title = element_text(size = label_size, face = "bold", hjust = 0, color = "black"),
        axis.text.x = element_text(size = 7.5, face = "bold", color = "black",
                                   margin = margin(t = 5, r = 0, b = 0, l = 0))
      ) +
      labs(title = "c")
    
    panel_d_final <- panel_plots[["D"]] +
      theme(
        plot.margin = margin(5, 5, 5, 5),
        plot.title = element_text(size = label_size, face = "bold", hjust = 0, color = "black"),
        axis.text.x = element_text(size = 7.5, face = "bold", color = "black",
                                   margin = margin(t = 5, r = 0, b = 0, l = 0))
      ) +
      labs(title = "d")
    
    panel_e_final <- panel_plots[["E"]] + 
      theme(
        plot.margin = margin(5, 5, 5, 5),
        plot.title = element_text(size = label_size, face = "bold", hjust = 0, color = "black"),
        axis.text.x = element_text(size = 7.5, face = "bold", color = "black",
                                   margin = margin(t = 5, r = 0, b = 0, l = 0))
      ) +
      labs(title = "e")
    
    panel_f_final <- panel_plots[["F"]] + 
      theme(
        plot.margin = margin(5, 5, 5, 5),
        plot.title = element_text(size = label_size, face = "bold", hjust = 0, color = "black")
      ) +
      labs(title = "f")
    
    # Create inner plot grid
    inner_plot <- plot_grid(
      panel_a_final, panel_b_final,
      panel_c_final, panel_d_final,
      panel_e_final, panel_f_final,
      nrow = 3, ncol = 2,
      rel_widths = c(1, 1),
      rel_heights = c(0.8, 0.8, 0.8),
      align = "hv",
      axis = "tblr"
    )
    
    # Add outer margins
    complete_figure <- plot_grid(
      inner_plot,
      nrow = 1, ncol = 1
    ) + 
      theme(plot.margin = margin(t = 60, r = 40, b = 60, l = 20, unit = "pt"))
    
    # Save combined figure
    output_file <- file.path(output_directory, "Figure5_Complete_eQTL_Analysis.pdf")
    ggsave(
      filename = output_file,
      plot = complete_figure,
      width = combined_width,
      height = combined_height,
      dpi = 600,
      device = "pdf",
      units = "in",
      bg = "white"
    )
    
    output_file <- file.path(output_directory, "Figure5_Complete_eQTL_Analysis.png")
    ggsave(
      filename = output_file,
      plot = complete_figure,
      width = combined_width,
      height = combined_height,
      dpi = 600,
      device = "png",
      units = "in",
      bg = "white"
    )
    
    
    results$complete_figure <- complete_figure
    results$combined_panels <- list(A = panel_a_final, B = panel_b_final, C = panel_c_final, 
                                    D = panel_d_final, E = panel_e_final, F = panel_f_final)
    
    cat(paste("Saved combined figure:", output_file, "\n"))
  }
  
  return(results)
}


# ==============================================================================
# Statistical significance testing for proportions
# ==============================================================================

perform_comprehensive_statistical_tests <- function(all_regulatory_loops, output_dir) {
  method_order <- names(method_colors)
  
  # ============================================================================
  # 1. Strong Effect Rate - Proportion Tests
  # ============================================================================
  strong_effect_data <- all_regulatory_loops %>%
    mutate(
      method = factor(method, levels = method_order),
      strong_effect = abs(eqtl_slope) > 0.4
    ) %>%
    filter(method %in% method_order, !is.na(method)) %>%
    group_by(method) %>%
    summarise(
      n_strong = sum(strong_effect, na.rm = TRUE),
      n_total = n(),
      prop_strong = n_strong / n_total,
      .groups = 'drop'
    )
  
  cat("\n=== Strong Effect Rate Statistical Tests ===\n")
  strong_test_results <- data.frame()
  
  for(i in 1:(nrow(strong_effect_data)-1)) {
    for(j in (i+1):nrow(strong_effect_data)) {
      method1 <- strong_effect_data$method[i]
      method2 <- strong_effect_data$method[j]
      
      # Create contingency table
      contingency <- matrix(c(
        strong_effect_data$n_strong[i], strong_effect_data$n_total[i] - strong_effect_data$n_strong[i],
        strong_effect_data$n_strong[j], strong_effect_data$n_total[j] - strong_effect_data$n_strong[j]
      ), nrow = 2, byrow = TRUE)
      
      # Chi-square test
      test_result <- chisq.test(contingency, correct = FALSE)
      
      # Two-proportion z-test
      prop_test <- prop.test(
        c(strong_effect_data$n_strong[i], strong_effect_data$n_strong[j]),
        c(strong_effect_data$n_total[i], strong_effect_data$n_total[j]),
        correct = FALSE
      )
      
      cat(sprintf("%s vs %s:\n", method1, method2))
      cat(sprintf("  Proportions: %.3f vs %.3f\n", 
                  strong_effect_data$prop_strong[i], 
                  strong_effect_data$prop_strong[j]))
      cat(sprintf("  Chi-square p-value: %.4f\n", test_result$p.value))
      cat(sprintf("  Significant at α=0.05: %s\n\n", 
                  ifelse(test_result$p.value < 0.05, "YES", "NO")))
      
      strong_test_results <- rbind(strong_test_results, data.frame(
        comparison = paste(method1, "vs", method2),
        prop1 = strong_effect_data$prop_strong[i],
        prop2 = strong_effect_data$prop_strong[j],
        chi_sq_pval = test_result$p.value,
        prop_test_pval = prop_test$p.value,
        significant_0.05 = test_result$p.value < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # ============================================================================
  # 2. High Significance Rate - Proportion Tests
  # ============================================================================
  high_sig_data <- all_regulatory_loops %>%
    mutate(
      method = factor(method, levels = method_order),
      high_sig = eqtl_qval < 1e-5
    ) %>%
    filter(method %in% method_order, !is.na(method)) %>%
    group_by(method) %>%
    summarise(
      n_high_sig = sum(high_sig, na.rm = TRUE),
      n_total = n(),
      prop_high_sig = n_high_sig / n_total,
      .groups = 'drop'
    )
  
  cat("\n=== High Significance Rate Statistical Tests ===\n")
  sig_test_results <- data.frame()
  
  for(i in 1:(nrow(high_sig_data)-1)) {
    for(j in (i+1):nrow(high_sig_data)) {
      method1 <- high_sig_data$method[i]
      method2 <- high_sig_data$method[j]
      
      contingency <- matrix(c(
        high_sig_data$n_high_sig[i], high_sig_data$n_total[i] - high_sig_data$n_high_sig[i],
        high_sig_data$n_high_sig[j], high_sig_data$n_total[j] - high_sig_data$n_high_sig[j]
      ), nrow = 2, byrow = TRUE)
      
      test_result <- chisq.test(contingency, correct = FALSE)
      
      prop_test <- prop.test(
        c(high_sig_data$n_high_sig[i], high_sig_data$n_high_sig[j]),
        c(high_sig_data$n_total[i], high_sig_data$n_total[j]),
        correct = FALSE
      )
      
      cat(sprintf("%s vs %s:\n", method1, method2))
      cat(sprintf("  Proportions: %.3f vs %.3f\n", 
                  high_sig_data$prop_high_sig[i], 
                  high_sig_data$prop_high_sig[j]))
      cat(sprintf("  Chi-square p-value: %.4f\n", test_result$p.value))
      cat(sprintf("  Significant at α=0.05: %s\n\n", 
                  ifelse(test_result$p.value < 0.05, "YES", "NO")))
      
      sig_test_results <- rbind(sig_test_results, data.frame(
        comparison = paste(method1, "vs", method2),
        prop1 = high_sig_data$prop_high_sig[i],
        prop2 = high_sig_data$prop_high_sig[j],
        chi_sq_pval = test_result$p.value,
        prop_test_pval = prop_test$p.value,
        significant_0.05 = test_result$p.value < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # ============================================================================
  # 3. Effect Size Distribution - Wilcoxon Rank-Sum Tests
  # ============================================================================
  cat("\n=== Effect Size Distribution Statistical Tests (Wilcoxon) ===\n")
  effect_size_test_results <- data.frame()
  
  methods <- unique(all_regulatory_loops$method[all_regulatory_loops$method %in% method_order])
  
  for(i in 1:(length(methods)-1)) {
    for(j in (i+1):length(methods)) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      # Extract effect sizes for each method
      effects1 <- all_regulatory_loops %>%
        filter(method == method1, !is.na(eqtl_slope)) %>%
        pull(eqtl_slope) %>%
        abs()
      
      effects2 <- all_regulatory_loops %>%
        filter(method == method2, !is.na(eqtl_slope)) %>%
        pull(eqtl_slope) %>%
        abs()
      
      # Wilcoxon rank-sum test (Mann-Whitney U test)
      wilcox_result <- wilcox.test(effects1, effects2, exact = FALSE)
      
      # Calculate medians
      median1 <- median(effects1, na.rm = TRUE)
      median2 <- median(effects2, na.rm = TRUE)
      
      cat(sprintf("%s vs %s:\n", method1, method2))
      cat(sprintf("  Medians: %.3f vs %.3f\n", median1, median2))
      cat(sprintf("  Wilcoxon p-value: %.4f\n", wilcox_result$p.value))
      cat(sprintf("  Significant at α=0.05: %s\n\n", 
                  ifelse(wilcox_result$p.value < 0.05, "YES", "NO")))
      
      effect_size_test_results <- rbind(effect_size_test_results, data.frame(
        comparison = paste(method1, "vs", method2),
        median1 = median1,
        median2 = median2,
        wilcox_pval = wilcox_result$p.value,
        significant_0.05 = wilcox_result$p.value < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # ============================================================================
  # 4. TSS Distance Distribution - Wilcoxon Rank-Sum Tests
  # ============================================================================
  cat("\n=== TSS Distance Distribution Statistical Tests (Wilcoxon) ===\n")
  tss_test_results <- data.frame()
  
  for(i in 1:(length(methods)-1)) {
    for(j in (i+1):length(methods)) {
      method1 <- methods[i]
      method2 <- methods[j]
      
      # Extract TSS distances for each method
      tss1 <- all_regulatory_loops %>%
        filter(method == method1, !is.na(tss_distance), abs(tss_distance) > -1) %>%
        pull(tss_distance) %>%
        abs()
      
      tss2 <- all_regulatory_loops %>%
        filter(method == method2, !is.na(tss_distance), abs(tss_distance) > -1) %>%
        pull(tss_distance) %>%
        abs()
      
      # Wilcoxon rank-sum test
      wilcox_result <- wilcox.test(tss1, tss2, exact = FALSE)
      
      # Calculate medians (in kb)
      median1_kb <- median(tss1, na.rm = TRUE) / 1000
      median2_kb <- median(tss2, na.rm = TRUE) / 1000
      
      cat(sprintf("%s vs %s:\n", method1, method2))
      cat(sprintf("  Medians (kb): %.1f vs %.1f\n", median1_kb, median2_kb))
      cat(sprintf("  Wilcoxon p-value: %.4f\n", wilcox_result$p.value))
      cat(sprintf("  Significant at α=0.05: %s\n\n", 
                  ifelse(wilcox_result$p.value < 0.05, "YES", "NO")))
      
      tss_test_results <- rbind(tss_test_results, data.frame(
        comparison = paste(method1, "vs", method2),
        median1_kb = median1_kb,
        median2_kb = median2_kb,
        wilcox_pval = wilcox_result$p.value,
        significant_0.05 = wilcox_result$p.value < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # ============================================================================
  # Save all results
  # ============================================================================
  write.csv(strong_test_results, 
            file.path(output_dir, "Figure5_strong_effect_statistical_tests.csv"), 
            row.names = FALSE)
  write.csv(sig_test_results, 
            file.path(output_dir, "Figure5_high_sig_statistical_tests.csv"), 
            row.names = FALSE)
  write.csv(effect_size_test_results, 
            file.path(output_dir, "Figure5_effect_size_statistical_tests.csv"), 
            row.names = FALSE)
  write.csv(tss_test_results, 
            file.path(output_dir, "Figure5_tss_distance_statistical_tests.csv"), 
            row.names = FALSE)
  
  # Summary of significant findings
  cat("\n=== SUMMARY OF SIGNIFICANT FINDINGS (α=0.05) ===\n")
  cat(sprintf("Strong effect rate: %d of %d comparisons significant\n", 
              sum(strong_test_results$significant_0.05), 
              nrow(strong_test_results)))
  cat(sprintf("High significance rate: %d of %d comparisons significant\n", 
              sum(sig_test_results$significant_0.05), 
              nrow(sig_test_results)))
  cat(sprintf("Effect size distribution: %d of %d comparisons significant\n", 
              sum(effect_size_test_results$significant_0.05), 
              nrow(effect_size_test_results)))
  cat(sprintf("TSS distance distribution: %d of %d comparisons significant\n", 
              sum(tss_test_results$significant_0.05), 
              nrow(tss_test_results)))
  
  return(list(
    strong_effects = strong_test_results,
    high_significance = sig_test_results,
    effect_size = effect_size_test_results,
    tss_distance = tss_test_results
  ))
}

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = 48, face = "bold", hjust = 0, color = "black"),
      axis.title = element_text(size = base_size + 2, face = "bold", color = "black"),
      axis.text = element_text(size = base_size - 1, face = "bold", color = "black"),
      axis.line = element_line(color = "black", size = 1),
      axis.ticks = element_line(color = "black", size = 1),
      
      # Adjust spacing for x-axis text (distance from axis line)
      axis.text.x = element_text(
        size = base_size - 1, 
        face = "bold", 
        color = "black",
        margin = margin(t = 10, r = 0, b = 0, l = 0)  # t = top margin (distance from axis)
      ),
      
      # Adjust spacing for x-axis title (distance from axis text)
      axis.title.x = element_text(
        size = base_size + 2, 
        face = "bold", 
        color = "black",
        margin = margin(t = 15, r = 0, b = 0, l = 0)  # t = top margin (distance from axis.text)
      ),
      
      # Same for y-axis if needed
      axis.text.y = element_text(
        size = base_size - 1, 
        face = "bold", 
        color = "black",
        margin = margin(t = 0, r = 10, b = 0, l = 0)  # r = right margin
      ),
      
      axis.title.y = element_text(
        size = base_size + 2, 
        face = "bold", 
        color = "black",
        margin = margin(t = 0, r = 15, b = 0, l = 0)  # r = right margin
      ),
      
      legend.title = element_text(size = base_size + 1, face = "bold", color = "black"),
      legend.text = element_text(size = base_size - 1, face = "bold", color = "black"),
      legend.key.size = unit(0.5, "cm"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.margin = margin(t = 20, r = 20, b = 40, l = 20, unit = "pt")
    )
}

# ==============================================================================
# Function to prepare all data for Figure 5 generation
# ==============================================================================
prepare_all_data <- function(method_files, eqtl_file) {
  # Initialize list to store loop data for each method
  loop_data <- list()
  
  # Load and process loop data for each method
  for (method in names(method_files)) {
    file_path <- method_files[method]
    cat(sprintf("Loading loop data for %s from %s\n", method, file_path))
    
    # Read loop data using the provided read_loop_data function
    loops <- read_loop_data(file_path, method)
    
    # Check if data was successfully loaded
    if (is.null(loops) || nrow(loops) == 0) {
      warning(sprintf("No valid loop data found for %s in %s", method, file_path))
      next
    }
    
    loop_data[[method]] <- loops
  }
  
  # Check if any loop data was loaded successfully
  if (length(loop_data) == 0) {
    stop("No valid loop data loaded for any method")
  }
  
  # Load eQTL data
  cat(sprintf("Loading eQTL data from %s\n", eqtl_file))
  eqtl_data <- load_eqtl_data(eqtl_file)
  
  # Check if eQTL data was loaded successfully
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    stop(sprintf("No valid eQTL data found in %s", eqtl_file))
  }
  
  # Create GenomicRanges objects for eQTL variants and genes
  eqtl_ranges <- create_eqtl_ranges(eqtl_data)
  
  # Initialize data frame to store all regulatory loops
  all_regulatory_loops <- data.frame()
  
  # Process each method to identify regulatory loops
  for (method in names(loop_data)) {
    cat(sprintf("Processing regulatory loops for %s\n", method))
    
    # Calculate overlaps between loops and eQTL variants
    eqtl_overlaps <- calculate_loop_eqtl_overlap(loop_data[[method]], 
                                                 eqtl_ranges$variants, 
                                                 method)
    
    # Identify regulatory loops where the other anchor overlaps with a gene
    regulatory_loops <- identify_regulatory_loops(loop_data[[method]], 
                                                  eqtl_overlaps, 
                                                  eqtl_ranges$genes, 
                                                  method)
    
    # Append to the combined regulatory loops data frame
    if (nrow(regulatory_loops) > 0) {
      all_regulatory_loops <- rbind(all_regulatory_loops, regulatory_loops)
    } else {
      warning(sprintf("No regulatory loops identified for %s", method))
    }
  }
  
  # Check if any regulatory loops were identified
  if (nrow(all_regulatory_loops) == 0) {
    warning("No regulatory loops identified across all methods")
  }
  
  # Return the prepared data as a list
  return(list(
    loop_data = loop_data,
    eqtl_ranges = eqtl_ranges,
    all_regulatory_loops = all_regulatory_loops
  ))
}


save_data_files <- function(prepared_data, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save loop data for each method
  for (method in names(prepared_data$loop_data)) {
    write.csv(prepared_data$loop_data[[method]], 
              file.path(output_dir, paste0("loop_data_", method, ".csv")), 
              row.names = FALSE)
    cat(sprintf("Saved loop data for %s\n", method))
  }
  
  # Save eQTL variants and genes
  write.csv(as.data.frame(prepared_data$eqtl_ranges$variants), 
            file.path(output_dir, "eqtl_variants.csv"), 
            row.names = FALSE)
  write.csv(as.data.frame(prepared_data$eqtl_ranges$genes), 
            file.path(output_dir, "eqtl_genes.csv"), 
            row.names = FALSE)
  
  # Save all regulatory loops
  write.csv(prepared_data$all_regulatory_loops, 
            file.path(output_dir, "all_regulatory_loops.csv"), 
            row.names = FALSE)
  
  cat("All data files saved to", output_dir, "\n")
}

# Main execution
prepared_data <- prepare_all_data(method_files, EQTL_FILE)
save_data_files(prepared_data, OUTPUT_DIR)

# NEW: Perform statistical tests
stat_tests <- perform_comprehensive_statistical_tests(prepared_data$all_regulatory_loops, OUTPUT_DIR)

results <- generate_figure5(prepared_data, OUTPUT_DIR, mode = "both")
