#!/usr/bin/env Rscript
# Clean HiChIP Recovery Script
# Generates: Overlap Matrix, Recovery Curves, Total Loops Barplot

suppressMessages({
    library(ggplot2)
    library(scales)
    library(dplyr)
    library(optparse)
    library(reshape2)
    library(cowplot)
    library(data.table)
})

# Command line options
option_list <- list(
    make_option(c("--input_dir"), type = "character", default = "comparison_output", 
                help = "Input directory containing recovery data files"),
    make_option(c("--output_dir"), type = "character", default = "publication_plots", 
                help = "Output directory for plots"),
    make_option(c("--file_pattern"), type = "character", default = "*_HiCCUPS.txt", 
                help = "File pattern to match recovery data files"),
    make_option(c("--ref_label"), type = "character", default = "HiCCUPS", 
                help = "Reference dataset label"),
    make_option(c("--reference_loop"), type = "character", default = NULL, 
                help = "Specific gold standard file name (e.g., K562_HiC_Gold.txt)"),
    make_option(c("--plot_title"), type = "character", 
                default = "HiChIP Loop Recovery", 
                help = "Main title for plots"),
    make_option(c("--output_name"), type = "character", default = "Recovery", 
                help = "Output file name prefix"),
    make_option(c("--plot_width"), type = "double", default = 12, 
                help = "Plot width in inches"),
    make_option(c("--plot_height"), type = "double", default = 8, 
                help = "Plot height in inches"),
    make_option(c("--save_png"), action = "store_true", default = FALSE, 
                help = "Save PNG versions"),
    make_option(c("--line_size"), type = "double", default = 1.5, 
                help = "Line thickness")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Color scheme for methods
get_method_colors <- function(method_names) {
    colors <- character(length(method_names))
    names(colors) <- method_names
    
    for (method in method_names) {
        if (grepl("sintHiChIP.*local", method, ignore.case = TRUE)) {
            colors[method] <- "#228B22"  # Forest Green
        } else if (grepl("sintHiChIP.*global", method, ignore.case = TRUE)) {
            colors[method] <- "#DC143C"  # Crimson Red
        } else if (grepl("HiC-DC\\+", method, ignore.case = TRUE)) {
            colors[method] <- "#FFD700"  # Gold
        } else if (grepl("MAPS", method, ignore.case = TRUE)) {
            colors[method] <- "#4169E1"  # Royal Blue
        } else if (grepl("hichipper", method, ignore.case = TRUE)) {
            colors[method] <- "#8B4513"  # Saddle Brown
        } else if (grepl("FitHiChIP", method, ignore.case = TRUE)) {
            colors[method] <- "#9932CC"  # Dark Orchid
        } else if (grepl("cLoops2", method, ignore.case = TRUE)) {
            colors[method] <- "#FF8C00"  # Dark Orange
        } else if (grepl("cLoops", method, ignore.case = TRUE)) {
            colors[method] <- "#708090"  # Slate Gray
        } else {
            colors[method] <- "#000000"  # Black
        }
    }
    return(colors)
}

# Extract and standardize method names
extract_method_name <- function(filename, ref_label = "HiCCUPS") {
    method <- gsub("Recovery_plot_data_", "", filename)
    method <- gsub(paste0("_", ref_label, ".txt"), "", method)
    method <- gsub("HiCDCPlus|HiCDC|HiC_DC_Plus", "HiC-DC+", method, ignore.case = TRUE)
    method <- gsub("sintHiChIP\\.local", "sintHiChIP (local)", method, ignore.case = TRUE)
    method <- gsub("sintHiChIP\\.global", "sintHiChIP (global)", method, ignore.case = TRUE)
    return(method)
}

# Format large numbers
format_numbers <- function(x) {
    ifelse(x >= 1000000, paste0(round(x/1000000, 1), "M"),
           ifelse(x >= 1000, paste0(round(x/1000, 1), "K"), as.character(x)))
}

# Find gold standard file
count_gold_loops <- function(input_dir, reference_loop = NULL) {
    if (!is.null(reference_loop)) {
        possible_paths <- c(
            file.path("Gold", reference_loop),
            file.path(input_dir, reference_loop),
            file.path(input_dir, "Gold", reference_loop),
            reference_loop
        )
        
        for (path in possible_paths) {
            if (file.exists(path)) {
                cat(sprintf("Using specified Gold standard file: %s\n", path))
                gold_data <- read.table(path, header = TRUE, sep = "\t")
                cat(sprintf("Successfully read %d gold standard loops\n", nrow(gold_data)))
                return(nrow(gold_data))
            }
        }
        
        stop(sprintf("Specified gold standard file not found: %s", reference_loop))
    }
    
    possible_patterns <- c("Gold.txt", "*Gold.txt", "*_Gold.txt")
    possible_dirs <- c(input_dir, "Gold", file.path(input_dir, "Gold"))
    
    for (dir in possible_dirs) {
        if (!dir.exists(dir)) next
        
        for (pattern in possible_patterns) {
            files <- list.files(path = dir, pattern = glob2rx(pattern), full.names = TRUE)
            
            if (length(files) > 0) {
                if (length(files) > 1) {
                    priority_files <- files[grepl("K562.*HiC", basename(files), ignore.case = TRUE)]
                    if (length(priority_files) > 0) {
                        files <- priority_files[1]
                    } else {
                        files <- files[1]
                    }
                }
                
                file_path <- files[1]
                cat(sprintf("Auto-detected Gold standard file: %s\n", file_path))
                gold_data <- read.table(file_path, header = TRUE, sep = "\t")
                cat(sprintf("Successfully read %d gold standard loops\n", nrow(gold_data)))
                return(nrow(gold_data))
            }
        }
    }
    
    cat("No Gold standard file found. Available files in Gold directory:\n")
    if (dir.exists("Gold")) {
        gold_files <- list.files("Gold", full.names = FALSE)
        for (file in gold_files) {
            cat(sprintf("  - %s\n", file))
        }
    }
    
    stop("Could not find Gold standard file. Use --reference_loop to specify the exact filename")
}

# Read original loop prediction files for spatial overlap calculation
read_loop_files <- function(input_dir) {
    loop_files <- list.files(path = input_dir, pattern = "Sorted_.*_qval_count_sorted\\.txt$", full.names = TRUE)
    
    if (length(loop_files) == 0) {
        cat("Warning: No loop prediction files found. Using recovery-based overlap estimation.\n")
        return(NULL)
    }
    
    # Column mapping for different file formats
    col_mapping <- list(
        "Chr1" = "chr1", "CHR1" = "chr1", "seqnames1" = "chr1", "V1" = "chr1",
        "Chr2" = "chr2", "CHR2" = "chr2", "seqnames2" = "chr2", "V4" = "chr2",
        "Start1" = "start1", "START1" = "start1", "start1" = "start1", "V2" = "start1",
        "Start2" = "start2", "START2" = "start2", "start2" = "start2", "V5" = "start2",
        "End1" = "end1", "END1" = "end1", "end1" = "end1", "V3" = "end1", 
        "End2" = "end2", "END2" = "end2", "end2" = "end2", "V6" = "end2",
        "ContactCount" = "count", "PETCount" = "count", "PETs" = "count", 
        "Contactcount" = "count", "counts" = "count", "V7" = "count",
        "PETCount.1" = "count2",
        "Q-Value" = "qvalue", "FDR" = "qvalue", "qvalue" = "qvalue", 
        "Q-value" = "qvalue", "pvalue" = "qvalue", "V8" = "qvalue"
    )
    
    loop_data <- list()
    
    for (file in loop_files) {
        method_name <- gsub("Sorted_", "", basename(file))
        method_name <- gsub("_qval_count_sorted\\.txt", "", method_name)
        method_name <- extract_method_name(paste0(method_name, "_HiCCUPS.txt"), "HiCCUPS")
        
        loops <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        
        # Apply column mapping
        original_names <- colnames(loops)
        for (i in 1:length(original_names)) {
            if (original_names[i] %in% names(col_mapping)) {
                colnames(loops)[i] <- col_mapping[[original_names[i]]]
            }
        }
        
        required_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
        if (all(required_cols %in% colnames(loops))) {
            loops_clean <- loops[, required_cols]
            loop_data[[method_name]] <- loops_clean
            cat(sprintf("Loaded %d loops for %s\n", nrow(loops_clean), method_name))
        } else {
            cat(sprintf("Warning: Required coordinate columns not found in %s\n", file))
            cat(sprintf("Available columns: %s\n", paste(colnames(loops), collapse = ", ")))
        }
    }
    
    return(loop_data)
}

# Fast spatial overlap calculation using data.table non-equi join
calculate_spatial_overlap_fast <- function(loops1, loops2, tolerance = 5000) {
    if (nrow(loops1) == 0 || nrow(loops2) == 0) return(0)
    
    # Convert to data.table and add expanded coordinate ranges
    dt1 <- data.table(loops1)
    dt1[, `:=`(
        start1_min = start1 - tolerance, start1_max = start1 + tolerance,
        end1_min = end1 - tolerance, end1_max = end1 + tolerance,
        start2_min = start2 - tolerance, start2_max = start2 + tolerance,
        end2_min = end2 - tolerance, end2_max = end2 + tolerance,
        id1 = .I
    )]
    
    dt2 <- data.table(loops2)
    dt2[, id2 := .I]
    
    # Perform non-equi join to find overlapping loops
    overlaps <- dt2[dt1, 
        .(id1, id2), 
        on = .(chr1, chr2,
               start1 >= start1_min, start1 <= start1_max,
               end1 >= end1_min, end1 <= end1_max,
               start2 >= start2_min, start2 <= start2_max,
               end2 >= end2_min, end2 <= end2_max),
        nomatch = NULL,
        allow.cartesian = TRUE
    ]
    
    # Count unique loops from loops1 that have at least one overlap
    unique_overlaps <- length(unique(overlaps$id1))
    
    return(unique_overlaps)
}

# Read all recovery data
read_recovery_data <- function(input_dir, file_pattern, ref_label, reference_loop = NULL) {
    files <- list.files(path = input_dir, pattern = file_pattern, full.names = FALSE)
    
    if (length(files) == 0) {
        stop(sprintf("No files found matching pattern '%s' in '%s'", file_pattern, input_dir))
    }
    
    combined_data <- data.frame()
    method_stats <- data.frame()
    total_ref_loops <- count_gold_loops(input_dir, reference_loop)
    
    for (file in files) {
        method <- extract_method_name(file, ref_label)
        full_path <- file.path(input_dir, file)
        data <- read.table(full_path, header = TRUE, sep = "\t")
        data$method <- method
        combined_data <- rbind(combined_data, data)
        
        max_idx <- which.max(data$frac)
        final_idx <- nrow(data)
        
        method_stats <- rbind(method_stats, data.frame(
            method = method,
            max_loops = data$loopcnt[max_idx],
            max_recovery = data$frac[max_idx],
            final_loops = data$loopcnt[final_idx],
            final_recovery = data$frac[final_idx]
        ))
    }
    
    return(list(
        combined_data = combined_data,
        method_stats = method_stats,
        total_ref_loops = total_ref_loops
    ))
}

# Plot 1: Recovery Curves
# Plot 1: Recovery Curves with improved spacing
create_recovery_plot <- function(combined_data, method_stats, total_ref_loops, opt) {
    # Sort methods by final recovery performance
    method_stats <- method_stats[order(-method_stats$final_recovery), ]
    colors <- get_method_colors(method_stats$method)
    
    # Create legend labels with recovery percentages
    legend_labels <- sapply(1:nrow(method_stats), function(i) {
        method <- method_stats$method[i]
        recovery <- method_stats$final_recovery[i]
        sprintf("%s (%.1f%%)", method, recovery * 100)
    })
    names(legend_labels) <- method_stats$method
    
    # Set factor levels for consistent ordering
    combined_data$method <- factor(combined_data$method, levels = method_stats$method)
    
    ggplot(combined_data, aes(x = loopcnt, y = frac, color = method)) +
        geom_line(linewidth = opt$line_size, alpha = 0.85) +
        scale_color_manual(values = colors, labels = legend_labels) +
        labs(
            title = opt$plot_title,
            x = "Number of predicted loops",
            y = sprintf("Recovery fraction (N = %s)", format(total_ref_loops, big.mark = ",")),
            color = "Method (final recovery)"  # Use sentence case
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 25)),
            axis.title = element_text(size = 16, face = "bold"),
            # Add spacing between axis labels and axis lines
            axis.text.x = element_text(size = 14, face = "bold", color = "black", 
                                     margin = margin(t = 8)),
            axis.text.y = element_text(size = 14, face = "bold", color = "black",
                                     margin = margin(r = 8)),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 12, face = "bold"),  # Reduced font size for legend
            # Improve legend spacing
            legend.key.height = unit(1.3, "lines"),     # Increase vertical spacing between legend items
            legend.key.width = unit(1.5, "lines"),      # Increase legend key width
            legend.spacing.y = unit(0.3, "cm"),         # Add vertical spacing between legend entries
            # Grid and border styling
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey95", linewidth = 0.5),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            # Legend positioning and styling
            legend.position = c(1, 0),
            legend.justification = c(1, 0),
            legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            legend.box.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            legend.margin = margin(8, 8, 10, 8)  # Increased margins around legend box
        ) +
        scale_x_continuous(labels = function(x) format_numbers(x)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        guides(color = guide_legend(override.aes = list(linewidth = 2)))
}

# Plot 2: Overlap Matrix using real spatial overlaps
create_overlap_plot <- function(method_stats, combined_data, opt, total_ref_loops, loop_data = NULL) {
    # Keep original method names for data lookup
    original_methods <- method_stats$method
    
    # Create display labels with line breaks for sintHiChIP methods
    display_methods <- gsub("sintHiChIP \\(global\\)", "sintHiChIP\n(global)", original_methods)
    display_methods <- gsub("sintHiChIP \\(local\\)", "sintHiChIP\n(local)", display_methods)
    
    n <- length(original_methods)
    
    # Initialize matrices using display names but track original names for data access
    overlap_matrix <- matrix(0, n, n, dimnames = list(display_methods, display_methods))
    diag(overlap_matrix) <- 1.0
    
    # Store raw overlap counts
    overlap_counts <- matrix(0, n, n, dimnames = list(display_methods, display_methods))
    
    if (!is.null(loop_data) && length(loop_data) > 0) {
        cat("Calculating spatial overlaps using data.table non-equi join (5K tolerance)...\n")
        
        for (i in 1:n) {
            for (j in 1:n) {
                if (i != j) {
                    # Use original method names for data lookup
                    method_i <- original_methods[i]
                    method_j <- original_methods[j]
                    
                    if (method_i %in% names(loop_data) && method_j %in% names(loop_data)) {
                        loops_i <- loop_data[[method_i]]
                        loops_j <- loop_data[[method_j]]
                        
                        overlap_count <- calculate_spatial_overlap_fast(loops_i, loops_j, tolerance = 5000)
                        overlap_fraction <- overlap_count / nrow(loops_i)
                        
                        overlap_matrix[i, j] <- overlap_fraction
                        overlap_counts[i, j] <- overlap_count
                        
                        cat(sprintf("  %s vs %s: %d/%d = %.3f\n", 
                                  method_i, method_j, overlap_count, nrow(loops_i), overlap_fraction))
                    } else {
                        cat(sprintf("  Warning: Missing loop data for %s or %s\n", method_i, method_j))
                    }
                }
            }
        }
        # Set diagonal counts to the total loops for each method
        for (i in 1:n) {
            # Use original method names for data lookup
            method_i <- original_methods[i]
            if (method_i %in% names(loop_data)) {
                overlap_counts[i, i] <- nrow(loop_data[[method_i]])
            } else {
                # Find matching method in method_stats using original name
                method_idx <- which(method_stats$method == method_i)
                if (length(method_idx) > 0) {
                    overlap_counts[i, i] <- method_stats$final_loops[method_idx]
                } else {
                    cat(sprintf("  Warning: Cannot find final_loops for method %s\n", method_i))
                    overlap_counts[i, i] <- 0
                }
            }
        }
    } else {
        cat("Warning: No loop data available. Using recovery-based overlap estimation...\n")
        
        for (i in 1:n) {
            for (j in 1:n) {
                if (i != j) {
                    # Use original method names for data lookup
                    method_i <- original_methods[i]
                    method_j <- original_methods[j]
                    
                    data_i <- combined_data[combined_data$method == method_i, ]
                    data_j <- combined_data[combined_data$method == method_j, ]
                    
                    # Find matching method in method_stats using original name
                    method_idx_i <- which(method_stats$method == method_i)
                    if (length(method_idx_i) > 0) {
                        loops_i_final <- method_stats$final_loops[method_idx_i]
                    } else {
                        cat(sprintf("  Warning: Cannot find final_loops for method %s\n", method_i))
                        next
                    }
                    
                    recovery_j_at_i_loops <- approx(data_j$loopcnt, data_j$frac, 
                                                   xout = loops_i_final, rule = 2)$y
                    
                    overlap_frac <- recovery_j_at_i_loops
                    
                    if (is.na(overlap_frac) || !is.finite(overlap_frac) || overlap_frac < 0) {
                        overlap_frac <- 0
                    } else if (overlap_frac > 1) {
                        overlap_frac <- 1
                    }
                    
                    overlap_matrix[i, j] <- overlap_frac
                }
            }
            # Set diagonal counts to final_loops for recovery-based case
            method_idx_i <- which(method_stats$method == original_methods[i])
            if (length(method_idx_i) > 0) {
                overlap_counts[i, i] <- method_stats$final_loops[method_idx_i]
            } else {
                cat(sprintf("  Warning: Cannot find final_loops for method %s\n", original_methods[i]))
                overlap_counts[i, i] <- 0
            }
        }
    }
    
    # Debug: Check if overlap_matrix has valid data
    if (all(overlap_matrix == 0 & diag(overlap_matrix) == 0)) {
        cat("Warning: Overlap matrix contains only zeros (no overlaps calculated). Check input data.\n")
    }
    
    overlap_long <- reshape2::melt(overlap_matrix)
    overlap_counts_long <- reshape2::melt(overlap_counts)
    
    overlap_long$label <- ifelse(
        overlap_long$Var1 == overlap_long$Var2,
        sprintf("%s\n(100%%)", format(overlap_counts_long$value[overlap_long$Var1 == overlap_long$Var2], big.mark = ",")),
        sprintf("%s\n(%.0f%%)", format(overlap_counts_long$value, big.mark = ","), overlap_long$value * 100)
    )
    
    if (grepl("from", opt$plot_title, ignore.case = TRUE)) {
        overlap_title <- trimws(sub(".*from\\s+", "", opt$plot_title, ignore.case = TRUE))
    } else {
        overlap_title <- "Method Overlap Matrix"
    }
    
    ggplot(overlap_long, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = label), 
                  color = ifelse(overlap_long$value > 0.6, "white", "black"), 
                  fontface = "bold", size = 4) +
        scale_fill_gradientn(
            colors = c("#E6F3FF", "#FEE090", "#D73027"),
            name = "Overlap\nfraction",
            labels = scales::percent_format(accuracy = 1),
            values = scales::rescale(c(0, 0.5, 1))
        ) +
        labs(
            title = overlap_title,
            x = "Method", y = "Method"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 25)),
            axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
            axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
            # FIX: Use 90-degree rotation for x-axis labels to prevent overlap
            axis.text.x = element_text(size = 14, face = "bold", color = "black", 
                          angle = 90, hjust = 1, vjust = 0.5,
                          margin = margin(t = 8)),  
            axis.text.y = element_text(size = 14, face = "bold", color = "black",
                          margin = margin(r = 8)), 
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 14, face = "bold"),
            panel.grid = element_blank(),
            aspect.ratio = 1,
            # Increase bottom margin for rotated labels
            plot.margin = margin(t = 15, r = 15, b = 50, l = 15)
        ) +
        coord_fixed() +
        scale_x_discrete() +
        scale_y_discrete()
}

# Plot 3: Total Loops Barplot (horizontal)
create_loops_barplot <- function(method_stats, opt) {
    method_stats_sorted <- method_stats[order(-method_stats$final_loops), ]
    colors <- get_method_colors(method_stats_sorted$method)
    
    if (grepl("from", opt$plot_title, ignore.case = TRUE)) {
        barplot_title <- trimws(sub(".*from\\s+", "", opt$plot_title, ignore.case = TRUE))
    } else {
        barplot_title <- "Total Loop Predictions"
    }
    
    ggplot(method_stats_sorted, aes(x = reorder(method, final_loops), 
                                   y = final_loops, fill = method)) +
        geom_col(alpha = 0.8, width = 0.7, color = "black", linewidth = 0.8) + 
        # geom_text(aes(label = format_numbers(final_loops)), 
        geom_text(aes(label = format(final_loops, big.mark = ",")), 
                 hjust = -0.1, size = 6, fontface = "bold", color = "black") +  
        coord_flip() +
        scale_fill_manual(values = colors) +
        labs(
            title = barplot_title,
            x = "Method",
            y = "Number of predicted loops"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(t = 20, b = 30)),  
            axis.title = element_text(size = 18, face = "bold"),  
            # Add spacing to axis labels to prevent crowding
            axis.text.x = element_text(size = 16, face = "bold", color = "black", 
                          margin = margin(t = 15)),  
            axis.text.y = element_text(size = 16, face = "bold", color = "black",
                          margin = margin(r = 15)),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey85", linewidth = 0.8),  
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),  
            legend.position = "none",
            plot.margin = margin(20, 20, 20, 20)  
        ) +
        # scale_y_continuous(labels = function(x) format_numbers(x),
        scale_y_continuous(labels = function(x) format(x, big.mark = ","),
                          expand = expansion(mult = c(0, 0.15)))
}
# Main execution
cat("=== HiChIP Recovery Script ===\n")

if (!dir.exists(opt$input_dir)) {
    stop(sprintf("Input directory does not exist: %s", opt$input_dir))
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading recovery data...\n")
results <- read_recovery_data(opt$input_dir, opt$file_pattern, opt$ref_label, opt$reference_loop)
combined_data <- results$combined_data
method_stats <- results$method_stats
total_ref_loops <- results$total_ref_loops

cat(sprintf("Loaded data for %d methods, %s reference loops\n", 
           nrow(method_stats), format(total_ref_loops, big.mark = ",")))

cat("Reading original loop files for overlap calculation...\n")
loop_data <- read_loop_files(opt$input_dir)

cat("Creating plots...\n")
plot_recovery <- create_recovery_plot(combined_data, method_stats, total_ref_loops, opt)
plot_overlap <- create_overlap_plot(method_stats, combined_data, opt, total_ref_loops, loop_data)
plot_loops <- create_loops_barplot(method_stats, opt)

cat("Saving individual plots...\n")

recovery_file <- file.path(opt$output_dir, paste0(opt$output_name, "_Recovery_Curves.pdf"))
ggsave(recovery_file, plot_recovery, width = 10, height = 8, dpi = 600)
cat(sprintf("Saved: %s\n", basename(recovery_file)))

overlap_file <- file.path(opt$output_dir, paste0(opt$output_name, "_Overlap_Matrix.pdf"))
ggsave(overlap_file, plot_overlap, width = 10, height = 10, dpi = 600)
cat(sprintf("Saved: %s\n", basename(overlap_file)))

loops_file <- file.path(opt$output_dir, paste0(opt$output_name, "_Total_Loops_Barplot.pdf"))
ggsave(loops_file, plot_loops, width = 14, height = 10, dpi = 600)
cat(sprintf("Saved: %s\n", basename(loops_file)))

if (opt$save_png) {
    ggsave(gsub("\\.pdf$", ".png", recovery_file), plot_recovery, 
           width = 10, height = 8, dpi = 300)
    ggsave(gsub("\\.pdf$", ".png", overlap_file), plot_overlap, 
           width = opt$plot_width, height = opt$plot_height, dpi = 300)
    ggsave(gsub("\\.pdf$", ".png", loops_file), plot_loops, 
           width = opt$plot_width, height = opt$plot_height, dpi = 300)
}

cat("\n=== Method Performance Summary ===\n")
ranking <- method_stats %>%
    arrange(desc(final_recovery)) %>%
    mutate(
        rank = row_number(),
        recovery_pct = sprintf("%.1f%%", final_recovery * 100),
        loops_fmt = format(final_loops, big.mark = ",")
    ) %>%
    select(rank, method, loops_fmt, recovery_pct)

names(ranking) <- c("Rank", "Method", "Total Loops", "Final Recovery")
print(ranking)

cat("\n=== Files Generated ===\n")
cat("Individual plots:\n")
cat(sprintf("  - %s\n", basename(recovery_file)))
cat(sprintf("  - %s\n", basename(overlap_file)))
cat(sprintf("  - %s\n", basename(loops_file)))

cat("\n=== Script Complete ===\n")