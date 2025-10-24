#!/usr/bin/env Rscript
# Complete HiChIP Loop Recovery Analysis with Multi-method Comparison
# Author: Bioinformatics Analysis Pipeline
# Description: Comprehensive script for analyzing HiChIP loop recovery rates
#              comparing multiple prediction methods against reference datasets

# Load required libraries
suppressMessages({
    library(GenomicRanges)
    library(data.table)
    library(optparse)
    library(ggplot2)
    library(scales)
    library(dplyr)
    library(RColorBrewer)
})

# Set global options
options(scipen = 10)  # Prevent scientific notation
options(datatable.fread.datatable = FALSE)  # Return data.frame instead of data.table

# ==============================================================================
# CORE FUNCTIONS FOR DATA PROCESSING
# ==============================================================================

#' Read contact data with various format handling
#' 
#' @param inpfile Input file path
#' @param headerInp Logical, whether file has header
#' @param chrNUM Logical, whether chromosomes are numeric (1,2,3) vs named (chr1,chr2,chr3)
#' @param chrSET Vector of chromosome column indices
#' @param mid Logical, whether file contains midpoints instead of intervals
#' @param midCOL Vector of midpoint column indices
#' @param binSize Bin size for converting midpoints to intervals
#' @return Data frame with processed contact data
ReadContactData <- function(inpfile, headerInp = TRUE, chrNUM = FALSE, 
                           chrSET = c(1, 4), mid = FALSE, midCOL = c(2), 
                           binSize = 5000) {
    
    # Read input file
    outDF <- data.table::fread(inpfile, header = headerInp)
    
    # Convert midpoints to intervals if specified
    if (mid == TRUE && length(midCOL) > 0 && binSize > 0) {
        n <- ncol(outDF)
        ModDF <- NULL
        
        for (i in 1:n) {
            if (i %in% midCOL) {
                # Convert midpoint to start and end coordinates
                s1 <- outDF[, i] - (binSize / 2)
                e1 <- outDF[, i] + (binSize / 2)
                
                if (is.null(ModDF)) {
                    ModDF <- cbind.data.frame(s1, e1)
                } else {
                    ModDF <- cbind.data.frame(ModDF, s1, e1)
                }
            } else {
                # Keep original column
                if (is.null(ModDF)) {
                    ModDF <- outDF[, i]
                } else {
                    ModDF <- cbind.data.frame(ModDF, outDF[, i])
                }
            }
        }
        outDF <- ModDF
    }
    
    # Convert numeric chromosomes to named format
    if (chrNUM == TRUE) {
        for (i in 1:length(chrSET)) {
            colno <- chrSET[i]
            outDF[, colno] <- paste0('chr', as.character(outDF[, colno]))
        }
    }
    
    return(outDF)
}

#' Compute overlap between two loop datasets
#' 
#' @param Inpdata1 First dataset (query)
#' @param Inpdata2 Second dataset (reference)
#' @param boundary Boundary offset for overlap detection
#' @param offset Additional offset parameter
#' @param uniqov Logical, whether to return unique overlaps
#' @param IDX Logical, whether to return indices only
#' @return List containing overlap statistics and data frames
OverlapLoop <- function(Inpdata1, Inpdata2, boundary = 1, offset = 0, 
                       uniqov = TRUE, IDX = FALSE) {
    
    # Compute overlap between first intervals (columns 1-3)
    ov1 <- as.data.frame(findOverlaps(
        GRanges(Inpdata1[, 1], IRanges(Inpdata1[, 2] + boundary - offset, 
                                      Inpdata1[, 3] - boundary + offset)),
        GRanges(Inpdata2[, 1], IRanges(Inpdata2[, 2] + boundary - offset, 
                                      Inpdata2[, 3] - boundary + offset))
    ))
    
    # Compute overlap between second intervals (columns 4-6)
    ov2 <- as.data.frame(findOverlaps(
        GRanges(Inpdata1[, 4], IRanges(Inpdata1[, 5] + boundary - offset, 
                                      Inpdata1[, 6] - boundary + offset)),
        GRanges(Inpdata2[, 4], IRanges(Inpdata2[, 5] + boundary - offset, 
                                      Inpdata2[, 6] - boundary + offset))
    ))
    
    # Find common overlaps (both intervals must overlap)
    overlap_uniq_mat <- ov1[unique(which(paste(ov1[, 1], ov1[, 2], sep = ".") %in% 
                                        paste(ov2[, 1], ov2[, 2], sep = "."))), ]
    
    # Handle reverse orientation (swap interval positions)
    ov1A <- as.data.frame(findOverlaps(
        GRanges(Inpdata1[, 4], IRanges(Inpdata1[, 5] + boundary - offset, 
                                      Inpdata1[, 6] - boundary + offset)),
        GRanges(Inpdata2[, 1], IRanges(Inpdata2[, 2] + boundary - offset, 
                                      Inpdata2[, 3] - boundary + offset))
    ))
    
    ov2A <- as.data.frame(findOverlaps(
        GRanges(Inpdata1[, 1], IRanges(Inpdata1[, 2] + boundary - offset, 
                                      Inpdata1[, 3] - boundary + offset)),
        GRanges(Inpdata2[, 4], IRanges(Inpdata2[, 5] + boundary - offset, 
                                      Inpdata2[, 6] - boundary + offset))
    ))
    
    overlap_uniq_mat_1 <- ov1A[unique(which(paste(ov1A[, 1], ov1A[, 2], sep = ".") %in% 
                                           paste(ov2A[, 1], ov2A[, 2], sep = "."))), ]
    
    # Combine overlaps from both orientations
    if (uniqov == TRUE) {
        ov_idx_file1 <- unique(c(overlap_uniq_mat[, 1], overlap_uniq_mat_1[, 1]))
        ov_idx_file2 <- unique(c(overlap_uniq_mat[, 2], overlap_uniq_mat_1[, 2]))
    } else {
        ov_idx_file1 <- c(overlap_uniq_mat[, 1], overlap_uniq_mat_1[, 1])
        ov_idx_file2 <- c(overlap_uniq_mat[, 2], overlap_uniq_mat_1[, 2])
    }
    
    # Find non-overlapping indices
    nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
    nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)
    
    # Return results
    if (IDX == FALSE) {
        newList <- list(
            A_AND_B = ov_idx_file1, 
            B_AND_A = ov_idx_file2, 
            A_MINUS_B = nonov_idx_file1, 
            B_MINUS_A = nonov_idx_file2,
            A_AND_B.df = Inpdata1[ov_idx_file1, ], 
            B_AND_A.df = Inpdata2[ov_idx_file2, ], 
            A_MINUS_B.df = Inpdata1[nonov_idx_file1, ], 
            B_MINUS_A.df = Inpdata2[nonov_idx_file2, ]
        )
    } else {
        newList <- list(
            A_AND_B = ov_idx_file1, 
            B_AND_A = ov_idx_file2, 
            A_MINUS_B = nonov_idx_file1, 
            B_MINUS_A = nonov_idx_file2
        )
    }
    
    return(newList)
}

# ==============================================================================
# RECOVERY ANALYSIS FUNCTIONS
# ==============================================================================

#' Generate pseudo q-values from contact counts for hichipper data
#' 
#' @param contact_counts Vector of contact count values
#' @return Vector of pseudo q-values (higher counts get lower q-values)
GeneratePseudoQvalues <- function(contact_counts) {
    # Rank contact counts (highest count gets rank 1)
    contact_ranks <- rank(-contact_counts, ties.method = "first")
    
    # Convert ranks to pseudo q-values (rank 1 gets smallest q-value)
    # Scale to range [1e-16, 1] to avoid zero q-values
    max_rank <- max(contact_ranks)
    pseudo_qvals <- (contact_ranks - 1) / (max_rank - 1)
    pseudo_qvals <- pseudo_qvals * (1 - 1e-16) + 1e-16
    
    return(pseudo_qvals)
}

#' Perform step-wise recovery analysis
#' 
#' @param opt Options list containing analysis parameters
#' @param method_label Method label for identifying hichipper data
#' @return List containing recovery data and sorted interactions
PerformRecoveryAnalysis <- function(opt, method_label = NULL) {
    
    cat("\n=== Starting Recovery Analysis ===\n")
    
    # Read reference contact data
    cat("Reading reference file:", opt$RefFile, "\n")
    RefData <- ReadContactData(
        opt$RefFile, 
        headerInp = opt$headerRef, 
        chrNUM = opt$chrRef, 
        chrSET = c(1, 4), 
        mid = opt$midRef, 
        midCOL = c(2, 4), 
        binSize = opt$binsizeRef
    )
    
    # Read input contact data
    cat("Reading input file:", opt$InpFile, "\n")
    InpData <- ReadContactData(
        opt$InpFile, 
        headerInp = opt$headerInp, 
        chrNUM = opt$chrInp, 
        chrSET = c(1, 4), 
        mid = opt$midInp, 
        midCOL = c(2, 4), 
        binSize = opt$binsizeInp
    )
    
    # Determine q-value column
    if (opt$QcolInp == 0) {
        QcolInp <- ncol(InpData)  # Last column
    } else {
        QcolInp <- opt$QcolInp
        if (opt$midInp) {
            QcolInp <- QcolInp + 2  # Adjust for added columns
        }
    }
    
    # Determine contact count column
    if (opt$CountcolInp == 0) {
        CountcolInp <- QcolInp - 1  # Column before q-value
    } else {
        CountcolInp <- opt$CountcolInp
        if (opt$midInp) {
            CountcolInp <- CountcolInp + 2  # Adjust for added columns
        }
    }
    
    # Check if this is hichipper data and handle pseudo q-value generation
    is_hichipper <- !is.null(method_label) && grepl("hichipper", method_label, ignore.case = TRUE)
    
    if (is_hichipper) {
        cat("Detected hichipper data - generating pseudo q-values from contact counts\n")
        
        # Validate contact count column
        if (CountcolInp <= 0 || CountcolInp > ncol(InpData)) {
            stop("Invalid contact count column index for hichipper data. Check --CountcolInp parameter.\n")
        }
        
        # Generate pseudo q-values from contact counts
        contact_counts <- InpData[, CountcolInp]
        pseudo_qvals <- GeneratePseudoQvalues(contact_counts)
        
        # Add pseudo q-values as new column
        InpData <- cbind(InpData, pseudo_qvals)
        QcolInp <- ncol(InpData)  # Update q-value column index
        
        cat(sprintf("Generated pseudo q-values: range [%.2e, %.2e]\n", 
                   min(pseudo_qvals), max(pseudo_qvals)))
    } else {
        # Validate column indices for non-hichipper data
        if (CountcolInp <= 0 || CountcolInp > ncol(InpData)) {
            stop("Invalid contact count column index. Check --CountcolInp parameter.\n")
        }
        if (QcolInp <= 0 || QcolInp > ncol(InpData)) {
            stop("Invalid q-value column index. Check --QcolInp parameter.\n")
        }
    }
    
    cat(sprintf("Using column %d for contact counts and column %d for q-values\n", 
                CountcolInp, QcolInp))
    
    # Filter data to required columns
    RefData <- RefData[, 1:6]
    InpData <- InpData[, c(1:6, CountcolInp, QcolInp)]
    
    cat(sprintf("Reference data: %d interactions\n", nrow(RefData)))
    cat(sprintf("Input data: %d interactions\n", nrow(InpData)))
    
    # Sort input data: first by q-value (ascending), then by contact count (descending)
    cat("Sorting interactions by q-value and contact count...\n")
    InpData <- InpData[order(InpData[, 8], -InpData[, 7]), ]
    
    # Prepare data for overlap analysis
    InpData_for_analysis <- InpData[, c(1:6, 8)]  # Remove contact count column
    
    # Perform step-wise overlap analysis
    cat("Starting step-wise overlap analysis...\n")
    stepsize <- 1000
    finalDF <- data.frame()
    
    for (startidx in seq(1, nrow(InpData_for_analysis), by = stepsize)) {
        endidx <- min((startidx + stepsize - 1), nrow(InpData_for_analysis))
        
        # Get subset of interactions
        loopdata <- InpData_for_analysis[1:endidx, ]
        
        # Compute overlap with reference
        ov <- OverlapLoop(loopdata, RefData, boundary = 1, offset = opt$offset)
        
        # Calculate recovery fraction
        fracval <- (length(ov$B_AND_A) * 1.0) / nrow(RefData)
        
        # Store results
        currDF <- data.frame(loopcnt = endidx, frac = fracval)
        finalDF <- rbind(finalDF, currDF)
        
        if (endidx %% 10000 == 0 || endidx == nrow(InpData_for_analysis)) {
            cat(sprintf("Processed %d interactions, recovery rate: %.4f\n", 
                       endidx, fracval))
        }
    }
    
    cat("Recovery analysis completed!\n")
    
    return(list(recovery_data = finalDF, sorted_interactions = InpData))
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create publication-quality recovery plot for single method
#' 
#' @param recovery_data Data frame with recovery analysis results
#' @param method_label Label for the method
#' @param plot_title Title for the plot
#' @param output_file Output file path
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @param line_color Color for the line
#' @param log_scale Whether to use log scale for x-axis
#' @return ggplot object
CreateSingleRecoveryPlot <- function(recovery_data, method_label = "Method",
                                   plot_title = "Recovery of reference loops",
                                   output_file = "recovery_plot.pdf",
                                   plot_width = 8, plot_height = 6,
                                   line_color = "#0072B2", log_scale = FALSE) {
    
    # Get maximum values for annotation
    max_idx <- which.max(recovery_data$loopcnt)
    max_count <- recovery_data$loopcnt[max_idx]
    max_frac <- recovery_data$frac[max_idx]
    
    # Create plot
    p <- ggplot(recovery_data, aes(x = loopcnt, y = frac)) +
        geom_line(color = line_color, size = 1.2, alpha = 0.8) +
        geom_point(data = recovery_data[max_idx, ], 
                  aes(x = loopcnt, y = frac), 
                  color = line_color, size = 3, alpha = 0.9) +
        labs(
            title = plot_title,
            x = "Number of loops predicted",
            y = "Fraction of reference loops",
            subtitle = sprintf("%s (N = %s)", method_label, format(max_count, big.mark = ","))
        ) +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey90", size = 0.5),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            plot.margin = margin(20, 20, 20, 20)
        )
    
    # Set axis scales
    if (log_scale) {
        p <- p + scale_x_log10(labels = comma_format())
    } else {
        p <- p + scale_x_continuous(labels = comma_format())
    }
    
    p <- p + scale_y_continuous(labels = percent_format(accuracy = 1))
    
    # Save plot
    ggsave(output_file, p, width = plot_width, height = plot_height, dpi = 300)
    cat(sprintf("Single method plot saved to: %s\n", output_file))
    
    return(p)
}

#' Create multi-method comparison recovery plot
#' 
#' @param recovery_files Vector of recovery data file paths
#' @param method_labels Vector of method labels
#' @param method_colors Vector of colors for each method
#' @param plot_title Plot title
#' @param output_file Output file path
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @param log_scale Whether to use log scale for x-axis
#' @param show_counts Whether to show sample counts in legend
#' @return List containing plot object and statistics
CreateMultiMethodRecoveryPlot <- function(recovery_files, method_labels, 
                                        method_colors = NULL,
                                        plot_title = "Recovery of reference loops",
                                        output_file = "multi_method_recovery.pdf",
                                        plot_width = 12, plot_height = 8,
                                        log_scale = FALSE, show_counts = TRUE) {
    
    # Validate inputs
    if (length(recovery_files) != length(method_labels)) {
        stop("Number of files must match number of method labels")
    }
    
    # Read and combine all recovery data
    combined_data <- data.frame()
    method_stats <- data.frame()
    
    for (i in 1:length(recovery_files)) {
        if (!file.exists(recovery_files[i])) {
            warning(sprintf("File does not exist: %s", recovery_files[i]))
            next
        }
        
        # Read recovery data
        data <- read.table(recovery_files[i], header = TRUE, sep = "\t")
        data$method <- method_labels[i]
        combined_data <- rbind(combined_data, data)
        
        # Calculate statistics for this method
        max_idx <- which.max(data$loopcnt)
        method_stats <- rbind(method_stats, 
                            data.frame(
                                method = method_labels[i],
                                max_count = data$loopcnt[max_idx],
                                max_frac = data$frac[max_idx],
                                final_frac = data$frac[nrow(data)]
                            ))
    }
    
    if (nrow(combined_data) == 0) {
        stop("No valid data files found")
    }
    
    # Set default colors if not provided
    if (is.null(method_colors)) {
        # Use color-blind friendly palette
        method_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                          "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999")
        method_colors <- method_colors[1:length(method_labels)]
    }
    
    # Order methods by performance (max recovery rate, then max count)
    method_stats <- method_stats[order(-method_stats$final_frac, -method_stats$max_count), ]
    method_order <- method_stats$method
    
    # Create legend labels
    if (show_counts) {
        legend_labels <- sapply(1:nrow(method_stats), function(i) {
            method <- method_stats$method[i]
            count <- method_stats$max_count[i]
            sprintf("%s (N = %s)", method, format(count, big.mark = ","))
        })
        names(legend_labels) <- method_order
    } else {
        legend_labels <- method_order
        names(legend_labels) <- method_order
    }
    
    # Set factor levels for proper ordering
    combined_data$method <- factor(combined_data$method, levels = method_order)
    method_stats$method <- factor(method_stats$method, levels = method_order)
    
    # Create the plot
    p <- ggplot(combined_data, aes(x = loopcnt, y = frac, color = method)) +
        geom_line(size = 1.2, alpha = 0.85) +
        geom_point(data = method_stats, 
                  aes(x = max_count, y = max_frac, color = method), 
                  size = 3, alpha = 0.9) +
        scale_color_manual(values = method_colors, labels = legend_labels) +
        labs(
            title = plot_title,
            x = "Number of loops predicted",
            y = "Fraction of reference loops",
            color = "Method"
        ) +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_line(color = "grey90", size = 0.5),
            panel.grid.major.y = element_line(color = "grey90", size = 0.5),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "black", size = 0.5),
            legend.key = element_rect(fill = "white"),
            legend.key.width = unit(1.5, "cm"),
            plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 14),
            legend.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 13),
            plot.margin = margin(25, 25, 25, 25)
        ) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    # Set axis scales
    if (log_scale) {
        p <- p + scale_x_log10(
            labels = comma_format(),
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            minor_breaks = scales::trans_breaks("log10", function(x) 10^x)
        )
    } else {
        p <- p + scale_x_continuous(
            labels = comma_format(),
            breaks = scales::pretty_breaks(n = 8)
        )
    }
    
    p <- p + scale_y_continuous(
        labels = percent_format(accuracy = 1),
        breaks = scales::pretty_breaks(n = 6),
        limits = c(0, max(combined_data$frac) * 1.05)
    )
    
    # Save plot
    ggsave(output_file, p, width = plot_width, height = plot_height, dpi = 300)
    cat(sprintf("Multi-method comparison plot saved to: %s\n", output_file))
    
    # Print performance summary
    cat("\n=== Method Performance Summary ===\n")
    for (i in 1:nrow(method_stats)) {
        cat(sprintf("%-20s: Max predictions = %8s, Max recovery = %5.1f%%, Final recovery = %5.1f%%\n", 
                   method_stats$method[i], 
                   format(method_stats$max_count[i], big.mark = ","),
                   method_stats$max_frac[i] * 100,
                   method_stats$final_frac[i] * 100))
    }
    
    return(list(plot = p, stats = method_stats, data = combined_data))
}

# ==============================================================================
# COMMAND LINE INTERFACE
# ==============================================================================

# Define command line options
option_list <- list(
    # Input files
    make_option(c("--RefFile"), type = "character", default = NULL, 
                help = "Reference loop file (ground truth interactions). Mandatory."),
    make_option(c("--InpFile"), type = "character", default = NULL, 
                help = "Input HiChIP loop file for analysis. Mandatory for single analysis."),
    make_option(c("--MultiFiles"), type = "character", default = NULL, 
                help = "Comma-separated list of input files for multi-method comparison."),
    make_option(c("--MultiLabels"), type = "character", default = NULL, 
                help = "Comma-separated list of method labels for multi-method comparison."),
    
    # File format options
    make_option(c("--headerRef"), action = "store_true", default = FALSE, 
                help = "Reference file has header. Default: FALSE"),
    make_option(c("--headerInp"), action = "store_true", default = FALSE, 
                help = "Input file has header. Default: FALSE"),
    make_option(c("--chrRef"), action = "store_true", default = FALSE, 
                help = "Reference file uses numeric chromosomes (1,2,3). Default: FALSE"),
    make_option(c("--chrInp"), action = "store_true", default = FALSE, 
                help = "Input file uses numeric chromosomes (1,2,3). Default: FALSE"),
    make_option(c("--midRef"), action = "store_true", default = FALSE, 
                help = "Reference file contains midpoints instead of intervals. Default: FALSE"),
    make_option(c("--midInp"), action = "store_true", default = FALSE, 
                help = "Input file contains midpoints instead of intervals. Default: FALSE"),
    
    # Column specifications
    make_option(c("--QcolInp"), type = "integer", default = 0, 
                help = "Q-value column number (0 = last column). Default: 0"),
    make_option(c("--CountcolInp"), type = "integer", default = 0, 
                help = "Contact count column number (0 = auto-detect). Default: 0"),
    
    # Analysis parameters
    make_option(c("--binsizeRef"), type = "integer", default = 5000, 
                help = "Bin size for reference file (bp). Default: 5000"),
    make_option(c("--binsizeInp"), type = "integer", default = 5000, 
                help = "Bin size for input file (bp). Default: 5000"),
    make_option(c("--offset"), type = "integer", default = 5000, 
                help = "Overlap tolerance (bp). Default: 5000"),
    
    # Output options
    make_option(c("--OutDir"), type = "character", default = "./output", 
                help = "Output directory. Default: ./output"),
    make_option(c("--LabelRef"), type = "character", default = "Reference", 
                help = "Reference dataset label. Default: Reference"),
    make_option(c("--LabelInp"), type = "character", default = "Input", 
                help = "Input dataset label. Default: Input"),
    
    # Plotting options
    make_option(c("--plotOutput"), action = "store_true", default = FALSE, 
                help = "Generate plots. Default: FALSE"),
    make_option(c("--plotTitle"), type = "character", default = "Recovery Analysis", 
                help = "Plot title. Default: Recovery Analysis"),
    make_option(c("--plotWidth"), type = "double", default = 10, 
                help = "Plot width (inches). Default: 10"),
    make_option(c("--plotHeight"), type = "double", default = 8, 
                help = "Plot height (inches). Default: 8"),
    make_option(c("--logScale"), action = "store_true", default = FALSE, 
                help = "Use log scale for x-axis. Default: FALSE"),
    make_option(c("--multiMethod"), action = "store_true", default = FALSE, 
                help = "Perform multi-method comparison. Default: FALSE")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

# Validate required parameters
if (is.null(opt$RefFile)) {
    cat("ERROR: Reference file is required. Use --RefFile parameter.\n")
    print_help(opt_parser)
    quit(status = 1)
}

if (!opt$multiMethod && is.null(opt$InpFile)) {
    cat("ERROR: Input file is required for single analysis. Use --InpFile parameter.\n")
    print_help(opt_parser)
    quit(status = 1)
}

if (opt$multiMethod && (is.null(opt$MultiFiles) || is.null(opt$MultiLabels))) {
    cat("ERROR: Multi-method analysis requires --MultiFiles and --MultiLabels parameters.\n")
    print_help(opt_parser)
    quit(status = 1)
}

# Create output directory
dir.create(opt$OutDir, recursive = TRUE, showWarnings = FALSE)

# Execute analysis based on mode
if (opt$multiMethod) {
    # ========== MULTI-METHOD COMPARISON MODE ==========
    cat("\n=== Multi-Method Recovery Analysis ===\n")
    
    # Parse file lists
    input_files <- strsplit(opt$MultiFiles, ",")[[1]]
    method_labels <- strsplit(opt$MultiLabels, ",")[[1]]
    
    # Trim whitespace
    input_files <- trimws(input_files)
    method_labels <- trimws(method_labels)
    
    if (length(input_files) != length(method_labels)) {
        stop("Number of input files must match number of method labels")
    }
    
    cat(sprintf("Analyzing %d methods:\n", length(method_labels)))
    for (i in 1:length(method_labels)) {
        cat(sprintf("  %d. %s: %s\n", i, method_labels[i], input_files[i]))
    }
    
    # Process each method
    recovery_files <- character(length(input_files))
    
    for (i in 1:length(input_files)) {
        cat(sprintf("\n--- Processing method %d/%d: %s ---\n", i, length(input_files), method_labels[i]))
        
        # Create temporary opt object for this method
        temp_opt <- opt
        temp_opt$InpFile <- input_files[i]
        temp_opt$LabelInp <- method_labels[i]
        
        # Perform recovery analysis with method label for hichipper detection
        results <- PerformRecoveryAnalysis(temp_opt, method_labels[i])
        
        # Save recovery data
        recovery_file <- file.path(opt$OutDir, sprintf("Recovery_%s_%s_offset_%d.txt", 
                                                      method_labels[i], opt$LabelRef, opt$offset))
        write.table(results$recovery_data, recovery_file, 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        recovery_files[i] <- recovery_file
        
        # Save sorted interactions
        sorted_file <- file.path(opt$OutDir, sprintf("Sorted_%s_qval_count_sorted.txt", method_labels[i]))
        write.table(results$sorted_interactions, sorted_file, 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        
        # cat(sprintf("Recovery data saved: %s\n", recovery_file))
        # cat(sprintf("Sorted data saved: %s\n", sorted_file))
    }
    
    # Generate multi-method comparison plot
    if (opt$plotOutput) {
        cat("\n--- Generating Multi-Method Comparison Plot ---\n")
        
        # Define color palette for methods
        method_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                          "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999")
        method_colors <- method_colors[1:length(method_labels)]
        
        # Create comparison plot
        plot_file <- file.path(opt$OutDir, sprintf("Recovery_Comparison_%s.pdf", opt$LabelRef))
        
        comparison_results <- CreateMultiMethodRecoveryPlot(
            recovery_files = recovery_files,
            method_labels = method_labels,
            method_colors = method_colors,
            plot_title = opt$plotTitle,
            output_file = plot_file,
            plot_width = opt$plotWidth,
            plot_height = opt$plotHeight,
            log_scale = opt$logScale,
            show_counts = TRUE
        )
        
        # Save comparison statistics
        stats_file <- file.path(opt$OutDir, sprintf("Method_Comparison_Stats_%s.txt", opt$LabelRef))
        write.table(comparison_results$stats, stats_file, 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        
        # cat(sprintf("Comparison plot saved: %s\n", plot_file))
        # cat(sprintf("Comparison statistics saved: %s\n", stats_file))
    }
    
} else {
    # ========== SINGLE METHOD ANALYSIS MODE ==========
    # cat("\n=== Single Method Recovery Analysis ===\n")
    # cat(sprintf("Reference file: %s\n", opt$RefFile))
    # cat(sprintf("Input file: %s\n", opt$InpFile))
    
    # Perform recovery analysis with method label for hichipper detection
    results <- PerformRecoveryAnalysis(opt, opt$LabelInp)
    recovery_data <- results$recovery_data
    sorted_interactions <- results$sorted_interactions
    
    # Save recovery analysis results
    recovery_output_file <- file.path(opt$OutDir, sprintf("Recovery_%s_%s_offset_%d_qval_count_sorted.txt", 
                                                         opt$LabelInp, opt$LabelRef, opt$offset))
    write.table(recovery_data, recovery_output_file, 
               row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    # Save sorted interaction data
    sorted_output_file <- file.path(opt$OutDir, sprintf("Sorted_%s_qval_count_sorted.txt", opt$LabelInp))
    write.table(sorted_interactions, sorted_output_file, 
               row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    # cat(sprintf("Recovery analysis saved: %s\n", recovery_output_file))
    # cat(sprintf("Sorted interactions saved: %s\n", sorted_output_file))
    
    # Generate single method plot if requested
    if (opt$plotOutput) {
        cat("\n--- Generating Single Method Plot ---\n")
        
        plot_output_file <- file.path(opt$OutDir, sprintf("Recovery_plot_%s_%s.pdf", 
                                                         opt$LabelInp, opt$LabelRef))
        
        single_plot <- CreateSingleRecoveryPlot(
            recovery_data = recovery_data,
            method_label = opt$LabelInp,
            plot_title = opt$plotTitle,
            output_file = plot_output_file,
            plot_width = opt$plotWidth,
            plot_height = opt$plotHeight,
            log_scale = opt$logScale
        )
        
        # Save plot data
        plot_data_file <- file.path(opt$OutDir, sprintf("Recovery_plot_data_%s_%s.txt", 
                                                        opt$LabelInp, opt$LabelRef))
        write.table(recovery_data, plot_data_file, 
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        
        # cat(sprintf("Single method plot saved: %s\n", plot_output_file))
        # cat(sprintf("Plot data saved: %s\n", plot_data_file))
    }
}

# ==============================================================================
# ANALYSIS COMPLETION
# ==============================================================================




# Rscript hichip_recovery.R \
#   --RefFile GM12878_HiCCUPS_loops.txt \
#   --multiMethod \
#   --MultiFiles "FitHiChIP_L.txt,FitHiChIP_S.txt,hichipper.txt,MAPS.txt" \
#   --MultiLabels "FitHiChIP (L),FitHiChIP (S),hichipper,MAPS" \
#   --OutDir comparison_output \
#   --LabelRef "HiCCUPS" \
#   --plotOutput \
#   --plotTitle "Recovery of GM12878 Hi-C HiCCUPS loops from H3K27ac HiChIP data" \
#   --plotWidth 12 \
#   --plotHeight 8