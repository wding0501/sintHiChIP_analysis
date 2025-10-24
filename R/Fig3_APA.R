suppressMessages({
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(optparse)
  library(RColorBrewer)
  library(cowplot)
  library(grid)
})

# Command line options
option_list = list(
  make_option(c("--plotWidth"), type="integer", default=6),
  make_option(c("--plotHeight"), type="integer", default=6),
  make_option(c("--fontSize"), type="integer", default=14),
  make_option(c("--window"), type="integer", default=100000),
  make_option(c("--useLog"), type="logical", default=FALSE),
  make_option(c("--usePercentile"), type="logical", default=TRUE),
  make_option(c("--percentileCut"), type="numeric", default=0.02),
  make_option(c("--colorScheme"), type="character", default="blue_red"),
  make_option(c("--debug"), type="logical", default=FALSE),
  make_option(c("--cellLine"), type="character", default="GM12878")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

PlotWidth <- opt$plotWidth
PlotHeight <- opt$plotHeight
FontSize <- opt$fontSize
WindowSize <- opt$window
use_log <- opt$useLog
color_scheme <- opt$colorScheme
cellLine <- opt$cellLine

# Function to extract and format method names
extract_method_name <- function(filename, ref_label = "HiCCUPS") {
  method <- gsub("Recovery_plot_data_", "", filename)
  method <- gsub(paste0("_", ref_label, ".txt"), "", method)
  method <- gsub("HiCDCPlus|HiCDC|HiC_DC_Plus", "HiC-DC+", method, ignore.case = TRUE)
  method <- gsub("sintHiChIP\\.local", "sintHiChIP (local)", method, ignore.case = TRUE)
  method <- gsub("sintHiChIP\\.global", "sintHiChIP (global)", method, ignore.case = TRUE)
  return(method)
}

# Global storage for all raw matrices
all_raw_matrices <- list()

# Global variables to store transformation parameters
global_log_min <- NULL
global_log_max <- NULL

# Function to calculate global parameters once for all methods
calculate_global_parameters <- function() {
  percentile_cut <- opt$percentileCut
  lower_percentile <- percentile_cut
  upper_percentile <- 1 - percentile_cut
  
  all_log_values <- unlist(lapply(all_raw_matrices, function(m) as.vector(log2(m + 1))))
  global_log_min <<- quantile(all_log_values, lower_percentile, na.rm = TRUE)
  global_log_max <<- quantile(all_log_values, upper_percentile, na.rm = TRUE)
  
  cat("Global log2 transformation parameters calculated:\n")
  cat("Using percentiles: ", round(lower_percentile*100, 1), "% - ", round(upper_percentile*100, 1), "%\n")
  cat("Log values range: [", round(global_log_min, 3), ",", round(global_log_max, 3), "]\n")
}

# Log2 transformation with FIXED global scaling
standardize_matrix_log <- function(matrix_data) {
  log_matrix <- log2(matrix_data + 1)
  
  scaled_log <- (log_matrix - global_log_min) / (global_log_max - global_log_min)
  scaled_log[scaled_log < 0] <- 0
  scaled_log[scaled_log > 1] <- 1
  
  return(scaled_log)
}

# Function to read APA matrix files
read_apa_matrix <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)

  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) == 0) return(NULL)

  matrix_data <- NULL

  if (grepl("^\\s*\\[", lines[1]) && grepl("\\]\\s*$", lines[1])) {
    matrix_rows <- list()
    for (i in 1:length(lines)) {
      line <- lines[i]
      clean_line <- gsub("^\\s*\\[|\\]\\s*$", "", line)
      if (nchar(trimws(clean_line)) == 0) next
      values <- strsplit(clean_line, ",")[[1]]
      numeric_values <- as.numeric(trimws(values))
      if (any(is.na(numeric_values))) next
      matrix_rows[[length(matrix_rows) + 1]] <- numeric_values
    }
    if (length(matrix_rows) > 0) {
      row_lengths <- sapply(matrix_rows, length)
      if (length(unique(row_lengths)) == 1) {
        matrix_data <- do.call(rbind, matrix_rows)
      }
    }
  }

  if (is.null(matrix_data)) {
    data <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE,
                      stringsAsFactors = FALSE, comment.char = "",
                      blank.lines.skip = TRUE)
    matrix_data <- as.matrix(data)
  }

  if (!is.numeric(matrix_data)) {
    matrix_data <- apply(matrix_data, c(1,2), function(x) {
      num_val <- suppressWarnings(as.numeric(as.character(x)))
      ifelse(is.na(num_val), 0, num_val)
    })
  }

  if (nrow(matrix_data) < 3 || ncol(matrix_data) < 3) return(NULL)
  return(matrix_data)
}

# Function to get loop count from .tmp files
get_loop_count <- function(method_name) {
  tmp_file <- paste0(method_name, "_150k_1Mb.txt")

  if (file.exists(tmp_file)) {
    lines <- readLines(tmp_file, warn = FALSE)
    non_empty_lines <- lines[nchar(trimws(lines)) > 0]

    first_line <- non_empty_lines[1]
    if (grepl("chr|start|end|#", first_line, ignore.case = TRUE)) {
      loop_count <- length(non_empty_lines) - 1
    } else {
      loop_count <- length(non_empty_lines)
    }
    return(loop_count)
  }
  return(NA)
}

# Function to calculate APA regions and scores using RAW COUNTS
calculate_apa_regions <- function(matrix_data) {
  n <- nrow(matrix_data)
  m <- ncol(matrix_data)

  center_row <- as.integer((n - 1) / 2) + 1
  center_col <- as.integer((m - 1) / 2) + 1
  center_value <- matrix_data[center_row, center_col]

  ll_size <- max(1, floor(min(n, m) / 4))
  ll_rows <- (n - ll_size + 1):n
  ll_cols <- 1:ll_size
  ll_rows <- ll_rows[ll_rows >= 1 & ll_rows <= n]
  ll_cols <- ll_cols[ll_cols >= 1 & ll_cols <= m]

  if (length(ll_rows) > 0 && length(ll_cols) > 0) {
    lower_left_mean <- mean(matrix_data[ll_rows, ll_cols], na.rm = TRUE)
  } else {
    lower_left_mean <- mean(c(matrix_data[n, 1], matrix_data[n-1, 1],
                             matrix_data[n, 2], matrix_data[n-1, 2]), na.rm = TRUE)
  }

  p2ll_score <- if(lower_left_mean > 0) center_value / lower_left_mean else center_value

  row_size <- floor(n / 3)
  col_size <- floor(m / 3)

  ul_rows <- 1:row_size
  ul_cols <- 1:col_size
  ul_mean <- mean(matrix_data[ul_rows, ul_cols], na.rm = TRUE)
  apa_ul <- if(ul_mean > 0) center_value / ul_mean else center_value

  ur_rows <- 1:row_size
  ur_cols <- (2 * col_size + 1):m
  ur_mean <- mean(matrix_data[ur_rows, ur_cols], na.rm = TRUE)
  apa_ur <- if(ur_mean > 0) center_value / ur_mean else center_value

  lr_rows <- (2 * row_size + 1):n
  lr_cols <- (2 * col_size + 1):m
  lr_mean <- mean(matrix_data[lr_rows, lr_cols], na.rm = TRUE)
  apa_lr <- if(lr_mean > 0) center_value / lr_mean else center_value

  ul_rows_r <- 1:(center_row - 1)
  ul_cols_r <- 1:(center_col - 1)
  ur_rows_r <- 1:(center_row - 1)
  ur_cols_r <- (center_col + 1):m
  ll_rows_r <- (center_row + 1):n
  ll_cols_r <- 1:(center_col - 1)
  lr_rows_r <- (center_row + 1):n
  lr_cols_r <- (center_col + 1):m

  quadrant_means <- c()
  if (length(ul_rows_r) > 0 && length(ul_cols_r) > 0) {
    quadrant_means <- c(quadrant_means, mean(matrix_data[ul_rows_r, ul_cols_r], na.rm = TRUE))
  }
  if (length(ur_rows_r) > 0 && length(ur_cols_r) > 0) {
    quadrant_means <- c(quadrant_means, mean(matrix_data[ur_rows_r, ur_cols_r], na.rm = TRUE))
  }
  if (length(ll_rows_r) > 0 && length(ll_cols_r) > 0) {
    quadrant_means <- c(quadrant_means, mean(matrix_data[ll_rows_r, ll_cols_r], na.rm = TRUE))
  }
  if (length(lr_rows_r) > 0 && length(lr_cols_r) > 0) {
    quadrant_means <- c(quadrant_means, mean(matrix_data[lr_rows_r, lr_cols_r], na.rm = TRUE))
  }

  if (length(quadrant_means) > 0) {
    four_quadrant_mean <- mean(quadrant_means)
    central_ratio <- if(four_quadrant_mean > 0) center_value / four_quadrant_mean else center_value
  } else {
    central_ratio <- center_value
  }

  return(list(
    apa_scores = list(ul = apa_ul, ur = apa_ur, ll = p2ll_score, lr = apa_lr, central_ratio = central_ratio),
    p2ll = p2ll_score,
    raw_matrix = matrix_data
  ))
}

# Enhanced color palette for log2 transformation visualization
create_log_palette <- function(n_colors = 100) {
  colorRampPalette(c(
    "#08306B",
    "#08519C",
    "#2171B5",
    "#4292C6",
    "#6BAED6",
    "#C6DBEF",
    "#FFFFFF",
    "#FDD0A2",
    "#FDAE6B",
    "#FD8D3C",
    "#E6550D",
    "#A63603",
    "#7F0000"
  ))(n_colors)
}

# Function to create APA plot with log2 transformation
create_apa_plot <- function(matrix_data, method_name, resolution, apa_analysis) {
  
  plot_matrix <- standardize_matrix_log(matrix_data)
  color_limits <- c(0, 1)
  transform_label <- "Log2 Normalized"
  colors <- create_log_palette()
  
  x_coords <- seq(-10, 10, length.out = ncol(plot_matrix))
  y_coords <- seq(10, -10, length.out = nrow(plot_matrix))

  melted_data <- expand.grid(X = x_coords, Y = y_coords)
  melted_data$Value <- as.vector(t(plot_matrix))

  if (grepl("sintHiChIP", method_name, ignore.case = TRUE)) {
    main_title <- paste0("* ", method_name, " *\nAPA: ", round(apa_analysis$p2ll, 2), 
                         " | R: ", round(apa_analysis$apa_scores$central_ratio, 2))
    title_color <- "#FF6600"
    border_color <- "gray50"
    border_size <- 1
  } else {
    main_title <- paste0(method_name, "\nAPA: ", round(apa_analysis$p2ll, 2), 
                         " | R: ", round(apa_analysis$apa_scores$central_ratio, 2))
    title_color <- "black"
    border_color <- "gray50"
    border_size <- 1
  }

  p <- ggplot(melted_data, aes(x = X, y = Y, fill = Value)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      limits = color_limits,
      name = transform_label,
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 1,
        barheight = 8
      )
    ) +
    labs(
      title = main_title,
      x = paste0("Distance from D anchor (", resolution, ")"),
      y = paste0("Distance from U anchor (", resolution, ")")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = FontSize - 2, face = "bold", color = title_color),
      axis.title = element_text(size = FontSize - 3, face = "bold"),
      axis.text = element_text(size = FontSize - 4),
      legend.title = element_text(size = FontSize - 3, face = "bold"),
      legend.text = element_text(size = FontSize - 4),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(color = border_color, fill = NA, size = border_size),
      aspect.ratio = 1
    ) +
    coord_fixed() +
    scale_x_continuous(limits = c(-10.5, 10.5), breaks = seq(-10, 10, by = 5)) +
    scale_y_continuous(limits = c(-10.5, 10.5), breaks = seq(-10, 10, by = 5), labels = function(x) { -x })

  region_size <- 20 / 3

  ul_xmin <- -10; ul_xmax <- -10 + region_size; ul_ymin <- 10 - region_size; ul_ymax <- 10
  ur_xmin <- 10 - region_size; ur_xmax <- 10; ur_ymin <- 10 - region_size; ur_ymax <- 10
  ll_xmin <- -10; ll_xmax <- -10 + region_size; ll_ymin <- -10; ll_ymax <- -10 + region_size
  lr_xmin <- 10 - region_size; lr_xmax <- 10; lr_ymin <- -10; lr_ymax <- -10 + region_size

  p <- p +
    annotate("rect", xmin = ul_xmin, xmax = ul_xmax, ymin = ul_ymin, ymax = ul_ymax,
             color = "gray60", fill = NA, linewidth = 0.8, linetype = "dashed") +
    annotate("rect", xmin = ur_xmin, xmax = ur_xmax, ymin = ur_ymin, ymax = ur_ymax,
             color = "gray60", fill = NA, linewidth = 0.8, linetype = "dashed") +
    annotate("rect", xmin = ll_xmin, xmax = ll_xmax, ymin = ll_ymin, ymax = ll_ymax,
             color = "gray60", fill = NA, linewidth = 0.8, linetype = "dashed") +
    annotate("rect", xmin = lr_xmin, xmax = lr_xmax, ymin = lr_ymin, ymax = lr_ymax,
             color = "gray60", fill = NA, linewidth = 0.8, linetype = "dashed")

  p <- p +
    annotate("text", x = (ul_xmin + ul_xmax)/2, y = (ul_ymin + ul_ymax)/2,
             label = paste0("APA: ", round(apa_analysis$apa_scores$ul, 2)),
             color = "black", size = 3, fontface = "bold") +
    annotate("text", x = (ur_xmin + ur_xmax)/2, y = (ur_ymin + ur_ymax)/2,
             label = paste0("APA: ", round(apa_analysis$apa_scores$ur, 2)),
             color = "black", size = 3, fontface = "bold") +
    annotate("text", x = (ll_xmin + ll_xmax)/2, y = (ll_ymin + ll_ymax)/2,
             label = paste0("APA: ", round(apa_analysis$apa_scores$ll, 2)),
             color = "black", size = 3, fontface = "bold") +
    annotate("text", x = (lr_xmin + lr_xmax)/2, y = (lr_ymin + lr_ymax)/2,
             label = paste0("APA: ", round(apa_analysis$apa_scores$lr, 2)),
             color = "black", size = 3, fontface = "bold") +
    annotate("text", x = 0, y = 0,
             label = paste0("R: ", round(apa_analysis$apa_scores$central_ratio, 2)),
             color = "black", size = 3, fontface = "bold")

  return(p)
}

# Function to create section header
create_section_header <- function(resolution, loop_count, cell_line) {
  if (!is.na(loop_count) && loop_count > 0) {
    return(paste0(cell_line, " H3K27ac: top ", loop_count, " HiChIP loops"))
  } else {
    return(paste0(cell_line, " H3K27ac: loops unavailable"))
  }
}

# Function to create and save plot grid (2x4 layout)
create_plot_grid <- function(plots_list, resolution_name, loop_count, cell_line) {
  if (length(plots_list) == 0) {
    cat("No plots available for", resolution_name, "\n")
    return(NULL)
  }

  cat("Creating", resolution_name, "plot with log2 transformation...\n")

  plot_grid <- plots_list
  while (length(plot_grid) < 8) {
    plot_grid <- c(plot_grid, list(ggplot() + theme_void()))
  }
  plot_grid <- plot_grid[1:8]

  temp_plot <- plots_list[[1]] + theme(legend.position = "right")
  shared_legend <- cowplot::get_legend(temp_plot)

  plot_grid <- lapply(plot_grid, function(p) p + theme(legend.position = "none"))

  grid_2x4 <- do.call(grid.arrange, c(plot_grid, ncol = 4, nrow = 2))

  section_header <- create_section_header(resolution_name, loop_count, cell_line)

  final_plot <- grid.arrange(
    grid.arrange(
      textGrob(section_header, gp = gpar(fontsize = 18, fontface = "bold")),
      grid_2x4,
      ncol = 1,
      heights = c(0.08, 0.92)
    ),
    shared_legend,
    ncol = 2,
    widths = c(0.88, 0.12)
  )

  filename <- paste0("APA_log2_", resolution_name, ".pdf")
  ggsave(filename, final_plot, width = 16, height = 10, dpi = 300, device = "pdf", limitsize = FALSE)
  cat("SUCCESS: Saved", filename, "\n")

  return(final_plot)
}

# Main analysis starts here
apa_dirs <- list.dirs(".", recursive = FALSE, full.names = FALSE)
apa_dirs <- apa_dirs[grepl("_APA_Results$", apa_dirs)]

if (length(apa_dirs) == 0) {
  cat("ERROR: No APA result directories found!\n")
  quit(status=1)
}

results <- data.frame(
  Method = character(),
  Resolution = character(),
  P2LL_Score = numeric(),
  Central_R = numeric(),
  Loop_Count = numeric(),
  stringsAsFactors = FALSE
)

all_plot_data <- list()

unified_loop_count <- NA
for (dir in apa_dirs) {
  raw_method_name <- gsub("_APA_Results$", "", dir)
  method_name <- extract_method_name(raw_method_name)
  loop_count <- get_loop_count(raw_method_name)
  if (!is.na(loop_count)) {
    unified_loop_count <- loop_count
    break
  }
}

cat("=== FIRST PASS: Loading all matrices and calculating global parameters ===\n")
for (dir in apa_dirs) {
  raw_method_name <- gsub("_APA_Results$", "", dir)
  method_name <- extract_method_name(raw_method_name)
  resolutions <- c("5000", "10000", "25000")

  for (resolution in resolutions) {
    apa_file <- file.path(dir, resolution, "gw", "APA.txt")
    if (file.exists(apa_file)) {
      apa_matrix <- read_apa_matrix(apa_file)
      if (!is.null(apa_matrix)) {
        res_label <- ifelse(resolution == "5000", "5kb", 
                           ifelse(resolution == "10000", "10kb", "25kb"))
        matrix_key <- paste0(method_name, "_", res_label)
        all_raw_matrices[[matrix_key]] <- apa_matrix
        cat("Loaded matrix for", matrix_key, "\n")
      }
    }
  }
}

cat("Total matrices loaded:", length(all_raw_matrices), "\n")

calculate_global_parameters()
cat("Global log2 parameters fixed - all plots will use identical color scales\n\n")

cat("=== SECOND PASS: Processing APA analysis ===\n")
for (dir in apa_dirs) {
  raw_method_name <- gsub("_APA_Results$", "", dir)
  method_name <- extract_method_name(raw_method_name)
  resolutions <- c("5000", "10000", "25000")

  cat("Processing method:", method_name, "\n")

  for (resolution in resolutions) {
    apa_file <- file.path(dir, resolution, "gw", "APA.txt")
    if (!file.exists(apa_file)) {
      cat("  APA.txt not found for", resolution, "\n")
      next
    }

    apa_matrix <- read_apa_matrix(apa_file)
    if (is.null(apa_matrix)) {
      cat("  Failed to read APA matrix for", resolution, "\n")
      next
    }

    apa_analysis <- calculate_apa_regions(apa_matrix)
    if (is.null(apa_analysis)) {
      cat("  Failed to analyze APA data for", resolution, "\n")
      next
    }

    res_label <- ifelse(resolution == "5000", "5kb", 
                       ifelse(resolution == "10000", "10kb", "25kb"))

    results <- rbind(results, data.frame(
      Method = method_name,
      Resolution = res_label,
      P2LL_Score = apa_analysis$p2ll,
      Central_R = apa_analysis$apa_scores$central_ratio,
      Loop_Count = ifelse(is.na(unified_loop_count), 0, unified_loop_count),
      stringsAsFactors = FALSE
    ))

    all_plot_data[[paste0(method_name, "_", res_label)]] <- list(
      matrix_data = apa_matrix,
      method_name = method_name,
      resolution = res_label,
      apa_analysis = apa_analysis
    )

    cat("   SUCCESS:", res_label, "| P2LL:", round(apa_analysis$p2ll, 3),
        "| R:", round(apa_analysis$apa_scores$central_ratio, 3), "\n")
  }
}

if (nrow(results) == 0) {
  cat("ERROR: No results generated!\n")
  quit(status=1)
}

results_5kb <- results[results$Resolution == "5kb", ]
results_10kb <- results[results$Resolution == "10kb", ]
results_25kb <- results[results$Resolution == "25kb", ]
results_5kb <- results_5kb[order(results_5kb$P2LL_Score, decreasing = TRUE), ]
results_10kb <- results_10kb[order(results_10kb$P2LL_Score, decreasing = TRUE), ]
results_25kb <- results_25kb[order(results_25kb$P2LL_Score, decreasing = TRUE), ]

cat("\n=== Generating plots with log2 transformation ===\n")

plots_5kb <- list()
plots_10kb <- list()
plots_25kb <- list()

for (i in 1:nrow(results_5kb)) {
  method_name <- results_5kb$Method[i]
  plot_data <- all_plot_data[[paste0(method_name, "_5kb")]]
  if (!is.null(plot_data)) {
    plot_title <- paste0(method_name, "\nAPA: ", round(plot_data$apa_analysis$p2ll, 2),
                         " | R: ", round(plot_data$apa_analysis$apa_scores$central_ratio, 2))

    p <- create_apa_plot(plot_data$matrix_data, method_name, "5kb", plot_data$apa_analysis) +
      labs(title = plot_title)
    plots_5kb[[i]] <- p
  }
}

for (i in 1:nrow(results_10kb)) {
  method_name <- results_10kb$Method[i]
  plot_data <- all_plot_data[[paste0(method_name, "_10kb")]]
  if (!is.null(plot_data)) {
    plot_title <- paste0(method_name, "\nAPA: ", round(plot_data$apa_analysis$p2ll, 2),
                         " | R: ", round(plot_data$apa_analysis$apa_scores$central_ratio, 2))

    p <- create_apa_plot(plot_data$matrix_data, method_name, "10kb", plot_data$apa_analysis) +
      labs(title = plot_title)
    plots_10kb[[i]] <- p
  }
}

for (i in 1:nrow(results_25kb)) {
  method_name <- results_25kb$Method[i]
  plot_data <- all_plot_data[[paste0(method_name, "_25kb")]]
  if (!is.null(plot_data)) {
    plot_title <- paste0(method_name, "\nAPA: ", round(plot_data$apa_analysis$p2ll, 2),
                         " | R: ", round(plot_data$apa_analysis$apa_scores$central_ratio, 2))

    p <- create_apa_plot(plot_data$matrix_data, method_name, "25kb", plot_data$apa_analysis) +
      labs(title = plot_title)
    plots_25kb[[i]] <- p
  }
}

create_plot_grid(plots_5kb, "5kb", unified_loop_count, cellLine)
create_plot_grid(plots_10kb, "10kb", unified_loop_count, cellLine)
create_plot_grid(plots_25kb, "25kb", unified_loop_count, cellLine)

results_filename <- "APA_Results.txt"
write.table(results, results_filename, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Creating method comparison summary ===\n")
summary_comparison <- results_5kb[, c("Method", "P2LL_Score", "Central_R")]
summary_comparison$Rank <- 1:nrow(summary_comparison)
summary_comparison$Performance_Level <- ifelse(summary_comparison$P2LL_Score >= 1.5, "Excellent",
                                              ifelse(summary_comparison$P2LL_Score >= 1.2, "Good",
                                                    ifelse(summary_comparison$P2LL_Score >= 1.0, "Fair", "Poor")))

summary_comparison$sintHiChIP_Method <- grepl("sintHiChIP", summary_comparison$Method, ignore.case = TRUE)

summary_filename <- "APA_Method_Ranking_5kb.txt"
write.table(summary_comparison, summary_filename, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n", rep("=", 80), "\n")
cat("=== APA ANALYSIS COMPLETE ===\n")
cat(rep("=", 80), "\n")

cat("\nGenerated files:\n")
cat("- APA_log2_5kb.pdf (log2 transformation, 5kb resolution)\n")
cat("- APA_log2_10kb.pdf (log2 transformation, 10kb resolution)\n")
cat("- APA_log2_25kb.pdf (log2 transformation, 25kb resolution)\n")
cat("- ", results_filename, " (comprehensive results)\n")
cat("- ", summary_filename, " (method ranking)\n")

cat("\nMethodology:\n")
cat("Biological P2LL and R scores calculated from RAW COUNTS\n")
cat("Visualization uses log2(value + 1) transformation\n")
cat("Global normalization to [0,1] using ", round(opt$percentileCut*100, 1), "%-", round((1-opt$percentileCut)*100, 1), "% percentiles\n")
cat("Unified color scaling ensures fair visual comparison\n")
cat("sintHiChIP methods highlighted with asterisks and gold borders\n")

cat("\nSummary Statistics:\n")
cat("- Cell line used:", cellLine, "\n")
cat("- Total loop count:", ifelse(is.na(unified_loop_count), "unavailable", unified_loop_count), "\n")
cat("- Methods analyzed (5kb):", nrow(results_5kb), "\n")
cat("- Methods analyzed (10kb):", nrow(results_10kb), "\n")
cat("- Methods analyzed (25kb):", nrow(results_25kb), "\n")

cat("\nTop 3 Performers (5kb resolution):\n")
for (i in 1:min(3, nrow(summary_comparison))) {
  method_info <- summary_comparison[i, ]
  star <- ifelse(method_info$sintHiChIP_Method, "* ", "")
  cat(sprintf("%d. %s%s - APA: %.3f, R: %.3f (%s)\n", 
              i, star, method_info$Method, method_info$P2LL_Score, 
              method_info$Central_R, method_info$Performance_Level))
}

cat(rep("=", 80), "\n")