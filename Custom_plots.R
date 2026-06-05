#!/usr/bin/env Rscript

######## Packages ##########

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(require(gridExtra))
suppressMessages(library(ggtext))
suppressMessages(library(tidyr))
suppressMessages(library(ggbeeswarm))

# Suppress warnings globally
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat("\nUsage:\n")
  cat("  Rscript Custom_plots.R -y <ylimit> -t <Y|N> [-i input_dir] [-m <peaks|bckg>] [-c <HEX1,HEX2,HEX3>] [-d <Y|N>]\n\n")
  cat("Parameters:\n")
  cat("  -i   Input directory (optional, default = current working directory)\n")
  cat("  -y   Numeric value for the y-axis limit (required)\n")
  cat("  -t   Apply signed sqrt transformation: Y or N (required)\n")
  cat("  -m   Mode: peaks (default) or bckg, case-insensitive (optional)\n")
  cat("  -c   3 comma-separated HEX colors for Positives,Negatives,Targets (optional)\n")
  cat("       Example: #FFA500,#B0B0B0,#FF0000\n")
  cat("  -d   Show dots in the plot: Y (default) or N (optional)\n\n")
  cat("Examples:\n")
  cat("  Rscript Custom_plots.R -y 10 -t Y\n")
  cat("  Rscript Custom_plots.R i /path/to/data -y 10 -t N -m bckg\n")
  cat("  Rscript Custom_plots.R -y 10 -t Y -c \"#FFA500,#B0B0B0,#FF0000\"\n")
  cat("  Rscript Custom_plots.R i /path/to/data -y 10 -t Y -m peaks -d N\n")
  cat("  Rscript Custom_plots.R -y 10 -t Y -c \"#FFA500,#B0B0B0,#FF0000\" -m bckg -d Y\n\n")
  cat("Output file naming convention:\n")
  cat("  outfile_<type>_<params>_<suffix>.<ext>\n\n")
  cat("  where <params> is a compact summary of selected options:\n")
  cat("    trX        : transform flag (trY or trN)\n")
  cat("    ylN        : y-axis limit (e.g., yl10)\n")
  cat("    peaks/bckg : normalization mode\n")
  cat("    dtX        : dots flag (dtY or dtN)\n")
  cat("    colors     :\n")
  cat("                 - cldefault → default colors\n")
  cat("                 - otherwise HEX codes joined with '-' (without #)\n\n")
  
  cat("Example:\n")
  cat("  outfile_peaks_trY_yl10_peaks_dtN_FFA500-B0B0B0-FF0000_sampleA.pdf\n\n")
}

fail <- function(msg) {
  cat("Error:", msg, "\n\n")
  print_help()
  quit(status = 1)
}

# Optional --help / -h
if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  print_help()
  quit(status = 0)
}

# All flags must be followed by a value, so total number of args must be even
if (length(args) %% 2 != 0) {
  
  # Detect likely unquoted HEX input
  if (any(args == "-c")) {
    fail("Parameter -c requires a value. If using HEX colors, wrap them in quotes (example: -c \"#FFA500,#B0B0B0,#FF0000\").")
  }
  
  fail("Each parameter flag must be followed by a value (example: -y 10).")
}

# Allowed flags
allowed_flags <- c("-y", "-t", "-m", "-c", "-d", "-i")

# Check that odd-positioned arguments are valid flags
input_flags <- args[seq(1, length(args), by = 2)]
input_values <- args[seq(2, length(args), by = 2)]

if (!all(input_flags %in% allowed_flags)) {
  bad_flags <- unique(input_flags[!(input_flags %in% allowed_flags)])
  fail(paste0(
    "Unknown parameter(s): ",
    paste(bad_flags, collapse = ", "),
    " (allowed: -y, -t, -m, -c, -d, -i)."
  ))
}

# Check duplicates
dup_flags <- unique(input_flags[duplicated(input_flags)])
if (length(dup_flags) > 0) {
  fail(paste0(
    "Duplicated parameter(s): ",
    paste(dup_flags, collapse = ", "),
    ". Each parameter must be provided only once."
  ))
}

# Convert to named list
param_list <- as.list(input_values)
names(param_list) <- input_flags

# Required params
if (is.null(param_list[["-y"]])) {
  fail("Missing required parameter -y (example: -y 10).")
}
if (is.null(param_list[["-t"]])) {
  fail("Missing required parameter -t (example: -t Y).")
}

# Parse -y
ylimit <- suppressWarnings(as.numeric(param_list[["-y"]]))
if (is.na(ylimit)) {
  fail("Parameter -y must be numeric (example: -y 10).")
}

# Parse -t
transform_flag <- toupper(param_list[["-t"]])
if (!(transform_flag %in% c("Y", "N"))) {
  fail("Parameter -t must be 'Y' or 'N' (example: -t Y).")
}

# Default colors
group_colors <- c("Positives" = "orange", "Negatives" = "grey", "Targets" = "red")

# Parse -c
if (!is.null(param_list[["-c"]])) {
  color_input <- param_list[["-c"]]
  color_values <- strsplit(color_input, ",", fixed = TRUE)[[1]]
  
  # Trim possible surrounding spaces
  color_values <- trimws(color_values)
  
  # Add # if missing
  color_values <- ifelse(grepl("^#", color_values), color_values, paste0("#", color_values))
  
  if (length(color_values) != 3) {
    fail("Parameter -c must contain exactly 3 comma-separated HEX colors (example: -c \"#FFA500,#B0B0B0,#FF0000\" or -c \"FFA500,B0B0B0,FF0000\").")
  }
  
  if (!all(grepl("^#([A-Fa-f0-9]{6})$", color_values))) {
    fail("Colors in -c must be valid 6-digit HEX codes (example: -c \"#FFA500,#B0B0B0,#FF0000\" or -c \"FFA500,B0B0B0,FF0000\").")
  }
  
  # Normalize to uppercase for consistency
  color_values <- toupper(color_values)
  
  group_colors <- setNames(color_values, c("Positives", "Negatives", "Targets"))
}

# Parse -m
mode <- "peaks"
if (!is.null(param_list[["-m"]])) {
  mode_input <- tolower(param_list[["-m"]])
  
  if (mode_input == "peaks") {
    mode <- "peaks"
  } else if (mode_input == "bckg") {
    mode <- "bckg"
  } else {
    fail("Parameter -m must be 'peaks' or 'bckg' (example: -m bckg).")
  }
}

# Parse -d
dots_flag <- "Y"
if (!is.null(param_list[["-d"]])) {
  dots_flag <- toupper(param_list[["-d"]])
  
  if (!(dots_flag %in% c("Y", "N"))) {
    fail("Parameter -d must be 'Y' or 'N' (example: -d N).")
  }
}

if (is.null(param_list[["-i"]])) {
  current_dir <- normalizePath(getwd())
  cat("No -i provided, using current working directory:", current_dir, "\n")
} else {
  current_dir <- param_list[["-i"]]
}

# Expand ~
current_dir <- path.expand(current_dir)

# Remove trailing slash
current_dir <- sub("/+$", "", current_dir)

# Check existence
if (!dir.exists(current_dir)) {
  fail(paste0("Input directory does not exist: ", current_dir))
}

# Check readability
if (file.access(current_dir, 4) != 0) {
  fail(paste0("Input directory is not readable: ", current_dir))
}

# Normalize path (absolute path, resolves symlinks)
current_dir <- normalizePath(current_dir)

# Debug / confirmation output
cat("Setting y-axis limit to:", ylimit, "\n")
cat("Apply signed sqrt transformation:", transform_flag, "\n")
cat("Mode:", mode, "\n")
cat("Show dots:", dots_flag, "\n")
cat("Group colors:\n")
print(group_colors)
cat("Using input directory:", current_dir, "\n")

########################################
##  INPUT FILE DISCOVERY AND PARSING  ##
########################################

input_dir <- current_dir
output_dir <- current_dir

lambda_flag <- "0"
suffix <- ""

fail <- function(msg) {
  stop(msg, call. = FALSE)
}

find_unique_file <- function(path, pattern, label, exclude_pattern = NULL) {
  
  matches <- list.files(path = path, pattern = pattern, full.names = TRUE)
  
  if (!is.null(exclude_pattern)) {
    matches <- matches[!grepl(exclude_pattern, basename(matches))]
  }
  
  # No match → stop
  if (length(matches) == 0) {
    stop(
      paste0(
        "No file found for ", label, " in directory:\n  ", path,
        "\nPattern used: ", pattern
      ),
      call. = FALSE
    )
  }
  
  # Multiple matches → stop with detailed message
  if (length(matches) > 1) {
    stop(
      paste0(
        "Ambiguous input for ", label, ": multiple files match.\n",
        "This likely indicates inconsistent or confusing file naming.\n\n",
        "Matching files:\n  - ",
        paste(basename(matches), collapse = "\n  - "),
        "\n\nPlease ensure only one file matches pattern:\n  ", pattern
      ),
      call. = FALSE
    )
  }
  
  cat("Found", label, "file:", matches[1], "\n")
  return(matches[1])
}

extract_suffix_from_top_peaks <- function(file_name, prefix) {
  sub(paste0("^", prefix, "(.*)\\.tsv$"), "\\1", basename(file_name))
}

########################################
##  TOP POSITIVES PEAKS               ##
########################################

top_positives_peaks_file <- find_unique_file(
  path = input_dir,
  pattern = "^top_positives_peaks(.*)\\.tsv$",
  label = "top positives peaks"
)

suffix_raw <- extract_suffix_from_top_peaks(
  file_name = top_positives_peaks_file,
  prefix = "top_positives_peaks"
)

cat("Detected suffix:", suffix_raw, "\n")

########################################
##  TOP NEGATIVES PEAKS               ##
########################################

top_negatives_peaks_file <- find_unique_file(
  path = input_dir,
  pattern = "^top_negatives_peaks(.*)\\.tsv$",
  label = "top negatives peaks"
)

########################################
##  TARGET WINDOWS FILE               ##
########################################

target_windows_pattern <- "^(filtered_)?[^[:space:]]+_win_count(_lambda_corrected)?(.*?)\\.tsv$"

target_windows_file <- find_unique_file(
  path = input_dir,
  pattern = target_windows_pattern,
  label = "target windows",
  exclude_pattern = "^(filtered_)?(positives|negatives)_"
)

target_windows_name <- basename(target_windows_file)

target_match_info <- regexec(target_windows_pattern, target_windows_name)
target_match_groups <- regmatches(target_windows_name, target_match_info)[[1]]

if (length(target_match_groups) == 0) {
  fail(paste0("Could not parse target windows filename: ", target_windows_name))
}

target_filtered_prefix <- if (length(target_match_groups) >= 2) target_match_groups[2] else ""
target_lambda_part     <- if (length(target_match_groups) >= 3) target_match_groups[3] else ""
target_extra_suffix    <- if (length(target_match_groups) >= 4) target_match_groups[4] else ""

if (!is.na(target_lambda_part) && target_lambda_part == "_lambda_corrected") {
  lambda_flag <- "1"
}

cat("Target windows filtered prefix:", target_filtered_prefix, "\n")
cat("Target windows lambda part:", target_lambda_part, "\n")
cat("Target windows extra suffix:", target_extra_suffix, "\n")

########################################
##  POSITIVES WINDOWS FILE            ##
########################################

positives_windows_file <- find_unique_file(
  path = input_dir,
  pattern = "^positives_win_count(_lambda_corrected_filtered)?(.*)?\\.tsv$",
  label = "positives windows"
)

########################################
##  NEGATIVES WINDOWS FILE            ##
########################################

negatives_windows_file <- find_unique_file(
  path = input_dir,
  pattern = "^negatives_win_count(_lambda_corrected_filtered)?(.*)?\\.tsv$",
  label = "negatives windows"
)

########################################
##  READ INPUT TABLES                 ##
########################################

pos_peaks <- read.table(top_positives_peaks_file, header = FALSE, sep = "\t", quote = "\"")
neg_peaks <- read.table(top_negatives_peaks_file, header = FALSE, sep = "\t", quote = "\"")
target    <- read.table(target_windows_file, header = FALSE, sep = "\t", quote = "\"")
pos_win   <- read.table(positives_windows_file, header = FALSE, sep = "\t", quote = "\"")
neg_win   <- read.table(negatives_windows_file, header = FALSE, sep = "\t", quote = "\"")

########################################
##  FINAL USER VARIABLES              ##
########################################

outdir <- output_dir
lambda <- lambda_flag

suffix <- suffix_raw
suffix <- sub("^_", "", suffix)

cat("Output directory:", outdir, "\n")
cat("Lambda corrected:", lambda, "\n")
cat("Final suffix:", suffix, "\n")

#######################
##  CALCULATE RATIO  ##
#######################
ratio_pos_peaks <- pos_peaks[,9] / pos_peaks[,8]

ratio_neg_peaks <- neg_peaks[,9] / neg_peaks[,8]

if(lambda == 0){
  target <- target[target[,8] != 0, ]
  ratio_target <- target[,9] / target[,8]
  
  pos_win <- pos_win[pos_win[,8] != 0, ]
  ratio_pos_win <- pos_win[,9] / pos_win[,8]
  
  neg_win <- neg_win[neg_win[,8] != 0, ]
  ratio_neg_win <- neg_win[,9] / neg_win[,8]
  
} else {
  target <- target[target[,11] != 0, ]
  ratio_target <- target[,9] / target[,11]
  
  pos_win <- pos_win[pos_win[,11] != 0, ]
  ratio_pos_win <- pos_win[,9] / pos_win[,11]
  
  neg_win <- neg_win[neg_win[,11] != 0, ]
  ratio_neg_win <- neg_win[,9] / neg_win[,11]
}

########################
##  CALCULATE MEDIAN  ##
########################
median_neg_peaks <- median(ratio_neg_peaks)
median_neg_win   <- median(ratio_neg_win)

median_pos_peaks <- median(ratio_pos_peaks)
median_pos_win   <- median(ratio_pos_win)

########################################
##  SELECT REFERENCE BASED ON MODE    ##
########################################
if (mode == "peaks") {
  ref_median_peaks <- median_pos_peaks
  ref_median_win   <- median_pos_win
} else if (mode == "bckg") {
  ref_median_peaks <- median_neg_peaks
  ref_median_win   <- median_neg_win
} else {
  stop("Invalid mode: must be 'peaks' or 'bckg'")
}

##############################
##  CALCULATE HOST POS&NEG  ##
##############################
pos_distribution      <- ratio_pos_peaks / ref_median_peaks
neg_distribution      <- ratio_neg_peaks / ref_median_peaks

pos_distribution_win  <- ratio_pos_win / ref_median_win
neg_distribution_win  <- ratio_neg_win / ref_median_win

########################
##  CALCULATE TARGET  ##
########################
target_distribution      <- ratio_target / ref_median_peaks
target_distribution_win  <- ratio_target / ref_median_win


######################
##  BUILD DATA PLOT ##
######################
pos_peaks_data <- data.frame(group = "Positives", value = pos_distribution)
neg_peaks_data <- data.frame(group = "Negatives", value = neg_distribution)
target_data <- data.frame(group = "Targets", value = target_distribution)

plot.data_peaks <- rbind(pos_peaks_data, neg_peaks_data, target_data)
stat_box_data_peaks <- function(y, upper_limit = max(plot.data_peaks$value) * 1.3) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'mean =', round(mean(y), 2), '\n',
                    'median =', round(median(y), 2), '\n')
    )
  )
}

pos_win_data <- data.frame(group = "Positives", value = pos_distribution_win)
neg_win_data <- data.frame(group = "Negatives", value = neg_distribution_win)
target_win_data <- data.frame(group = "Targets", value = target_distribution_win)

plot.data_win <- rbind(pos_win_data, neg_win_data, target_win_data)
stat_box_data_win <- function(y, upper_limit = max(plot.data_win$value) * 1.3) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'mean =', round(mean(y), 2), '\n',
                    'median =', round(median(y), 2), '\n')
    )
  )
}

############
##  PLOT  ##
############
create_chip_signal_plot <- function(plot.data, plot_title, x_label) {
  
  # Expected globals:
  # - transform_flag : "Y" or "N"
  # - dots_flag      : "Y" or "N"
  # - ylimit         : numeric
  # - group_colors   : named vector with Positives / Negatives / Targets
  
  required_groups <- c("Positives", "Negatives", "Targets")
  
  # Check required columns
  if (!all(c("group", "value") %in% colnames(plot.data))) {
    stop("plot.data must contain at least the columns 'group' and 'value'.")
  }
  
  # Force group order
  plot.data <- plot.data %>%
    dplyr::filter(!is.na(group), !is.na(value)) %>%
    dplyr::mutate(
      group = factor(as.character(group), levels = required_groups)
    ) %>%
    dplyr::filter(!is.na(group))
  
  if (nrow(plot.data) == 0) {
    stop("plot.data is empty after filtering valid groups and non-missing values.")
  }
  
  # Validate flags
  transform_flag_local <- toupper(transform_flag)
  dots_flag_local <- toupper(dots_flag)
  
  if (!(transform_flag_local %in% c("Y", "N"))) {
    stop("transform_flag must be 'Y' or 'N'.")
  }
  if (!(dots_flag_local %in% c("Y", "N"))) {
    stop("dots_flag must be 'Y' or 'N'.")
  }
  
  # Validate ylimit
  if (!is.numeric(ylimit) || length(ylimit) != 1 || is.na(ylimit) || ylimit <= 0) {
    stop("ylimit must be a single numeric value > 0.")
  }
  
  # Apply optional transformation
  if (transform_flag_local == "Y") {
    plot.data <- plot.data %>%
      dplyr::mutate(value = sign(value) * sqrt(abs(value)))
    ylegend <- "Signed Sqrt Transformed Relative Enrichment"
    transformtag <- "sqrt_transformed"
  } else {
    message("Skipping signed sqrt transformation")
    ylegend <- "Relative Enrichment"
    transformtag <- "raw"
  }
  
  # Compute stats
  stats <- plot.data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      count = dplyr::n(),
      .groups = "drop"
    )
  
  # Ensure every group exists in stats, even if absent from data
  stats <- dplyr::left_join(
    tibble::tibble(group = factor(required_groups, levels = required_groups)),
    stats,
    by = "group"
  ) %>%
    dplyr::mutate(count = ifelse(is.na(count), 0, count))
  
  # Optional: keep this only if you still need it later
  cap_threshold <- plot.data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(threshold = stats::quantile(value, 0.999, na.rm = TRUE), .groups = "drop") %>%
    dplyr::summarise(max_threshold = max(threshold, na.rm = TRUE)) %>%
    dplyr::pull(max_threshold)
  
  # Custom y-axis labels
  custom_left_labels <- function(y) {
    pretty_breaks <- pretty(c(0, ylimit))
    
    sapply(y, function(val) {
      if (is.na(val)) {
        return("")
      }
      
      if (val %in% pretty_breaks || isTRUE(all.equal(val, ylimit))) {
        return(as.character(signif(val, 4)))
      }
      
      med_idx <- which(abs(stats$median - val) < .Machine$double.eps^0.5)
      med_idx <- med_idx[!is.na(med_idx)]
      
      if (length(med_idx) > 0) {
        grp <- as.character(stats$group[med_idx[1]])
        return(
          paste0(
            "<b><span style='color:",
            group_colors[grp],
            "; font-size:15px;'>",
            round(val, 2),
            "</span></b>"
          )
        )
      }
      
      ""
    })
  }
  
  # Counts for x-axis labels
  n_pos <- stats$count[stats$group == "Positives"]
  n_neg <- stats$count[stats$group == "Negatives"]
  n_tar <- stats$count[stats$group == "Targets"]
  
  xlabels_expr <- c(
    "Positives" = bquote(atop("Positives", bold(n == .(scales::comma(n_pos))))),
    "Negatives" = bquote(atop("Background", bold(n == .(scales::comma(n_neg))))),
    "Targets"   = bquote(atop("Targets", bold(n == .(scales::comma(n_tar)))))
  )
  
  # Base plot
  plot_peaks <- ggplot(
    plot.data,
    aes(
      x = factor(group, levels = required_groups),
      y = value,
      fill = group
    )
  ) +
    geom_violin(
      trim = FALSE,
      alpha = 0.6,
      scale = if (dots_flag_local == "Y") "width" else "area"
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_x_discrete(
      limits = required_groups,
      labels = xlabels_expr
    ) +
    scale_y_continuous(
      name = ylegend,
      breaks = function(limits) {
        base_breaks <- pretty(c(0, ylimit))
        medians <- stats$median[!is.na(stats$median)]
        unique(c(0, base_breaks, medians, ylimit))
      },
      labels = custom_left_labels,
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      axis.title.y.left = element_text(size = 16, face = "bold"),
      axis.title.y.right = element_text(size = 16, face = "bold"),
      axis.text.y.left = ggtext::element_markdown(size = 12),
      axis.text.y.right = ggtext::element_markdown(size = 12),
      axis.text.x = element_text(size = 12),
      legend.position = "none",
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.4),
      axis.line.y = element_line(color = "black", linewidth = 0.4),
      axis.ticks.x = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y = element_line(color = "black", linewidth = 0.4)
    ) +
    labs(
      title = plot_title,
      x = x_label,
      y = NULL,
      caption = "Dashed lines represent the median signal."
    ) +
    geom_hline(
      data = stats %>% dplyr::filter(!is.na(median)),
      aes(yintercept = median, color = group),
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    coord_cartesian(ylim = c(0, ylimit))
  
  # Add dots only if requested
  if (dots_flag_local == "Y") {
    plot_peaks <- plot_peaks +
      ggbeeswarm::geom_quasirandom(alpha = 0.6, size = 1.5)
  }
  
  return(plot_peaks)
}

############################
##  BUILD COMPACT TAGS    ##
############################

# transform flag
transformtag <- paste0("tr", toupper(transform_flag))

# y-limit
ylimit_tag <- paste0("yl", format(ylimit, scientific = FALSE, trim = TRUE))

# mode (now explicit)
mode_tag <- tolower(mode)  # "peaks" or "bckg"

# dots
dots_tag <- paste0("dt", toupper(dots_flag))

# colors
default_group_colors <- c("Positives" = "orange", "Negatives" = "grey", "Targets" = "red")

if (identical(unname(group_colors), unname(default_group_colors))) {
  color_tag <- "cldefault"
} else {
  # remove "#" and join
  clean_hex <- gsub("#", "", group_colors)
  color_tag <- paste(clean_hex, collapse = "-")
}

# final param tag
param_tag <- paste(transformtag, ylimit_tag, mode_tag, dots_tag, color_tag, sep = "_")

############################
##  OPTIONAL: SAFE SUFFIX ##
############################
safe_suffix <- gsub("[^A-Za-z0-9._-]+", "_", suffix)

############################
##  OUTPUT FILE NAMES     ##
############################

# Plot files
out_peaks_name <- paste0("outfile_peaks_", param_tag, "_", safe_suffix, ".pdf")
out_peaks <- file.path(outdir, out_peaks_name)

out_win_name <- paste0("outfile_win_", param_tag, "_", safe_suffix, ".pdf")
out_win <- file.path(outdir, out_win_name)

# Data files
out_data_peaks_name <- paste0("outfile_peaks_", param_tag, "_", safe_suffix, ".tsv")
out_data_peaks <- file.path(outdir, out_data_peaks_name)

out_data_win_name <- paste0("outfile_win_", param_tag, "_", safe_suffix, ".tsv")
out_data_win <- file.path(outdir, out_data_win_name)

############################
##  PLOT TITLES           ##
############################
out_title_peaks <- paste("ChIP Signal Distribution Across Genome Peaks", suffix)
out_title_win   <- paste("ChIP Signal Distribution Across Genome Bins", suffix)

subset_for_plot_if_needed <- function(plot.data, dots_flag, max_points = 50000, seed = 1, plot_label = "plot") {
  
  dots_flag <- toupper(dots_flag)
  
  if (!(dots_flag %in% c("Y", "N"))) {
    stop("dots_flag must be 'Y' or 'N'.")
  }
  
  n_total <- nrow(plot.data)
  
  if (dots_flag == "N") {
    return(plot.data)
  }
  
  if (n_total <= max_points) {
    return(plot.data)
  }
  
  if (!("group" %in% colnames(plot.data))) {
    stop("plot.data must contain a 'group' column for proportional subsetting.")
  }
  
  set.seed(seed)
  
  group_counts <- table(plot.data$group)
  group_props <- group_counts / sum(group_counts)
  group_targets <- pmax(1, round(max_points * group_props))
  
  diff_n <- max_points - sum(group_targets)
  if (diff_n != 0) {
    ord <- order(group_targets, decreasing = (diff_n > 0))
    for (i in seq_len(abs(diff_n))) {
      j <- ord[((i - 1) %% length(ord)) + 1]
      group_targets[j] <- group_targets[j] + sign(diff_n)
    }
  }
  
  split_df <- split(plot.data, plot.data$group)
  
  sampled_list <- mapply(
    FUN = function(df, n_keep) {
      n_keep <- min(nrow(df), n_keep)
      df[sample.int(nrow(df), size = n_keep, replace = FALSE), , drop = FALSE]
    },
    df = split_df,
    n_keep = as.numeric(group_targets[names(split_df)]),
    SIMPLIFY = FALSE
  )
  
  plot_subset <- do.call(rbind, sampled_list)
  rownames(plot_subset) <- NULL
  
  subset_fraction <- nrow(plot_subset) / n_total
  
  warning(
    paste0(
      "Too many points to draw for ", plot_label, ". ",
      "Using a random plotting subset of ",
      format(nrow(plot_subset), big.mark = ","),
      " / ",
      format(n_total, big.mark = ","),
      " rows (",
      sprintf("%.2f", 100 * subset_fraction),
      "% of total data). ",
      "The full dataset is still written to output tables."
    ),
    call. = FALSE
  )
  
  plot_subset
}

max_plot_points <- 50000

############################
##  SAVE PEAKS PLOT       ##
############################
plot.data_peaks_for_plot <- subset_for_plot_if_needed(
  plot.data = plot.data_peaks,
  dots_flag = dots_flag,
  max_points = max_plot_points,
  seed = 1,
  plot_label = "peaks plot"
)

g_peaks <- create_chip_signal_plot(plot.data_peaks_for_plot, out_title_peaks, "Genome Regions")

ggsave(
  filename = out_peaks,
  plot = g_peaks,
  width = 35,
  height = 20,
  units = "cm",
  device = grDevices::cairo_pdf
)

############################
##  SAVE WINDOWS PLOT     ##
############################
plot.data_win_for_plot <- subset_for_plot_if_needed(
  plot.data = plot.data_win,
  dots_flag = dots_flag,
  max_points = max_plot_points,
  seed = 1,
  plot_label = "windows plot"
)

g_win <- create_chip_signal_plot(plot.data_win_for_plot, out_title_win, "Genome Bins")

ggsave(
  filename = out_win,
  plot = g_win,
  width = 35,
  height = 20,
  units = "cm",
  device = grDevices::cairo_pdf
)


############################
##  WRITE PEAKS DATA      ##
############################
myls <- list(
  as.data.frame(pos_peaks_data),
  as.data.frame(neg_peaks_data),
  as.data.frame(target_data)
)

max.rows <- max(
  nrow(as.data.frame(pos_peaks_data)),
  nrow(as.data.frame(neg_peaks_data)),
  nrow(as.data.frame(target_data))
)

new_myls <- lapply(myls, function(x) x[seq_len(max.rows), , drop = FALSE])
value_columns <- lapply(new_myls, function(x) x$value)
new_myls <- as.data.frame(do.call(cbind, value_columns))
colnames(new_myls) <- c("Positives", "Negatives", "Targets")

write.table(new_myls, file = out_data_peaks, sep = "\t", row.names = FALSE, quote = FALSE)

############################
##  WRITE WINDOWS DATA    ##
############################
myls <- list(
  as.data.frame(pos_win_data),
  as.data.frame(neg_win_data),
  as.data.frame(target_win_data)
)

max.rows <- max(
  nrow(as.data.frame(pos_win_data)),
  nrow(as.data.frame(neg_win_data)),
  nrow(as.data.frame(target_win_data))
)

new_myls <- lapply(myls, function(x) x[seq_len(max.rows), , drop = FALSE])
value_columns <- lapply(new_myls, function(x) x$value)
new_myls <- as.data.frame(do.call(cbind, value_columns))
colnames(new_myls) <- c("Positives", "Negatives", "Targets")

write.table(new_myls, file = out_data_win, sep = "\t", row.names = FALSE, quote = FALSE)
