#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggtext))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(tibble))
suppressMessages(library(scales))

args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat("\nUsage:\n")
  cat("  Standalone mode (automatic file discovery):\n")
  cat("    Rscript Plots_ViPeaQ.R [-i <input_dir>] [-t <Y|N>] [-y <ylimit>] [-m <peaks|bckg>] [-c <HEX1,HEX2,HEX3>] [-d <Y|N>]\n\n")
  
  cat("  Explicit-input mode:\n")
  cat("    Rscript Plots_ViPeaQ.R [-t <Y|N>] [-y <ylimit>] -P <top_pos.tsv> -N <top_neg.tsv> -T <target.tsv> -A <pos_win.tsv> -B <neg_win.tsv> -o <outdir> -l <0|1> -s <suffix> [-m <peaks|bckg>] [-c <HEX1,HEX2,HEX3>] [-d <Y|N>]\n\n")
  
  cat("General parameters:\n")
  cat("  -t   Apply signed sqrt transformation: Y or N (optional, default = Y)\n")
  cat("  -y   Numeric value for the upper y-axis limit (optional)\n")
  cat("       If omitted, the plot uses automatic y-limits\n")
  cat("  -m   Mode: peaks (default) or bckg (optional, case-insensitive)\n")
  cat("  -c   3 comma-separated HEX colors for Positives,Negatives,Targets (optional)\n")
  cat("       Examples: FFA500,B0B0B0,FF0000   or   #FFA500,#B0B0B0,#FF0000\n")
  cat("  -d   Show dots in the plot: N (default) or Y (optional)\n")
  cat("  -i   Input directory for automatic file discovery (optional, default = current working directory)\n\n")
  
  cat("Explicit-input mode parameters:\n")
  cat("  -P   top_positives_peaks file\n")
  cat("  -N   top_negatives_peaks file\n")
  cat("  -T   target windows file\n")
  cat("  -A   positives windows file\n")
  cat("  -B   negatives windows file\n")
  cat("  -o   output directory\n")
  cat("  -l   lambda flag: 0 or 1\n")
  cat("  -s   suffix\n\n")
  
  cat("Notes:\n")
  cat("  - If none of -P/-N/-T/-A/-B/-o/-l/-s are provided, the script runs in standalone discovery mode.\n")
  cat("  - If any of -P/-N/-T/-A/-B/-o/-l/-s are provided, then all of them must be provided.\n")
  cat("  - With the new defaults, standalone mode can be run with no required parameters.\n\n")
  
  cat("Output file naming convention:\n")
  cat("  outfile_<type>_<params>_<suffix>.<ext>\n\n")
  cat("  where <params> is a compact summary of selected options:\n")
  cat("    trX          : transform flag (trY or trN)\n")
  cat("    ylN / ylauto : y-axis limit (e.g. yl10) or automatic scaling if -y is omitted\n")
  cat("    peaks/bckg   : normalization mode\n")
  cat("    dtX          : dots flag (dtY or dtN)\n")
  cat("    colors       :\n")
  cat("                   - cldefault -> default colors\n")
  cat("                   - otherwise HEX codes joined with '-' (without #)\n\n")
  
  cat("Examples:\n")
  cat("  Rscript Plots_ViPeaQ.R\n")
  cat("  Rscript Plots_ViPeaQ.R -i /path/to/data -m bckg -d N\n")
  cat("  Rscript Plots_ViPeaQ.R -y 10 -t N -c \"FFA500,B0B0B0,FF0000\"\n")
  cat("  Rscript Plots_ViPeaQ.R -P top_pos.tsv -N top_neg.tsv -T target.tsv -A pos_win.tsv -B neg_win.tsv -o ./out -l 1 -s sampleA\n\n")
}

fail <- function(msg) {
  cat(paste0("ERROR: ", msg ,"\n",print_help(), "\n"), file = stderr())
  quit(status = 1)
}

file_must_exist <- function(path, label) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    fail(paste0("Missing path for ", label, "."))
  }
  path <- path.expand(path)
  if (!file.exists(path)) {
    fail(paste0("File does not exist for ", label, ":\n  ", path))
  }
  if (dir.exists(path)) {
    fail(paste0("Expected a file for ", label, " but got a directory:\n  ", path))
  }
  if (file.access(path, 4) != 0) {
    fail(paste0("File is not readable for ", label, ":\n  ", path))
  }
  normalizePath(path)
}

dir_must_exist <- function(path, label = "directory") {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    fail(paste0("Missing path for ", label, "."))
  }
  path <- path.expand(path)
  path <- sub("/+$", "", path)
  if (!dir.exists(path)) {
    fail(paste0("Directory does not exist for ", label, ":\n  ", path))
  }
  if (file.access(path, 4) != 0) {
    fail(paste0("Directory is not readable for ", label, ":\n  ", path))
  }
  normalizePath(path)
}

find_unique_file <- function(path, pattern, label, exclude_pattern = NULL) {
  matches <- list.files(path = path, pattern = pattern, full.names = TRUE)
  
  if (!is.null(exclude_pattern)) {
    matches <- matches[!grepl(exclude_pattern, basename(matches))]
  }
  
  if (length(matches) == 0) {
    fail(
      paste0(
        "No file found for ", label, " in directory:\n  ", path,
        "\nPattern used: ", pattern
      )
    )
  }
  
  if (length(matches) > 1) {
    fail(
      paste0(
        "Ambiguous input for ", label, ": multiple files match.\n",
        "This likely indicates inconsistent or confusing file naming.\n\n",
        "Matching files:\n  - ",
        paste(basename(matches), collapse = "\n  - "),
        "\n\nPlease ensure only one file matches pattern:\n  ", pattern
      )
    )
  }
  
  cat("Found", label, "file:", matches[1], "\n")
  normalizePath(matches[1])
}

extract_suffix_from_top_peaks <- function(file_name, prefix) {
  sub(paste0("^", prefix, "(.*)\\.tsv$"), "\\1", basename(file_name))
}

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
  
  # keep only valid, non-missing groups
  plot.data <- plot.data[!is.na(plot.data$group), , drop = FALSE]
  plot.data$group <- as.character(plot.data$group)
  
  if (nrow(plot.data) == 0) {
    stop("plot.data is empty after removing rows with missing group.")
  }
  
  set.seed(seed)
  
  group_counts <- table(plot.data$group)
  group_names <- names(group_counts)
  group_props <- as.numeric(group_counts) / sum(group_counts)
  names(group_props) <- group_names
  
  group_targets <- pmax(1, round(max_points * group_props))
  names(group_targets) <- group_names
  
  # correct rounding drift
  diff_n <- max_points - sum(group_targets)
  if (diff_n != 0) {
    ord <- order(group_targets, decreasing = (diff_n > 0))
    for (i in seq_len(abs(diff_n))) {
      j <- ord[((i - 1) %% length(ord)) + 1]
      group_targets[j] <- group_targets[j] + sign(diff_n)
    }
  }
  
  split_df <- split(plot.data, plot.data$group, drop = TRUE)
  
  sampled_list <- lapply(names(split_df), function(grp) {
    df <- split_df[[grp]]
    n_keep <- group_targets[[grp]]
    
    if (is.null(n_keep) || is.na(n_keep) || n_keep <= 0) {
      stop(paste0("Internal error while subsetting plot data: invalid subset size for group '", grp, "'."))
    }
    
    n_keep <- min(nrow(df), n_keep)
    df[sample.int(nrow(df), size = n_keep, replace = FALSE), , drop = FALSE]
  })
  
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

create_chip_signal_plot <- function(plot.data, plot_title, x_label) {
  required_groups <- c("Positives", "Negatives", "Targets")
  
  if (!all(c("group", "value") %in% colnames(plot.data))) {
    stop("plot.data must contain at least the columns 'group' and 'value'.")
  }
  
  plot.data <- plot.data %>%
    dplyr::filter(!is.na(group), !is.na(value)) %>%
    dplyr::mutate(group = factor(as.character(group), levels = required_groups)) %>%
    dplyr::filter(!is.na(group))
  
  if (nrow(plot.data) == 0) {
    stop("plot.data is empty after filtering valid groups and non-missing values.")
  }
  
  transform_flag_local <- toupper(transform_flag)
  dots_flag_local <- toupper(dots_flag)
  
  if (!(transform_flag_local %in% c("Y", "N"))) {
    stop("transform_flag must be 'Y' or 'N'.")
  }
  if (!(dots_flag_local %in% c("Y", "N"))) {
    stop("dots_flag must be 'Y' or 'N'.")
  }
  
  if (!is.null(ylimit)) {
    if (!is.numeric(ylimit) || length(ylimit) != 1 || is.na(ylimit) || ylimit <= 0) {
      stop("ylimit must be NULL or a single numeric value > 0.")
    }
  }
  
  if (transform_flag_local == "Y") {
    plot.data <- plot.data %>%
      dplyr::mutate(value = sign(value) * sqrt(abs(value)))
    ylegend <- "Signed Sqrt Transformed Relative Enrichment"
  } else {
    message("Skipping signed sqrt transformation")
    ylegend <- "Relative Enrichment"
  }
  
  stats <- plot.data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      count = dplyr::n(),
      .groups = "drop"
    )
  
  stats <- dplyr::left_join(
    tibble::tibble(group = factor(required_groups, levels = required_groups)),
    stats,
    by = "group"
  ) %>%
    dplyr::mutate(count = ifelse(is.na(count), 0, count))
  
  custom_left_labels <- function(y) {
    sapply(y, function(value) {
      if (is.na(value)) {
        return("")
      } else if (value %in% stats$median) {
        group <- stats$group[which(stats$median == value)[1]]
        return(
          paste0(
            "<b><span style='color:",
            group_colors[group],
            "; font-size:15px;'>",
            round(value, 2),
            "</span></b>"
          )
        )
      } else {
        return(as.character(value))
      }
    })
  }
  
  n_pos <- stats$count[stats$group == "Positives"]
  n_neg <- stats$count[stats$group == "Negatives"]
  n_tar <- stats$count[stats$group == "Targets"]
  
  xlabels_expr <- c(
    "Positives" = bquote(atop("Positives", bold(n == .(scales::comma(n_pos))))),
    "Negatives" = bquote(atop("Background", bold(n == .(scales::comma(n_neg))))),
    "Targets"   = bquote(atop("Targets", bold(n == .(scales::comma(n_tar)))))
  )
  
  y_expand <- ggplot2::expansion(mult = c(0, 0.05))
  
  plot_peaks <- ggplot(
    plot.data,
    aes(
      x = factor(group, levels = required_groups),
      y = value,
      fill = group
    )
  ) +
    geom_violin(
      trim = TRUE,
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
        base_breaks <- pretty(limits)
        medians <- stats$median[!is.na(stats$median)]
        unique(sort(c(base_breaks, medians)))
      },
      labels = custom_left_labels,
      expand = y_expand
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
    scale_color_manual(values = group_colors, drop = FALSE)
  
  if (dots_flag_local == "Y") {
    plot_peaks <- plot_peaks +
      ggbeeswarm::geom_quasirandom(alpha = 0.6, size = 1.5)
  }
  
  if (!is.null(ylimit)) {
    plot_peaks <- plot_peaks + coord_cartesian(ylim = c(0, ylimit))
  }
  
  return(plot_peaks)
}

parse_common_options <- function(param_list) {
  transform_flag_local <- "Y"
  if (!is.null(param_list[["-t"]])) {
    transform_flag_local <- toupper(param_list[["-t"]])
    if (!(transform_flag_local %in% c("Y", "N"))) {
      fail("Parameter -t must be 'Y' or 'N' (example: -t Y).")
    }
  }
  
  ylimit_local <- NULL
  if (!is.null(param_list[["-y"]])) {
    ylimit_local <- suppressWarnings(as.numeric(param_list[["-y"]]))
    if (is.na(ylimit_local) || length(ylimit_local) != 1 || ylimit_local <= 0) {
      fail("Parameter -y must be a single numeric value > 0 (example: -y 10).")
    }
  }
  
  mode_local <- "peaks"
  if (!is.null(param_list[["-m"]])) {
    mode_input <- tolower(param_list[["-m"]])
    if (mode_input == "peaks") {
      mode_local <- "peaks"
    } else if (mode_input == "bckg") {
      mode_local <- "bckg"
    } else {
      fail("Parameter -m must be 'peaks' or 'bckg' (example: -m bckg).")
    }
  }
  
  dots_flag_local <- "N"
  if (!is.null(param_list[["-d"]])) {
    dots_flag_local <- toupper(param_list[["-d"]])
    if (!(dots_flag_local %in% c("Y", "N"))) {
      fail("Parameter -d must be 'Y' or 'N' (example: -d N).")
    }
  }
  
  group_colors_local <- c("Positives" = "orange", "Negatives" = "grey", "Targets" = "red")
  if (!is.null(param_list[["-c"]])) {
    color_input <- param_list[["-c"]]
    color_values <- strsplit(color_input, ",", fixed = TRUE)[[1]]
    color_values <- trimws(color_values)
    color_values <- ifelse(grepl("^#", color_values), color_values, paste0("#", color_values))
    
    if (length(color_values) != 3) {
      fail("Parameter -c must contain exactly 3 comma-separated HEX colors (example: -c \"FFA500,B0B0B0,FF0000\").")
    }
    
    if (!all(grepl("^#([A-Fa-f0-9]{6})$", color_values))) {
      fail("Colors in -c must be valid 6-digit HEX codes (example: -c \"FFA500,B0B0B0,FF0000\").")
    }
    
    color_values <- toupper(color_values)
    group_colors_local <- setNames(color_values, c("Positives", "Negatives", "Targets"))
  }
  
  list(
    transform_flag = transform_flag_local,
    ylimit = ylimit_local,
    mode = mode_local,
    dots_flag = dots_flag_local,
    group_colors = group_colors_local
  )
}

run_standalone_mode <- function(param_list) {
  current_dir_local <- if (is.null(param_list[["-i"]])) getwd() else param_list[["-i"]]
  current_dir_local <- dir_must_exist(current_dir_local, "input directory")
  
  input_dir <- current_dir_local
  output_dir <- current_dir_local
  lambda_flag <- "0"
  
  cat("Input mode: automatic file discovery\n")
  cat("Using input directory:", input_dir, "\n")
  
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
  
  top_negatives_peaks_file <- find_unique_file(
    path = input_dir,
    pattern = "^top_negatives_peaks(.*)\\.tsv$",
    label = "top negatives peaks"
  )
  
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
  
  positives_windows_file <- find_unique_file(
    path = input_dir,
    pattern = "^positives_win_count(_lambda_corrected_filtered)?(.*)?\\.tsv$",
    label = "positives windows"
  )
  
  negatives_windows_file <- find_unique_file(
    path = input_dir,
    pattern = "^negatives_win_count(_lambda_corrected_filtered)?(.*)?\\.tsv$",
    label = "negatives windows"
  )
  
  suffix_local <- sub("^_", "", suffix_raw)
  
  list(
    top_positives_peaks_file = top_positives_peaks_file,
    top_negatives_peaks_file = top_negatives_peaks_file,
    target_windows_file = target_windows_file,
    positives_windows_file = positives_windows_file,
    negatives_windows_file = negatives_windows_file,
    output_dir = output_dir,
    lambda_flag = lambda_flag,
    suffix = suffix_local,
    input_dir = input_dir
  )
}

run_explicit_mode <- function(param_list) {
  explicit_file_flags <- c("-P", "-N", "-T", "-A", "-B", "-o", "-l", "-s")
  missing_explicit <- explicit_file_flags[!(explicit_file_flags %in% names(param_list))]
  
  if (length(missing_explicit) > 0) {
    fail(paste0(
      "Explicit-input mode was detected because some of -P/-N/-T/-A/-B/-o/-l/-s were provided, ",
      "but the following required parameter(s) are missing: ",
      paste(missing_explicit, collapse = ", "),
      "."
    ))
  }
  
  top_positives_peaks_file <- file_must_exist(param_list[["-P"]], "top positives peaks")
  top_negatives_peaks_file <- file_must_exist(param_list[["-N"]], "top negatives peaks")
  target_windows_file      <- file_must_exist(param_list[["-T"]], "target windows")
  positives_windows_file   <- file_must_exist(param_list[["-A"]], "positives windows")
  negatives_windows_file   <- file_must_exist(param_list[["-B"]], "negatives windows")
  output_dir               <- dir_must_exist(param_list[["-o"]], "output directory")
  
  lambda_flag <- as.character(param_list[["-l"]])
  if (!(lambda_flag %in% c("0", "1"))) {
    fail("Parameter -l must be '0' or '1' (example: -l 1).")
  }
  
  suffix_local <- as.character(param_list[["-s"]])
  suffix_local <- sub("^_", "", suffix_local)
  
  cat("Input mode: explicit file parameters\n")
  cat("Top positives peaks:", top_positives_peaks_file, "\n")
  cat("Top negatives peaks:", top_negatives_peaks_file, "\n")
  cat("Target windows:", target_windows_file, "\n")
  cat("Positives windows:", positives_windows_file, "\n")
  cat("Negatives windows:", negatives_windows_file, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Lambda corrected:", lambda_flag, "\n")
  cat("Final suffix:", suffix_local, "\n")
  
  list(
    top_positives_peaks_file = top_positives_peaks_file,
    top_negatives_peaks_file = top_negatives_peaks_file,
    target_windows_file = target_windows_file,
    positives_windows_file = positives_windows_file,
    negatives_windows_file = negatives_windows_file,
    output_dir = output_dir,
    lambda_flag = lambda_flag,
    suffix = suffix_local,
    input_dir = NULL
  )
}

if (any(args %in% c("-h", "--help"))) {
  print_help()
  quit(status = 0)
}

allowed_flags <- c("-i", "-y", "-t", "-m", "-c", "-d", "-P", "-N", "-T", "-A", "-B", "-o", "-l", "-s")

if (length(args) > 0 && length(args) %% 2 != 0) {
  if (any(args == "-c")) {
    fail("Parameter -c requires a value. If using HEX colors, wrap them in quotes (example: -c \"FFA500,B0B0B0,FF0000\").")
  }
  fail("Each parameter flag must be followed by a value (example: -y 10).")
}

if (length(args) == 0) {
  param_list <- list()
} else {
  input_flags  <- args[seq(1, length(args), by = 2)]
  input_values <- args[seq(2, length(args), by = 2)]
  
  if (!all(input_flags %in% allowed_flags)) {
    bad_flags <- unique(input_flags[!(input_flags %in% allowed_flags)])
    fail(paste0(
      "Unknown parameter(s): ",
      paste(bad_flags, collapse = ", "),
      " (allowed: ", paste(allowed_flags, collapse = ", "), ")."
    ))
  }
  
  dup_flags <- unique(input_flags[duplicated(input_flags)])
  if (length(dup_flags) > 0) {
    fail(paste0(
      "Duplicated parameter(s): ",
      paste(dup_flags, collapse = ", "),
      ". Each parameter must be provided only once."
    ))
  }
  
  param_list <- as.list(input_values)
  names(param_list) <- input_flags
}

common_opts <- parse_common_options(param_list)
transform_flag <- common_opts$transform_flag
ylimit         <- common_opts$ylimit
mode           <- common_opts$mode
dots_flag      <- common_opts$dots_flag
group_colors   <- common_opts$group_colors

explicit_file_flags <- c("-P", "-N", "-T", "-A", "-B", "-o", "-l", "-s")
use_explicit_inputs <- sum(explicit_file_flags %in% names(param_list)) > 0

inputs <- if (use_explicit_inputs) run_explicit_mode(param_list) else run_standalone_mode(param_list)

top_positives_peaks_file <- inputs$top_positives_peaks_file
top_negatives_peaks_file <- inputs$top_negatives_peaks_file
target_windows_file      <- inputs$target_windows_file
positives_windows_file   <- inputs$positives_windows_file
negatives_windows_file   <- inputs$negatives_windows_file
outdir                   <- inputs$output_dir
lambda                   <- as.integer(inputs$lambda_flag)
suffix                   <- inputs$suffix
current_dir              <- if (is.null(inputs$input_dir)) outdir else inputs$input_dir

cat("Output directory:", outdir, "\n")
cat("Lambda corrected:", lambda, "\n")
cat("Final suffix:", suffix, "\n")
cat("Mode:", mode, "\n")
cat("Dots:", dots_flag, "\n")
cat("Transform:", transform_flag, "\n")
cat("Y-limit:", if (is.null(ylimit)) "auto" else ylimit, "\n")
cat("Group colors:\n")
print(group_colors)

########################################
##  READ INPUT TABLES                 ##
########################################
pos_peaks <- read.table(top_positives_peaks_file, header = FALSE, sep = "\t", quote = "\"")
neg_peaks <- read.table(top_negatives_peaks_file, header = FALSE, sep = "\t", quote = "\"")
target    <- read.table(target_windows_file,      header = FALSE, sep = "\t", quote = "\"")
pos_win   <- read.table(positives_windows_file,   header = FALSE, sep = "\t", quote = "\"")
neg_win   <- read.table(negatives_windows_file,   header = FALSE, sep = "\t", quote = "\"")

########################################
##  CALCULATE RATIO                   ##
########################################
ratio_pos_peaks <- pos_peaks[, 9] / pos_peaks[, 8]
ratio_neg_peaks <- neg_peaks[, 9] / neg_peaks[, 8]

if (lambda == 0) {
  target  <- target[target[, 8]  != 0, , drop = FALSE]
  pos_win <- pos_win[pos_win[, 8] != 0, , drop = FALSE]
  neg_win <- neg_win[neg_win[, 8] != 0, , drop = FALSE]
  
  ratio_target  <- target[, 9]  / target[, 8]
  ratio_pos_win <- pos_win[, 9] / pos_win[, 8]
  ratio_neg_win <- neg_win[, 9] / neg_win[, 8]
} else {
  target  <- target[target[, 11]  != 0, , drop = FALSE]
  pos_win <- pos_win[pos_win[, 11] != 0, , drop = FALSE]
  neg_win <- neg_win[neg_win[, 11] != 0, , drop = FALSE]
  
  ratio_target  <- target[, 9]  / target[, 11]
  ratio_pos_win <- pos_win[, 9] / pos_win[, 11]
  ratio_neg_win <- neg_win[, 9] / neg_win[, 11]
}

########################################
##  CALCULATE MEDIAN                  ##
########################################
median_neg_peaks <- median(ratio_neg_peaks, na.rm = TRUE)
median_neg_win   <- median(ratio_neg_win, na.rm = TRUE)
median_pos_peaks <- median(ratio_pos_peaks, na.rm = TRUE)
median_pos_win   <- median(ratio_pos_win, na.rm = TRUE)

########################################
##  SAFETY CHECK FOR BCKG MODE        ##
########################################
# if (mode == "bckg"&& (median_neg_peaks == 0 || median_neg_win == 0)) {
#   warning(
#     paste0(
#       "Background/control signal is not noisy enough to run mode = 'bckg' ",
#       "\nMedian of background distribution computed from the peaks = ", median_neg_peaks,
#       "\nMedian of background distribution computed from the bins = ", median_neg_win, ".\n",
#       "Please use mode = 'peaks' instead.\n",
#       "This is often the case in CUT&Tag data."
#     ),
#     call. = FALSE
#   )
#   quit(status = 1)
# }
if (mode == "bckg" && (median_neg_peaks == 0 || median_neg_win == 0)) {
  fail(
    paste0(
      "Background/control signal is not noisy enough to run mode = 'bckg' ",
      "\nMedian of background distribution computed from the peaks = ", median_neg_peaks,
      "\nMedian of background distribution computed from the bins = ", median_neg_win, ".\n",
      "Please use mode = 'peaks' instead. ",
      "This is often the case in CUT&Tag data."
    )
  )
}

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
  fail("Invalid mode: must be 'peaks' or 'bckg'.")
}

########################################
##  CALCULATE DISTRIBUTIONS           ##
########################################
pos_distribution     <- ratio_pos_peaks / ref_median_peaks
neg_distribution     <- ratio_neg_peaks / ref_median_peaks
target_distribution  <- ratio_target    / ref_median_peaks

pos_distribution_win    <- ratio_pos_win / ref_median_win
neg_distribution_win    <- ratio_neg_win / ref_median_win
target_distribution_win <- ratio_target  / ref_median_win

########################################
##  PREPARE PLOT DATA                 ##
########################################
pos_peaks_data <- data.frame(group = "Positives", value = pos_distribution)
neg_peaks_data <- data.frame(group = "Negatives", value = neg_distribution)
target_data    <- data.frame(group = "Targets",   value = target_distribution)

pos_win_data    <- data.frame(group = "Positives", value = pos_distribution_win)
neg_win_data    <- data.frame(group = "Negatives", value = neg_distribution_win)
target_win_data <- data.frame(group = "Targets",   value = target_distribution_win)

plot.data_peaks <- dplyr::bind_rows(pos_peaks_data, neg_peaks_data, target_data)
plot.data_win   <- dplyr::bind_rows(pos_win_data, neg_win_data, target_win_data)

########################################
##  BUILD OUTPUT TAGS                 ##
########################################
transformtag <- paste0("tr", toupper(transform_flag))
ylimit_tag <- if (is.null(ylimit)) "ylauto" else paste0("yl", format(ylimit, scientific = FALSE, trim = TRUE))
mode_tag <- tolower(mode)
dots_tag <- paste0("dt", toupper(dots_flag))

default_group_colors <- c("Positives" = "orange", "Negatives" = "grey", "Targets" = "red")
if (identical(unname(group_colors), unname(default_group_colors))) {
  color_tag <- "cldefault"
} else {
  clean_hex <- gsub("#", "", group_colors)
  color_tag <- paste(clean_hex, collapse = "-")
}

param_tag <- paste(transformtag, ylimit_tag, mode_tag, dots_tag, color_tag, sep = "_")
safe_suffix <- gsub("[^A-Za-z0-9._-]+", "_", suffix)

out_peaks_name <- paste0("outfile_peaks_", param_tag, "_", safe_suffix, ".pdf")
out_peaks <- file.path(outdir, out_peaks_name)

out_win_name <- paste0("outfile_win_", param_tag, "_", safe_suffix, ".pdf")
out_win <- file.path(outdir, out_win_name)

out_data_peaks_name <- paste0("outfile_peaks_", mode_tag, "_", safe_suffix, ".tsv")
out_data_peaks <- file.path(outdir, out_data_peaks_name)

out_data_win_name <- paste0("outfile_win_", mode_tag, "_", safe_suffix, ".tsv")
out_data_win <- file.path(outdir, out_data_win_name)

out_title_peaks <- paste("ChIP Signal Distribution Across Genome Peaks", suffix)
out_title_win   <- paste("ChIP Signal Distribution Across Genome Bins", suffix)

########################################
##  SAVE PLOTS                        ##
########################################
max_plot_points <- 50000

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

########################################
##  WRITE FULL OUTPUT TABLES          ##
########################################
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