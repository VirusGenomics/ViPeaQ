#!/usr/bin/env Rscript

# Suppress warnings globally
options(warn = -1)

args = commandArgs(trailingOnly=TRUE)

######## Packages ##########

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(require(gridExtra))
suppressMessages(library(ggtext))
suppressMessages(library(ggbreak))

######## Load data ########
if (length(args)<8) {
  stop("At least 8 arguments must be supplied", call.=FALSE)
}

pos_peaks <- read.table(args[1], header=F, sep="\t", quote="\"")

neg_peaks <- read.table(args[2], header=F, sep="\t", quote="\"")

target <- read.table(args[3], header=F, sep="\t", quote="\"")

pos_win <- read.table(args[4], header=F, sep="\t", quote="\"")

neg_win <- read.table(args[5], header=F, sep="\t", quote="\"")

outdir <- args[6]

lambda <- args[7]

suffix <- args[8]
suffix <- substring(suffix, 2)

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

median_neg_win <- median(ratio_neg_win)

##############################
##  CALCULATE HOST POS&NEG  ##
##############################
pos_distribution <- ratio_pos_peaks / median_neg_peaks
neg_distribution <- ratio_neg_peaks / median_neg_peaks

pos_distribution_win <- ratio_pos_win / median_neg_win
neg_distribution_win <- ratio_neg_win / median_neg_win

########################
##  CALCULATE TARGET  ##
########################
target_distribution <- ratio_target / median_neg_peaks
target_distribution_win <- ratio_target / median_neg_win

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
target_data <- data.frame(group = "Targets", value = target_distribution_win)

plot.data_win <- rbind(pos_win_data, neg_win_data, target_data)
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
  
  plot.data <- plot.data %>%
    mutate(value = sign(value) * sqrt(abs(value)))
  
  # Compute stats for the plot
  stats <- plot.data %>%
    group_by(group) %>%
    summarise(
      mean = mean(value),
      median = median(value),
      count = n()
    ) %>%
    ungroup()
  
  # Calculate the 99.9th percentile threshold for each group and select the maximum
  cap_threshold <- plot.data %>%
    group_by(group) %>%
    summarise(threshold = quantile(value, 0.999)) %>%
    summarise(max_threshold = max(threshold)) %>%
    pull(max_threshold)
  
  # Cap values exceeding the threshold
  # plot.data <- plot.data %>%
  #   mutate(value = ifelse(value > cap_threshold, cap_threshold, value))
  
  # Define custom labels for the left and right axes
  custom_left_labels <- function(y) {
    sapply(y, function(value) {
      if (value %in% pretty(range(plot.data$value))) {
        as.character(value)
      } else if (value %in% stats$mean) {
        group <- stats$group[which(stats$mean == value)]
        paste0("<span style='color:", group_colors[group], "; font-size:12px;'>", round(value, 2), "</span>")
      } else {
        ""
      }
    })
  }
  
  custom_right_labels <- function(y) {
    sapply(y, function(value) {
      if (value %in% stats$median) {
        group <- stats$group[which(stats$median == value)]
        paste0("<span style='color:", group_colors[group], "; font-size:12px;'>", round(value, 2), "</span>")
      } else {
        ""
      }
    })
  }
  
  # Define group colors
  group_colors <- c("Positives" = "grey", "Negatives" = "orange", "Targets" = "red")
  
  # Generate the plot
  plot_peaks <- ggplot(plot.data, aes(x = group, y = value, fill = group)) + 
    geom_violin(trim = TRUE, alpha = 0.6) + 
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(limits = c("Positives", "Negatives", "Targets")) +
    scale_y_continuous(
      name = "Mean",
      breaks = sort(unique(c(pretty(range(plot.data$value)), stats$mean))),
      labels = custom_left_labels,
      sec.axis = sec_axis(~., name = "Median", breaks = stats$median, labels = custom_right_labels)
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      axis.title.y.left = element_text(size = 16, face = "bold"),
      axis.title.y.right = element_text(size = 16, face = "bold"),
      axis.text.y.left = element_markdown(size = 12),
      axis.text.y.right = element_markdown(size = 12),
      axis.text.x = element_text(size = 12),
      legend.position = "none",
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = plot_title,
      x = x_label,
      y = NULL,
      #caption = paste("Values capped at the 99.9th percentile (", round(cap_threshold, 2), ").", sep = "", "\nDashed lines represent the mean signal; solid lines show the median.")
      caption = paste("Dashed lines represent the mean signal; solid lines show the median.")
    ) +
    
    # Add dashed lines for the mean
    geom_hline(data = stats, aes(yintercept = mean, color = group), linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
    geom_hline(data = stats, aes(yintercept = median, color = group), linetype = "solid", linewidth = 0.5, alpha = 0.5) +
    coord_cartesian(ylim = c(0, cap_threshold)) +
    geom_point(
      data = stats,
      aes(x = group, y = max(cap_threshold) * 1, size = count, color = group),
      shape = 21, fill = "white", stroke = 1.5, alpha = 0.9
    ) +
    geom_text(
      data = stats,
      aes(x = group, y = max(cap_threshold) * 1, label = scales::comma(count)),
      size = 2.5, fontface = "bold", color = "black"
    ) +
    scale_size(range = c(6, 12), guide = "none") +
    scale_color_manual(values = group_colors)
  
  return(plot_peaks)
}

# plot
out_peaks_name <- paste("outfile_peaks",suffix,".pdf",sep="")
out_peaks <- paste(outdir,out_peaks_name,sep="/")
out_title_peaks <- paste("ChiP Signal Distribution Across Genome Peaks",suffix,sep=" ")
g_peaks <- grid.arrange(create_chip_signal_plot(plot.data_peaks, out_title_peaks, "Genome Regions"), nrow=1, ncol=1)
ggsave(out_peaks, g_peaks, width = 35, height = 20, units = "cm")
dev.off()

out_win_name <- paste("outfile_win",suffix,".pdf",sep="")
out_win <- paste(outdir,out_win_name,sep="/")
out_title_win <- paste("ChiP Signal Distribution Across Genome Bins",suffix,sep=" ")
g_win <- grid.arrange(create_chip_signal_plot(plot.data_win, out_title_win, "Genome Bins"), nrow=1, ncol=1)
ggsave(out_win, g_win, width = 35, height = 20, units = "cm")
dev.off()

out_data_peaks_name <- paste("outfile_peaks",suffix,".tsv",sep="")
out_data_peaks <- paste(outdir,out_data_peaks_name,sep="/")
myls <- list(as.data.frame(pos_peaks_data), as.data.frame(neg_peaks_data), as.data.frame(target_data))
max.rows <- max(length(pos_peaks_data$value), length(neg_peaks_data$value), length(target_data$value))
new_myls <- lapply(myls, function(x) { x[1:max.rows,] })
value_columns <- lapply(new_myls, function(x) x$value)
new_myls <- as.data.frame(do.call(cbind, value_columns))
colnames(new_myls) <- c("Positives", "Negatives", "Targets")
write.table(new_myls, file = out_data_peaks, sep = "\t", row.names = FALSE)

out_data_win_name <- paste("outfile_win",suffix,".tsv",sep="")
out_data_win <- paste(outdir,out_data_win_name,sep="/")
myls <- list(as.data.frame(pos_win_data), as.data.frame(neg_win_data), as.data.frame(target_data))
max.rows <- max(length(pos_win_data$value), length(neg_win_data$value), length(target_data$value))
new_myls <- lapply(myls, function(x) { x[1:max.rows,] })
value_columns <- lapply(new_myls, function(x) x$value)
new_myls <- as.data.frame(do.call(cbind, value_columns))
colnames(new_myls) <- c("Positives", "Negatives", "Targets")
write.table(new_myls, file = out_data_win, sep = "\t", row.names = FALSE)
