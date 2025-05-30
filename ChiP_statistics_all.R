# Suppress warnings globally
options(warn = -1)

args = commandArgs(trailingOnly=TRUE)
# args = c()

######## Packages ##########

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(require(gridExtra))
suppressMessages(library(ggtext))
suppressMessages(library(tidyr))

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
  
  # Expand stats so each median appears for every category
  expanded_stats <- stats %>%
    crossing(group_plot = c("Positives", "Negatives", "Targets"))  # Ensure all medians appear on all violins
  
  # Calculate the 99.9th percentile threshold for each group and select the maximum
  cap_threshold <- plot.data %>%
    group_by(group) %>%
    summarise(threshold = quantile(value, 0.999)) %>%
    summarise(max_threshold = max(threshold)) %>%
    pull(max_threshold)
  
  # Define custom labels for the left and right axes
  # custom_left_labels <- function(y) {
  #   sapply(y, function(value) {
  #     if (value %in% pretty(range(plot.data$value))) {
  #       as.character(value)
  #     } else if (value %in% stats$median) {
  #       group <- stats$group[which(stats$median == value)]
  #       return(paste0("<b><span style='color:", group_colors[group], "; font-size:15px;'>", round(value, 2), "</span></b>"))
  #     } else {
  #       ""
  #     }
  #   })
  # }
  
  custom_left_labels <- function(y) {
    sapply(y, function(value) {
      if (is.na(value)) {  
        return("")  # Skip NA values to prevent errors
      } else if (value %in% pretty(range(plot.data$value, na.rm = TRUE))) {  
        return(as.character(value))  # Regular tick marks and ylimit
      } else if (value %in% stats$median) {
        group <- stats$group[which(stats$median == value)]
        return(paste0("<b><span style='color:", group_colors[group], "; font-size:15px;'>", round(value, 2), "</span></b>"))
      } else {
        return("")  # Hide everything else
      }
    })
  }
  
  
  # Define group colors
  group_colors <- c("Positives" = "orange", "Negatives" = "grey", "Targets" = "red")
  
  # Generate the plot
  plot_peaks <- ggplot(plot.data, aes(x = factor(group, levels = c("Positives", "Negatives", "Targets")), y = value, fill = group)) + 
    geom_violin(trim = TRUE, alpha = 0.6) + 
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(limits = c("Positives", "Negatives", "Targets"),
                     labels = c(
                       "Positives" = paste0("Positives<br><b>n=", scales::comma(stats$count[stats$group == "Positives"]), "</b>"),
                       "Negatives" = paste0("Background<br><b>n=", scales::comma(stats$count[stats$group == "Negatives"]), "</b>"),
                       "Targets" = paste0("Targets<br><b>n=", scales::comma(stats$count[stats$group == "Targets"]), "</b>")
                     )) +
    scale_y_continuous(
      name = "Signed Sqrt Transformed Relative Enrichment",
      breaks = sort(unique(c(0, pretty(range(plot.data$value)), stats$median))),
      labels = custom_left_labels,
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      axis.title.y.left = element_text(size = 16, face = "bold"),
      axis.title.y.right = element_text(size = 16, face = "bold"),
      axis.text.y.left = element_markdown(size = 12),
      axis.text.y.right = element_markdown(size = 12),
      axis.text.x = element_markdown(size = 12),
      legend.position = "none",
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      # 🔹 Ensure Axes Are Visible 🔹
      axis.line.x = element_line(color = "black", size = 0.4),  # X-axis line
      axis.line.y = element_line(color = "black", size = 0.4),  # Y-axis line
      axis.ticks.x = element_line(color = "black", size = 0.4),  # X-axis ticks
      axis.ticks.y = element_line(color = "black", size = 0.4)   # Y-axis ticks
    ) +
    labs(
      title = plot_title,
      x = x_label,
      y = NULL,
      caption = paste("Dashed lines represent the median signal.")
    ) +
    
    geom_hline(data = stats, aes(yintercept = median, color = group), linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
    # Add horizontal dashed median lines across all violin plots
    # geom_segment(
    #   data = expanded_stats, 
    #   aes(
    #     x = as.numeric(factor(group_plot, levels = c("Positives", "Negatives", "Targets"))) - 0.3,  
    #     xend = as.numeric(factor(group_plot, levels = c("Positives", "Negatives", "Targets"))) + 0.3,  
    #     y = median,        
    #     yend = median, 
    #     color = group
    #   ), 
    #   inherit.aes = FALSE,  
    #   linetype = "dashed", 
    #   size = 1, 
    #   alpha = 0.8
    # ) +
    
    # Ensure correct color mapping for dashed lines
    scale_color_manual(values = group_colors) + 
    
    coord_cartesian(ylim = c(0, cap_threshold * 1.1)) + 
    
    # geom_text(
    #   data = stats,
    #   aes(x = group, y = min(plot.data$value) - 0.5, label = paste0("n=", scales::comma(count))),
    #   size = 3.5, fontface = "bold", color = "black"
    # ) +
    scale_size(range = c(6, 12), guide = "none")  
  
  return(plot_peaks)
}


# plot
out_peaks_name <- paste("outfile_peaks",suffix,".pdf",sep="")
out_peaks <- paste(outdir,out_peaks_name,sep="/")
out_title_peaks <- paste("ChIP Signal Distribution Across Genome Peaks",suffix,sep=" ")
g_peaks <- create_chip_signal_plot(plot.data_peaks, out_title_peaks, "Genome Regions")
ggsave(filename = out_peaks, plot = g_peaks, width = 35, height = 20, units = "cm")

out_win_name <- paste("outfile_win",suffix,".pdf",sep="")
out_win <- paste(outdir,out_win_name,sep="/")
out_title_win <- paste("ChIP Signal Distribution Across Genome Bins",suffix,sep=" ")
g_win <- create_chip_signal_plot(plot.data_win, out_title_win, "Genome Bins")
ggsave(filename = out_win, plot = g_win, width = 35, height = 20, units = "cm")


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
