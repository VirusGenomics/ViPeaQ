#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

######## Packages ##########

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(require(gridExtra))

######## Load data ########
# test if there is at least five argument: if not, return an error
if (length(args)<7) {
  stop("At least 6 arguments must be supplied", call.=FALSE)
}

args[1] <- "/home/robitaillea/test/out/top_positives_peaks.tsv"
args[2] <- "/home/robitaillea/test/out/top_negatives_peaks.tsv"
args[3] <- "/home/robitaillea/test/out/HQ404500_win_count_lambda_corrected.tsv"
args[4] <- "/home/robitaillea/test/out/positives_win_count_lambda_corrected_filtered.tsv"
args[5] <- "/home/robitaillea/test/out/negatives_win_count_lambda_corrected_filtered.tsv"
args[6] <- "/home/robitaillea/test/out/"
args[7] <- 1

pos_peaks <- read.table(args[1], header=F, sep="\t", quote="\"")

neg_peaks <- read.table(args[2], header=F, sep="\t", quote="\"")

target <- read.table(args[3], header=F, sep="\t", quote="\"")

pos_win <- read.table(args[4], header=F, sep="\t", quote="\"")

neg_win <- read.table(args[5], header=F, sep="\t", quote="\"")

outdir <- args[6]

lambda <- args[7]

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

means_peaks <- plot.data_peaks %>%
  group_by(group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

plot_peaks <- ggplot(plot.data_peaks, aes(x=group, y=value, fill=group)) + 
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values=c("grey","orange","red")) +
  scale_x_discrete(limits=c("Positives", "Negatives", "Targets")) +
  # stat_boxplot(geom ='errorbar',width = 0.025) +
  theme_classic() +
  theme(legend.position="none") +
  stat_summary(fun = "mean", geom="point", shape=23, size=2) +
  stat_summary(fun.data = stat_box_data_peaks, geom = "text", hjust = 0.5, vjust = 0.9) +
  geom_boxplot(width=0.05) +
  labs(x="Random 200 (neg sites)", y = "ChiP Signal") +
  geom_hline(data = means_peaks, aes(yintercept = mean_value, color = group), linetype = "dashed", linewidth = 0.25) +
  scale_color_manual(values = c("darkgrey", "orange", "red"))

means_win <- plot.data_win %>%
  group_by(group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

plot_win <- ggplot(plot.data_win, aes(x=group, y=value, fill=group)) + 
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values=c("grey","orange","red")) +
  scale_x_discrete(limits=c("Positives", "Negatives", "Targets")) +
  # stat_boxplot(geom ='errorbar',width = 0.025) +
  theme_classic() +
  theme(legend.position="none") +
  stat_summary(fun = "mean", geom="point", shape=23, size=2) +
  stat_summary(fun.data = stat_box_data_win, geom = "text", hjust = 0.5, vjust = 0.9) +
  geom_boxplot(width=0.05) +
  labs(x="Genome bins (neg sites)", y = "ChiP Signal") +
  geom_hline(data = means_win, aes(yintercept = mean_value, color = group), linetype = "dashed", linewidth = 0.25) +
  scale_color_manual(values = c("darkgrey", "orange", "red"))

# plot
out_peaks <- paste(outdir,"outfile_peaks.pdf",sep="/")
g_peaks <- grid.arrange(plot_peaks, nrow=1, ncol=1)
ggsave(out_peaks, g_peaks, width = 35, height = 20, units = "cm")
dev.off()

out_win <- paste(outdir,"outfile_win.pdf",sep="/")
g_win <- grid.arrange(plot_win, nrow=1, ncol=1)
ggsave(out_win, g_win, width = 35, height = 20, units = "cm")
dev.off()


out_data_peaks <- paste(outdir,"outfile_peaks.tsv",sep="/")
myls <- list(as.data.frame(pos_peaks_data), as.data.frame(neg_peaks_data), as.data.frame(target_data))
max.rows <- max(length(pos_peaks_data), length(neg_peaks_data), length(target_data))
new_myls <- lapply(myls, function(x) { x[1:max.rows,] })
new_myls <- as.data.frame(do.call(cbind, new_myls))
colnames(new_myls) <- c("Positives", "Negatives", "Targets")
write.table(new_myls, file = out_data_peaks, sep = "\t", row.names = FALSE)

out_data_win <- paste(outdir,"outfile_win.tsv",sep="/")
myls <- list(as.data.frame(pos_win_data), as.data.frame(neg_win_data), as.data.frame(target_data))
max.rows <- max(length(pos_win_data), length(neg_win_data), length(target_data))
new_myls <- lapply(myls, function(x) { x[1:max.rows,] })
new_myls <- as.data.frame(do.call(cbind, new_myls))
colnames(new_myls) <- c("Positives", "Negatives", "Targets")
write.table(new_myls, file = out_data_win, sep = "\t", row.names = FALSE)
