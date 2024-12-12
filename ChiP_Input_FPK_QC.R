#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Suppress warnings globally
options(warn = -1)

######## Packages ##########

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# Rscript ${BASEDIR}/ChiP_Input_FPK_QC.R ${outdir}/genome_win_count_lambda_corrected.tsv 1 ${outdir} genome_FPK.pdf
# Rscript ${BASEDIR}/ChiP_Input_FPK_QC.R ${outdir}/genome_win_count.tsv 0 ${outdir} genome_FPK.pdf

######## Load data ########
# test if there is at least five argument: if not, return an error
if (length(args)<4) {
  stop("At least four arguments must be supplied", call.=FALSE)
}

# data <- read.table("/home/robitaillea/test2/out/genome_win_count.tsv", header=T, sep="\t", quote="\"")
# lambda <- 0
# outdir <- "/mnt/R7525/robitaillea/test/out"
# outfile <- "genome_FPK.pdf"

data <- read.table(args[1], header=T, sep="\t", quote="\"")
lambda <- args[2]
outdir <- args[3]
outfile <- args[4]

if (lambda>0){
  fpk <- data$FPKInputCorrected
} else {
  fpk <- data$FPKInput
}

fpk_min10 <- fpk[fpk >= 10]
q <- quantile(fpk_min10, prob=c(.1,.5,.9), type=1)

fpk_min10_q <- fpk_min10[fpk_min10 > q[[1]] & fpk_min10 < q[[3]]]
# q[[1]]
# q[[2]]
# q[[3]]

data_plot <- data.frame(group = "FPK", value = fpk)
data_cov10_plot <- data.frame(group = "FPK_mincov10", value = fpk_min10)
data_filter_plot <- data.frame(group = "FPK_mincov10_percentiles", value = fpk_min10_q)

plot.data <- rbind(data_plot, data_cov10_plot, data_filter_plot)

stat_box_data <- function(y, upper_limit = q[[2]] * 2) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'mean =', round(mean(y), 2), '\n',
                    'median =', round(median(y), 2), '\n')
    )
  )
}

plot <- ggplot(plot.data, aes(x=group, y=value, fill=group)) +
  theme_classic() +
  theme(legend.position="none") +
  stat_summary(fun.y = "mean", geom="point", shape=23, size=2) +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9) +
  geom_boxplot(width=0.25,outlier.shape = NA) +
  ylim(0,q[[2]]*2) +
  labs(x="Input windows distributions", y = "FPK values")

# plot

out <- paste(outdir,outfile,sep="/")

pdf(file=out)
suppressMessages(print(plot))
garbage <- dev.off()
