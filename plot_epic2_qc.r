#!/usr/bin/env Rscript

# Suppress warnings globally
options(warn = -1)

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(reshape2)
library(scales)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--peak_files"), type="character", default=NULL, help="Comma-separated list of peak files.", metavar="path"),
                    make_option(c("-s", "--sample_ids"), type="character", default=NULL, help="Comma-separated list of sample ids associated with peak files. Must be unique and in same order as peaks files input.", metavar="string"),
                    make_option(c("-o", "--outdir"), type="character", default='./', help="Output directory", metavar="path"),
                    make_option(c("-p", "--outprefix"), type="character", default='epic2_peakqc', help="Output prefix", metavar="string"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$peak_files)){
  print_help(opt_parser)
  stop("At least one peak file must be supplied", call.=FALSE)
}
if (is.null(opt$sample_ids)){
  print_help(opt_parser)
  stop("Please provide sample ids associated with peak files.", call.=FALSE)
}

if (file.exists(opt$outdir) == FALSE) {
  dir.create(opt$outdir,recursive=TRUE)
}

# PeakFiles <- unlist(strsplit("/mnt/R7525/robitaillea/normalization_scripts/20221130_hap1_prc1_prc2_ko2_chip/epic2_w200_hap1_wt_parental_h3k9me3.txt",","))
# SampleIDs <- unlist(strsplit("hap1_wt_parental_h3k9me3_w200",","))

PeakFiles <- unlist(strsplit(opt$peak_files,","))
SampleIDs <- unlist(strsplit(opt$sample_ids,","))
if (length(PeakFiles) != length(SampleIDs)) {
  print_help(opt_parser)
  stop("Number of sample ids must equal number of homer annotated files.", call.=FALSE)
}

################################################
################################################
## READ IN DATA                               ##
################################################
################################################

plot.dat <- data.frame()
summary.dat <- data.frame()
for (idx in 1:length(PeakFiles)) {
  
  # idx <- 1
  
  sampleid = SampleIDs[idx]
  isAnnot <- TRUE
  header <- c("#Chromosome","Start","End","PValue","Score", "Strand", "ChIPCount", "InputCount","FDR","log2FoldChange")
  fsplit <- unlist(strsplit(basename(PeakFiles[idx]), split='.',fixed=TRUE))
  if (fsplit[length(fsplit)] == 'bed') {
    isAnnot <- FALSE
  }
  if(isAnnot == FALSE){
    print_help(opt_parser)
    stop("Please provide the .txt epic2 peak files.", call.=FALSE)
  }
  peaks <- read.table(PeakFiles[idx], sep="\t", header=TRUE)
  colnames(peaks) <- header
  
  ## GET SUMMARY STATISTICS
  #peaks.dat <- peaks[,c('fold','-log10(qvalue)','-log10(pvalue)')]
  # peaks.dat <- peaks[,c('PValue','Score','FDR','log2FoldChange')]
  peaks.dat <- peaks[,c('Score','log2FoldChange')]
  peaks.dat$length <- (peaks$End - peaks$Start)
  peaks.dat$"TransformedPValue" <- -log10(peaks$PValue+1e-323)
  peaks.dat$"PValue" <- peaks$PValue
  peaks.dat$"TransformedFDR" <- -log10(peaks$FDR+1e-323)
  peaks.dat$"FDR" <- peaks$FDR
  
  for (cname in colnames(peaks.dat)) {
    sdat <- summary(peaks.dat[,cname])
    sdat["num_peaks"] <- nrow(peaks.dat)
    sdat["measure"] <- cname
    sdat["sample"] <- sampleid
    sdat <- t(data.frame(x=matrix(sdat),row.names=names(sdat)))
    summary.dat <- rbind(summary.dat,sdat)
  }
  colnames(peaks.dat) <- c('Score','log2FoldChange','length','TransformedPValue','PValue','TransformedFDR','FDR')
  peaks.dat$name <- rep(sampleid,nrow(peaks.dat))
  plot.dat <- rbind(plot.dat,peaks.dat)
}
plot.dat$name <- factor(plot.dat$name, levels=sort(unique(as.character(plot.dat$name))))

SummaryFile <- file.path(opt$outdir,paste(opt$outprefix,".summary.txt",sep=""))
write.table(summary.dat,file=SummaryFile,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


################################################
################################################
## PLOTS                                      ##
################################################
################################################

## RETURNS VIOLIN PLOT OBJECT
violin.plot <- function(plot.dat,x,y,ylab,title,log) {
  
  plot  <- ggplot(plot.dat, aes_string(x=x, y=y)) +
    geom_violin(aes_string(colour=x,fill=x), alpha = 0.3) +
    geom_boxplot(width=0.1) +
    xlab("") +
    ylab(ylab) +
    ggtitle(title) +
    theme(legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(colour="black"),
          axis.text.x= element_text(colour="black",face="bold"),
          axis.line.x = element_line(size = 1, colour = "black", linetype = "solid"),
          axis.line.y = element_line(size = 1, colour = "black", linetype = "solid"))
  if (log == 10) {
    plot <- plot + scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
  }
  if (log == 2) {
    plot <- plot + scale_y_continuous(trans='log2',breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))
  }
  return(plot)
}

############################

PlotFile <- file.path(opt$outdir,paste(opt$outprefix,".plots.pdf",sep=""))
pdf(PlotFile,height=6,width=3*length(unique(plot.dat$name)))

## PEAK COUNT PLOT
peak.count.dat <- as.data.frame(table(plot.dat$name))
colnames(peak.count.dat) <- c("name","count")
plot  <- ggplot(peak.count.dat, aes(x=name, y=count)) +
  geom_bar(stat="identity",aes(colour=name,fill=name), position = "dodge", width = 0.8, alpha = 0.3) +
  xlab("") +
  ylab("Number of peaks") +
  ggtitle("Peak count") +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour="black"),
        axis.text.x= element_text(colour="black",face="bold"),
        axis.line.x = element_line(size = 1, colour = "black", linetype = "solid"),
        axis.line.y = element_line(size = 1, colour = "black", linetype = "solid")) +
  geom_text(aes(label = count, x = name, y = count), position = position_dodge(width = 0.8), vjust = -0.6)
print(plot)

## VIOLIN PLOTS
print(violin.plot(plot.dat=plot.dat,x="name",y="length",ylab=expression(log[10]*" peak length"),title="Peak length distribution",log=10))
print(violin.plot(plot.dat=plot.dat,x="name",y="TransformedPValue",ylab=expression(-log[10]*" PValue +1e-323"),title="PValue distribution",log=-1))
print(violin.plot(plot.dat=plot.dat,x="name",y="Score",ylab=expression(log[10]*" Score"),title="Score distribution",log=10))
print(violin.plot(plot.dat=plot.dat,x="name",y="TransformedFDR",ylab=expression(-log[10]*" FDR +1e-323"),title="FDR distribution",log=-1))
print(violin.plot(plot.dat=plot.dat,x="name",y="log2FoldChange",ylab=expression(log[2]*" fold-enrichment"),title="Log2 Fold-change distribution",log=2))


garbage <- dev.off()

################################################
################################################
################################################
################################################
