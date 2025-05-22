# Suppress warnings globally
options(warn = -1)

args = commandArgs(trailingOnly=TRUE)
# args = c()

######## Packages ##########

suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbio))
suppressMessages(library(dplyr))

# ---- File paths ----

bed_file <- args[1]

gff_file <- args[2]

suffix <- args[3]

outdir <- args[4]

# ---- Read BED file as GRanges ----
bed_df <- read.table(bed_file, header = FALSE, sep = "\t",
                     col.names = c("seqid", "start", "end", "coverage"))

bed_df$midpoint <- (bed_df$start + bed_df$end) / 2

# Filter a single sequence if needed
seqid <- unique(bed_df$seqid)[1]
bed_df_filtered <- bed_df[bed_df$seqid == seqid, ]

# ---- Read GFF3 file and filter ----
gff <- import(gff_file, format = "gff3")
gff_filtered <- gff[seqnames(gff) == seqid & gff$type == "gene"]
gff_df <- as.data.frame(gff_filtered)

# Then extract a label column
gff_df$label <- ifelse(!is.na(gff_df$gene), gff_df$gene,
                       ifelse(!is.na(gff_df$Name), gff_df$Name,
                              ifelse(!is.na(gff_df$description), gff_df$description,
                                     ifelse(!is.na(gff_df$ID), gff_df$ID, NA))))

# ---- Greedy layer assignment to avoid overlap ----
gff_df <- gff_df[order(gff_df$start), ]
layers <- list()
assign_layer <- function(start, end) {
  for (i in seq_along(layers)) {
    if (start > layers[[i]]) {
      layers[[i]] <<- end
      return(i)
    }
  }
  layers[[length(layers) + 1]] <<- end
  return(length(layers))
}
gff_df$layer <- mapply(assign_layer, gff_df$start, gff_df$end)


# Position annotations in negative space
gff_df$track_y <- -0.25 * gff_df$layer

# Calculate y-limits
ymax <- max(bed_df_filtered$coverage, na.rm = TRUE)
ymin <- min(gff_df$track_y, na.rm = TRUE) - 0.2

# Compute y-axis breaks (simpler grid)
y_breaks <- pretty(c(0, ymax), n = 4)
y_breaks <- y_breaks[y_breaks >= 0]

p <- ggplot(bed_df_filtered, aes(x = midpoint, y = coverage)) +
  # Raw and smoothed coverage
  geom_line(alpha = 0.4, color = "gray40") +
  geom_smooth(se = FALSE, span = 0.5, color = "steelblue") +
  
  # Gene annotation layers below zero
  geom_segment(data = gff_df, aes(x = start, xend = end, y = track_y, yend = track_y),
               color = "darkred", size = 1.2, inherit.aes = FALSE) +
  geom_text(data = gff_df, aes(x = (start + end)/2, y = track_y - 0.05, label = label),
            size = 2.5, angle = 45, hjust = 1, inherit.aes = FALSE) +
  
  # Manual y-grid lines only above 0
  geom_hline(yintercept = y_breaks, color = "grey85", linewidth = 0.4) +
  
  # Axis scaling
  scale_y_continuous(
    limits = c(ymin, ymax * 1.1),
    breaks = y_breaks,
    expand = c(0, 0)
  ) +
  scale_x_continuous(expand = c(0.002, 0)) +  # Removes gap at x=0
  
  # Labels
  labs(
    x = "Genomic Position",
    y = "Normalized Coverage (FPKM)"
  ) +
  
  # Theme configuration
  theme_minimal() +  # Set base theme FIRST
  theme(
    # ✅ Remove x-axis line & vertical gridlines
    axis.line.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    # ✅ Keep subtle grey y-axis line
    axis.line.y = element_line(color = "grey60"),
    
    # ✅ Clean y-axis (left only)
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    
    # ✅ Suppress y-grid entirely (we add it manually above 0)
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # ✅ Tight margin to eliminate whitespace on left
    plot.margin = margin(t = 5, r = 10, b = 5, l = 0)
  )

outfilename <- paste("coverage",suffix,".pdf",sep="")
outfile <- paste(outdir,outfilename,sep="/")

ggsave(outfile, plot = p,
       width = 14, height = 6, units = "in", dpi = 300, device = cairo_pdf)
