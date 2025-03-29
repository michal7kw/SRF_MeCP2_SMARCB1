#!/usr/bin/env Rscript

# Script for annotating ChIP-seq peaks and generating pie charts
# This script uses ChIPseeker to annotate peaks based on genomic features
# and generates pie charts showing the distribution of peaks

# Load required libraries
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(GenomicFeatures)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

# Set working directory
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1")

# Create output directory for annotations
dir.create("results/annotations", showWarnings = FALSE, recursive = TRUE)
dir.create("results/annotations/plots", showWarnings = FALSE, recursive = TRUE)

# Define sample names
samples <- c("BG1", "BG2", "BG3", "BM3", "BG_consistent")

# Define consistent color palette for genomic features
# This ensures the same colors are used across all plots
feature_colors <- c(
  "Promoter (<=1kb)" = "#E41A1C",     # Bright red
  "Promoter (1-2kb)" = "#984EA3",     # Purple
  "Promoter (2-3kb)" = "#4DAF4A",     # Green
  "5' UTR" = "#377EB8",               # Blue
  "3' UTR" = "#FF7F00",               # Orange
  "1st Exon" = "#FFFF33",             # Yellow
  "Other Exon" = "#A65628",           # Brown
  "1st Intron" = "#F781BF",           # Pink
  "Other Intron" = "#999999",         # Gray
  "Downstream (<=300)" = "#66C2A5",   # Teal
  "Distal Intergenic" = "#FC8D62",    # Coral
  "Intergenic" = "#8DA0CB"            # Light blue
)

# Import custom annotation from Gencode
message("Importing Gencode annotation...")
gencode_gtf <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1/data/gencode.vM10.basic.annotation.gtf.gz"
txdb_gencode <- makeTxDbFromGFF(gencode_gtf, format="gtf")

# Function to process each sample
process_sample <- function(sample) {
  message(paste0("Processing sample: ", sample))
  
  # Read peak file
  peak_file <- paste0("results/peaks/", sample, "_peaks.narrowPeak")
  peaks <- readPeakFile(peak_file)
  
  # Add chromosome lengths
  seqlevelsStyle(peaks) <- "UCSC"
  
  # Annotate peaks using ChIPseeker
  peakAnno <- annotatePeak(peaks, 
                          tssRegion=c(-3000, 3000),
                          TxDb=txdb_gencode,
                          annoDb="org.Mm.eg.db")
  
  # Save annotation results
  anno_df <- as.data.frame(peakAnno)
  write.csv(anno_df, paste0("results/annotations/", sample, "_peak_annotation.csv"), row.names = FALSE)
  
  # Generate pie chart for genomic feature distribution with consistent colors
  pdf(paste0("results/annotations/plots/", sample, "_genomic_feature_pie.pdf"), width=8, height=8)
  pie1 <- plotAnnoPie(peakAnno, col=feature_colors)
  title(main=sample, line=0)
  print(pie1)
  dev.off()
  
  # Generate bar plot for genomic feature distribution with consistent colors
  pdf(paste0("results/annotations/plots/", sample, "_genomic_feature_bar.pdf"), width=10, height=6)
  bar1 <- plotAnnoBar(peakAnno) + 
    ggtitle(sample) + 
    scale_fill_manual(values=feature_colors)
  print(bar1)
  dev.off()
  
  # Generate distribution of peaks relative to TSS
  pdf(paste0("results/annotations/plots/", sample, "_TSS_distribution.pdf"), width=10, height=6)
  tss1 <- plotDistToTSS(peakAnno) + ggtitle(sample)
  print(tss1)
  dev.off()
  
  return(peakAnno)
}

# Process all samples
anno_list <- lapply(samples, process_sample)
names(anno_list) <- samples

# Create a combined plot with all samples
message("Creating combined plots...")

# Combined pie charts
pdf("results/annotations/plots/combined_genomic_feature_pie.pdf", width=12, height=10)
# Create a layout for multiple plots
par(mfrow=c(3,2))
for(sample in names(anno_list)) {
  plotAnnoPie(anno_list[[sample]], col=feature_colors)
  title(main=sample, line=0)
}
dev.off()

# Combined bar plots
pdf("results/annotations/plots/combined_genomic_feature_bar.pdf", width=12, height=10)
bar_plots <- lapply(names(anno_list), function(sample) {
  p <- plotAnnoBar(anno_list[[sample]])
  p <- p + ggtitle(sample) + scale_fill_manual(values=feature_colors)
  return(p)
})
grid.arrange(grobs=bar_plots, ncol=1)
dev.off()

# Combined TSS distribution
pdf("results/annotations/plots/combined_TSS_distribution.pdf", width=12, height=10)
tss_plots <- lapply(names(anno_list), function(sample) {
  p <- plotDistToTSS(anno_list[[sample]])
  p <- p + ggtitle(sample)
  return(p)
})
grid.arrange(grobs=tss_plots, ncol=2)
dev.off()

# Generate a summary table with peak counts per feature for all samples
summary_table <- data.frame(Sample = character(),
                           Promoter = integer(),
                           Exon = integer(),
                           Intron = integer(),
                           Downstream = integer(),
                           Intergenic = integer(),
                           Total = integer(),
                           stringsAsFactors = FALSE)

for (sample in names(anno_list)) {
  anno_summary <- as.data.frame(anno_list[[sample]])
  
  # Count peaks by annotation
  counts <- table(anno_summary$annotation)
  
  # Create summary row
  row <- data.frame(
    Sample = sample,
    Promoter = sum(grepl("Promoter", names(counts)) * counts),
    Exon = sum(grepl("Exon", names(counts)) * counts),
    Intron = sum(grepl("Intron", names(counts)) * counts),
    Downstream = sum(grepl("Downstream", names(counts)) * counts),
    Intergenic = sum(grepl("Intergenic", names(counts)) * counts),
    Total = nrow(anno_summary),
    stringsAsFactors = FALSE
  )
  
  summary_table <- rbind(summary_table, row)
}

# Save summary table
write.csv(summary_table, "results/annotations/peak_annotation_summary.csv", row.names = FALSE)

message("Peak annotation completed successfully!") 