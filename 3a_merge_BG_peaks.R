#!/usr/bin/env Rscript

# Script to merge BG peaks and identify consistent peaks across replicates
# Requires peaks present in at least 2 replicates

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
})

# Set working directory
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1")

# Function to read narrowPeak files
read_narrowPeak <- function(file_path) {
  peaks <- read.table(file_path, 
                     col.names = c("chr", "start", "end", "name", "score", 
                                 "strand", "signalValue", "pValue", "qValue", "peak"))
  return(peaks)
}

# Read BG peak files
BG1_peaks <- read_narrowPeak("results/peaks/BG1_peaks.narrowPeak")
BG2_peaks <- read_narrowPeak("results/peaks/BG2_peaks.narrowPeak")
BG3_peaks <- read_narrowPeak("results/peaks/BG3_peaks.narrowPeak")

# Convert to GRanges objects
BG1_gr <- GRanges(seqnames = BG1_peaks$chr,
                  ranges = IRanges(start = BG1_peaks$start, end = BG1_peaks$end),
                  score = BG1_peaks$score,
                  signalValue = BG1_peaks$signalValue,
                  pValue = BG1_peaks$pValue,
                  qValue = BG1_peaks$qValue,
                  peak = BG1_peaks$peak)

BG2_gr <- GRanges(seqnames = BG2_peaks$chr,
                  ranges = IRanges(start = BG2_peaks$start, end = BG2_peaks$end),
                  score = BG2_peaks$score,
                  signalValue = BG2_peaks$signalValue,
                  pValue = BG2_peaks$pValue,
                  qValue = BG2_peaks$qValue,
                  peak = BG2_peaks$peak)

BG3_gr <- GRanges(seqnames = BG3_peaks$chr,
                  ranges = IRanges(start = BG3_peaks$start, end = BG3_peaks$end),
                  score = BG3_peaks$score,
                  signalValue = BG3_peaks$signalValue,
                  pValue = BG3_peaks$pValue,
                  qValue = BG3_peaks$qValue,
                  peak = BG3_peaks$peak)

# Find overlapping peaks
all_peaks <- c(BG1_gr, BG2_gr, BG3_gr)
hits <- findOverlaps(all_peaks, all_peaks)
overlaps <- data.frame(query = queryHits(hits), 
                      subject = subjectHits(hits))

# Create groups of overlapping peaks
peak_groups <- list()
processed_peaks <- c()

for (i in 1:nrow(overlaps)) {
  if (!(overlaps$query[i] %in% processed_peaks)) {
    # Find all peaks that overlap with current peak
    group <- unique(c(
      overlaps$query[overlaps$subject == overlaps$query[i]],
      overlaps$subject[overlaps$query == overlaps$query[i]]
    ))
    
    # Add group if it contains peaks from at least 2 different samples
    peak_sources <- ifelse(group <= length(BG1_gr), "BG1",
                          ifelse(group <= length(BG1_gr) + length(BG2_gr), "BG2", "BG3"))
    if (length(unique(peak_sources)) >= 2) {
      peak_groups[[length(peak_groups) + 1]] <- group
    }
    
    processed_peaks <- c(processed_peaks, group)
  }
}

# Merge overlapping peaks and calculate mean scores
consistent_peaks <- GRanges()
for (group in peak_groups) {
  peaks_in_group <- all_peaks[group]
  merged_peak <- reduce(peaks_in_group)
  
  # Calculate mean scores for the merged peak
  mcols(merged_peak)$score <- mean(peaks_in_group$score)
  mcols(merged_peak)$signalValue <- mean(peaks_in_group$signalValue)
  mcols(merged_peak)$pValue <- mean(peaks_in_group$pValue)
  mcols(merged_peak)$qValue <- mean(peaks_in_group$qValue)
  mcols(merged_peak)$peak <- mean(peaks_in_group$peak)
  
  consistent_peaks <- c(consistent_peaks, merged_peak)
}

# Convert back to data frame
consistent_peaks_df <- data.frame(
  chr = seqnames(consistent_peaks),
  start = start(consistent_peaks),
  end = end(consistent_peaks),
  name = paste0("peak_", 1:length(consistent_peaks)),
  score = mcols(consistent_peaks)$score,
  strand = ".",
  signalValue = mcols(consistent_peaks)$signalValue,
  pValue = mcols(consistent_peaks)$pValue,
  qValue = mcols(consistent_peaks)$qValue,
  peak = mcols(consistent_peaks)$peak
)

# Save consistent peaks
write.table(consistent_peaks_df,
            file = "results/peaks/BG_consistent_peaks.narrowPeak",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Print summary
message(paste0("Number of consistent peaks found: ", nrow(consistent_peaks_df)))
message("Results saved to: results/peaks/BG_consistent_peaks.narrowPeak") 