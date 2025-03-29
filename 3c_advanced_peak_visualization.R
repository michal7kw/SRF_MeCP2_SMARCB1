#!/usr/bin/env Rscript

# Advanced script for visualizing ChIP-seq peak annotations
# This script creates customized pie charts and additional visualizations

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(reshape2)
  library(RColorBrewer)
  library(ChIPseeker)
  library(GenomicFeatures)
})

# Set working directory
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1")

# Create output directory
dir.create("results/annotations/advanced_plots", showWarnings = FALSE, recursive = TRUE)

# Define sample names
samples <- c("BG1", "BG2", "BG3", "BM3", "BG_consistent")

# Define custom colors for pie charts (similar to the example image)
custom_colors <- c(
  "Promoter" = "#F8E16C",       # Light yellow
  "5' UTR" = "#F8E16C",         # Light yellow (same as Promoter)
  "3' UTR" = "#F8E16C",         # Light yellow (same as Promoter)
  "Exon" = "#D4F0F7",           # Light blue
  "Intron" = "#2A9DF4",         # Medium blue
  "Downstream" = "#FF8C42",     # Orange
  "Intergenic" = "#E0BBE4"       # Light purple
)

# Function to simplify annotation categories
simplify_annotation <- function(annotation) {
  # Map detailed annotations to simplified categories
  simplified <- case_when(
    grepl("Promoter", annotation) ~ "Promoter",
    grepl("5' UTR", annotation) ~ "5' UTR",
    grepl("3' UTR", annotation) ~ "3' UTR",
    grepl("Exon", annotation) ~ "Exon",
    grepl("Intron", annotation) ~ "Intron",
    grepl("Downstream", annotation) ~ "Downstream",
    grepl("Intergenic", annotation) ~ "Intergenic",
    TRUE ~ "Other"
  )
  return(simplified)
}

# Function to create custom pie chart
create_custom_pie <- function(data, title) {
  # Calculate percentages
  data$percentage <- data$count / sum(data$count) * 100
  
  # Add label with percentages
  data$label <- paste0(data$category, "\n(", round(data$percentage, 1), "%)")
  
  # Create pie chart
  p <- ggplot(data, aes(x = "", y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = custom_colors[data$category]) +
    theme_void() +
    theme(legend.position = "right") +
    labs(title = title, fill = "Genomic Feature") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              position = position_stack(vjust = 0.5))
  
  return(p)
}

# Process each sample
all_data <- data.frame()

for (sample in samples) {
  # Read annotation data
  anno_file <- paste0("results/annotations/", sample, "_peak_annotation.csv")
  
  # Check if file exists
  if (!file.exists(anno_file)) {
    message(paste0("Warning: Annotation file for ", sample, " not found. Run 2_peak_annotation.R first."))
    next
  }
  
  # Read and process data
  anno_data <- read.csv(anno_file)
  
  # Check if annotation column exists
  if (!"annotation" %in% colnames(anno_data)) {
    message(paste0("Warning: 'annotation' column not found in file for ", sample, ". Skipping."))
    next
  }
  
  # Simplify annotation categories
  anno_data$simplified_annotation <- simplify_annotation(anno_data$annotation)
  
  # Count peaks by simplified category - using more explicit approach
  category_counts <- data.frame(
    category = names(table(anno_data$simplified_annotation)),
    count = as.numeric(table(anno_data$simplified_annotation)),
    sample = sample
  )
  
  # Add sample information
  category_counts$sample <- sample
  
  # Append to all data
  all_data <- rbind(all_data, category_counts)
  
  # Create custom pie chart for this sample
  pie_chart <- create_custom_pie(category_counts, paste0("Peak Distribution - ", sample))
  
  # Save pie chart
  ggsave(paste0("results/annotations/advanced_plots/", sample, "_custom_pie.pdf"), 
         pie_chart, width = 8, height = 8)
  ggsave(paste0("results/annotations/advanced_plots/", sample, "_custom_pie.png"), 
         pie_chart, width = 8, height = 8, dpi = 300)
}

# Create combined visualization for all samples
if (nrow(all_data) > 0) {
  # Create a faceted pie chart
  all_samples_plot <- ggplot(all_data, aes(x = "", y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = custom_colors) +
    facet_wrap(~sample, ncol = 2) +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(title = "Peak Distribution Across Samples", fill = "Genomic Feature")
  
  # Save combined plot
  ggsave("results/annotations/advanced_plots/all_samples_pie.pdf", 
         all_samples_plot, width = 10, height = 10)
  ggsave("results/annotations/advanced_plots/all_samples_pie.png", 
         all_samples_plot, width = 10, height = 10, dpi = 300)
  
  # Create a stacked bar chart for comparison
  stacked_bar <- ggplot(all_data, aes(x = sample, y = count, fill = category)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = custom_colors) +
    theme_minimal() +
    labs(title = "Relative Peak Distribution Across Samples", 
         y = "Proportion", 
         x = "Sample",
         fill = "Genomic Feature")
  
  # Save stacked bar chart
  ggsave("results/annotations/advanced_plots/all_samples_stacked_bar.pdf", 
         stacked_bar, width = 10, height = 6)
  ggsave("results/annotations/advanced_plots/all_samples_stacked_bar.png", 
         stacked_bar, width = 10, height = 6, dpi = 300)
  
  # Create a heatmap of peak distribution
  # Reshape data for heatmap
  heatmap_data <- all_data %>%
    group_by(sample, category) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup() %>%
    dplyr::select(sample, category, percentage) %>%
    reshape2::dcast(category ~ sample, value.var = "percentage")
  
  # Convert to matrix for heatmap
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- heatmap_data$category
  
  # Save heatmap data
  write.csv(heatmap_data, "results/annotations/advanced_plots/peak_distribution_heatmap_data.csv")
  
  # Print summary
  message("Advanced visualizations completed successfully!")
  message(paste0("Results saved to: results/annotations/advanced_plots/"))
} else {
  message("No annotation data found. Please run 2_peak_annotation.R first.")
} 