#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(data.table)
    library(tidyverse)
})

# Set working directory
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1")

# Function to read and process matrix data
process_matrix <- function(matrix_file) {
    # Read the matrix using data.table for better performance with large files
    mat <- fread(cmd = paste("zcat", matrix_file), skip = 1)
    
    # Get gene names from the first column
    genes <- mat$V1
    
    # Calculate number of bins based on matrix dimensions
    n_bins <- (ncol(mat) - 6) / 2  # Subtract metadata columns and divide by 2 (BG and BM)
    
    # Extract BG and BM signal matrices
    bg_cols <- 7:(7 + n_bins - 1)
    bm_cols <- (7 + n_bins):(7 + 2*n_bins - 1)
    
    bg_matrix <- as.matrix(mat[, ..bg_cols])
    bm_matrix <- as.matrix(mat[, ..bm_cols])
    
    # Calculate log2 fold change between BM and BG
    # Add small constant to avoid log(0)
    epsilon <- 1e-5
    log2fc <- log2((bm_matrix + epsilon) / (bg_matrix + epsilon))
    
    # Calculate mean signal for sorting
    mean_signal <- rowMeans(bg_matrix + bm_matrix)
    
    # Sort genes by mean signal
    sort_idx <- order(mean_signal, decreasing = TRUE)
    
    # Return sorted matrices and gene names
    list(
        genes = genes[sort_idx],
        bg_matrix = bg_matrix[sort_idx,],
        bm_matrix = bm_matrix[sort_idx,],
        log2fc = log2fc[sort_idx,]
    )
}

# Function to create and save heatmap
create_heatmap <- function(data, output_file, title) {
    # Set up color mapping for signal strength
    max_signal <- max(c(data$bg_matrix, data$bm_matrix))
    if(max_signal == 0) max_signal <- 1  # Prevent zero maximum
    
    signal_colors <- colorRamp2(
        c(0, max_signal/2, max_signal),
        c("#FFFFFF", "#FFF7BC", "#D95F02")
    )
    
    # Set up color mapping for log2 fold change
    # Use symmetric scale around 0
    max_abs_fc <- max(abs(data$log2fc), na.rm = TRUE)
    if(is.infinite(max_abs_fc) || is.na(max_abs_fc)) max_abs_fc <- 1
    
    fc_colors <- colorRamp2(
        c(-max_abs_fc, 0, max_abs_fc),
        c("#313695", "#FFFFFF", "#A50026")
    )
    
    # Create heatmaps
    ht_list = Heatmap(data$bg_matrix,
                      name = "BG CPM",
                      col = signal_colors,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = TRUE,
                      row_names_gp = gpar(fontsize = 6),
                      column_title = "BG (mean of 3)",
                      row_names_max_width = unit(10, "cm")) +
              Heatmap(data$bm_matrix,
                      name = "BM CPM",
                      col = signal_colors,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = FALSE,
                      column_title = "BM") +
              Heatmap(data$log2fc,
                      name = "log2(BM/BG)",
                      col = fc_colors,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      show_row_names = FALSE,
                      column_title = "log2 Fold Change")
    
    # Save heatmap to PDF
    pdf(output_file, width = 16, height = 20)
    draw(ht_list, 
         column_title = title, 
         column_title_gp = gpar(fontsize = 14, fontface = "bold"))
    dev.off()
}

# Process each gene list
gene_lists <- c("enriched_down_regulated", 
                "enriched_not_disregulated", 
                "enriched_up_regulated")

for (list_name in gene_lists) {
    message(paste("Processing", list_name))
    
    # Read and process matrix
    matrix_file <- file.path("results/metaprofiles", paste0(list_name, "_matrix.gz"))
    data <- process_matrix(matrix_file)
    
    # Create heatmap
    output_file <- file.path("results/metaprofiles", paste0(list_name, "_heatmap.pdf"))
    create_heatmap(data, 
                  output_file, 
                  paste("SMARCB1 binding around TSS -", gsub("_", " ", list_name)))
} 