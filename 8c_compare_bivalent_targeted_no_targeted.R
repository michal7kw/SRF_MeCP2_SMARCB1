#!/usr/bin/env Rscript

#' Compare SMARCB1 binding profiles between bivalent targeted and non-targeted genes
#' 
#' Input files:
#' - results/bigwig/BG1_CPM.bw: Control replicate 1 bigWig file
#' - results/bigwig/BG2_CPM.bw: Control replicate 2 bigWig file 
#' - results/bigwig/BG3_CPM.bw: Control replicate 3 bigWig file
#' - results/bigwig/BM3_CPM.bw: SMARCB1 ChIP-seq bigWig file
#' - data/gencode.vM10.annotation.gtf.gz: Gene annotations
#' - Gene_lists/bivalent/expressed_targeted_bivalent_NPCs_1000.csv: List of bivalent targeted genes
#' - Gene_lists/bivalent/expressed_not_targeted_bivalent_NPCs_1000.csv: List of bivalent non-targeted genes
#'
#' Output files:
#' - results/metaprofiles_comparison_R/combined_bivalent_targeted_vs_no_targeted_profile.pdf: 
#'   Combined plot showing binding profiles and fold changes
#' - Console output: Summary statistics for signal levels and fold changes
#'
#' The script compares SMARCB1 binding profiles around TSS regions between
#' bivalent targeted and non-targeted genes, generating metaprofiles and fold change plots.

# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(data.table)
library(parallel)
library(patchwork)  # For combining plots

# Import bigWig files and calculate mean signal
import_bigwig <- function(file_paths, bin_size = 10) {
  # Use parallel processing for multiple files
  if (length(file_paths) > 1) {
    n_cores <- max(1, parallel::detectCores() %/% 2)
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl))
    
    clusterEvalQ(cl, {
      library(rtracklayer)
      library(GenomicRanges)
    })
    
    bw_list <- parLapply(cl, file_paths, function(f) {
      tryCatch({
        cat(sprintf("\nProcessing file: %s\n", basename(f)))
        bw <- BigWigFile(f)
        import.bw(bw, as = "RleList")
      }, error = function(e) {
        stop(sprintf("Error processing file %s: %s", basename(f), e$message))
      })
    })
    
    if (length(bw_list) > 1) {
      seqinfo <- seqinfo(bw_list[[1]])
      bw_list <- lapply(bw_list, function(x) {
        seqinfo(x) <- seqinfo
        x
      })
      mean_coverage <- Reduce("+", bw_list) / length(bw_list)
    } else {
      mean_coverage <- bw_list[[1]]
    }
  } else {
    mean_coverage <- import.bw(file_paths[[1]], as = "RleList")
  }
  return(mean_coverage)
}

# Create TSS regions from gene list
create_tss_regions <- function(gene_list_file, gtf_file, upstream = 2500, downstream = 2500) {
  genes <- fread(gene_list_file, header = FALSE)$V1
  gtf <- rtracklayer::import(gtf_file)
  gene_records <- gtf[gtf$type == "gene" & gtf$gene_name %in% genes]
  
  tss_regions <- GRanges(
    seqnames = seqnames(gene_records),
    ranges = IRanges(
      start = ifelse(strand(gene_records) == "+",
                    start(gene_records) - upstream,
                    end(gene_records) - downstream),
      end = ifelse(strand(gene_records) == "+",
                  start(gene_records) + downstream,
                  end(gene_records) + upstream)
    ),
    strand = strand(gene_records),
    gene_name = gene_records$gene_name
  )
  sort(tss_regions)
}

# Calculate metaprofile matrix
calculate_metaprofile <- function(coverage, regions, bins = 100) {
  profile_matrix <- matrix(0, nrow = length(regions), ncol = bins)
  available_chrs <- names(coverage)
  chunk_size <- 1000
  n_chunks <- ceiling(length(regions) / chunk_size)
  
  # Print debug info
  cat("\nProcessing", length(regions), "regions\n")
  cat("Available chromosomes:", paste(available_chrs, collapse=", "), "\n")
  
  for (chunk in seq_len(n_chunks)) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, length(regions))
    chunk_regions <- regions[start_idx:end_idx]
    
    for (chr in unique(seqnames(chunk_regions))) {
      if (!as.character(chr) %in% available_chrs) {
        warning(sprintf("Chromosome %s not found in coverage data, skipping...", chr))
        next
      }
      
      chr_regions <- chunk_regions[seqnames(chunk_regions) == chr]
      if (length(chr_regions) == 0) next
      
      # Get chromosome length from coverage
      chr_length <- length(coverage[[as.character(chr)]])
      
      # Check if any regions are outside chromosome bounds
      out_of_bounds <- which(end(ranges(chr_regions)) > chr_length | start(ranges(chr_regions)) < 1)
      if (length(out_of_bounds) > 0) {
        warning(sprintf("Found %d regions outside chromosome bounds for %s", 
                      length(out_of_bounds), chr))
        # Filter out out-of-bounds regions
        chr_regions <- chr_regions[-out_of_bounds]
        if (length(chr_regions) == 0) next
      }
      
      cat(sprintf("Processing %d regions on chromosome %s (length: %d)\n", 
                length(chr_regions), chr, chr_length))
      
      chunk_scores <- tryCatch({
        Views(coverage[[as.character(chr)]], ranges(chr_regions))
      }, error = function(e) {
        warning(sprintf("Error processing chromosome %s: %s", chr, e$message))
        return(NULL)
      })
      
      if (is.null(chunk_scores)) next
      
      chunk_matrix <- vapply(chunk_scores, function(x) {
        x_len <- length(x)
        if (x_len == 0) return(numeric(bins))
        if (x_len == bins) return(as.numeric(x))
        approx(seq_len(x_len), as.numeric(x), n = bins)$y
      }, numeric(bins))
      
      chr_indices <- which(seqnames(chunk_regions) == chr)
      
      # Add dimension checks and debugging
      tryCatch({
        if (length(chr_indices) > 0) {
          target_indices <- start_idx + chr_indices - 1
          if (length(target_indices) != ncol(t(chunk_matrix))) {
            warning(sprintf("Dimension mismatch: target_indices length (%d) != chunk_matrix cols (%d)", 
                          length(target_indices), ncol(t(chunk_matrix))))
            # Ensure dimensions match by taking only what we can fit
            n_cols <- min(length(target_indices), ncol(t(chunk_matrix)))
            profile_matrix[target_indices[1:n_cols], ] <- t(chunk_matrix)[1:n_cols, ]
          } else {
            profile_matrix[target_indices, ] <- t(chunk_matrix)
          }
        }
      }, error = function(e) {
        warning(sprintf("Error assigning chunk matrix: %s", e$message))
      })
    }
  }
  return(profile_matrix)
}

# Calculate fold change between conditions
calculate_fold_change <- function(bg_profile, bm_profile) {
  bg_mean <- colMeans(bg_profile)
  bm_mean <- colMeans(bm_profile)
  log2(bm_mean / bg_mean)
}

# Main execution
main <- function() {
  # Set paths
  bg_files <- c("results/bigwig/BG1_CPM.bw",
                "results/bigwig/BG2_CPM.bw",
                "results/bigwig/BG3_CPM.bw")
  bm_file <- "results/bigwig/BM3_CPM.bw"
  gtf_file <- "data/gencode.vM10.annotation.gtf.gz"
  
  # Define gene lists
  bivalent_targeted_genes <- "Gene_lists/bivalent/expressed_targeted_bivalent_NPCs_1000.csv"
  bivalent_no_targeted_genes <- "Gene_lists/bivalent/expressed_not_targeted_bivalent_NPCs_1000.csv"
  
  # Create output directory
  dir.create("results/metaprofiles_comparison_R", showWarnings = FALSE)
  
  # Import bigWig data
  cat("Importing BG files...\n")
  bg_data <- import_bigwig(bg_files)
  cat("Importing BM file...\n")
  bm_data <- import_bigwig(list(bm_file))
  
  # Process gene lists
  cat("Processing bivalent targeted genes...\n")
  bivalent_targeted_regions <- create_tss_regions(bivalent_targeted_genes, gtf_file)
  cat("Processing bivalent no targeted genes...\n")
  bivalent_no_targeted_regions <- create_tss_regions(bivalent_no_targeted_genes, gtf_file)
  
  # Calculate profiles
  bivalent_targeted_bg_profile <- calculate_metaprofile(bg_data, bivalent_targeted_regions)
  bivalent_targeted_bm_profile <- calculate_metaprofile(bm_data, bivalent_targeted_regions)
  bivalent_no_targeted_bg_profile <- calculate_metaprofile(bg_data, bivalent_no_targeted_regions)
  bivalent_no_targeted_bm_profile <- calculate_metaprofile(bm_data, bivalent_no_targeted_regions)
  
  # Calculate means and standard errors
  x_pos <- seq(-2500, 2500, length.out = ncol(bivalent_targeted_bg_profile))
  
  # Create plotting data for individual profiles
  create_plot_data <- function(bg_profile, bm_profile, group_name) {
    bg_mean <- colMeans(bg_profile)
    bm_mean <- colMeans(bm_profile)
    bg_se <- apply(bg_profile, 2, function(x) sd(x)/sqrt(length(x)))
    bm_se <- apply(bm_profile, 2, function(x) sd(x)/sqrt(length(x)))
    
    data.frame(
      Position = rep(x_pos, 2),
      Signal = c(bg_mean, bm_mean),
      SE = c(bg_se, bm_se),
      Condition = rep(c("BG", "BM"), each = length(x_pos)),
      Group = group_name
    )
  }
  
  plot_data <- rbind(
    create_plot_data(bivalent_targeted_bg_profile, bivalent_targeted_bm_profile, "Bivalent Targeted"),
    create_plot_data(bivalent_no_targeted_bg_profile, bivalent_no_targeted_bm_profile, "Bivalent Not Targeted")
  )
  
  # Calculate fold changes
  bivalent_targeted_fc <- calculate_fold_change(bivalent_targeted_bg_profile, bivalent_targeted_bm_profile)
  bivalent_no_targeted_fc <- calculate_fold_change(bivalent_no_targeted_bg_profile, bivalent_no_targeted_bm_profile)
  
  fc_data <- data.frame(
    Position = rep(x_pos, 2),
    FoldChange = c(bivalent_targeted_fc, bivalent_no_targeted_fc),
    Group = rep(c("Bivalent Targeted", "Bivalent Not Targeted"), each = length(x_pos))
  )
  
  # Create main profile plot
  p1 <- ggplot(plot_data, aes(x = Position, y = Signal, color = Condition, linetype = Group)) +
    geom_line() +
    geom_ribbon(aes(ymin = Signal - SE, ymax = Signal + SE, fill = Condition), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("BG" = "#1f77b4", "BM" = "#d62728")) +
    scale_fill_manual(values = c("BG" = "#1f77b4", "BM" = "#d62728")) +
    labs(title = "SMARCB1 binding around TSS (Bivalent Genes)",
         x = "Distance from TSS (bp)",
         y = "Average RPKM") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Create fold change plot
  p2 <- ggplot(fc_data, aes(x = Position, y = FoldChange, color = Group)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Bivalent Targeted" = "#ff7f0e", "Bivalent Not Targeted" = "#2ca02c")) +
    labs(title = "Log2 Fold Change (BM/BG)",
         x = "Distance from TSS (bp)",
         y = "Log2 Fold Change") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Combine plots
  combined_plot <- p1 / p2
  
  # Save plots
  ggsave(
    filename = "results/metaprofiles_comparison_R/combined_bivalent_targeted_vs_no_targeted_profile.pdf",
    plot = combined_plot,
    width = 10,
    height = 12
  )
  
  # Calculate and print summary statistics
  cat("\nSummary Statistics:\n")
  cat("\nBivalent Targeted Genes:\n")
  cat("Mean BG signal:", mean(colMeans(bivalent_targeted_bg_profile)), "\n")
  cat("Mean BM signal:", mean(colMeans(bivalent_targeted_bm_profile)), "\n")
  cat("Mean log2 fold change:", mean(bivalent_targeted_fc), "\n")
  
  cat("\nBivalent Not Targeted Genes:\n")
  cat("Mean BG signal:", mean(colMeans(bivalent_no_targeted_bg_profile)), "\n")
  cat("Mean BM signal:", mean(colMeans(bivalent_no_targeted_bm_profile)), "\n")
  cat("Mean log2 fold change:", mean(bivalent_no_targeted_fc), "\n")
}

# Run the main function
if (!interactive()) {
  main()
} 