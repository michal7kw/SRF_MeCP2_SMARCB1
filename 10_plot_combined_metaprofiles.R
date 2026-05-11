#!/usr/bin/env Rscript

#' Combined Mecp2 and SMARCB1 binding profile analysis
#' 
#' This script generates metaprofiles comparing:
#' - Endogenous Mecp2 (NSCM)
#' - Exogenous Mecp2 (NSCv)
#' - SMARCB1 control (BG)
#' - SMARCB1 ChIP (BM)
#'
#' Input files:
#' - SMARCB1 control: results/bigwig/BG{1,2,3}_CPM.bw
#' - SMARCB1 ChIP: results/bigwig/BM3_CPM.bw
#' - Mecp2 endo: iterative_alternative/results_1b/bigwig/NSCM{1,2,3}.bw
#' - Mecp2 exo: iterative_alternative/results_1b/bigwig/NSCv{1,2,3}.bw
#' 
#' Output files:
#' - results/metaprofiles_comparison_R/combined_profiles.pdf

library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(data.table)
library(parallel)
library(patchwork)

# Import and average bigWig files
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
create_tss_regions <- function(gtf_file, upstream = 2500, downstream = 2500) {
  gtf <- rtracklayer::import(gtf_file)
  gene_records <- gtf[gtf$type == "gene"]
  
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

#' Calculate metaprofile matrix for genomic regions
#' 
#' This function calculates a matrix of signal values across genomic regions,
#' normalizing each region to a fixed number of bins. The calculation is done
#' in chunks to manage memory usage efficiently.
#'
#' @param coverage RleList containing the signal values across chromosomes
#' @param regions GRanges object containing the genomic regions to profile
#' @param bins Integer specifying number of bins to divide each region into (default 100)
#' @return Matrix where each row represents a region and each column represents a bin
calculate_metaprofile <- function(coverage, regions, bins = 100) {
  # Initialize matrix to store binned signal values for all regions
  profile_matrix <- matrix(0, nrow = length(regions), ncol = bins)
  available_chrs <- names(coverage)
  
  # Process regions in chunks to avoid memory issues
  chunk_size <- 1000  # Number of regions to process at once
  n_chunks <- ceiling(length(regions) / chunk_size)
  
  # Iterate through chunks of regions
  for (chunk in seq_len(n_chunks)) {
    # Calculate indices for current chunk
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, length(regions))
    chunk_regions <- regions[start_idx:end_idx]
    
    # Process each chromosome separately
    for (chr in unique(seqnames(chunk_regions))) {
      # Skip chromosomes not present in coverage data
      if (!as.character(chr) %in% available_chrs) next
      
      # Get regions for current chromosome
      chr_regions <- chunk_regions[seqnames(chunk_regions) == chr]
      if (length(chr_regions) == 0) next
      
      # Filter out regions extending beyond chromosome boundaries
      chr_length <- length(coverage[[as.character(chr)]])
      valid_regions <- end(ranges(chr_regions)) <= chr_length & start(ranges(chr_regions)) >= 1
      chr_regions <- chr_regions[valid_regions]
      if (length(chr_regions) == 0) next
      
      # Extract signal values for valid regions
      chunk_scores <- Views(coverage[[as.character(chr)]], ranges(chr_regions))
      
      # Process each region's signal values
      for (i in seq_along(chunk_scores)) {
        x <- chunk_scores[[i]]
        x_len <- length(x)
        
        # Handle different cases for binning:
        if (x_len == 0) {
          # No data available - fill with zeros
          profile <- numeric(bins)
        } else if (x_len == bins) {
          # Data length matches desired bins - use as is
          profile <- as.numeric(x)
        } else {
          # Interpolate data to match desired number of bins
          profile <- approx(seq_len(x_len), as.numeric(x), n = bins)$y
        }
        
        # Store binned profile in result matrix
        region_idx <- start_idx + i - 1
        if (region_idx <= nrow(profile_matrix)) {
          profile_matrix[region_idx, ] <- profile
        }
      }
    }
  }
  
  return(profile_matrix)
}

# Create plotting data
create_plot_data <- function(profile_matrix, group_name, x_pos) {
  mean_signal <- colMeans(profile_matrix)
  se <- apply(profile_matrix, 2, function(x) sd(x)/sqrt(length(x)))
  
  data.frame(
    Position = x_pos,
    Signal = mean_signal,
    SE = se,
    Group = group_name
  )
}

# Main execution
main <- function() {
  # Set paths
  bg_files <- Sys.glob("results/bigwig/BG*_CPM.bw")
  nscm_files <- Sys.glob("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_1b/bigwig/NSCM*.bw")
  gtf_file <- "data/gencode.vM10.annotation.gtf.gz"
  
  # Read target and non-target genes
  target_genes <- read.csv("Gene_lists/targets/all_targets_final.csv", header = FALSE)$V1
  no_target_genes <- read.csv("Gene_lists/targets/all_no_targets_final.csv", header = FALSE)$V1
  
  # Create output directory
  output_dir <- "results/metaprofiles_comparison_R"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Import data
  cat("Importing bigWig files...\n")
  bg_data <- import_bigwig(bg_files)
  nscm_data <- import_bigwig(nscm_files)
  
  # Create TSS regions for both gene sets
  cat("Creating TSS regions...\n")
  all_tss_regions <- create_tss_regions(gtf_file)
  
  # Filter TSS regions for target and non-target genes
  target_tss <- all_tss_regions[all_tss_regions$gene_name %in% target_genes]
  no_target_tss <- all_tss_regions[all_tss_regions$gene_name %in% no_target_genes]
  
  # Calculate profiles for both gene sets
  cat("Calculating profiles...\n")
  # Target genes
  bg_profile_targets <- calculate_metaprofile(bg_data, target_tss)
  nscm_profile_targets <- calculate_metaprofile(nscm_data, target_tss)
  
  # Non-target genes
  bg_profile_no_targets <- calculate_metaprofile(bg_data, no_target_tss)
  nscm_profile_no_targets <- calculate_metaprofile(nscm_data, no_target_tss)
  
  # Set up x-axis positions
  x_pos <- seq(-2500, 2500, length.out = ncol(bg_profile_targets))
  
  # Create plotting data for both sets
  plot_data_targets <- rbind(
    create_plot_data(bg_profile_targets, "SMARCB1 (BG)", x_pos),
    create_plot_data(nscm_profile_targets, "Mecp2 endo (NSCM)", x_pos)
  )
  plot_data_targets$Set <- "Target Genes"
  
  plot_data_no_targets <- rbind(
    create_plot_data(bg_profile_no_targets, "SMARCB1 (BG)", x_pos),
    create_plot_data(nscm_profile_no_targets, "Mecp2 endo (NSCM)", x_pos)
  )
  plot_data_no_targets$Set <- "Non-target Genes"
  
  # Combine all data
  all_plot_data <- rbind(plot_data_targets, plot_data_no_targets)
  
  # Find global y-axis limits for consistent scaling
  y_min <- min(all_plot_data$Signal - all_plot_data$SE)
  y_max <- max(all_plot_data$Signal + all_plot_data$SE)
  
  # Create plot function with independent scales for each protein
  create_scaled_plot <- function(data, title) {
    # Split data by protein type
    smarcb1_data <- data[data$Group == "SMARCB1 (BG)",]
    mecp2_data <- data[data$Group == "Mecp2 endo (NSCM)",]
    
    # Calculate scaling factor for secondary axis
    smarcb1_range <- range(smarcb1_data$Signal)
    mecp2_range <- range(mecp2_data$Signal)
    
    # Create plot with two y-axes
    p <- ggplot() +
      # SMARCB1 data on primary y-axis
      geom_line(data = smarcb1_data, 
                aes(x = Position, y = Signal, color = "SMARCB1 (BG)"),
                size = 1) +
      geom_ribbon(data = smarcb1_data,
                 aes(x = Position, 
                     ymin = Signal - SE, 
                     ymax = Signal + SE, 
                     fill = "SMARCB1 (BG)"), 
                 alpha = 0.2) +
      # Mecp2 data on secondary y-axis with transformed scale
      geom_line(data = mecp2_data,
                aes(x = Position, y = Signal, color = "Mecp2 endo (NSCM)"),
                size = 1) +
      geom_ribbon(data = mecp2_data,
                 aes(x = Position, 
                     ymin = Signal - SE, 
                     ymax = Signal + SE,
                     fill = "Mecp2 endo (NSCM)"),
                 alpha = 0.2) +
      scale_y_continuous(
        name = "SMARCB1 signal",
        sec.axis = sec_axis(~., name = "Mecp2 signal")
      ) +
      scale_color_manual(values = c(
        "SMARCB1 (BG)" = "#1f77b4",
        "Mecp2 endo (NSCM)" = "#2ca02c"
      )) +
      scale_fill_manual(values = c(
        "SMARCB1 (BG)" = "#1f77b4",
        "Mecp2 endo (NSCM)" = "#2ca02c"
      )) +
      labs(title = title,
           x = "Distance from TSS (bp)") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "right",
        axis.title.y.left = element_text(color = "#1f77b4"),
        axis.title.y.right = element_text(color = "#2ca02c"),
        axis.text.y.left = element_text(color = "#1f77b4"),
        axis.text.y.right = element_text(color = "#2ca02c")
      )
    
    return(p)
  }
  
  # Create both plots
  p1 <- create_scaled_plot(plot_data_targets, "Target Genes - Protein Binding Profiles")
  p2 <- create_scaled_plot(plot_data_no_targets, "Non-target Genes - Protein Binding Profiles")
  
  # Combine plots using patchwork
  combined_plot <- p1 / p2
  
  # Save plots
  ggsave(
    filename = file.path(output_dir, "combined_profiles_targets_vs_nonTargets.pdf"),
    plot = combined_plot,
    width = 10,
    height = 12
  )
  
  # Print summary statistics
  cat("\nSummary Statistics:\n")
  for (set in c("Target Genes", "Non-target Genes")) {
    cat(sprintf("\n=== %s ===\n", set))
    set_data <- all_plot_data[all_plot_data$Set == set, ]
    for (group in unique(set_data$Group)) {
      group_data <- set_data[set_data$Group == group, ]
      cat(sprintf("\n%s:\n", group))
      cat("Mean signal:", mean(group_data$Signal), "\n")
      cat("Max signal:", max(group_data$Signal), "\n")
      cat("Position of max signal:", group_data$Position[which.max(group_data$Signal)], "bp\n")
    }
  }
}

# Run the main function
if (!interactive()) {
  main()
} 