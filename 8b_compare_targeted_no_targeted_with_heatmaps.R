#!/usr/bin/env Rscript

# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(data.table)
library(parallel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

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
create_tss_regions <- function(gene_list_file, gtf_file, upstream = 5000, downstream = 5000) {
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
calculate_metaprofile <- function(coverage, regions, bins = 100, no_chunking = FALSE) {
  profile_matrix <- matrix(0, nrow = length(regions), ncol = bins)
  available_chrs <- names(coverage)
  
  # Print debug info
  cat("\nProcessing", length(regions), "regions\n")
  cat("Available chromosomes:", paste(available_chrs, collapse=", "), "\n")
  
  if (no_chunking) {
    # Process all regions at once without chunking
    for (chr in unique(seqnames(regions))) {
      if (!as.character(chr) %in% available_chrs) {
        warning(sprintf("Chromosome %s not found in coverage data, skipping...", chr))
        next
      }
      
      chr_regions <- regions[seqnames(regions) == chr]
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
        # Use mean for binning to preserve signal intensity
        values <- as.numeric(x)
        bin_size <- x_len / bins
        result <- numeric(bins)
        for (i in 1:bins) {
          start_idx <- floor((i-1) * bin_size) + 1
          end_idx <- floor(i * bin_size)
          result[i] <- mean(values[start_idx:end_idx])
        }
        return(result)
      }, numeric(bins))
      
      chr_indices <- which(seqnames(regions) == chr)
      
      # Add dimension checks and debugging
      tryCatch({
        if (length(chr_indices) > 0) {
          if (length(chr_indices) != ncol(t(chunk_matrix))) {
            warning(sprintf("Dimension mismatch: chr_indices length (%d) != chunk_matrix cols (%d)", 
                          length(chr_indices), ncol(t(chunk_matrix))))
            # Ensure dimensions match by taking only what we can fit
            n_cols <- min(length(chr_indices), ncol(t(chunk_matrix)))
            profile_matrix[chr_indices[1:n_cols], ] <- t(chunk_matrix)[1:n_cols, ]
          } else {
            profile_matrix[chr_indices, ] <- t(chunk_matrix)
          }
        }
      }, error = function(e) {
        warning(sprintf("Error assigning chunk matrix: %s", e$message))
      })
    }
  } else {
    # Process in chunks (original implementation)
    chunk_size <- 10000
    n_chunks <- ceiling(length(regions) / chunk_size)
    
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
          # Use mean for binning to preserve signal intensity
          values <- as.numeric(x)
          bin_size <- x_len / bins
          result <- numeric(bins)
          for (i in 1:bins) {
            start_idx <- floor((i-1) * bin_size) + 1
            end_idx <- floor(i * bin_size)
            result[i] <- mean(values[start_idx:end_idx])
          }
          return(result)
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
  }
  return(profile_matrix)
}

# Modify the create_heatmap function to fix the legend issue
create_heatmap <- function(chip_matrix, nscm_matrix, title) {
  # Check if matrices are empty or all zero
  if (nrow(chip_matrix) == 0 || nrow(nscm_matrix) == 0) {
    stop("Empty matrix provided to create_heatmap")
  }
  
  # Calculate mean profiles and correlation
  chip_profile <- colMeans(chip_matrix)
  nscm_profile <- colMeans(nscm_matrix)
  correlation <- cor(chip_profile, nscm_profile)
  
  # Create profile annotations with improved styling
  ha_profiles <- HeatmapAnnotation(
    "Signal Profile" = anno_lines(
      cbind(chip_profile, nscm_profile),
      gp = gpar(col = c("#E41A1C", "#377EB8"), lwd = 2),
      height = unit(5, "cm")
    ),
    which = "column",
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 12)
  )
  
  # Create the main heatmaps with improved styling
  chip_hm <- Heatmap(chip_matrix,
                     name = "ChIP Signal",
                     column_title = "TSS ± 5kb",
                     cluster_rows = FALSE,
                     show_row_names = FALSE,
                     use_raster = TRUE,
                     raster_quality = 2,
                     top_annotation = ha_profiles,
                     col = colorRamp2(c(-4, 0, 4),
                                    c("blue", "white", "red")),
                     column_title_gp = gpar(fontsize = 14),
                     heatmap_legend_param = list(
                       title_gp = gpar(fontsize = 12),
                       labels_gp = gpar(fontsize = 10),
                       title_position = "leftcenter-rot"
                     ))
  
  nscm_hm <- Heatmap(nscm_matrix,
                     name = "NSCM Signal",
                     column_title = "TSS ± 5kb",
                     cluster_rows = FALSE,
                     show_row_names = FALSE,
                     use_raster = TRUE,
                     raster_quality = 2,
                     col = colorRamp2(c(-4, 0, 4),
                                    c("blue", "white", "red")),
                     column_title_gp = gpar(fontsize = 14),
                     heatmap_legend_param = list(
                       title_gp = gpar(fontsize = 12),
                       labels_gp = gpar(fontsize = 10),
                       title_position = "leftcenter-rot"
                     ))
  
  # Create a legend for the profile lines
  lgd = Legend(
    labels = c("ChIP-seq", "NSCM"),
    title = "Signal Type",
    type = "lines",
    legend_gp = gpar(col = c("#E41A1C", "#377EB8"), lwd = 2)
  )
  
  # Add correlation information to the title
  title <- sprintf("%s\nPearson correlation: %.2f", title, correlation)
  
  # Return the heatmap with additional legend
  ht_list <- chip_hm + nscm_hm
  draw(ht_list, annotation_legend_list = list(lgd))
}

# Modify the create_heatmap_matrix function to improve sorting and scaling
create_heatmap_matrix <- function(coverage, regions, bins = 100) {
  # Check if we have valid regions
  if (length(regions) == 0) {
    stop("No regions provided to create_heatmap_matrix")
  }
  
  profile_matrix <- calculate_metaprofile(coverage, regions, bins, no_chunking = TRUE)
  profile_matrix[!is.finite(profile_matrix)] <- 0
  
  # Calculate TSS-centered signal for sorting
  center_idx <- bins/2
  window_size <- 10  # bins around TSS to consider
  tss_signal <- rowMeans(profile_matrix[, (center_idx-window_size):(center_idx+window_size)])
  
  # Sort by TSS signal strength
  sorted_indices <- order(tss_signal, decreasing = TRUE)
  profile_matrix <- profile_matrix[sorted_indices, ]
  
  # Scale the matrix
  scaled_matrix <- t(scale(t(profile_matrix)))
  scaled_matrix[scaled_matrix > 4] <- 4
  scaled_matrix[scaled_matrix < -4] <- -4
  
  return(list(
    matrix = scaled_matrix,
    indices = sorted_indices
  ))
}

# Modify the create_comparison_heatmap function to handle different row counts
create_comparison_heatmap <- function(chip_targeted, nscm_targeted, 
                                    chip_nontargeted, nscm_nontargeted) {
  # Create two separate heatmap lists for targeted and non-targeted genes
  
  # Calculate correlations
  cor_targeted <- cor(colMeans(chip_targeted$matrix), colMeans(nscm_targeted$matrix))
  cor_nontargeted <- cor(colMeans(chip_nontargeted$matrix), colMeans(nscm_nontargeted$matrix))
  
  # Create profile annotations
  create_annotation <- function(chip_mat, nscm_mat, name) {
    chip_mean <- colMeans(chip_mat$matrix)
    nscm_mean <- colMeans(nscm_mat$matrix)
    
    # Scale means to 0-1 range
    scale_to_01 <- function(x) {
      (x - min(x)) / (max(x) - min(x))
    }
    chip_mean_scaled <- scale_to_01(chip_mean)
    nscm_mean_scaled <- scale_to_01(nscm_mean)
    
    HeatmapAnnotation(
      "Signal" = anno_lines(
        cbind(chip_mean_scaled, nscm_mean_scaled),
        gp = gpar(col = c("#E41A1C", "#377EB8"), lwd = 2),
        height = unit(3, "cm")
      ),
      show_annotation_name = FALSE
    )
  }
  
  # Create legends
  signal_lgd <- Legend(
    labels = c("ChIP-seq", "NSCM"),
    title = "Signal Type",
    type = "lines",
    legend_gp = gpar(col = c("#E41A1C", "#377EB8"), lwd = 2)
  )
  
  # Save targeted heatmaps
  pdf("results/metaprofiles_comparison_R/targeted_heatmaps.pdf", 
      width = 8, height = 12)
  
  # Draw targeted heatmaps
  ht_targeted <- Heatmap(chip_targeted$matrix,
                        name = "ChIP Signal",
                        column_title = "Targeted Genes (TSS ± 5kb)",
                        column_title_gp = gpar(fontsize = 14),
                        top_annotation = create_annotation(chip_targeted, nscm_targeted, "targeted"),
                        cluster_rows = FALSE,
                        show_row_names = FALSE,
                        use_raster = TRUE,
                        raster_quality = 2,
                        col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))) +
    
    Heatmap(nscm_targeted$matrix,
            name = "NSCM Signal",
            column_title = sprintf("Correlation: %.2f", cor_targeted),
            cluster_rows = FALSE,
            show_row_names = FALSE,
            use_raster = TRUE,
            raster_quality = 2,
            col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")))
  
  draw(ht_targeted, 
       annotation_legend_list = list(signal_lgd),
       padding = unit(c(2, 2, 2, 2), "cm"))
  
  dev.off()
  
  # Save non-targeted heatmaps
  pdf("results/metaprofiles_comparison_R/non_targeted_heatmaps.pdf", 
      width = 8, height = 12)
  
  # Draw non-targeted heatmaps
  ht_nontargeted <- Heatmap(chip_nontargeted$matrix,
                           name = "ChIP Signal",
                           column_title = "Non-targeted Genes (TSS ± 5kb)",
                           column_title_gp = gpar(fontsize = 14),
                           top_annotation = create_annotation(chip_nontargeted, nscm_nontargeted, "non-targeted"),
                           cluster_rows = FALSE,
                           show_row_names = FALSE,
                           use_raster = TRUE,
                           raster_quality = 2,
                           col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))) +
    
    Heatmap(nscm_nontargeted$matrix,
            name = "NSCM Signal",
            column_title = sprintf("Correlation: %.2f", cor_nontargeted),
            cluster_rows = FALSE,
            show_row_names = FALSE,
            use_raster = TRUE,
            raster_quality = 2,
            col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")))
  
  draw(ht_nontargeted, 
       annotation_legend_list = list(signal_lgd),
       padding = unit(c(2, 2, 2, 2), "cm"))
  
  dev.off()
  
  # Print summary statistics
  cat("\nSummary Statistics:\n")
  cat(sprintf("Targeted genes (n=%d):\n", nrow(chip_targeted$matrix)))
  cat(sprintf("  Correlation: %.3f\n", cor_targeted))
  cat(sprintf("Non-targeted genes (n=%d):\n", nrow(chip_nontargeted$matrix)))
  cat(sprintf("  Correlation: %.3f\n", cor_nontargeted))
}

# After the main function definition, add this helper function
check_cached_data <- function(file_path) {
  if (file.exists(file_path)) {
    cat(sprintf("Loading cached data from %s\n", file_path))
    return(readRDS(file_path))
  }
  return(NULL)
}

# Modify the main function to include better error handling
main <- function() {
  # Create cache directory
  cache_dir <- "results/metaprofiles_comparison_R/cache"
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Set paths
  bg_files <- c("results/bigwig/BG1_CPM.bw",
                "results/bigwig/BG2_CPM.bw",
                "results/bigwig/BG3_CPM.bw")
  bm_files <- c("results/bigwig/BM3_CPM.bw")
  nscm_files <- list.files("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_1b/bigwig/",
                          pattern = "NSCM.*\\.bw$", 
                          full.names = TRUE)
  
  gtf_file <- "data/gencode.vM10.basic.annotation.gtf.gz"
  targeted_genes <- "Gene_lists/targets/all_targets_final.csv"
  no_targeted_genes <- "Gene_lists/targets/all_no_targets_mm10.csv"
  
  # Check for cached coverage data
  chip_data <- check_cached_data(file.path(cache_dir, "chip_coverage.rds"))
  if (is.null(chip_data)) {
    cat("Importing BG and BM files...\n")
    all_chip_files <- c(bg_files, bm_files)
    chip_data <- import_bigwig(all_chip_files)
    saveRDS(chip_data, file.path(cache_dir, "chip_coverage.rds"))
  }

  nscm_data <- check_cached_data(file.path(cache_dir, "nscm_coverage.rds"))
  if (is.null(nscm_data)) {
    cat("Importing NSCM files...\n")
    nscm_data <- import_bigwig(nscm_files)
    saveRDS(nscm_data, file.path(cache_dir, "nscm_coverage.rds"))
  }

  # Check for cached regions
  regions_cache <- check_cached_data(file.path(cache_dir, "regions.rds"))
  if (is.null(regions_cache)) {
    cat("Processing regions...\n")
    targeted_regions <- create_tss_regions(targeted_genes, gtf_file)
    no_targeted_regions <- create_tss_regions(no_targeted_genes, gtf_file)
    regions_cache <- list(
      targeted = targeted_regions,
      no_targeted = no_targeted_regions
    )
    saveRDS(regions_cache, file.path(cache_dir, "regions.rds"))
  } else {
    targeted_regions <- regions_cache$targeted
    no_targeted_regions <- regions_cache$no_targeted
  }

  # Check for cached profiles
  profiles_cache <- check_cached_data(file.path(cache_dir, "profiles.rds"))
  if (is.null(profiles_cache)) {
    cat("Calculating profiles...\n")
    targeted_chip_profile <- calculate_metaprofile(chip_data, targeted_regions)
    no_targeted_chip_profile <- calculate_metaprofile(chip_data, no_targeted_regions)
    targeted_nscm_profile <- calculate_metaprofile(nscm_data, targeted_regions)
    no_targeted_nscm_profile <- calculate_metaprofile(nscm_data, no_targeted_regions)
    
    profiles_cache <- list(
      targeted_chip = targeted_chip_profile,
      no_targeted_chip = no_targeted_chip_profile,
      targeted_nscm = targeted_nscm_profile,
      no_targeted_nscm = no_targeted_nscm_profile
    )
    saveRDS(profiles_cache, file.path(cache_dir, "profiles.rds"))
  } else {
    targeted_chip_profile <- profiles_cache$targeted_chip
    no_targeted_chip_profile <- profiles_cache$no_targeted_chip
    targeted_nscm_profile <- profiles_cache$targeted_nscm
    no_targeted_nscm_profile <- profiles_cache$no_targeted_nscm
  }

  # Create plotting data for metaprofiles
  create_plot_data <- function(profile_matrix, group_name) {
    mean_profile <- colMeans(profile_matrix)
    se_profile <- apply(profile_matrix, 2, function(x) sd(x)/sqrt(length(x)))
    
    data.frame(
      Position = seq(-2500, 2500, length.out = 100),
      Signal = mean_profile,
      SE = se_profile,
      Group = group_name
    )
  }
  
  # Create plot data
  plot_data_targeted <- create_plot_data(targeted_chip_profile, "Targeted")
  plot_data_no_targeted <- create_plot_data(no_targeted_chip_profile, "No Targeted")
  
  # Create metaprofile plots
  create_profile_plot <- function(plot_data, title) {
    ggplot(plot_data, aes(x = Position, y = Signal)) +
      geom_line(color = "#1f77b4") +
      geom_ribbon(aes(ymin = Signal - SE, ymax = Signal + SE), 
                 alpha = 0.2, fill = "#1f77b4") +
      labs(title = title,
           x = "Distance from TSS (bp)",
           y = "Average Signal") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  p1 <- create_profile_plot(plot_data_targeted, "Targeted Genes")
  p2 <- create_profile_plot(plot_data_no_targeted, "Non-targeted Genes")
  
  # Create heatmap data
  targeted_chip_data <- create_heatmap_matrix(chip_data, targeted_regions)
  targeted_nscm_data <- create_heatmap_matrix(nscm_data, targeted_regions)
  no_targeted_chip_data <- create_heatmap_matrix(chip_data, no_targeted_regions)
  no_targeted_nscm_data <- create_heatmap_matrix(nscm_data, no_targeted_regions)
  
  # Create output directory
  dir.create("results/metaprofiles_comparison_R", showWarnings = FALSE)
  
  # Save combined visualization
  pdf("results/metaprofiles_comparison_R/combined_comparison.pdf", 
      width = 16, height = 12)
  
  create_comparison_heatmap(
    targeted_chip_data, targeted_nscm_data,
    no_targeted_chip_data, no_targeted_nscm_data
  )
  
  dev.off()
  
  # Save as PNG with higher resolution
  png("results/metaprofiles_comparison_R/combined_comparison.png",
      width = 1600, height = 1200, res = 150)
  
  create_comparison_heatmap(
    targeted_chip_data, targeted_nscm_data,
    no_targeted_chip_data, no_targeted_nscm_data
  )
  
  dev.off()
  
  # Save heatmap matrices in cache as well
  heatmap_matrices <- list(
    targeted_chip = targeted_chip_data$matrix,
    targeted_nscm = targeted_nscm_data$matrix,
    no_targeted_chip = no_targeted_chip_data$matrix,
    no_targeted_nscm = no_targeted_nscm_data$matrix
  )
  saveRDS(heatmap_matrices, file.path(cache_dir, "heatmap_matrices.rds"))
}

# Add a cleanup function to remove cache if needed
cleanup_cache <- function() {
  cache_dir <- "results/metaprofiles_comparison_R/cache"
  if (dir.exists(cache_dir)) {
    unlink(cache_dir, recursive = TRUE)
    cat("Cache directory cleaned up\n")
  }
}

# Run the main function
if (!interactive()) {
  # Add argument parsing for cache control
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0 && args[1] == "--clean") {
    cleanup_cache()
  }
  main()
} 