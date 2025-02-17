#!/usr/bin/env Rscript

#' Combined SMARCB1 binding profile comparisons
#' 
#' This script combines four different comparisons:
#' 1. Bivalent vs Non-bivalent genes
#' 2. Targeted vs Non-targeted genes
#' 3. Bivalent Targeted vs Bivalent Non-targeted genes
#' 4. Non-bivalent Targeted vs Non-bivalent Non-targeted genes
#' 
#' Input files:
#' - results/bigwig/BG1_CPM.bw: Control replicate 1 bigWig file
#' - results/bigwig/BG2_CPM.bw: Control replicate 2 bigWig file 
#' - results/bigwig/BG3_CPM.bw: Control replicate 3 bigWig file
#' - results/bigwig/BM3_CPM.bw: SMARCB1 ChIP-seq bigWig file
#' - data/gencode.vM10.annotation.gtf.gz: Gene annotations
#' - Gene_lists/bivalent/expressed_targeted_bivalent_NPCs_1000.csv
#' - Gene_lists/bivalent/expressed_targeted_non_bivalent_NPCs_1000.csv
#' - Gene_lists/bivalent/expressed_not_targeted_bivalent_NPCs_1000.csv
#' - Gene_lists/bivalent/expressed_not_targeted_non_bivalent_NPCs_1000.csv
#' - Gene_lists/targets/high_expression_targets2_1000.0.csv
#' - Gene_lists/targets/high_expression_no_targets_1000.0.csv
#'
#' Output files:
#' - results/metaprofiles_comparison_R_combined/combined_all_comparisons.pdf: 
#'   Combined plot showing all comparisons
#' - Console output: Summary statistics for all comparisons

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

# Create plotting data for individual profiles
create_plot_data <- function(bg_profile, bm_profile, group_name, x_pos) {
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

# Create comparison plot
create_comparison_plot <- function(plot_data, fc_data, title, subtitle = NULL, group1_n, group2_n) {
  # Add sample counts to subtitle
  if (!is.null(subtitle)) {
    subtitle <- sprintf("%s\n(%s: n=%d, %s: n=%d)", 
                       subtitle,
                       unique(plot_data$Group)[1], group1_n,
                       unique(plot_data$Group)[2], group2_n)
  } else {
    subtitle <- sprintf("%s: n=%d, %s: n=%d", 
                       unique(plot_data$Group)[1], group1_n,
                       unique(plot_data$Group)[2], group2_n)
  }
  
  # Create main profile plot
  p1 <- ggplot(plot_data, aes(x = Position, y = Signal, color = Condition, linetype = Group)) +
    geom_line() +
    geom_ribbon(aes(ymin = Signal - SE, ymax = Signal + SE, fill = Condition), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("BG" = "#1f77b4", "BM" = "#d62728")) +
    scale_fill_manual(values = c("BG" = "#1f77b4", "BM" = "#d62728")) +
    labs(title = title,
         subtitle = subtitle,
         x = "Distance from TSS (bp)",
         y = "Average RPKM") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # Create fold change plot
  p2 <- ggplot(fc_data, aes(x = Position, y = FoldChange, color = Group)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Bivalent" = "#ff7f0e", "Non-bivalent" = "#2ca02c",
                                 "Targeted" = "#ff7f0e", "No Targeted" = "#2ca02c",
                                 "Bivalent Targeted" = "#ff7f0e", "Bivalent Not Targeted" = "#2ca02c",
                                 "Non-bivalent Targeted" = "#ff7f0e", "Non-bivalent Not Targeted" = "#2ca02c")) +
    labs(title = "Log2 Fold Change (BM/BG)",
         x = "Distance from TSS (bp)",
         y = "Log2 Fold Change") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1 / p2
}

# Print summary statistics
print_summary_stats <- function(name1, bg_profile1, bm_profile1, fc1,
                              name2, bg_profile2, bm_profile2, fc2) {
  cat("\nSummary Statistics:\n")
  cat(sprintf("\n%s:\n", name1))
  cat("Mean BG signal:", mean(colMeans(bg_profile1)), "\n")
  cat("Mean BM signal:", mean(colMeans(bm_profile1)), "\n")
  cat("Mean log2 fold change:", mean(fc1), "\n")
  
  cat(sprintf("\n%s:\n", name2))
  cat("Mean BG signal:", mean(colMeans(bg_profile2)), "\n")
  cat("Mean BM signal:", mean(colMeans(bm_profile2)), "\n")
  cat("Mean log2 fold change:", mean(fc2), "\n")
}

# Main execution
main <- function() {
  input_dir_bw <- "results/bigwig"
  input_dir_gene_lists <- "Gene_lists"
  output_dir <- "results/metaprofiles_comparison_R_combined"
  
  # Set paths
  bg_files <- c(file.path(input_dir_bw, "BG1_CPM.bw"),
                file.path(input_dir_bw, "BG2_CPM.bw"),
                file.path(input_dir_bw, "BG3_CPM.bw"))
  bm_file <- file.path(input_dir_bw, "BM3_CPM.bw")
  gtf_file <- "data/gencode.vM10.annotation.gtf.gz"
  
  files_suffix <- "_all.csv"
  # files_suffix <- "_1000.csv"

  # Define gene lists
  bivalent_genes <- paste0(file.path(input_dir_gene_lists, "bivalent/expressed_targeted_bivalent_NPCs"), files_suffix)
  non_bivalent_genes <- paste0(file.path(input_dir_gene_lists, "bivalent/expressed_targeted_non_bivalent_NPCs"), files_suffix)

  # targeted_genes <- paste0(file.path(input_dir_gene_lists, "targets/high_expression_targets2"), files_suffix)
  # no_targeted_genes <- paste0(file.path(input_dir_gene_lists, "targets/high_expression_no_targets"), files_suffix)

  targeted_genes <- file.path(input_dir_gene_lists, "targets/all_targets_final.csv")
  no_targeted_genes <- file.path(input_dir_gene_lists, "targets/all_no_targets_final.csv")
  
  bivalent_targeted_genes <- paste0(file.path(input_dir_gene_lists, "bivalent/expressed_targeted_bivalent_NPCs"), files_suffix)
  bivalent_no_targeted_genes <- paste0(file.path(input_dir_gene_lists, "bivalent/expressed_not_targeted_bivalent_NPCs"), files_suffix)
  nonbivalent_targeted_genes <- paste0(file.path(input_dir_gene_lists, "bivalent/expressed_targeted_non_bivalent_NPCs"), files_suffix)
  nonbivalent_no_targeted_genes <- paste0(file.path(input_dir_gene_lists, "bivalent/expressed_not_targeted_non_bivalent_NPCs"), files_suffix)
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE)
  
  # Import bigWig data
  cat("Importing BG files...\n")
  bg_data <- import_bigwig(bg_files)
  cat("Importing BM file...\n")
  bm_data <- import_bigwig(list(bm_file))
  
  # Process all gene lists
  cat("Processing gene lists...\n")
  regions_list <- list(
    bivalent = create_tss_regions(bivalent_genes, gtf_file),
    non_bivalent = create_tss_regions(non_bivalent_genes, gtf_file),
    targeted = create_tss_regions(targeted_genes, gtf_file),
    no_targeted = create_tss_regions(no_targeted_genes, gtf_file),
    bivalent_targeted = create_tss_regions(bivalent_targeted_genes, gtf_file),
    bivalent_no_targeted = create_tss_regions(bivalent_no_targeted_genes, gtf_file),
    nonbivalent_targeted = create_tss_regions(nonbivalent_targeted_genes, gtf_file),
    nonbivalent_no_targeted = create_tss_regions(nonbivalent_no_targeted_genes, gtf_file)
  )
  
  # Calculate all profiles
  cat("Calculating profiles...\n")
  profiles <- list()
  for (name in names(regions_list)) {
    profiles[[paste0(name, "_bg")]] <- calculate_metaprofile(bg_data, regions_list[[name]])
    profiles[[paste0(name, "_bm")]] <- calculate_metaprofile(bm_data, regions_list[[name]])
  }
  
  # Set up x-axis positions
  x_pos <- seq(-2500, 2500, length.out = ncol(profiles[[1]]))
  
  # Calculate sample sizes
  sample_sizes <- list()
  for (name in names(regions_list)) {
    sample_sizes[[name]] <- length(regions_list[[name]])
  }
  
  # Create all comparison plots
  cat("Creating plots...\n")
  
  # 1. Bivalent vs Non-bivalent comparison
  plot_data1 <- rbind(
    create_plot_data(profiles$bivalent_bg, profiles$bivalent_bm, "Bivalent", x_pos),
    create_plot_data(profiles$non_bivalent_bg, profiles$non_bivalent_bm, "Non-bivalent", x_pos)
  )
  
  fc_data1 <- data.frame(
    Position = rep(x_pos, 2),
    FoldChange = c(
      calculate_fold_change(profiles$bivalent_bg, profiles$bivalent_bm),
      calculate_fold_change(profiles$non_bivalent_bg, profiles$non_bivalent_bm)
    ),
    Group = rep(c("Bivalent", "Non-bivalent"), each = length(x_pos))
  )
  
  # 2. Targeted vs Non-targeted comparison
  plot_data2 <- rbind(
    create_plot_data(profiles$targeted_bg, profiles$targeted_bm, "Targeted", x_pos),
    create_plot_data(profiles$no_targeted_bg, profiles$no_targeted_bm, "No Targeted", x_pos)
  )
  
  fc_data2 <- data.frame(
    Position = rep(x_pos, 2),
    FoldChange = c(
      calculate_fold_change(profiles$targeted_bg, profiles$targeted_bm),
      calculate_fold_change(profiles$no_targeted_bg, profiles$no_targeted_bm)
    ),
    Group = rep(c("Targeted", "No Targeted"), each = length(x_pos))
  )
  
  # 3. Bivalent Targeted vs Non-targeted comparison
  plot_data3 <- rbind(
    create_plot_data(profiles$bivalent_targeted_bg, profiles$bivalent_targeted_bm, "Bivalent Targeted", x_pos),
    create_plot_data(profiles$bivalent_no_targeted_bg, profiles$bivalent_no_targeted_bm, "Bivalent Not Targeted", x_pos)
  )
  
  fc_data3 <- data.frame(
    Position = rep(x_pos, 2),
    FoldChange = c(
      calculate_fold_change(profiles$bivalent_targeted_bg, profiles$bivalent_targeted_bm),
      calculate_fold_change(profiles$bivalent_no_targeted_bg, profiles$bivalent_no_targeted_bm)
    ),
    Group = rep(c("Bivalent Targeted", "Bivalent Not Targeted"), each = length(x_pos))
  )
  
  # 4. Non-bivalent Targeted vs Non-targeted comparison
  plot_data4 <- rbind(
    create_plot_data(profiles$nonbivalent_targeted_bg, profiles$nonbivalent_targeted_bm, "Non-bivalent Targeted", x_pos),
    create_plot_data(profiles$nonbivalent_no_targeted_bg, profiles$nonbivalent_no_targeted_bm, "Non-bivalent Not Targeted", x_pos)
  )
  
  fc_data4 <- data.frame(
    Position = rep(x_pos, 2),
    FoldChange = c(
      calculate_fold_change(profiles$nonbivalent_targeted_bg, profiles$nonbivalent_targeted_bm),
      calculate_fold_change(profiles$nonbivalent_no_targeted_bg, profiles$nonbivalent_no_targeted_bm)
    ),
    Group = rep(c("Non-bivalent Targeted", "Non-bivalent Not Targeted"), each = length(x_pos))
  )
  
  # Create combined plots with sample sizes
  p1 <- create_comparison_plot(plot_data1, fc_data1, "SMARCB1 binding around TSS", 
                             "Bivalent vs Non-bivalent",
                             sample_sizes$bivalent, sample_sizes$non_bivalent)
  
  p2 <- create_comparison_plot(plot_data2, fc_data2, "SMARCB1 binding around TSS", 
                             "Targeted vs Non-targeted",
                             sample_sizes$targeted, sample_sizes$no_targeted)
  
  p3 <- create_comparison_plot(plot_data3, fc_data3, "SMARCB1 binding around TSS", 
                             "Bivalent Targeted vs Non-targeted",
                             sample_sizes$bivalent_targeted, sample_sizes$bivalent_no_targeted)
  
  p4 <- create_comparison_plot(plot_data4, fc_data4, "SMARCB1 binding around TSS", 
                             "Non-bivalent Targeted vs Non-targeted",
                             sample_sizes$nonbivalent_targeted, sample_sizes$nonbivalent_no_targeted)
  
  # Combine all plots
  combined_plot <- (p1 | p2) / (p3 | p4)
  
  # Save combined plot
  ggsave(
    filename = file.path(output_dir, "combined_all_comparisons.pdf"),
    plot = combined_plot,
    width = 20,
    height = 24
  )
  
  # Print summary statistics for all comparisons with sample sizes
  cat("\n=== Bivalent vs Non-bivalent ===\n")
  cat(sprintf("Bivalent genes: n=%d\n", sample_sizes$bivalent))
  cat(sprintf("Non-bivalent genes: n=%d\n", sample_sizes$non_bivalent))
  print_summary_stats("Bivalent", profiles$bivalent_bg, profiles$bivalent_bm,
                     calculate_fold_change(profiles$bivalent_bg, profiles$bivalent_bm),
                     "Non-bivalent", profiles$non_bivalent_bg, profiles$non_bivalent_bm,
                     calculate_fold_change(profiles$non_bivalent_bg, profiles$non_bivalent_bm))
  
  cat("\n=== Targeted vs Non-targeted ===\n")
  cat(sprintf("Targeted genes: n=%d\n", sample_sizes$targeted))
  cat(sprintf("Non-targeted genes: n=%d\n", sample_sizes$no_targeted))
  print_summary_stats("Targeted", profiles$targeted_bg, profiles$targeted_bm,
                     calculate_fold_change(profiles$targeted_bg, profiles$targeted_bm),
                     "Non-targeted", profiles$no_targeted_bg, profiles$no_targeted_bm,
                     calculate_fold_change(profiles$no_targeted_bg, profiles$no_targeted_bm))
  
  cat("\n=== Bivalent Targeted vs Non-targeted ===\n")
  cat(sprintf("Bivalent targeted genes: n=%d\n", sample_sizes$bivalent_targeted))
  cat(sprintf("Bivalent non-targeted genes: n=%d\n", sample_sizes$bivalent_no_targeted))
  print_summary_stats("Bivalent Targeted", profiles$bivalent_targeted_bg, profiles$bivalent_targeted_bm,
                     calculate_fold_change(profiles$bivalent_targeted_bg, profiles$bivalent_targeted_bm),
                     "Bivalent Non-targeted", profiles$bivalent_no_targeted_bg, profiles$bivalent_no_targeted_bm,
                     calculate_fold_change(profiles$bivalent_no_targeted_bg, profiles$bivalent_no_targeted_bm))
  
  cat("\n=== Non-bivalent Targeted vs Non-targeted ===\n")
  cat(sprintf("Non-bivalent targeted genes: n=%d\n", sample_sizes$nonbivalent_targeted))
  cat(sprintf("Non-bivalent non-targeted genes: n=%d\n", sample_sizes$nonbivalent_no_targeted))
  print_summary_stats("Non-bivalent Targeted", profiles$nonbivalent_targeted_bg, profiles$nonbivalent_targeted_bm,
                     calculate_fold_change(profiles$nonbivalent_targeted_bg, profiles$nonbivalent_targeted_bm),
                     "Non-bivalent Non-targeted", profiles$nonbivalent_no_targeted_bg, profiles$nonbivalent_no_targeted_bm,
                     calculate_fold_change(profiles$nonbivalent_no_targeted_bg, profiles$nonbivalent_no_targeted_bm))
}

# Run the main function
if (!interactive()) {
  main()
} 