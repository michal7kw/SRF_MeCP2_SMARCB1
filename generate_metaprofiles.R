#!/usr/bin/env Rscript

# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(data.table)  # For faster data operations
library(parallel)    # For parallel processing

# Function to get bigWig binning info
get_bigwig_info <- function(bw_file) {
  # Import bigWig file
  bw <- BigWigFile(bw_file)
  
  # Get chromosome sizes
  seqinfo <- seqinfo(bw)
  
  # Create a GRanges object covering all chromosomes
  all_ranges <- GRanges(
    seqnames = seqnames(seqinfo),
    ranges = IRanges(start = 1, end = seqlengths(seqinfo))
  )
  
  # Get summary statistics for the entire genome
  stats <- summary(bw, which = all_ranges, size = 1000)  # using 1000bp windows for summary
  
  # Calculate statistics from non-empty regions
  valid_widths <- stats$width[!is.na(stats$mean)]
  
  list(
    seqinfo = seqinfo,
    mean_bin_size = mean(valid_widths),
    min_bin_size = min(valid_widths),
    max_bin_size = max(valid_widths),
    total_regions = length(valid_widths)
  )
}

# Function to import bigWig files and calculate mean signal
import_bigwig <- function(file_paths, bin_size = 10) {
  # Use parallel processing for multiple files
  if (length(file_paths) > 1) {
    # Determine number of cores (use half of available cores)
    n_cores <- max(1, parallel::detectCores() %/% 2)
    
    # Create cluster
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl))
    
    # Export required packages to cluster
    clusterEvalQ(cl, {
      library(rtracklayer)
      library(GenomicRanges)
    })
    
    # Process files in parallel
    bw_list <- parLapply(cl, file_paths, function(f) {
      tryCatch({
        cat(sprintf("\nProcessing file: %s\n", basename(f)))
        bw <- BigWigFile(f)
        import.bw(bw, as = "RleList")  # Import as RleList for memory efficiency
      }, error = function(e) {
        stop(sprintf("Error processing file %s: %s", basename(f), e$message))
      })
    })
    
    # Calculate mean coverage
    if (length(bw_list) > 1) {
      # Convert to common seqinfo
      seqinfo <- seqinfo(bw_list[[1]])
      bw_list <- lapply(bw_list, function(x) {
        seqinfo(x) <- seqinfo
        x
      })
      
      # Calculate mean coverage
      mean_coverage <- Reduce("+", bw_list) / length(bw_list)
    } else {
      mean_coverage <- bw_list[[1]]
    }
  } else {
    # Single file processing
    mean_coverage <- import.bw(file_paths[[1]], as = "RleList")
  }
  
  return(mean_coverage)
}

# Function to process gene lists and create TSS regions
create_tss_regions <- function(gene_list_file, gtf_file, upstream = 2500, downstream = 2500) {
  # Read gene list
  genes <- fread(gene_list_file, header = FALSE)$V1
  
  # Import GTF file
  gtf <- rtracklayer::import(gtf_file)
  
  # Filter for genes and matching gene names
  gene_records <- gtf[gtf$type == "gene" & gtf$gene_name %in% genes]
  
  # Create TSS regions based on strand
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
  
  # Sort regions for faster operations
  sort(tss_regions)
}

# Function to calculate metaprofile matrix
calculate_metaprofile <- function(coverage, regions, bins = 100) {
  # Initialize matrix
  profile_matrix <- matrix(0, nrow = length(regions), ncol = bins)
  
  # Get available chromosomes in coverage
  available_chrs <- names(coverage)
  
  # Process in chunks for memory efficiency
  chunk_size <- 1000
  n_chunks <- ceiling(length(regions) / chunk_size)
  
  for (chunk in seq_len(n_chunks)) {
    # Get chunk of regions
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, length(regions))
    chunk_regions <- regions[start_idx:end_idx]
    
    # Process each chromosome separately
    for (chr in unique(seqnames(chunk_regions))) {
      # Skip if chromosome not in coverage
      if (!as.character(chr) %in% available_chrs) {
        warning(sprintf("Chromosome %s not found in coverage data, skipping...", chr))
        next
      }
      
      # Get regions for this chromosome
      chr_regions <- chunk_regions[seqnames(chunk_regions) == chr]
      if (length(chr_regions) == 0) next
      
      # Calculate scores for chunk
      chunk_scores <- tryCatch({
        Views(coverage[[as.character(chr)]], ranges(chr_regions))
      }, error = function(e) {
        warning(sprintf("Error processing chromosome %s: %s", chr, e$message))
        return(NULL)
      })
      
      if (is.null(chunk_scores)) next
      
      # Convert to matrix and resize to desired number of bins
      chunk_matrix <- vapply(chunk_scores, function(x) {
        # Interpolate to desired number of bins
        x_len <- length(x)
        if (x_len == 0) return(numeric(bins))
        if (x_len == bins) return(as.numeric(x))
        # Use approx for interpolation
        approx(seq_len(x_len), as.numeric(x), n = bins)$y
      }, numeric(bins))
      
      # Get indices in the original chunk
      chr_indices <- which(seqnames(chunk_regions) == chr)
      # Store in profile matrix
      profile_matrix[start_idx + chr_indices - 1,] <- t(chunk_matrix)
    }
  }
  
  return(profile_matrix)
}

# Main execution
main <- function() {
  # Set paths
  bg_files <- c("results/bigwig/BG1_RPKM.bw",
                "results/bigwig/BG2_RPKM.bw",
                "results/bigwig/BG3_RPKM.bw")
  bm_file <- "results/bigwig/BM3_RPKM.bw"
  gtf_file <- "data/gencode.vM10.annotation.gtf.gz"
  gene_lists <- c("Gene_lists/by_regulation_type/enriched_down_regulated.csv",
                  "Gene_lists/by_regulation_type/enriched_not_disregulated.csv",
                  "Gene_lists/by_regulation_type/enriched_up_regulated.csv")
  
  # Create output directory
  dir.create("results/metaprofiles_R", showWarnings = FALSE)
  
  # Import bigWig data
  cat("Importing BG files...\n")
  bg_data <- import_bigwig(bg_files)
  cat("Importing BM file...\n")
  bm_data <- import_bigwig(list(bm_file))
  
  # Process each gene list
  for (gene_list in gene_lists) {
    cat(sprintf("\nProcessing gene list: %s\n", basename(gene_list)))
    
    # Get base name
    base_name <- tools::file_path_sans_ext(basename(gene_list))
    
    # Create TSS regions
    tss_regions <- create_tss_regions(gene_list, gtf_file)
    
    # Calculate profiles
    bg_profile <- calculate_metaprofile(bg_data, tss_regions)
    bm_profile <- calculate_metaprofile(bm_data, tss_regions)
    
    # Calculate mean and standard error
    bg_mean <- colMeans(bg_profile)
    bm_mean <- colMeans(bm_profile)
    bg_se <- apply(bg_profile, 2, function(x) sd(x)/sqrt(length(x)))
    bm_se <- apply(bm_profile, 2, function(x) sd(x)/sqrt(length(x)))
    
    # Create plotting data
    x_pos <- seq(-2500, 2500, length.out = length(bg_mean))
    plot_data <- data.frame(
      Position = rep(x_pos, 2),
      Signal = c(bg_mean, bm_mean),
      SE = c(bg_se, bm_se),
      Group = rep(c("BG", "BM"), each = length(x_pos))
    )
    
    # Create plot
    p <- ggplot(plot_data, aes(x = Position, y = Signal, color = Group)) +
      geom_line() +
      geom_ribbon(aes(ymin = Signal - SE, ymax = Signal + SE, fill = Group),
                  alpha = 0.2, color = NA) +
      scale_color_manual(values = c("BG" = "#1f77b4", "BM" = "#d62728")) +
      scale_fill_manual(values = c("BG" = "#1f77b4", "BM" = "#d62728")) +
      labs(title = paste("SMARCB1 binding around TSS -", base_name),
           x = "Distance from TSS (bp)",
           y = "Average RPKM") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Save plot
    ggsave(
      filename = file.path("results/metaprofiles_R", paste0(base_name, "_profile.pdf")),
      plot = p,
      width = 8,
      height = 6
    )
  }
}

# Run the main function
if (!interactive()) {
  main()
} 