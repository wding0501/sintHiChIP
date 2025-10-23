#' sintHiChIP Global Mode Interface
#'
#' This file provides global mode functions for sintHiChIP.
#'
#' @author sintHiChIP Development Team

#' Run sintHiChIP Global Mode
#'
#' This function processes HiC-Pro output in global mode with peak-to-all interactions.
#'
#' @param outdir Character string. Output directory for global mode results
#' @param hicpro_output Character string. Path to HiC-Pro output directory containing allValidPairs
#' @param sample_name Character string. Sample identifier for file naming
#' @param peaks Character string. Path to the peak file (e.g., from MACS2)
#' @param chr_size Character string. Path to the chromosome size file
#' @param build_matrix Character string. Path to the build matrix tool (from HiC-Pro)
#' @param bin_size Numeric. Bin size for matrix generation (default: 5000 bp)
#' @param normSiteFile Character string. Path to normalization restriction enzyme cut site density file
#' @param FDR Numeric. False Discovery Rate threshold (default: 0.01)
#' @param min_dist Numeric. Minimum distance threshold (default: 20000)
#' @param max_dist Numeric. Maximum distance threshold (default: 2000000)
#' @param keep_temp Logical. Keep temporary intermediate files (default: FALSE)
#' @return Invisible NULL. Results are written to output directory
#' @export
#' @examples
#' \dontrun{
#' sintHiChIP_global(
#'   outdir = "/home/user/hichip/global_output",
#'   hicpro_output = "/home/user/hicpro_results",
#'   sample_name = "sample1",
#'   peaks = "/home/user/hichip/sample_peaks.bed",
#'   chr_size = "/home/user/hichip/chrom_sizes.txt",
#'   build_matrix = "/home/user/hichip/build_matrix",
#'   normSiteFile = "/home/user/hichip/norm_sites.bed"
#' )
#' }
sintHiChIP_global <- function(outdir,
                              hicpro_output,
                              sample_name,
                              peaks,
                              chr_size,
                              build_matrix,
                              bin_size = 5000,
                              normSiteFile,
                              FDR = 0.01,
                              min_dist = 20000,
                              max_dist = 2000000,
                              keep_temp = FALSE) {
  
  # Delegate to the unified interface
  run_sintHiChIP(
    mode = "global",
    outdir = outdir,
    hicpro_output = hicpro_output,
    sample_name = sample_name,
    peaks = peaks,
    chr_size = chr_size,
    build_matrix = build_matrix,
    bin_size = bin_size,
    normSiteFile = normSiteFile,
    FDR = FDR,
    min_dist = min_dist,
    max_dist = max_dist,
    keep_temp = keep_temp
  )
}

#' Run sintHiChIP Global Mode with ValidPairs File
#'
#' Processes a single sample in sintHiChIP global mode from valid pairs to filtered interactions,
#' with statistical testing and WashU track generation. Sample name is derived from the
#' valid_pairs file prefix. Accepts only allValidPairs file input.
#'
#' @param valid_pairs Character string. Path to the allValidPairs file
#' @param peaks Character string. Path to the peak file (BED format)
#' @param chr_size Character string. Path to chromosome sizes file
#' @param build_matrix Character string. Path to build matrix tool
#' @param normSiteFile Character string. Path to the normalization site file
#' @param outdir Character string. Output directory
#' @param bin_size Numeric. Bin size for matrix generation (default: 5000)
#' @param FDR Numeric. False Discovery Rate threshold (default: 0.01)
#' @param min_dist Integer. Minimum interaction distance (default: 20000)
#' @param max_dist Integer. Maximum interaction distance (default: 2000000)
#' @param nbins Numeric. Number of bins for statistical modeling (default: 10)
#' @param keep_temp Logical. Keep temporary files (default: FALSE)
#' 
#' @return Invisible NULL
#' 
#' @examples
#' \dontrun{
#' sintHiChIP_global_single(
#'   valid_pairs = "/path/to/sample.allValidPairs",
#'   peaks = "/path/to/peaks.bed",
#'   chr_size = "/path/to/chrom_sizes.txt",
#'   build_matrix = "/path/to/build_matrix",
#'   normSiteFile = "/path/to/normsite.bed",
#'   outdir = "/path/to/output"
#' )
#' }
#' 
#' @export
sintHiChIP_global_single <- function(valid_pairs,
                                     peaks,
                                     chr_size,
                                     build_matrix,
                                     normSiteFile,
                                     outdir,
                                     bin_size = 5000,
                                     FDR = 0.01,
                                     min_dist = 20000L,
                                     max_dist = 2000000L,
                                     nbins = 10,
                                     keep_temp = FALSE) {
  
  # Validate inputs
  if (!file.exists(valid_pairs)) stop("Valid pairs file not found: ", valid_pairs)
  if (file.size(valid_pairs) == 0) stop("Valid pairs file is empty: ", valid_pairs)
  if (!file.exists(peaks)) stop("Peak file not found: ", peaks)
  if (!file.exists(chr_size)) stop("Chromosome sizes file not found: ", chr_size)
  if (!file.exists(build_matrix)) stop("Build matrix tool not found: ", build_matrix)
  if (!file.exists(normSiteFile)) stop("Normalization site file not found: ", normSiteFile)
  if (missing(outdir) || nchar(outdir) == 0) stop("Output directory is required")
  
  # Derive sample name from valid_pairs file prefix
  sample_name <- sub("\\.allValidPairs.*$", "", basename(valid_pairs))
  if (nchar(sample_name) == 0) stop("Could not derive sample name from valid_pairs file")
  
  # Validate distance parameters
  min_dist <- as.integer(min_dist)
  max_dist <- as.integer(max_dist)
  
  # Create output directory
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  cat("=== sintHiChIP Global Mode ===\n")
  cat("Sample:", sample_name, "\n")
  cat("Bin size:", bin_size, "| FDR:", FDR, "\n")
  cat("Distance filter:", min_dist, "-", max_dist, "bp\n")
  
  # Step 1: Fast matrix generation and peak filtering using pipeline
  cat("Step 1: Fast matrix generation and peak filtering...\n")
  result <- global_process_exec(
    valid_pairs_file = valid_pairs,
    sample_name = sample_name,
    peaks_file = peaks,
    chr_size_file = chr_size,
    build_matrix_tool = build_matrix,
    output_dir = outdir,
    bin_size = bin_size,
    min_dist = min_dist,
    max_dist = max_dist,
    target_chunks = 10,
    keep_temp = keep_temp
  )
  
  if (!file.exists(result)) {
    stop("Fast pipeline failed to generate output file")
  }
  
  # Step 2: Perform statistical significance testing
  cat("Step 2: Statistical significance testing...\n")
  sintHiChIP_sigloops(
    sname = sample_name,
    cwd = outdir,
    outdir = outdir,
    normSiteFile = normSiteFile,
    local = FALSE,  # Global mode flag
    FDR = FDR,
    nbins = nbins
  )
  
  # Step 3: Generate visualization tracks
  cat("Step 3: Generating WashU browser tracks...\n")
  sintHiChIP_make_washu(outdir, FDR, sample_name, mode = "global")
  
  # Step 4: Clean up temporary files if requested
  if (!keep_temp) {
    sintHiChIP_clean_temp(outdir)
  }
  
  cat("Global mode completed for", sample_name, "!\n")
  invisible(NULL)
}