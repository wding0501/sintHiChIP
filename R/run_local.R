#' sintHiChIP Local Mode Interface
#'
#' This file provides local mode functions for sintHiChIP.
#'
#' @author sintHiChIP Development Team

#' Run sintHiChIP Local Mode
#'
#' This function processes HiC-Pro output in local mode with peak-to-peak interactions.
#'
#' @param outdir Character string. Output directory for results
#' @param hicpro_output Character string. Path to HiC-Pro output directory containing allValidPairs
#' @param sample_name Character string. Sample identifier for file naming
#' @param peaks Character string. Path to peak file (BED format)
#' @param resfrags Character string. Path to restriction fragments file (BED format)
#' @param normSiteFile Character string. Path to normalization site density file
#' @param FDR Numeric. False Discovery Rate threshold (default: 0.01)
#' @param min_dist Numeric. Minimum interaction distance (default: 20000)
#' @param max_dist Numeric. Maximum interaction distance (default: 2000000)
#' @param half_length Numeric. Read extension length (default: 73)
#' @param no_merge Logical. Skip merging of overlapping anchors (default: FALSE)
#' @param max_anchor_width Numeric. Maximum allowed anchor width (default: 50000)
#' @param keep_temp Logical. Keep temporary intermediate files (default: FALSE)
#' @param nbins Numeric. Number of bins for statistical modeling (default: 10)
#' @param peak_pad Numeric. Peak padding in bp (default: 500)
#' @param merge_gap Numeric. Merge gap for bedtools merge (default: 500)
#' @param make_washu Logical. Create WashU/UCSC compatible files (default: TRUE)
#' @return Invisible NULL. Results are written to output directory
#' @export
#' @examples
#' \dontrun{
#' sintHiChIP_local(
#'   outdir = "/path/to/output",
#'   hicpro_output = "/path/to/hicpro_results",
#'   sample_name = "sample1",
#'   peaks = "/path/to/peaks.bed",
#'   resfrags = "/path/to/restriction_fragments.bed",
#'   normSiteFile = "/path/to/normsite.bed"
#' )
#' }
sintHiChIP_local <- function(outdir,
                             hicpro_output,
                             sample_name,
                             peaks,
                             resfrags,
                             normSiteFile,
                             FDR = 0.01,
                             min_dist = 20000,
                             max_dist = 2000000,
                             half_length = 73,
                             no_merge = FALSE,
                             max_anchor_width = 50000,
                             keep_temp = FALSE,
                             nbins = 10,
                             peak_pad = 500,
                             merge_gap = 500,
                             make_washu = TRUE) {
  
  # Delegate to the unified interface with all parameters
  run_sintHiChIP(
    mode = "local",
    outdir = outdir,
    hicpro_output = hicpro_output,
    sample_name = sample_name,
    peaks = peaks,
    resfrags = resfrags,
    normSiteFile = normSiteFile,
    FDR = FDR,
    min_dist = min_dist,
    max_dist = max_dist,
    half_length = half_length,
    no_merge = no_merge,
    max_anchor_width = max_anchor_width,
    keep_temp = keep_temp,
    nbins = nbins,
    peak_pad = peak_pad,
    merge_gap = merge_gap,
    make_washu = make_washu
  )
}

#' Run sintHiChIP Local Mode with ValidPairs File
#'
#' Processes a single sample in sintHiChIP local mode from valid pairs to filtered loop counts,
#' with statistical testing and WashU track generation. Sample name is derived from the
#' valid_pairs file prefix. Accepts only allValidPairs file input.
#'
#' @param valid_pairs Character string. Path to the allValidPairs file
#' @param peaks Character string. Path to the peak file (BED format)
#' @param resfrags Character string. Path to the restriction fragments file
#' @param normSiteFile Character string. Path to the normalization site file
#' @param outdir Character string. Output directory
#' @param FDR Numeric. False Discovery Rate threshold (default: 0.01)
#' @param min_dist Integer. Minimum interaction distance (default: 20000)
#' @param max_dist Integer. Maximum interaction distance (default: 2000000)
#' @param half_length Numeric. Read extension length (default: 73)
#' @param no_merge Logical. Skip anchor merging (default: FALSE)
#' @param max_anchor_width Numeric. Maximum anchor width (default: 50000)
#' @param nbins Numeric. Number of bins for statistical modeling (default: 10)
#' 
#' @return Invisible NULL
#' 
#' @examples
#' \dontrun{
#' sintHiChIP_local_single(
#'   valid_pairs = "/path/to/sample.allValidPairs",
#'   peaks = "/path/to/peaks.bed",
#'   resfrags = "/path/to/restriction_fragments.bed",
#'   normSiteFile = "/path/to/normsite.bed",
#'   outdir = "/path/to/output"
#' )
#' }
#' 
#' @export
sintHiChIP_local_single <- function(valid_pairs,
                                    peaks,
                                    resfrags,
                                    normSiteFile,
                                    outdir,
                                    FDR = 0.01,
                                    min_dist = 20000L,
                                    max_dist = 2000000L,
                                    half_length = 73,
                                    no_merge = FALSE,
                                    max_anchor_width = 50000,
                                    nbins = 10) {
  
  # Validate inputs
  if (!file.exists(valid_pairs)) stop("Valid pairs file not found: ", valid_pairs)
  if (file.size(valid_pairs) == 0) stop("Valid pairs file is empty: ", valid_pairs)
  if (!file.exists(peaks)) stop("Peak file not found: ", peaks)
  if (!file.exists(resfrags)) stop("Restriction fragments file not found: ", resfrags)
  if (file.size(resfrags) == 0) stop("Restriction fragments file is empty: ", resfrags)
  if (!file.exists(normSiteFile)) stop("Normalization site file not found: ", normSiteFile)
  if (missing(outdir) || nchar(outdir) == 0) stop("Output directory is required")
  
  # Derive sample name from valid_pairs file prefix
  sample_name <- sub("\\.allValidPairs.*$", "", basename(valid_pairs))
  if (nchar(sample_name) == 0) stop("Could not derive sample name from valid_pairs file")
  
  # Create temporary directory to mimic hicpro_output structure
  temp_dir <- file.path(outdir, "temp_hicpro")
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  temp_data_dir <- file.path(temp_dir, "hic_results", "data", sample_name)
  if (!dir.exists(temp_data_dir)) {
    dir.create(temp_data_dir, recursive = TRUE)
  }
  file.copy(valid_pairs, file.path(temp_data_dir, basename(valid_pairs)))
  
  # Validate distance parameters
  min_dist <- as.integer(min_dist)
  max_dist <- as.integer(max_dist)
  
  # Create output directory
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  cat("=== sintHiChIP Local Mode ===\n")
  cat("Sample:", sample_name, "\n")
  
  # Step 1: Core processing
  cat("Step 1: Core processing (allValidPairs to filtered interactions)...\n")
  bedpe_file <- local_process_exec(
    hicpro_output = temp_dir,
    sample_name = sample_name,
    peak_file = peaks,
    resfrags_file = resfrags,
    output_dir = outdir,
    min_dist = min_dist,
    max_dist = max_dist,
    half_length = half_length,
    peak_pad = 500,
    merge_gap = 500,
    no_merge = no_merge,
    max_anchor_width = max_anchor_width
  )
  
  # Verify output
  if (!file.exists(bedpe_file)) stop("Output file not created: ", bedpe_file)
  if (file.size(bedpe_file) == 0) {
    warning("Output file is empty: ", bedpe_file)
    cat("Local mode workflow completed with empty results\n")
    unlink(temp_dir, recursive = TRUE)
    return(invisible(NULL))
  }
  
  # Step 2: Statistical significance testing
  cat("Step 2: Running statistical significance testing...\n")
  sintHiChIP_sigloops(
    sname = sample_name,
    cwd = outdir,
    outdir = outdir,
    normSiteFile = normSiteFile,
    local = TRUE,
    FDR = FDR,
    nbins = nbins
  )
  
  # Step 3: Generate WashU genome browser tracks
  cat("Step 3: Generating WashU genome browser tracks...\n")
  sintHiChIP_make_washu(outdir, FDR, sample_name, mode = "local")
  
  # Clean up temporary directory
  unlink(temp_dir, recursive = TRUE)
  
  cat("Local mode completed for", sample_name, "!\n")
  invisible(NULL)
}