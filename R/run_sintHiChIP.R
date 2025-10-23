#' Main sintHiChIP Pipeline
#'
#' Main entry point for sintHiChIP processing, supporting both Local and Global modes.
#' This function properly delegates to the specialized functions with unified parameter handling.
#'
#' @param mode Character string. Processing mode: "local" or "global"
#' @param outdir Character string. Output directory for results
#' @param hicpro_output Character string. Path to HiC-Pro output directory (contains allValidPairs)
#' @param sample_name Character string. Sample identifier for file naming (optional, will auto-detect if NULL)
#' @param peaks Character string. Path to peak file (BED format)
#' @param resfrags Character string. Path to restriction fragments file (Local mode)
#' @param chr_size Character string. Path to chromosome sizes file (Global mode)
#' @param build_matrix Character string. Path to build_matrix tool (Global mode)
#' @param bin_size Numeric. Genomic bin size for matrix generation (Global mode, default: 5000)
#' @param normSiteFile Character string. Path to normalization site density file
#' @param FDR Numeric. False Discovery Rate threshold (default: 0.01)
#' @param min_dist Numeric. Minimum interaction distance in base pairs (default: 20000)
#' @param max_dist Numeric. Maximum interaction distance in base pairs (default: 2000000)
#' @param half_length Numeric. Read extension length (default: 73)
#' @param no_merge Logical. Skip anchor merging in Local mode (default: FALSE)
#' @param max_anchor_width Numeric. Maximum anchor width in Local mode (default: 50000)
#' @param keep_temp Logical. Keep temporary files (default: FALSE)
#' @param nbins Numeric. Number of bins for statistical modeling (default: 10)
#' @param peak_pad Numeric. Peak padding in bp (default: 500)
#' @param merge_gap Numeric. Merge gap for bedtools merge (default: 500)
#' @param make_washu Logical. Create WashU/UCSC compatible files (default: TRUE)
#' 
#' @return Invisible NULL. Results written to output directory
#' 
#' @examples
#' \dontrun{
#' # Local Mode
#' run_sintHiChIP(
#'   mode = "local",
#'   outdir = "/path/to/output",
#'   hicpro_output = "/path/to/hicpro_results",
#'   sample_name = "sample1",
#'   peaks = "/path/to/peaks.bed",
#'   resfrags = "/path/to/restriction_fragments.bed",
#'   normSiteFile = "/path/to/norm_sites.bed"
#' )
#' 
#' # Global Mode
#' run_sintHiChIP(
#'   mode = "global",
#'   outdir = "/path/to/output",
#'   hicpro_output = "/path/to/hicpro_results",
#'   sample_name = "sample1",
#'   peaks = "/path/to/peaks.bed",
#'   chr_size = "/path/to/chrom_sizes.txt",
#'   build_matrix = "/path/to/build_matrix",
#'   normSiteFile = "/path/to/norm_sites.bed"
#' )
#' }
#' 
#' @export
run_sintHiChIP <- function(mode,
                           outdir,
                           hicpro_output = NULL,
                           sample_name = NULL,
                           peaks,
                           resfrags = NULL,
                           chr_size = NULL,
                           build_matrix = NULL,
                           bin_size = 5000,
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
  
  # Validate mode
  mode <- match.arg(mode, choices = c("local", "global"))
  
  # Auto-detect sample name if not provided
  if (is.null(sample_name) && !is.null(hicpro_output)) {
    sample_name <- sintHiChIP_detect_sample(hicpro_output)
  }
  
  # Validate common parameters
  sintHiChIP_check_params(outdir, sample_name, peaks, normSiteFile, FDR, min_dist, max_dist)
  
  # Create output directory
  sintHiChIP_make_dir(outdir)
  
  cat("=== sintHiChIP Main Pipeline ===\n")
  cat("Mode:", mode, "\n")
  cat("Sample:", sample_name, "\n")
  cat("Output directory:", outdir, "\n")
  
  # Delegate to mode-specific functions
  if (mode == "local") {
    # Validate local mode specific parameters
    if (is.null(resfrags)) {
      stop("resfrags parameter is required for local mode")
    }
    if (is.null(hicpro_output)) {
      stop("hicpro_output parameter is required for local mode")
    }
    
    # Call local mode function
    sintHiChIP_local(
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
    
  } else {  # mode == "global"
    # Validate global mode specific parameters
    if (is.null(chr_size)) {
      stop("chr_size parameter is required for global mode")
    }
    if (is.null(build_matrix)) {
      stop("build_matrix parameter is required for global mode")
    }
    if (is.null(hicpro_output)) {
      stop("hicpro_output parameter is required for global mode")
    }
    
    # Call global mode function
    sintHiChIP_global(
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
  
  cat("=== sintHiChIP Pipeline Complete ===\n")
  invisible(NULL)
}