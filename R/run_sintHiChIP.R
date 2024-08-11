#' Run sintHiChIP Pipeline
#'
#' This function runs the sintHiChIP pipeline in either local or global mode.
#'
#' @param mode Character string, either "local" or "global".
#' @param outdir Character string. Output directory name for local mode; full path to output for global mode.
#' @param yaml Character string. Name of the yaml configuration file only (local mode).
#' @param peaks Character string. Peak options for local mode (default "EACH,ALL"), or path to input peaks file.
#' @param resfrags Character string. Path to restriction fragments file (local mode).
#' @param hicpro_output Character string. Path to HiC-Pro output directory (local mode).
#' @param keep_temp Logical. Keep temporary files (local mode).
#' @param valid_pairs Character string. Path to valid pairs file (global mode).
#' @param chr_size Character string. Path to chromosome size file (global mode).
#' @param build_matrix Character string. Path to build matrix tool (global mode).
#' @param bin_size Numeric. Bin size for matrix generation (global mode, default 5000).
#' @param normSiteFile Character string. Path to normalization restriction enzyme cut site density file.
#' @param FDR Numeric. False Discovery Rate threshold for significant HiChIP loops (default 0.01).
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' # Run in local mode
#' run_sintHiChIP(mode = "local",
#'                outdir = "output_local",
#'                yaml = "config.yaml",
#'                peaks = "EACH,ALL",
#'                resfrags = "/home/user/hichip/restriction_fragments.bed",
#'                hicpro_output = "/home/user/hichip/hicpro_output",
#'                normSiteFile = "/home/user/hichip/norm_sites.bed",
#'                FDR = 0.01)
#'
#' # Run in global mode
#' run_sintHiChIP(mode = "global",
#'                outdir = "/home/user/hichip/output_global",
#'                valid_pairs = "/home/user/hichip/valid_pairs.txt",
#'                chr_size = "/home/user/hichip/chrom_sizes.txt",
#'                build_matrix = "/home/user/hichip/build_matrix",
#'                peaks = "/home/user/hichip/peaks.bed",
#'                normSiteFile = "/home/user/hichip/norm_sites.bed",
#'                bin_size = 5000,
#'                FDR = 0.01)
#' }
run_sintHiChIP <- function(mode = "local", outdir = "output", yaml = NULL,
                           peaks = "EACH,ALL", resfrags = NULL, hicpro_output = NULL,
                           keep_temp = FALSE, valid_pairs = NULL, chr_size = NULL,
                           build_matrix = NULL, bin_size = 5000, normSiteFile = NULL,
                           FDR = 0.01) {
  if (mode == "local") {
    run_sintHiChIP_local(outdir = outdir, yaml = yaml, peaks = peaks, resfrags = resfrags,
              hicpro_output = hicpro_output, keep_temp = keep_temp,
              normSiteFile = normSiteFile, FDR = FDR)
  } else if (mode == "global") {
    run_sintHiChIP_global(outdir = outdir, valid_pairs = valid_pairs, chr_size = chr_size,
               build_matrix = build_matrix, bin_size = bin_size, peaks = peaks,
               normSiteFile = normSiteFile, FDR = FDR)
  } else {
    stop("Invalid mode. Please choose 'local' or 'global'.")
  }
  
  invisible(NULL)
}