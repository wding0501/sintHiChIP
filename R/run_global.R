#' Run Global Mode of sintHiChIP
#'
#' This function executes the global mode of the sintHiChIP pipeline, which processes
#' HiChIP data to identify significant chromatin interactions across the entire genome.
#'
#' @param outdir Character string. Full path to the output directory for global mode results.
#' @param valid_pairs Character string. Path to the valid pairs file from HiC-Pro.
#' @param chr_size Character string. Path to the chromosome size file.
#' @param build_matrix Character string. Path to the build matrix tool (typically from HiC-Pro).
#' @param bin_size Numeric. Bin size for matrix generation (default: 5000 bp).
#' @param peaks Character string. Path to the peak file (e.g., from MACS2).
#' @param normSiteFile Character string. Path to normalization restriction enzyme cut site density file.
#' @param FDR Numeric. False Discovery Rate threshold for significant interactions (default: 0.01).
#' @return Invisible NULL. The function generates output files in the specified directory.
#' @export
#' @examples
#' \dontrun{
#' run_sintHiChIP_global(
#'   outdir = "/home/user/hichip/global_output",
#'   valid_pairs = "/home/user/hichip/sample_valid_pairs.txt",
#'   chr_size = "/home/user/hichip/chrom_sizes.txt",
#'   build_matrix = "/home/user/hichip/build_matrix",
#'   bin_size = 5000,
#'   peaks = "/home/user/hichip/sample_peaks.bed",
#'   normSiteFile = "/home/user/hichip/norm_sites.bed",
#'   FDR = 0.01
#' )
#' }
run_sintHiChIP_global <- function(outdir, valid_pairs, chr_size, build_matrix, bin_size = 5000,
                       peaks, normSiteFile, FDR = 0.01) {
  prefix <- tools::file_path_sans_ext(basename(valid_pairs))
  
  # Step 1: Binning interactions and generating matrix
  create_hichip_matrix_global(valid_pairs, chr_size, build_matrix, outdir, bin_size, prefix = prefix)
  
  # Step 2: Generating interactions file with distance filtered
  filter_intra_global(outdir, bin_size, prefix)
  interaction_file <- file.path(outdir, paste0(prefix, ".all.intra.loop_counts.bedpe"))
  # Step 3: Peak overlapped interactions
  overlap_peaks_global(peaks, interaction_file, outdir, prefix)
  # Step 4: Normalization and significance testing
  identify_sig_global_loops(outdir, normSiteFile, FDR, prefix)
  
  # Step 5: Formatting results for WashU Epigenome Browser
  create_washu_global(outdir, FDR, prefix)
  
  cat("Global mode processing completed.\n")
  
  invisible(NULL)
}
#' Filter Interactions for Global Analysis
#'
#' This function filters interactions based on distance thresholds for global analysis in sintHiChIP.
#'
#' @param outdir Character string. Output directory.
#' @param bin_size Numeric. Bin size used for matrix generation.
#' @param prefix Character string. Prefix for output files.
#' @return Invisible NULL
#' @keywords internal
#' @examples
#' \dontrun{
#' filter_intra_global(
#'   outdir = "/home/user/hichip/global_output",
#'   bin_size = 5000,
#'   prefix = "sample1"
#' )
#' }
filter_intra_global <- function(outdir, bin_size, prefix) {
  BinIntervalFile <- file.path(outdir, paste0(prefix, "_abs.bed"))
  MatrixFile <- file.path(outdir, paste0(prefix, ".matrix"))
  Interaction_Initial_File <- file.path(outdir, paste0(prefix, ".interactions.initial.bed"))
  Interaction_File <- file.path(outdir, paste0(prefix, ".all.intra.loop_counts.bedpe"))
  
  filter_distance_global(BinIntervalFile, MatrixFile, Interaction_Initial_File, Interaction_File, 20000, 2000000)
}

#' Identify Significant Global Loops
#'
#' This function identifies significant loops in the global analysis of sintHiChIP.
#'
#' @param outdir Character string. Output directory.
#' @param normSiteFile Character string. Path to normalization site file.
#' @param FDR Numeric. False Discovery Rate threshold.
#' @param prefix Character string. Prefix for output files.
#' @return Invisible NULL
#' @keywords internal
#' @examples
#' \dontrun{
#' identify_sig_global_loops(
#'   outdir = "/home/user/hichip/global_output",
#'   normSiteFile = "/home/user/hichip/norm_sites.bed",
#'   FDR = 0.01,
#'   prefix = "sample1"
#' )
#' }
identify_sig_global_loops <- function(outdir, normSiteFile, FDR, prefix) {
  sintHiChIP_loops(prefix, outdir, outdir, normSiteFile, FALSE, FDR, 10)
}

#' Create Global WashU Genome Browser Tracks
#'
#' This function creates WashU Genome Browser tracks for the global mode of sintHiChIP.
#'
#' @param outdir Character string. Output directory.
#' @param FDR Numeric. False Discovery Rate threshold.
#' @param prefix Character string. Prefix for output files.
#' @return Invisible NULL
#' @keywords internal
#' @examples
#' \dontrun{
#' create_washu_global(
#'   outdir = "/home/user/hichip/global_output",
#'   FDR = 0.01,
#'   prefix = "sample1"
#' )
#' }
create_washu_global <- function(outdir, FDR, prefix) {
  cat("Generating WashU Genome Browser track for global mode.\n")
  
  input_file <- file.path(outdir, paste0(prefix, ".interaction.Q", FDR, ".txt"))
  output_file <- file.path(outdir, paste0(prefix, ".interaction.Q", FDR, ".washu.txt"))
  
  # Create WashU compatible file
  awk_command <- "awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\":\"$5\"-\"$6\",\"$7\"\\t\"(NR*2-1)\"\\t.\\n\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$1\":\"$2\"-\"$3\",\"$7\"\\t\"(NR*2)\"\\t.\"}'"
  full_command <- paste(awk_command, input_file, "| bedtools sort >", output_file)
  system(full_command)
  
  # Compress and index the output file
  system2("bgzip", args = c("-f", output_file))
  system2("tabix", args = c("-p", "bed", paste0(output_file, ".gz")))
  
  cat("Completed WashU track creation for sample:", prefix, "\n")
}