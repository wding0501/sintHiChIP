#' Run Local Mode of sintHiChIP
#'
#' This function executes the local mode of the sintHiChIP pipeline, which analyzes
#' HiChIP data focusing on specific genomic regions or peaks.
#'
#' @param outdir Character string. Output directory name (folder name only, created in the same directory as the yaml file).
#' @param yaml Character string. Name of the yaml configuration file.
#' @param peaks Character string. Peak options ("EACH", "ALL") or file path to custom peak file.
#' @param resfrags Character string. Path to restriction fragments file.
#' @param hicpro_output Character string. Path to HiC-Pro output directory.
#' @param keep_temp Logical. Whether to keep temporary files (default: FALSE).
#' @param normSiteFile Character string. Path to normalization restriction enzyme cut site density file.
#' @param FDR Numeric. False Discovery Rate threshold for significant interactions (default: 0.01).
#' @param ... Additional parameters to be passed to internal functions.
#' @return Invisible NULL. The function produces output files in the specified directory.
#' @export
#' @examples
#' \dontrun{
#' run_sintHiChIP_local(
#'   outdir = "local_analysis",
#'   yaml = "config.yaml",
#'   peaks = "EACH,ALL",
#'   resfrags = "/home/user/hichip/hg19_mboi.bed",
#'   hicpro_output = "/home/user/hichip/sample_HiCPro",
#'   normSiteFile = "/home/user/hichip/normsite_hg19_5000.bed",
#'   FDR = 0.01
#' )
#' }
run_sintHiChIP_local <- function(outdir, yaml = NULL, peaks = "EACH,ALL", resfrags = NULL,
                      hicpro_output = NULL, keep_temp = FALSE, normSiteFile = NULL,
                      FDR = 0.01, ...) {
  # Ensure outdir is just a folder name, not a path
  if (grepl("[/\\]", outdir)) {
    stop("outdir should be a folder name only, not a path")
  }
  
  # Get the directory of the yaml file
  yaml_dir <- dirname(yaml)
  
  # Construct the full path for outdir
  full_outdir <- file.path(yaml_dir, outdir)
  # Step 1: Generate peak-to-peak (p2p) interactions
  generate_p2p_local(yaml, peaks, resfrags, hicpro_output, outdir, keep_temp)

  # Step 2: Generate significant interactions
  compute_sig_local(full_outdir, normSiteFile, FDR, ...)

  # Step 3: Generate WashU Genome Browser tracks
  create_washu_local(full_outdir, FDR)

  cat("Local mode processing completed. Output saved in:", full_outdir, "\n")
}

#' Generate Peak-to-Peak Interactions for Local Mode
#'
#' This function generates peak-to-peak interactions for the local mode of sintHiChIP.
#'
#' @param yaml Character string. Path to output YAML file.
#' @param peaks Character string. Peak options or file path.
#' @param resfrags Character string. Path to restriction fragments file.
#' @param hicpro_output Character string. Path to HiC-Pro output directory.
#' @param outdir Character string. Output directory name (folder name only).
#' @param keep_temp Logical. Whether to keep temporary files.
#' @return Invisible NULL
#' @keywords internal
#' @examples
#' \dontrun{
#' generate_p2p_local(
#'   yaml = "/home/user/hichip/config.yaml",
#'   peaks = "EACH,ALL",
#'   resfrags = "/home/user/hichip/restriction_fragments.bed",
#'   hicpro_output = "/home/user/hichip/hicpro_output",
#'   outdir = "p2p_output",
#'   keep_temp = FALSE
#' )
#' }
generate_p2p_local <- function(yaml, peaks, resfrags, hicpro_output, outdir, keep_temp) {
  # Generate YAML configuration
  generate_hichipper_config(yaml, peaks, resfrags, hicpro_output)
  # Run hichipper
  execute_hichipper(yaml, outdir, keep_temp)
}

#' Compute Significant Local Interactions
#'
#' This function computes significant interactions for the local mode of sintHiChIP.
#'
#' @param outdir Character string. Output directory.
#' @param normSiteFile Character string. Path to normalization site file.
#' @param FDR Numeric. False Discovery Rate threshold.
#' @param nbins Integer. Number of bins for loop calling (default: 10).
#' @return Invisible NULL
#' @keywords internal
#' @examples
#' \dontrun{
#' compute_sig_local(
#'   outdir = "/home/user/hichip/local_output",
#'   normSiteFile = "/home/user/hichip/norm_sites.bed",
#'   FDR = 0.01,
#'   nbins = 10
#' )
#' }
compute_sig_local <- function(outdir, normSiteFile, FDR = 0.01, nbins = 10) {
  cat("Generating significant interactions for local mode.\n")
  
  # Find sample files and extract sample names
  file_list <- list.files(path = outdir, pattern = "filt.intra.loop_counts.bedpe$", full.names = TRUE)
  if (length(file_list) == 0) {
    stop("No .filt.intra.loop_counts.bedpe files found in the specified directory")
  }
  snames <- sub("\\.filt\\.intra\\.loop_counts\\.bedpe$", "", basename(file_list))
  
  # Process each sample
  for (sname in snames) {
    cat(paste("Processing sample:", sname, "\n"))
    
    # Call sintHiChIP_loops function
    sintHiChIP_loops(
      sname = sname,
      cwd = outdir,
      outdir = outdir,
      normSiteFile = normSiteFile,
      local = TRUE,
      FDR = FDR,
      nbins = nbins
    )
    
    cat(paste("Completed processing for sample:", sname, "\n"))
  }
  
  cat("Completed generating significant interactions for all samples.\n")
}

#' Create Local WashU Genome Browser Tracks
#'
#' This function creates WashU Genome Browser tracks for the local mode of sintHiChIP.
#'
#' @param outdir Character string. Output directory.
#' @param FDR Numeric. False Discovery Rate threshold.
#' @return Invisible NULL
#' @keywords internal
create_washu_local <- function(outdir, FDR) {
  cat("Generating WashU Genome Browser tracks for local mode.\n")
  
  # Find sample files
  file_list <- list.files(path = outdir, pattern = paste0("interaction.Q", FDR, ".txt$"), full.names = TRUE)
  
  for (input_file in file_list) {
    sname <- sub(paste0("\\.interaction\\.Q", FDR, "\\.txt$"), "", basename(input_file))
    output_file <- file.path(outdir, paste0(sname, ".interaction.Q", FDR, ".washu.txt"))
    
    # Create WashU compatible file
    awk_command <- "awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\":\"$5\"-\"$6\",\"$7\"\\t\"(NR*2-1)\"\\t.\\n\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$1\":\"$2\"-\"$3\",\"$7\"\\t\"(NR*2)\"\\t.\"}'"
    
    if (Sys.info()['sysname'] == "Darwin") {
      system2("sh", args = c("-c", paste(awk_command, input_file, "| sort -k1,1 -k2,2n >", output_file)))
    } else {
      system(paste(awk_command, input_file, "| bedtools sort >", output_file))
    }
    
    # Compress and index the output file
    system2("bgzip", args = c("-f", output_file))
    system2("tabix", args = c("-p", "bed", paste0(output_file, ".gz")))
    
    cat("Completed WashU track creation for sample:", sname, "\n")
  }
}