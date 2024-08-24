#' Create HiChIP Matrix for Global Analysis
#'
#' This function generates a HiChIP matrix using the specified parameters for global mode analysis.
#' It utilizes the build_matrix tool from HiC-Pro to create the interaction matrix.
#'
#' @param valid_pairs Character string. Path to valid pairs file from HiC-Pro.
#' @param chr_size Character string. Path to chromosome size file.
#' @param build_matrix Character string. Path to build_matrix executable from HiC-Pro utilities.
#' @param outdir Character string. Output directory for the generated matrix.
#' @param bin_size Numeric. Bin size for matrix generation (default: 5000 bp).
#' @param prefix Character string. Prefix for output files (default: NULL, derived from valid_pairs filename).
#' @return Invisible NULL. The function generates matrix files in the specified output directory.
#' @export
#' @examples
#' \dontrun{
#' create_hichip_matrix_global(
#'   valid_pairs = "/home/user/hichip/sample_valid_pairs.txt",
#'   chr_size = "/home/user/hichip/chrom_sizes.txt",
#'   build_matrix = "/home/user/hichip/build_matrix",
#'   outdir = "/home/user/hichip/matrix_output",
#'   bin_size = 5000,
#'   prefix = "sample_analysis"
#' )
#' }
create_hichip_matrix_global <- function(valid_pairs, chr_size, build_matrix, outdir, bin_size, prefix = NULL) {
  if (is.null(prefix)) {
    prefix <- tools::file_path_sans_ext(basename(valid_pairs))
  }
  
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  input_command <- if(grepl("\\.gz$", valid_pairs)) "zcat" else "cat"
  
  command <- sprintf(
    '%s %s | %s --binsize %d --chrsizes %s --ifile /dev/stdin --oprefix %s --matrix-format upper',
    input_command,
    valid_pairs,
    build_matrix,
    bin_size,
    chr_size,
    file.path(outdir, prefix)
  )
  
  if (Sys.info()['sysname'] == "Darwin") {
    system2("sh", args = c("-c", command))
  } else {
    system(command)
  }
  
  output_file <- file.path(outdir, paste0(prefix, ".matrix"))
  while (!file.exists(output_file) || file.size(output_file) == 0) {
    Sys.sleep(10)
  }
  
  invisible(NULL)
}