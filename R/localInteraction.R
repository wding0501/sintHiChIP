#' Generate hichipper Config
#'
#' This function generates a YAML configuration file for hichipper.
#'
#' @param output_file Character string. Name of the yaml configuration file only.
#' @param peaks Character string. Peaks options or path to peak file.
#' @param resfrags Character string. Paths to restriction fragments files.
#' @param hicpro_output Character string. Paths to HiC-Pro output directories.
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' generate_hichipper_config(
#'   output_file = "config.yaml",
#'   peaks = "EACH,ALL",
#'   resfrags = "/home/user/hichip/restriction_fragments.bed",
#'   hicpro_output = "/home/user/hichip/hicpro_output"
#' )
#' }
generate_hichipper_config <- function(output_file, peaks, resfrags, hicpro_output) {
  config <- list(
    peaks = list(peaks),
    resfrags = as.list(resfrags),
    hicpro_output = as.list(hicpro_output)
  )
  
  yaml_string <- yaml::as.yaml(config)
  yaml_string <- as.character(yaml_string)
  writeLines(yaml_string, output_file)
  invisible(NULL)
}

#' Execute hichipper
#'
#' This function runs the hichipper tool.
#'
#' @param yaml_file Character string. Path to YAML configuration file.
#' @param outdir Character string. Output directory.
#' @param keep_temp Logical. Keep temporary files.
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' execute_hichipper(
#'   yaml_file = "/home/user/hichip/config.yaml",
#'   outdir = "hichipper_output",
#'   keep_temp = FALSE
#' )
#' }
execute_hichipper <- function(yaml_file, outdir, keep_temp = FALSE) {
  cmd <- paste(
    "hichipper", yaml_file,
    "--out", outdir,
    "--skip-diffloop",
    "--basic-qc",
    if (keep_temp) "--keep-temp-files"
  )
  
  if (.Platform$OS.type == "windows") {
    # For Windows, we need an alternative approach or a Windows-compatible hichipper
    stop("hichipper is not directly available for Windows. Please use WSL or a Unix-like environment.")
  } else if (Sys.info()['sysname'] == "Darwin") {
    system2("hichipper", args = c(yaml_file, "--out", outdir, "--skip-diffloop", "--basic-qc", if (keep_temp) "--keep-temp-files"))
  } else {
    system(cmd)
  }
  invisible(NULL)
}