#' Common Utilities for sintHiChIP
#'
#' This file contains unified utility functions shared across different modes
#' to eliminate redundancy and ensure consistency.
#'
#' @author sintHiChIP Development Team

#' Find allValidPairs File
#' 
#' Function to locate allValidPairs file with comprehensive error handling.
#' This replaces all the mode-specific versions to ensure consistency.
#' 
#' @param hicpro_output Character string. HiC-Pro output directory
#' @param sample_name Character string. Sample name
#' 
#' @return Character string. Path to allValidPairs file
#' 
#' @examples
#' \dontrun{
#' # Find allValidPairs file
#' valid_pairs_path <- sintHiChIP_find_validpairs(
#'   hicpro_output = "/path/to/hicpro_results",
#'   sample_name = "sample1"
#' )
#' print(valid_pairs_path)
#' }
#' 
#' @export
sintHiChIP_find_validpairs <- function(hicpro_output, sample_name) {
  
  if (!dir.exists(hicpro_output)) {
    stop("HiC-Pro output directory does not exist: ", hicpro_output)
  }
  
  if (is.null(sample_name) || nchar(sample_name) == 0) {
    stop("sample_name cannot be empty")
  }
  
  # Define possible paths in order of preference
  possible_paths <- c(
    file.path(hicpro_output, "hic_results", "data", sample_name),
    file.path(hicpro_output, "hic_results", "data"),
    file.path(hicpro_output, "data", sample_name),
    file.path(hicpro_output, "data"),
    file.path(hicpro_output, sample_name),
    hicpro_output
  )
  
  # Define possible filename patterns in order of preference
  patterns <- c(
    paste0("^", sample_name, ".*allValidPairs.*"),  # Most specific pattern
    paste0("^", sample_name, ".*\\.allValidPairs"),
    paste0("^", sample_name, "_allValidPairs.*"),
    ".*allValidPairs.*"  # Most general pattern
  )
  
  found_files <- character()
  
  # Search through all combinations of paths and patterns
  for (path in possible_paths) {
    if (dir.exists(path)) {
      for (pattern in patterns) {
        files <- list.files(path, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
        if (length(files) > 0) {
          # Filter for valid files (exist and non-empty)
          valid_files <- files[file.exists(files) & file.size(files) > 0]
          if (length(valid_files) > 0) {
            found_files <- c(found_files, valid_files)
          }
        }
      }
    }
  }
  
  # Remove duplicates and validate results
  found_files <- unique(found_files)
  
  if (length(found_files) == 0) {
    stop("Could not find allValidPairs file for sample: ", sample_name,
         "\nSearched in directories: ", paste(possible_paths[file.exists(possible_paths)], collapse = ", "),
         "\nUsed patterns: ", paste(patterns, collapse = ", "),
         "\nPlease verify that HiC-Pro processing completed successfully and files are accessible.")
  } else if (length(found_files) == 1) {
    cat("Found allValidPairs file:", found_files[1], "\n")
    return(found_files[1])
  } else {
    # Multiple files found - select the most appropriate one
    # Prefer files that match sample name most closely
    exact_matches <- found_files[grepl(paste0(sample_name, ".*allValidPairs"), basename(found_files))]
    if (length(exact_matches) > 0) {
      selected_file <- exact_matches[1]
    } else {
      selected_file <- found_files[1]
    }
    
    warning("Multiple allValidPairs files found for sample: ", sample_name,
            "\nFound files: ", paste(found_files, collapse = ", "),
            "\nUsing: ", selected_file)
    return(selected_file)
  }
}

#' Create WashU Genome Browser Tracks
#' 
#' Internal helper function that generates visualization tracks
#' compatible with WashU Epigenome Browser for both local and global modes.
#' 
#' @param outdir Character string. Output directory
#' @param FDR Numeric. FDR threshold used for file naming
#' @param prefix Character string. Sample name prefix
#' @param mode Character string. Processing mode ("local" or "global")
#' 
#' @return Invisible NULL
#' 
#' @examples
#' \dontrun{
#' # Generate WashU tracks for local mode
#' sintHiChIP_make_washu(
#'   outdir = "/path/to/output",
#'   FDR = 0.01,
#'   prefix = "sample1",
#'   mode = "local"
#' )
#' 
#' # Generate WashU tracks for global mode
#' sintHiChIP_make_washu(
#'   outdir = "/path/to/output",
#'   FDR = 0.05,
#'   prefix = "sample1",
#'   mode = "global"
#' )
#' }
#' 
#' @export
sintHiChIP_make_washu <- function(outdir, FDR, prefix, mode = "local") {
  
  # Validate mode
  mode <- match.arg(mode, choices = c("local", "global"))
  
  cat("Generating WashU Genome Browser track for", mode, "mode...\n")
  
  # Check for interaction file with mode-specific naming convention
  input_file <- file.path(outdir, paste0(prefix, ".interaction.", mode, ".Q", FDR, ".txt"))
  output_file <- file.path(outdir, paste0(prefix, ".interaction.", mode, ".Q", FDR, ".washu.txt"))
  
  # Fallback to legacy naming convention if new one doesn't exist
  if (!file.exists(input_file)) {
    legacy_input_file <- file.path(outdir, paste0(prefix, ".interaction.Q", FDR, ".txt"))
    if (file.exists(legacy_input_file)) {
      input_file <- legacy_input_file
      output_file <- file.path(outdir, paste0(prefix, ".interaction.Q", FDR, ".washu.txt"))
      cat("Using legacy naming convention for interaction file\n")
    } else {
      cat("Warning: No interaction file found for FDR:", FDR, "\n")
      cat("Checked files:\n")
      cat("  -", file.path(outdir, paste0(prefix, ".interaction.", mode, ".Q", FDR, ".txt")), "\n")
      cat("  -", legacy_input_file, "\n")
      return(invisible(NULL))
    }
  }
  
  # Validate input file is not empty
  if (file.size(input_file) == 0) {
    cat("Warning: Input interaction file is empty:", basename(input_file), "\n")
    return(invisible(NULL))
  }
  
  # Create WashU compatible file using awk command
  # This converts the interaction format to WashU browser format
  awk_command <- "awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\":\"$5\"-\"$6\",\"$7\"\\t\"(NR*2-1)\"\\t.\\n\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$1\":\"$2\"-\"$3\",\"$7\"\\t\"(NR*2)\"\\t.\"}'"
  
  # Execute command based on operating system
  if (Sys.info()['sysname'] == "Darwin") {
    # macOS specific command execution
    full_command <- paste(awk_command, shQuote(input_file), "| sort -k1,1 -k2,2n >", shQuote(output_file))
    system2("sh", args = c("-c", full_command))
  } else {
    # Linux/Unix command execution with bedtools if available
    if (system("which bedtools", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
      full_command <- paste(awk_command, shQuote(input_file), "| bedtools sort >", shQuote(output_file))
    } else {
      full_command <- paste(awk_command, shQuote(input_file), "| sort -k1,1 -k2,2n >", shQuote(output_file))
    }
    system(full_command)
  }
  
  # Check if output file was created successfully
  if (file.exists(output_file) && file.size(output_file) > 0) {
    # Try compression but continue if tools are not available
    bgzip_result <- system2("bgzip", args = c("-f", output_file), stderr = FALSE, stdout = FALSE)
    if (bgzip_result == 0) {
      tabix_result <- system2("tabix", args = c("-p", "bed", paste0(output_file, ".gz")), 
                             stderr = FALSE, stdout = FALSE)
      if (tabix_result == 0) {
        cat("Created compressed WashU track:", paste0(basename(output_file), ".gz"), "\n")
      } else {
        cat("Created WashU track (indexing failed):", basename(output_file), "\n")
      }
    } else {
      cat("Created WashU track (compression failed):", basename(output_file), "\n")
    }
  } else {
    warning("Failed to create WashU track file: ", output_file)
  }
  
  invisible(NULL)
}

#' Validate Distance Parameters with Type Conversion
#' 
#' Ensures distance parameters are properly validated and converted to int64_t compatible format
#' to match C++ function expectations.
#' 
#' @param min_dist Numeric. Minimum distance value
#' @param max_dist Numeric. Maximum distance value
#' @param param_name Character string. Parameter name for error messages
#' 
#' @return List with validated min_dist and max_dist as integers
#' 
#' @examples
#' \dontrun{
#' # Validate distance parameters
#' distance_params <- sintHiChIP_check_distance(
#'   min_dist = 20000,
#'   max_dist = 2000000,
#'   param_name = "interaction_distance"
#' )
#' print(distance_params$min_dist)
#' print(distance_params$max_dist)
#' }
#' 
#' @export
sintHiChIP_check_distance <- function(min_dist, max_dist, param_name = "distance") {
  
  # Check for missing parameters
  if (missing(min_dist) || missing(max_dist)) {
    stop(param_name, " parameters (min_dist, max_dist) are required")
  }
  
  # Validate numeric types
  if (!is.numeric(min_dist) || !is.numeric(max_dist)) {
    stop(param_name, " parameters must be numeric values")
  }
  
  # Check for valid ranges
  if (min_dist < 0) {
    stop("min_dist must be non-negative, got: ", min_dist)
  }
  
  if (max_dist <= min_dist) {
    stop("max_dist must be greater than min_dist, got min_dist=", min_dist, ", max_dist=", max_dist)
  }
  
  # Convert to integers to match C++ int64_t expectations
  # Check for potential overflow (R integers are 32-bit)
  if (max_dist > .Machine$integer.max) {
    warning("max_dist (", max_dist, ") exceeds R integer max (", .Machine$integer.max, "). Using as-is.")
    min_dist_int <- as.numeric(min_dist)
    max_dist_int <- as.numeric(max_dist)
  } else {
    min_dist_int <- as.integer(min_dist)
    max_dist_int <- as.integer(max_dist)
  }
  
  return(list(
    min_dist = min_dist_int,
    max_dist = max_dist_int
  ))
}

#' Validate Common Input Parameters
#' 
#' Unified parameter validation function to ensure consistency across all modes.
#' 
#' @param outdir Character string. Output directory
#' @param sample_name Character string. Sample name
#' @param peaks Character string. Peak file path
#' @param normSiteFile Character string. Normalization site file path
#' @param FDR Numeric. FDR threshold
#' @param min_dist Numeric. Minimum distance
#' @param max_dist Numeric. Maximum distance
#' 
#' @return Invisible NULL (stops execution if validation fails)
#' 
#' @examples
#' \dontrun{
#' # Validate common parameters
#' sintHiChIP_check_params(
#'   outdir = "/path/to/output",
#'   sample_name = "sample1",
#'   peaks = "/path/to/peaks.bed",
#'   normSiteFile = "/path/to/norm_sites.bed",
#'   FDR = 0.01,
#'   min_dist = 20000,
#'   max_dist = 2000000
#' )
#' }
#' 
#' @export
sintHiChIP_check_params <- function(outdir, sample_name, peaks, normSiteFile, FDR, min_dist, max_dist) {
  
  # Validate output directory
  if (missing(outdir) || is.null(outdir) || nchar(outdir) == 0) {
    stop("outdir is required and cannot be empty")
  }
  
  # Validate sample name
  if (missing(sample_name) || is.null(sample_name) || nchar(sample_name) == 0) {
    stop("sample_name is required and cannot be empty")
  }
  
  # Validate peaks file
  if (missing(peaks) || is.null(peaks) || !file.exists(peaks)) {
    stop("Peak file not found: ", peaks)
  }
  if (file.size(peaks) == 0) {
    stop("Peak file is empty: ", peaks)
  }
  
  # Validate normalization site file
  if (missing(normSiteFile) || is.null(normSiteFile) || !file.exists(normSiteFile)) {
    stop("Normalization site file not found: ", normSiteFile)
  }
  if (file.size(normSiteFile) == 0) {
    stop("Normalization site file is empty: ", normSiteFile)
  }
  
  # Validate FDR
  if (!missing(FDR)) {
    if (!is.numeric(FDR) || FDR <= 0 || FDR >= 1) {
      stop("FDR must be a numeric value between 0 and 1, got: ", FDR)
    }
  }
  
  # Validate distance parameters if provided
  if (!missing(min_dist) && !missing(max_dist)) {
    sintHiChIP_check_distance(min_dist, max_dist)
  }
  
  invisible(NULL)
}

#' Create Directory with Error Handling
#' 
#' Unified directory creation function with consistent error handling.
#' 
#' @param dir_path Character string. Directory path to create
#' @param recursive Logical. Create parent directories if needed (default: TRUE)
#' 
#' @return Invisible NULL
#' 
#' @examples
#' \dontrun{
#' # Create a new directory
#' sintHiChIP_make_dir("/path/to/new/directory")
#' 
#' # Create directory without recursive creation
#' sintHiChIP_make_dir("/existing/path/new_dir", recursive = FALSE)
#' }
#' 
#' @export
sintHiChIP_make_dir <- function(dir_path, recursive = TRUE) {
  
  if (is.null(dir_path) || nchar(dir_path) == 0) {
    stop("Directory path cannot be empty")
  }
  
  if (!dir.exists(dir_path)) {
    success <- dir.create(dir_path, recursive = recursive)
    if (!success) {
      stop("Failed to create directory: ", dir_path)
    }
    cat("Created directory:", dir_path, "\n")
  }
  
  invisible(NULL)
}

#' Generate Standard File Paths
#' 
#' Centralized function to generate consistent file paths across modes.
#' This eliminates hardcoded path inconsistencies.
#' 
#' @param outdir Character string. Output directory
#' @param sample_name Character string. Sample name
#' @param file_type Character string. Type of file to generate path for
#' @param mode Character string. Processing mode ("local" or "global")
#' @param bin_size Numeric. Bin size for global mode (optional)
#' @param FDR Numeric. FDR threshold (optional)
#' 
#' @return Character string. Generated file path
#' 
#' @examples
#' \dontrun{
#' # Generate path for significant interactions file
#' sig_file_path <- sintHiChIP_file_path(
#'   outdir = "/path/to/output",
#'   sample_name = "sample1",
#'   file_type = "significant_interactions",
#'   mode = "local",
#'   FDR = 0.01
#' )
#' 
#' # Generate path for matrix file
#' matrix_path <- sintHiChIP_file_path(
#'   outdir = "/path/to/output",
#'   sample_name = "sample1",
#'   file_type = "matrix",
#'   mode = "global"
#' )
#' }
#' 
#' @export
sintHiChIP_file_path <- function(outdir, sample_name, file_type, mode = "local", bin_size = NULL, FDR = NULL) {
  
  # Validate inputs
  if (missing(outdir) || missing(sample_name) || missing(file_type)) {
    stop("outdir, sample_name, and file_type are required")
  }
  
  mode <- match.arg(mode, choices = c("local", "global"))
  
  # Generate paths based on file type
  file_path <- switch(file_type,
    # Matrix files (global mode)
    "matrix" = file.path(outdir, paste0(sample_name, ".matrix")),
    "bin_intervals" = {
      if (is.null(bin_size)) stop("bin_size required for bin_intervals file")
      file.path(outdir, paste0(sample_name,"_abs.bed"))
    },
    
    # Interaction files
    "initial_interactions" = file.path(outdir, paste0(sample_name, ".interactions.initial.bed")),
    "filtered_interactions" = file.path(outdir, paste0(sample_name, ".all.intra.loop_counts.bedpe")),
    "final_bedpe" = file.path(outdir, paste0(sample_name, ".filt.intra.loop_counts.bedpe")),
    
    # Results files
    "significant_interactions" = {
      if (is.null(FDR)) stop("FDR required for significant_interactions file")
      file.path(outdir, paste0(sample_name, ".interaction.", mode, ".Q", FDR, ".txt"))
    },
    "washu_track" = {
      if (is.null(FDR)) stop("FDR required for washu_track file")
      file.path(outdir, paste0(sample_name, ".interaction.", mode, ".Q", FDR, ".washu.txt"))
    },
    
    # Temporary files
    "temp_interactions" = file.path(outdir, paste0(sample_name, ".interactions.tmp")),
    "temp_loops" = file.path(outdir, paste0(sample_name, ".loops.tmp")),
    
    # Anchor files (local mode)
    "anchors" = file.path(outdir, paste0(sample_name, ".anchors.bed")),
    
    # Other output files
    "inter_bedpe" = file.path(outdir, paste0(sample_name, ".inter.bedpe")),
    "intra_bedpe" = file.path(outdir, paste0(sample_name, ".intra.bedpe")),
    
    # Default case
    stop("Unknown file_type: ", file_type)
  )
  
  return(file_path)
}

#' Clean Temporary Files
#' 
#' Function to clean up temporary files with consistent error handling.
#' 
#' @param outdir Character string. Output directory
#' @param pattern Character string. File pattern to match (default: "\\.tmp$")
#' @param verbose Logical. Print cleanup messages (default: TRUE)
#' @return Invisible NULL
#' @keywords internal
sintHiChIP_clean_temp <- function(outdir, pattern = "\\.tmp$", verbose = TRUE) {
  
  if (!dir.exists(outdir)) {
    if (verbose) cat("Output directory does not exist, skipping cleanup\n")
    return(invisible(NULL))
  }
  
  temp_files <- list.files(outdir, pattern = pattern, full.names = TRUE)
  
  if (length(temp_files) > 0) {
    # Single core file removal
    success <- file.remove(temp_files)
    removed_count <- sum(success)
    
    if (verbose) {
      if (removed_count > 0) {
        cat("Cleaned up", removed_count, "temporary files\n")
      }
      if (removed_count < length(temp_files)) {
        warning("Failed to remove some temporary files")
      }
    }
  } else {
    if (verbose) cat("No temporary files found for cleanup\n")
  }
  
  invisible(NULL)
}

#' Auto-detect Sample Name from HiC-Pro Output
#' 
#' Function to automatically detect sample name with better error handling.
#' This is the unified version used across all modes.
#' 
#' @param hicpro_output Character string. HiC-Pro output directory
#' @return Character string. Detected sample name
#' @export
sintHiChIP_detect_sample <- function(hicpro_output) {
  
  if (!dir.exists(hicpro_output)) {
    stop("HiC-Pro output directory does not exist: ", hicpro_output)
  }
  
  # Check common HiC-Pro directory structures
  data_paths <- c(
    file.path(hicpro_output, "hic_results", "data"),
    file.path(hicpro_output, "data"),
    hicpro_output
  )
  
  samples <- character()
  
  for (data_path in data_paths) {
    if (dir.exists(data_path)) {
      # Look for sample subdirectories containing allValidPairs files
      subdirs <- list.dirs(data_path, full.names = FALSE, recursive = FALSE)
      for (subdir in subdirs) {
        if (nchar(subdir) > 0) {  # Skip empty directory names
          subdir_path <- file.path(data_path, subdir)
          allvalidpairs_files <- list.files(subdir_path, pattern = "allValidPairs", full.names = TRUE)
          if (length(allvalidpairs_files) > 0) {
            # Verify files are not empty
            valid_files <- allvalidpairs_files[file.size(allvalidpairs_files) > 0]
            if (length(valid_files) > 0) {
              samples <- c(samples, subdir)
            }
          }
        }
      }
      
      # If no subdirectories with samples found, look for files directly in data_path
      if (length(samples) == 0) {
        files <- list.files(data_path, pattern = "allValidPairs", full.names = FALSE)
        for (file in files) {
          file_path <- file.path(data_path, file)
          if (file.size(file_path) > 0) {  # Check file is not empty
            # Extract sample name from filename
            sample_name <- gsub("_allValidPairs.*|\\.allValidPairs.*", "", file)
            if (nchar(sample_name) > 0) {
              samples <- c(samples, sample_name)
            }
          }
        }
      }
      
      # Stop searching if samples found
      if (length(samples) > 0) break
    }
  }
  
  # Remove duplicates and validate results
  samples <- unique(samples)
  
  if (length(samples) == 0) {
    stop("No valid allValidPairs files detected in HiC-Pro output directory: ", hicpro_output,
         "\nPlease check that HiC-Pro processing completed successfully.",
         "\nSearched in paths: ", paste(data_paths, collapse = ", "))
  } else if (length(samples) == 1) {
    cat("Auto-detected sample name:", samples[1], "\n")
    return(samples[1])
  } else {
    stop("Multiple samples detected: ", paste(samples, collapse = ", "), 
         "\nPlease specify sample_name manually to avoid ambiguity.")
  }
}


#' Create GRanges with data.table compatibility
#' 
#' Helper function that converts data.frame or data.table to GRanges object
#' while preserving additional metadata columns.
#' 
#' @param df Data.frame or data.table. Data frame or data.table with genomic coordinates
#' @param seqnames.field Character. Column name for chromosome. Default "seqnames"
#' @param start.field Character. Column name for start position. Default "start" 
#' @param end.field Character. Column name for end position. Default "end"
#' 
#' @return GRanges object with genomic coordinates and metadata
#' 
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom data.table is.data.table
#' 
#' @examples
#' \dontrun{
#' # Basic usage with default column names
#' df <- data.frame(
#'   seqnames = c("chr1", "chr2"), 
#'   start = c(1000, 2000), 
#'   end = c(1500, 2500),
#'   score = c(10, 20)
#' )
#' gr <- sintHiChIP_create_granges_fast(df)
#' 
#' # Custom column names
#' df2 <- data.frame(
#'   chromosome = c("chr1", "chr2"),
#'   pos_start = c(1000, 2000),
#'   pos_end = c(1500, 2500)
#' )
#' gr2 <- sintHiChIP_create_granges_fast(df2, 
#'                                       seqnames.field = "chromosome",
#'                                       start.field = "pos_start", 
#'                                       end.field = "pos_end")
#' 
#' # With data.table input
#' library(data.table)
#' dt <- data.table(seqnames = "chr1", start = 1000, end = 2000, value = 5)
#' gr3 <- sintHiChIP_create_granges_fast(dt)
#' }
#' 
#' @export
sintHiChIP_create_granges_fast <- function(df, seqnames.field = "seqnames", start.field = "start", end.field = "end") {
  # Ensure compatibility with both data.frame and data.table
  if (data.table::is.data.table(df)) {
    df <- as.data.frame(df)
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = df[[seqnames.field]],
    ranges = IRanges::IRanges(start = df[[start.field]], end = df[[end.field]])
  )
  
  other_columns <- setdiff(names(df), c(seqnames.field, start.field, end.field, "width", "strand"))
  if (length(other_columns) > 0) {
    S4Vectors::mcols(gr) <- df[other_columns] 
  }
  
  return(gr)
}