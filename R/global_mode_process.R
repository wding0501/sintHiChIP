#' Convert matrix to bedpe format for sintHiChIP global mode
#'
#' @param matrix_file Path to the matrix file
#' @param bed_file Path to the bed file containing bin coordinates
#' @param output_file Path to the output bedpe file
#' @param min_dist Minimum interaction distance in base pairs (default: 20000)
#' @param max_dist Maximum interaction distance in base pairs (default: 2000000)
#'
#' @return Number of interactions written to output file
#'
#' @examples
#' \dontrun{
#' matrix_to_bedpe_global(
#'   matrix_file = "sample_matrix.matrix",
#'   bed_file = "sample_matrix_abs.bed", 
#'   output_file = "sample_temp.bedpe",
#'   min_dist = 20000,
#'   max_dist = 2000000
#' )
#' }
matrix_to_bedpe_global <- function(matrix_file, bed_file, output_file, 
                                   min_dist = 20000, max_dist = 2000000) {
  
  cat("Converting matrix to bedpe format with distance filtering...\n")
  
  # Read bin coordinates
  bins <- fread(bed_file, header = FALSE, col.names = c("chr", "start", "end", "bin_id"))
  setkey(bins, bin_id)
  
  cat("Loaded", format(nrow(bins), big.mark = ","), "bins\n")
  
  # Read matrix data and convert
  matrix_data <- fread(matrix_file, header = FALSE, col.names = c("bin1", "bin2", "count"))
  cat("Loaded", format(nrow(matrix_data), big.mark = ","), "matrix entries\n")
  
  # Join with first bin coordinates
  result <- matrix_data[bins, on = c(bin1 = "bin_id"), nomatch = 0]
  setnames(result, c("chr", "start", "end"), c("chr1", "start1", "end1"))
  
  # Join with second bin coordinates
  result <- result[bins, on = c(bin2 = "bin_id"), nomatch = 0]
  setnames(result, c("chr", "start", "end"), c("chr2", "start2", "end2"))
  
  cat("After coordinate joining:", format(nrow(result), big.mark = ","), "entries\n")
  
  # Filter for intra-chromosomal interactions only
  result <- result[chr1 == chr2]
  cat("After intra-chromosomal filtering:", format(nrow(result), big.mark = ","), "entries\n")
  
  # Apply distance filtering
  if (min_dist > 0 || is.finite(max_dist)) {
    result[, distance := abs(start1 - start2)]
    result <- result[distance >= min_dist & distance <= max_dist]
    result[, distance := NULL]  # Remove temporary column
    
    cat("After distance filtering (", format(min_dist, big.mark = ","), "-", format(max_dist, big.mark = ","), "bp):", format(nrow(result), big.mark = ","), "entries\n")
  }
  
  # Prepare final bedpe format: chr1, start1, end1, chr2, start2, end2, ., ., count
  final_result <- result[, .(chr1, start1, end1, chr2, start2, end2, ".", ".", count)]
  
  # Write to file
  fwrite(final_result, output_file, col.names = FALSE, sep = "\t")
  
  cat("Matrix to bedpe conversion completed:", format(nrow(final_result), big.mark = ","), "interactions\n")
  return(nrow(final_result))
}

#' Peak overlap detection for sintHiChIP global mode processing
#'
#' @param bedpe_file Path to the input bedpe file
#' @param peaks_file Path to the peaks file in BED format
#' @param output_file Path to the output filtered bedpe file
#' @param target_chunks Number of chunks to divide the processing into (default: 10)
#'
#' @return Number of interactions with peak overlaps
#'
#' @examples
#' \dontrun{
#' sintHiChIP_peak_filter_global(
#'   bedpe_file = "sample_temp.bedpe",
#'   peaks_file = "peaks.bed",
#'   output_file = "sample.filt.intra.loop_counts.bedpe",
#'   target_chunks = 10
#' )
#' }
sintHiChIP_peak_filter_global <- function(bedpe_file, peaks_file, output_file, target_chunks = 10) {
  
  # Read peaks
  cat("Reading peaks file...\n")
  peaks_data <- fread(peaks_file, header = FALSE, col.names = c("chr", "start", "end"))
  cat("Loaded", format(nrow(peaks_data), big.mark = ","), "peaks\n")
  
  # Create GenomicRanges object for peaks
  peaks_gr <- GenomicRanges::GRanges(
    seqnames = peaks_data$chr,
    ranges = IRanges::IRanges(start = peaks_data$start, end = peaks_data$end)
  )
  
  # Get exact total lines using wc -l
  cat("Counting total interactions...\n")
  total_lines <- as.numeric(system(sprintf("wc -l < %s", shQuote(bedpe_file)), intern = TRUE))
  
  # Calculate chunk_size based on target chunks
  chunk_size <- ceiling(total_lines / target_chunks)
  
  cat(sprintf("Total interactions: %s\n", format(total_lines, big.mark = ",")))
  cat(sprintf("Target chunks: %d, calculated chunk_size: %s\n", 
              target_chunks, format(chunk_size, big.mark = ",")))
  
  # Initialize output file
  cat("", file = output_file)
  
  total_filtered <- 0
  skip_rows <- 0
  chunk_count <- 0
  
  repeat {
    chunk_count <- chunk_count + 1
    
    # Read current chunk
    chunk_data <- fread(bedpe_file,
                        skip = skip_rows,
                        nrows = chunk_size,
                        header = FALSE,
                        col.names = c("chr1", "start1", "end1", "chr2", "start2", "end2", "strand1", "strand2", "count"),
                        showProgress = FALSE)
    
    if (is.null(chunk_data) || nrow(chunk_data) == 0) break
    
    cat(sprintf("  Chunk %d/%d: %s interactions...", 
                chunk_count, target_chunks, format(nrow(chunk_data), big.mark = ",")))
    
    # Find peak overlaps for current chunk
    filtered_chunk <- find_peak_overlaps_global(chunk_data, peaks_gr)
    
    if (nrow(filtered_chunk) > 0) {
      # Write to output file (standard bedpe format)
      output_chunk <- filtered_chunk[, .(chr1, start1, end1, chr2, start2, end2, ".", count)]
      fwrite(output_chunk, output_file, append = TRUE, col.names = FALSE, sep = "\t")
      total_filtered <- total_filtered + nrow(filtered_chunk)
    }
    
    kept_pct <- ifelse(nrow(chunk_data) > 0, 100 * nrow(filtered_chunk) / nrow(chunk_data), 0)
    cat(sprintf(" kept %s (%.1f%%)\n", format(nrow(filtered_chunk), big.mark = ","), kept_pct))
    
    skip_rows <- skip_rows + chunk_size
    
    # Memory cleanup
    rm(chunk_data, filtered_chunk)
    invisible(gc())
    
    # If we've processed the expected chunks, check for remaining data
    if (chunk_count >= target_chunks) {
      remaining_lines <- total_lines - skip_rows
      if (remaining_lines > 0) {
        cat(sprintf("Processing remaining %s interactions in final chunk...\n", 
                    format(remaining_lines, big.mark = ",")))
        
        # Process remaining data
        final_chunk <- fread(bedpe_file,
                            skip = skip_rows,
                            header = FALSE,
                            col.names = c("chr1", "start1", "end1", "chr2", "start2", "end2", "strand1", "strand2", "count"),
                            showProgress = FALSE)
        
        if (!is.null(final_chunk) && nrow(final_chunk) > 0) {
          cat(sprintf("  Final chunk: %s interactions...", 
                      format(nrow(final_chunk), big.mark = ",")))
          
          filtered_final <- sintHiChIP::find_peak_overlaps_global(final_chunk, peaks_gr)
          
          if (nrow(filtered_final) > 0) {
            output_final <- filtered_final[, .(chr1, start1, end1, chr2, start2, end2, ".", count)]
            fwrite(output_final, output_file, append = TRUE, col.names = FALSE, sep = "\t")
            total_filtered <- total_filtered + nrow(filtered_final)
          }
          
          final_pct <- ifelse(nrow(final_chunk) > 0, 100 * nrow(filtered_final) / nrow(final_chunk), 0)
          cat(sprintf(" kept %s (%.1f%%)\n", format(nrow(filtered_final), big.mark = ","), final_pct))
        }
      }
      break
    }
  }
  
  cat("Peak overlap detection completed:", format(total_filtered, big.mark = ","), "interactions\n")
  return(total_filtered)
}

#' Find overlaps between HiChIP interactions and peak regions for Global Mode
#'
#' @param interactions Data table with interaction data
#' @param peaks_gr GenomicRanges object containing peak regions
#'
#' @return Data table with interactions that overlap with peaks
#'
#' @examples
#' \dontrun{
#' peaks_gr <- GenomicRanges::GRanges(
#'   seqnames = c("chr1", "chr2"),
#'   ranges = IRanges::IRanges(start = c(1000, 2000), end = c(1500, 2500))
#' )
#' interactions_dt <- data.table(
#'   chr1 = "chr1", start1 = 1200, end1 = 1300,
#'   chr2 = "chr1", start2 = 2200, end2 = 2300,
#'   count = 10
#' )
#' overlapped <- find_peak_overlaps_global(interactions_dt, peaks_gr)
#' }
find_peak_overlaps_global <- function(interactions, peaks_gr) {
  
  if (nrow(interactions) == 0) return(interactions[0, ])
  
  # Create GenomicRanges for left and right anchors
  left_gr <- GenomicRanges::GRanges(
    seqnames = interactions$chr1,
    ranges = IRanges::IRanges(start = interactions$start1, end = interactions$end1)
  )
  
  right_gr <- GenomicRanges::GRanges(
    seqnames = interactions$chr2, 
    ranges = IRanges::IRanges(start = interactions$start2, end = interactions$end2)
  )
  
  # Find overlaps with peaks
  left_overlaps <- GenomicRanges::findOverlaps(left_gr, peaks_gr)
  right_overlaps <- GenomicRanges::findOverlaps(right_gr, peaks_gr)
  
  # Get indices of interactions with at least one anchor overlapping peaks
  left_overlap_idx <- unique(S4Vectors::queryHits(left_overlaps))
  right_overlap_idx <- unique(S4Vectors::queryHits(right_overlaps))
  overlap_idx <- unique(c(left_overlap_idx, right_overlap_idx))
  
  # Return interactions with peak overlaps
  if (length(overlap_idx) > 0) {
    return(interactions[overlap_idx, ])
  } else {
    return(interactions[0, ])
  }
}

#' Generate Matrix Files from allValidPairs for Global Mode
#'
#' This function only handles matrix generation from allValidPairs files,
#' producing matrix and bed files that serve as input for matrix_to_bedpe_global.
#'
#' @param valid_pairs_file Character string. Path to allValidPairs file
#' @param chr_size_file Character string. Path to chromosome sizes file
#' @param build_matrix_tool Character string. Path to build_matrix executable
#' @param output_dir Character string. Output directory
#' @param sample_name Character string. Sample name prefix
#' @param bin_size Numeric. Bin size in bp (default: 5000)
#'
#' @return List containing matrix file path, bed file path, and processing statistics
#'
#' @examples
#' \dontrun{
#' result <- sintHiChIP_Pairs_Matrix_global(
#'   valid_pairs_file = "/path/to/sample.allValidPairs",
#'   chr_size_file = "/path/to/chrom_sizes.txt",
#'   build_matrix_tool = "/path/to/build_matrix",
#'   output_dir = "/path/to/output",
#'   sample_name = "sample1",
#'   bin_size = 5000
#' )
#' }
#'
#' @export
sintHiChIP_Pairs_Matrix_global <- function(valid_pairs_file, 
                                           chr_size_file, 
                                           build_matrix_tool, 
                                           output_dir, 
                                           sample_name, 
                                           bin_size = 5000) {
  
  cat("=== sintHiChIP Matrix Generation (Global Mode) ===\n")
  cat("Sample:", sample_name, "\n")
  cat("Input:", basename(valid_pairs_file), "\n")
  cat("Bin size:", bin_size, "bp\n")
  
  # Validate input files
  if (!file.exists(valid_pairs_file)) {
    stop("Valid pairs file not found: ", valid_pairs_file)
  }
  if (!file.exists(chr_size_file)) {
    stop("Chromosome sizes file not found: ", chr_size_file)
  }
  if (!file.exists(build_matrix_tool)) {
    stop("Build matrix tool not found: ", build_matrix_tool)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  start_time <- Sys.time()
  file_size_gb <- file.size(valid_pairs_file) / (1024^3)
  cat("Input file size:", round(file_size_gb, 1), "GB\n")
  
  # Generate matrix using build_matrix tool
  cat("\nGenerating matrix with build_matrix tool...\n")
  matrix_prefix <- file.path(output_dir, paste0(sample_name, "_matrix"))
  
  # Build command with proper parameters
  matrix_cmd <- sprintf("%s --binsize %d --chrsizes %s --ifile %s --oprefix %s --matrix-format upper",
                        shQuote(build_matrix_tool), 
                        bin_size, 
                        shQuote(chr_size_file), 
                        shQuote(valid_pairs_file), 
                        shQuote(matrix_prefix))
  
  cat("Matrix command:", matrix_cmd, "\n")
  matrix_result <- system(matrix_cmd)
  
  if (matrix_result != 0) {
    stop("Matrix building failed with error code: ", matrix_result)
  }
  
  processing_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Verify output files
  matrix_file <- paste0(matrix_prefix, ".matrix")
  bed_file <- paste0(matrix_prefix, "_abs.bed")
  
  if (!file.exists(matrix_file) || !file.exists(bed_file)) {
    stop("Matrix generation failed. Expected files not found: ", matrix_file, ", ", bed_file)
  }
  
  cat("Matrix generation completed in", round(processing_time, 1), "seconds\n")
  cat("Matrix file:", basename(matrix_file), "(", round(file.size(matrix_file)/(1024^2), 1), "MB)\n")
  cat("Bed file:", basename(bed_file), "(", round(file.size(bed_file)/(1024^2), 1), "MB)\n")
  
  # Return file paths for next step
  return(list(
    matrix_file = matrix_file,
    bed_file = bed_file,
    processing_time = processing_time,
    success = TRUE
  ))
}

#' Global Mode Processing Workflow
#'
#' Complete global mode processing workflow integrating three phases:
#' 1. sintHiChIP_Pairs_Matrix_global - Generate matrix from allValidPairs
#' 2. matrix_to_bedpe_global - Convert matrix to BEDPE format
#' 3. sintHiChIP_peak_filter_global - Filter interactions by peak overlaps
#'
#' @param valid_pairs_file Character string. Path to allValidPairs file
#' @param sample_name Character string. Sample name for identification
#' @param peaks_file Character string. Peak file path (BED format)
#' @param chr_size_file Character string. Path to chromosome sizes file
#' @param build_matrix_tool Character string. Path to build_matrix executable
#' @param output_dir Character string. Output directory
#' @param bin_size Numeric. Bin size in bp (default: 5000)
#' @param min_dist Integer. Minimum interaction distance (default: 20000)
#' @param max_dist Integer. Maximum interaction distance (default: 2000000)
#' @param target_chunks Integer. Target number of chunks for peak detection (default: 10)
#' @param keep_temp Logical. Keep temporary files (default: FALSE)
#' 
#' @return Character string. Path to final output file
#' 
#' @examples
#' \dontrun{
#' output_path <- global_process_exec(
#'   valid_pairs_file = "/path/to/sample.allValidPairs",
#'   sample_name = "sample1",
#'   peaks_file = "/path/to/peaks.bed",
#'   chr_size_file = "/path/to/chrom_sizes.txt",
#'   build_matrix_tool = "/path/to/build_matrix",
#'   output_dir = "/path/to/output"
#' )
#' }
#' 
#' @keywords internal
global_process_exec <- function(valid_pairs_file,
                                sample_name,
                                peaks_file,
                                chr_size_file, 
                                build_matrix_tool, 
                                output_dir,
                                bin_size = 5000,
                                min_dist = 20000L,
                                max_dist = 2000000L,
                                target_chunks = 10,
                                keep_temp = FALSE) {
  
  # Input validation
  if (missing(valid_pairs_file) || is.null(valid_pairs_file) || !file.exists(valid_pairs_file)) {
    stop("Valid pairs file not found: ", valid_pairs_file)
  }
  if (missing(sample_name) || is.null(sample_name) || nchar(sample_name) == 0) {
    stop("sample_name is required and cannot be empty")
  }
  if (missing(peaks_file) || is.null(peaks_file) || !file.exists(peaks_file)) {
    stop("Peak file not found: ", peaks_file)
  }
  if (missing(chr_size_file) || is.null(chr_size_file) || !file.exists(chr_size_file)) {
    stop("Chromosome sizes file not found: ", chr_size_file)
  }
  if (missing(build_matrix_tool) || is.null(build_matrix_tool) || !file.exists(build_matrix_tool)) {
    stop("Build matrix tool not found: ", build_matrix_tool)
  }
  if (missing(output_dir) || is.null(output_dir) || nchar(output_dir) == 0) {
    stop("output_dir is required and cannot be empty")
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("=== sintHiChIP Global Mode Processing ===\n")
  cat("Sample:", sample_name, "\n")
  cat("Bin size:", bin_size, "bp\n")
  cat("Distance filter:", min_dist, "-", max_dist, "bp\n")
  cat("Target chunks:", target_chunks, "\n")
  
  total_start <- Sys.time()
  
  # Phase 1: Generate matrix files using sintHiChIP_Pairs_Matrix_global
  cat("Phase 1: Generating matrix files...\n")
  phase1_start <- Sys.time()
  
  matrix_result <- sintHiChIP_Pairs_Matrix_global(
    valid_pairs_file = valid_pairs_file,
    chr_size_file = chr_size_file,
    build_matrix_tool = build_matrix_tool,
    output_dir = output_dir,
    sample_name = sample_name,
    bin_size = bin_size
  )
  
  if (!matrix_result$success) {
    stop("Phase 1 failed: Matrix generation unsuccessful")
  }
  
  phase1_time <- as.numeric(difftime(Sys.time(), phase1_start, units = "secs"))
  cat("Phase 1 completed in", round(phase1_time, 1), "seconds - matrix generated\n")
  
  # Phase 2: Convert matrix to BEDPE format using matrix_to_bedpe_global
  cat("Phase 2: Converting matrix to BEDPE format...\n")
  phase2_start <- Sys.time()
  
  temp_bedpe <- file.path(output_dir, paste0(sample_name, "_temp.bedpe"))
  
  interaction_count <- matrix_to_bedpe_global(
    matrix_file = matrix_result$matrix_file,
    bed_file = matrix_result$bed_file,
    output_file = temp_bedpe,
    min_dist = min_dist,
    max_dist = max_dist
  )
  
  phase2_time <- as.numeric(difftime(Sys.time(), phase2_start, units = "secs"))
  cat("Phase 2 completed in", round(phase2_time, 1), "seconds -", format(interaction_count, big.mark = ","), "interactions\n")
  
  # Phase 3: Filter interactions by peak overlaps using sintHiChIP_peak_filter_global
  cat("Phase 3: Filtering interactions by peak overlaps...\n")
  phase3_start <- Sys.time()
  
  final_output <- file.path(output_dir, paste0(sample_name, ".filt.intra.loop_counts.bedpe"))
  
  peak_filtered_count <- sintHiChIP_peak_filter_global(
    bedpe_file = temp_bedpe,
    peaks_file = peaks_file,
    output_file = final_output,
    target_chunks = target_chunks
  )
  
  phase3_time <- as.numeric(difftime(Sys.time(), phase3_start, units = "secs"))
  cat("Phase 3 completed in", round(phase3_time, 1), "seconds -", format(peak_filtered_count, big.mark = ","), "peak interactions\n")
  
  # Clean up temporary files if requested
  if (!keep_temp) {
    cat("\nCleaning up temporary files...\n")
    temp_files <- c(matrix_result$matrix_file, matrix_result$bed_file, temp_bedpe)
    existing_temp_files <- temp_files[file.exists(temp_files)]
    if (length(existing_temp_files) > 0) {
      unlink(existing_temp_files)
      cat("Removed", length(existing_temp_files), "temporary files\n")
    }
  }
  
  # Summary
  total_time <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  
  # Verify final output
  if (!file.exists(final_output)) {
    stop("Global processing pipeline failed to generate output file")
  }
  
  cat("\n=== Global Mode Processing Summary ===\n")
  cat("Phase 1 (matrix generation):", round(phase1_time, 1), "s\n")
  cat("Phase 2 (matrix to BEDPE):", round(phase2_time, 1), "s\n")
  cat("Phase 3 (peak filtering):", round(phase3_time, 1), "s\n")
  cat("Total time:", round(total_time, 1), "s (", round(total_time/60, 1), "minutes)\n")
  cat("Performance: Matrix-based processing for comprehensive interaction mapping\n")
  cat("Final interactions:", format(peak_filtered_count, big.mark = ","), "\n")
  cat("Output file:", basename(final_output), "\n")
  
  return(final_output)
}