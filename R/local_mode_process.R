#' sintHiChIP Local Mode Processing - Fixed Version
#' 
#' Core workflow: allValidPairs -> *.filt.intra.loop_counts.bedpe
#'
#' Essential processing steps:
#' 1. Process peaks to anchors  
#' 2. Combined sorting and overlap detection
#' 3. Generate filtered loop counts

#' Process Peaks to Anchors for Local Mode
#'
#' @param peak_file Character string. Peak file path (BED format)
#' @param resfrags_file Character string. Restriction fragments file path (unused)
#' @param no_merge Logical. Skip anchor merging (default: FALSE)
#' @param pad Numeric. Padding to add to peaks (default: 500)
#' @param merge_gap Numeric. Merge gap for bedtools merge (default: 500)
#' @param max_anchor_width Numeric. Maximum anchor width (default: 50000)
#' @param output_file Character string. Output file path (optional)
#' 
#' @return data.table with processed anchors
#' 
#' @keywords internal
process_anchors_local <- function(peak_file,
                                  resfrags_file = NULL,
                                  no_merge = FALSE,
                                  pad = 500,
                                  merge_gap = 500,
                                  max_anchor_width = 50000,
                                  output_file = NULL) {
  
  # Input validation
  if (missing(peak_file) || is.null(peak_file) || !file.exists(peak_file)) {
    stop("Peak file not found: ", peak_file)
  }
  if (file.size(peak_file) == 0) {
    stop("Peak file is empty: ", peak_file)
  }
  
  cat("=== Processing peaks to anchors ===\n")
  
  # Load and sort peaks
  if (endsWith(peak_file, ".gz")) {
    peaks <- data.table::fread(cmd = paste("zcat <", peak_file), header = FALSE)
  } else {
    peaks <- data.table::fread(peak_file, header = FALSE)
  }
  
  # Validate peaks file format (must have at least 3 columns)
  if (ncol(peaks) < 3) {
    stop("Peak file must have at least 3 columns (chr, start, end): ", peak_file)
  }
  
  # Take first 3 columns and apply padding
  peaks <- peaks[, 1:3]
  data.table::setnames(peaks, c("chr", "start_orig", "end_orig"))
  
  # Validate chromosome names are not empty
  if (any(is.na(peaks$chr)) || any(peaks$chr == "")) {
    stop("Peak file contains invalid chromosome names")
  }
  
  # Validate coordinates and apply padding
  peaks[, ':='(
    start = pmax(0L, as.integer(start_orig - pad)),
    end = as.integer(end_orig + pad),
    start_orig = NULL,
    end_orig = NULL
  )]
  
  # Remove invalid intervals
  invalid_intervals <- peaks[start >= end]
  if (nrow(invalid_intervals) > 0) {
    warning("Removed ", nrow(invalid_intervals), " invalid intervals (start >= end)")
    peaks <- peaks[start < end]
  }
  
  if (nrow(peaks) == 0) {
    stop("No valid peaks remaining after processing")
  }
  
  # Sort by chromosome and position
  data.table::setorder(peaks, chr, start, end)
  
  if (no_merge) {
    anchors_final <- peaks[, .(chr, start, end)]
  } else {
    # Merge overlapping intervals with merge_gap
    anchors_final <- peaks[, {
      if (.N == 1) {
        .(start = start, end = end)
      } else {
        merged_starts <- integer()
        merged_ends <- integer()
        current_start <- start[1]
        current_end <- end[1]
        
        for (i in 2:.N) {
          # Check if intervals should be merged (gap <= merge_gap)
          if (start[i] - current_end <= merge_gap) {
            current_end <- max(current_end, end[i])
          } else {
            merged_starts <- c(merged_starts, current_start)
            merged_ends <- c(merged_ends, current_end)
            current_start <- start[i]
            current_end <- end[i]
          }
        }
        
        merged_starts <- c(merged_starts, current_start)
        merged_ends <- c(merged_ends, current_end)
        .(start = merged_starts, end = merged_ends)
      }
    }, by = chr]
  }
  
  # Filter by maximum anchor width
  anchors_final[, width := end - start]
  original_count <- nrow(anchors_final)
  anchors_filtered <- anchors_final[width <= max_anchor_width]
  anchors_filtered[, width := NULL]
  
  if (nrow(anchors_filtered) == 0) {
    stop("No anchors remain after width filtering (max_anchor_width = ", max_anchor_width, ")")
  }
  
  if (nrow(anchors_filtered) < original_count) {
    cat("Filtered", original_count - nrow(anchors_filtered), "anchors exceeding max width of", max_anchor_width, "bp\n")
  }
  
  # Write output file if specified
  if (!is.null(output_file)) {
    data.table::fwrite(anchors_filtered[, .(chr, start, end)], 
                       output_file, sep = "\t", col.names = FALSE, scipen = 999)
  }
  
  return(anchors_filtered)
}

#' Fixed Combined Processing: Sort Interactions + Dual Overlap (Pure R Implementation)
#'
#' This function combines steps 2 and 3 of local processing for better performance
#' by processing allValidPairs in fixed chunks. Fixed memory management issues
#' that cause hanging at specific chunks.
#'
#' @param hicpro_output Character string. HiC-Pro output directory
#' @param sample_name Character string. Sample name
#' @param anchors data.table. Processed anchors from step 1
#' @param half_length Numeric. Read extension length
#' @param output_dir Character string. Output directory
#' @param chunk_num Integer. Number of chunks to process (default:
sort_and_overlap_local <- function(hicpro_output, 
                                   sample_name, 
                                   anchors, 
                                   half_length, 
                                   output_dir,
                                   chunk_num=10) {
  
  row_id<-NULL
  cat(sprintf("=== Conservative Processing: %d chunks with data.table overlap ===\n", chunk_num))
  
  # Find files
  valid_pairs_pattern <- file.path(hicpro_output, "hic_results", "data", sample_name, "*allValidPairs")
  valid_pairs_files <- Sys.glob(valid_pairs_pattern)
  
  if (length(valid_pairs_files) == 0) {
    valid_pairs_files <- Sys.glob(file.path(hicpro_output, "hic_results", "data", sample_name, paste0(sample_name, ".allValidPairs")))
  }
  
  if (length(valid_pairs_files) == 0) {
    stop("Cannot find allValidPairs file for sample: ", sample_name)
  }
  
  valid_pairs_file <- valid_pairs_files[1]
  cat("Processing file:", valid_pairs_file, "\n")
  
  # Count total lines
  if (grepl("\\.gz$", valid_pairs_file)) {
    total_lines <- as.numeric(system(sprintf("pigz -dc %s | wc -l", shQuote(valid_pairs_file)), intern = TRUE))
  } else {
    total_lines <- as.numeric(system(sprintf("wc -l < %s", shQuote(valid_pairs_file)), intern = TRUE))
  }
  
  chunk_size <- ceiling(total_lines / chunk_num)
  cat(sprintf("Total interactions: %s\n", format(total_lines, big.mark = ",")))
  cat(sprintf("Processing in %d chunks\n", chunk_num))
  
  # Prepare anchors
  anchors_dt <- data.table::as.data.table(anchors)
  data.table::setorder(anchors_dt, chr, start, end)
  anchors_dt[, anchor_id := .I]
  data.table::setkey(anchors_dt, chr, start, end)
  
  # Main processing loop
  result_list <- vector("list", chunk_num)
  for (chunk_id in 1:chunk_num) {
    chunk_start <- Sys.time()
    skip_lines <- (chunk_id - 1) * chunk_size
    
    if (chunk_id == chunk_num) {
      nrows_to_read <- -1
    } else {
      nrows_to_read <- chunk_size
    }
    
    cat(sprintf("Chunk %d/%d: ", chunk_id, chunk_num))
    
    # Aggressive memory cleanup before reading large chunks (starting from chunk 4)
    if (chunk_id >= 4) {
      invisible(gc(verbose = FALSE))
      if (chunk_id >= 6) {
        invisible(gc(verbose = FALSE))  # Double gc for later chunks
      }
    }
    
    # Read data
    chunk_data <- data.table::fread(
      valid_pairs_file,
      skip = skip_lines,
      nrows = nrows_to_read,
      header = FALSE,
      select = c(2, 3, 5, 6),
      showProgress = FALSE,
      verbose = FALSE
    )
    
    if (nrow(chunk_data) == 0) {
      cat("empty\n")
      next
    }
    
    data.table::setnames(chunk_data, c("chr1", "pos1", "chr2", "pos2"))
    
    # Data cleaning
    chunk_data <- chunk_data[pos1 > 0 & pos2 > 0 & !is.na(chr1) & !is.na(chr2)]
    
    if (nrow(chunk_data) == 0) {
      cat("no valid data\n")
      next
    }
    
    # Extend read length
    chunk_data[, ':='(
      start1 = pmax(1L, as.integer(pos1 - half_length)),
      end1 = as.integer(pos1 + half_length),
      start2 = pmax(1L, as.integer(pos2 - half_length)),
      end2 = as.integer(pos2 + half_length),
      row_id = .I
    )]
    
    # Perform overlap detection with data.table
    overlap_start <- Sys.time()
    
    # Left anchor
    left_query <- chunk_data[, .(chr = chr1, start = start1, end = end1, row_id)]
    data.table::setkey(left_query, chr, start, end)
    
    left_overlaps <- data.table::foverlaps(left_query, anchors_dt, type = "any", nomatch = 0L)
    
    # Right anchor
    right_query <- chunk_data[, .(chr = chr2, start = start2, end = end2, row_id)]
    data.table::setkey(right_query, chr, start, end)
    
    right_overlaps <- data.table::foverlaps(right_query, anchors_dt, type = "any", nomatch = 0L)
    
    overlap_time <- as.numeric(difftime(Sys.time(), overlap_start, units = "secs"))
    
    # Process results
    left_anchor_map <- rep(NA_integer_, nrow(chunk_data))
    right_anchor_map <- rep(NA_integer_, nrow(chunk_data))
    
    if (nrow(left_overlaps) > 0) {
      left_unique <- left_overlaps[, .SD[1], by = row_id]
      left_anchor_map[left_unique$row_id] <- left_unique$anchor_id
    }
    
    if (nrow(right_overlaps) > 0) {
      right_unique <- right_overlaps[, .SD[1], by = row_id]
      right_anchor_map[right_unique$row_id] <- right_unique$anchor_id
    }
    
    valid_mask <- !is.na(left_anchor_map) & !is.na(right_anchor_map)
    chunk_valid_count <- sum(valid_mask)
    
    if (chunk_valid_count > 0) {
      valid_indices <- which(valid_mask)
      
      chunk_result <- data.table::data.table(
        anchor1_chr = anchors_dt$chr[left_anchor_map[valid_indices]],
        anchor1_start = anchors_dt$start[left_anchor_map[valid_indices]],
        anchor1_end = anchors_dt$end[left_anchor_map[valid_indices]],
        anchor2_chr = anchors_dt$chr[right_anchor_map[valid_indices]],
        anchor2_start = anchors_dt$start[right_anchor_map[valid_indices]],
        anchor2_end = anchors_dt$end[right_anchor_map[valid_indices]]
      )
      
      result_list[[chunk_id]] <- chunk_result
    }
    
    chunk_time <- as.numeric(difftime(Sys.time(), chunk_start, units = "secs"))
    cat(sprintf("%.1fs (overlap %.1fs) -> %s valid\n", 
                chunk_time, overlap_time, format(chunk_valid_count, big.mark = ",")))
    
    # Enhanced memory cleanup after each chunk
    rm(chunk_data, left_query, right_query, left_overlaps, right_overlaps)
    if (exists("left_unique")) rm(left_unique)
    if (exists("right_unique")) rm(right_unique)
    if (exists("valid_indices")) rm(valid_indices)
    if (exists("chunk_result")) rm(chunk_result)
    if (exists("left_anchor_map")) rm(left_anchor_map)
    if (exists("right_anchor_map")) rm(right_anchor_map)
    if (exists("valid_mask")) rm(valid_mask)
    
    # More frequent garbage collection for later chunks
    if (chunk_id %% 2 == 0) {
      invisible(gc(verbose = FALSE))
    }
    
    # Additional memory management for problematic chunks (7-9 based on your data)
    if (chunk_id >= 7 && chunk_id <= 9) {
      invisible(gc(verbose = FALSE))
      Sys.sleep(0.05)  # Brief pause to allow system memory management
    }
  }
  
  # Combine results
  if (length(result_list) > 0) {
    result_list <- result_list[!sapply(result_list, is.null)]
    if (length(result_list) > 0) {
      final_result <- data.table::rbindlist(result_list)
    } else {
      final_result <- data.table::data.table(
        anchor1_chr = character(0), anchor1_start = integer(0), anchor1_end = integer(0),
        anchor2_chr = character(0), anchor2_start = integer(0), anchor2_end = integer(0)
      )
    }
  } else {
    final_result <- data.table::data.table(
      anchor1_chr = character(0), anchor1_start = integer(0), anchor1_end = integer(0),
      anchor2_chr = character(0), anchor2_start = integer(0), anchor2_end = integer(0)
    )
  }
  
  invisible(gc())
  return(final_result)
}

#' Generate Final Loop Counts with Standardized Paths for Local Mode
#' 
#' @param anchor_pairs data.table with valid anchor pairs OR character string path to anchor interactions file
#' @param output_dir Character string. Output directory
#' @param sample_name Character string. Sample name
#' @param min_dist Integer. Minimum interaction distance (int64_t compatible)
#' @param max_dist Integer. Maximum interaction distance (int64_t compatible)
#' 
#' @return Character string. Path to final output file
#' 
#' @keywords internal
final_loops_local <- function(anchor_pairs, output_dir, sample_name, min_dist, max_dist) {
  
  cat("=== Using R-based loop generation ===\n")
  
  # Handle case where anchor_pairs is a file path
  if (is.character(anchor_pairs)) {
    if (!file.exists(anchor_pairs) || file.size(anchor_pairs) == 0) {
      warning("No anchor interactions file provided or file is empty")
      output_file <- file.path(output_dir, paste0(sample_name, ".filt.intra.loop_counts.bedpe"))
      file.create(output_file)
      return(output_file)
    }
    
    # Read the file
    anchor_pairs <- data.table::fread(anchor_pairs, header = FALSE)
    data.table::setnames(anchor_pairs, c("anchor1_chr", "anchor1_start", "anchor1_end", 
                                         "anchor2_chr", "anchor2_start", "anchor2_end"))
  }
  
  # Validate inputs
  if (nrow(anchor_pairs) == 0) {
    warning("No anchor pairs provided for loop generation")
    # Create empty output file
    output_file <- file.path(output_dir, paste0(sample_name, ".filt.intra.loop_counts.bedpe"))
    file.create(output_file)
    return(output_file)
  }
  
  # Count occurrences of each unique interaction
  loop_counts <- anchor_pairs[, .N, by = .(anchor1_chr, anchor1_start, anchor1_end, 
                                           anchor2_chr, anchor2_start, anchor2_end)]
  
  # Add dot column for format consistency
  loop_counts[, dot := "."]
  data.table::setcolorder(loop_counts, c("anchor1_chr", "anchor1_start", "anchor1_end", 
                                         "anchor2_chr", "anchor2_start", "anchor2_end", 
                                         "dot", "N"))
  
  # Write main loop_counts file
  main_loop_file <- file.path(output_dir, paste0(sample_name, ".loop_counts.bedpe.tmp"))
  data.table::fwrite(loop_counts, main_loop_file, sep = "\t", col.names = FALSE)
  
  # Separate inter- and intra-chromosomal interactions
  inter_loops <- loop_counts[anchor1_chr != anchor2_chr]
  intra_loops <- loop_counts[anchor1_chr == anchor2_chr]
  
  # Write inter-chromosomal file
  inter_file <- file.path(output_dir, paste0(sample_name, ".inter.loop_counts.bedpe"))
  data.table::fwrite(inter_loops, inter_file, sep = "\t", col.names = FALSE)
  cat("Inter-chromosomal interactions:", nrow(inter_loops), "\n")
  
  # Write intra-chromosomal file
  intra_file <- file.path(output_dir, paste0(sample_name, ".intra.loop_counts.bedpe"))
  data.table::fwrite(intra_loops, intra_file, sep = "\t", col.names = FALSE)
  cat("Intra-chromosomal interactions:", nrow(intra_loops), "\n")
  
  # Apply distance filtering for final output
  if (nrow(intra_loops) > 0) {
    # Calculate distances and apply filters
    intra_loops[, distance := abs((anchor2_start + anchor2_end)/2 - (anchor1_start + anchor1_end)/2)]
    
    original_count <- nrow(intra_loops)
    filtered_loops <- intra_loops[
      anchor1_start != anchor2_start & distance >= min_dist & distance <= max_dist
    ]
    filtered_loops[, distance := NULL]
    
    cat("Distance filtering: ", original_count, " -> ", nrow(filtered_loops), " interactions\n")
    cat("Distance range: ", min_dist, " - ", max_dist, " bp\n")
  } else {
    filtered_loops <- data.table::data.table()
    cat("No intra-chromosomal interactions found\n")
  }
  
  # Write final filtered output
  output_file <- file.path(output_dir, paste0(sample_name, ".filt.intra.loop_counts.bedpe"))
  data.table::fwrite(filtered_loops, output_file, sep = "\t", col.names = FALSE)
  cat("Final filtered interactions:", nrow(filtered_loops), "\n")
  
  return(output_file)
}

#' sintHiChIP Local Mode Processing Workflow (Fixed Version)
#'
#' Local mode processing with 3 steps, fixing issues that cause hanging
#' during GenomicRanges operations on large datasets.
#'
#' @param hicpro_output Character string. HiC-Pro output directory path
#' @param sample_name Character string. Sample name for identification
#' @param peak_file Character string. Peak file path (BED format)
#' @param resfrags_file Character string. Restriction fragments file path (unused)
#' @param output_dir Character string. Output directory
#' @param min_dist Integer. Minimum interaction distance
#' @param max_dist Integer. Maximum interaction distance
#' @param half_length Numeric. Read extension length (default: 73)
#' @param peak_pad Numeric. Peak padding in bp (default: 500)
#' @param merge_gap Numeric. Merge gap for bedtools merge (default: 500)
#' @param no_merge Logical. Skip anchor merging (default: FALSE)
#' @param max_anchor_width Numeric. Maximum anchor width (default: 50000)
#' 
#' @return Character string. Path to generated output file
#' 
#' @examples
#' \dontrun{
#' output_path <- local_process_exec(
#'   hicpro_output = "/path/to/hicpro",
#'   sample_name = "sample1",
#'   peak_file = "/path/to/peaks.bed",
#'   output_dir = "/path/to/output",
#'   min_dist = 20000,
#'   max_dist = 2000000
#' )
#' }
#' 
#' @keywords internal
local_process_exec <- function(hicpro_output,
                               sample_name,
                               peak_file,
                               resfrags_file = NULL,
                               output_dir,
                               min_dist = 20000L,
                               max_dist = 2000000L,
                               half_length = 73,
                               peak_pad = 500,
                               merge_gap = 500,
                               no_merge = FALSE,
                               max_anchor_width = 50000) {
  
  # Input validation
  if (missing(hicpro_output) || is.null(hicpro_output) || !dir.exists(hicpro_output)) {
    stop("HiC-Pro output directory not found: ", hicpro_output)
  }
  if (missing(sample_name) || is.null(sample_name) || nchar(sample_name) == 0) {
    stop("sample_name is required and cannot be empty")
  }
  if (missing(peak_file) || is.null(peak_file) || !file.exists(peak_file)) {
    stop("Peak file not found: ", peak_file)
  }
  if (missing(output_dir) || is.null(output_dir) || nchar(output_dir) == 0) {
    stop("output_dir is required and cannot be empty")
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("=== sintHiChIP Local Mode Processing ===\n")
  cat("Sample:", sample_name, "\n")
  cat("Distance filter:", min_dist, "-", max_dist, "bp\n")
  
  total_start <- Sys.time()
  
  # Phase 1: Process peaks to anchors
  cat("Phase 1: Processing peaks to anchors...\n")
  step1_start <- Sys.time()
  
  anchor_file <- file.path(output_dir, paste0(sample_name, ".anchors.bed"))
  anchors <- process_anchors_local(
    peak_file = peak_file,
    resfrags_file = resfrags_file,
    no_merge = no_merge,
    pad = peak_pad,
    merge_gap = merge_gap,
    max_anchor_width = max_anchor_width,
    output_file = anchor_file
  )
  
  step1_time <- as.numeric(difftime(Sys.time(), step1_start, units = "secs"))
  cat("Phase 1 completed in", round(step1_time, 1), "seconds -", nrow(anchors), "anchors\n")
  
  # Phase 2: Combined sorting and overlap detection
  cat("Phase 2: Combined sorting and overlap detection...\n")
  step2_start <- Sys.time()
  
  anchor_pairs <- sort_and_overlap_local(
    hicpro_output = hicpro_output,
    sample_name = sample_name,
    anchors = anchors,
    half_length = half_length,
    output_dir = output_dir
  )
  
  step2_time <- as.numeric(difftime(Sys.time(), step2_start, units = "secs"))
  cat("Phase 2 completed in", round(step2_time, 1), "seconds\n")
  
  # Phase 3: Generate final loop counts and apply distance filtering
  cat("Phase 3: Generating final loop counts with distance filtering...\n")
  step3_start <- Sys.time()
  
  output_file <- final_loops_local(
    anchor_pairs = anchor_pairs, 
    output_dir = output_dir, 
    sample_name = sample_name, 
    min_dist = min_dist, 
    max_dist = max_dist
  )
  
  step3_time <- as.numeric(difftime(Sys.time(), step3_start, units = "secs"))
  cat("Phase 3 completed in", round(step3_time, 1), "seconds\n")
  
  # Summary
  total_time <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  
  cat("\n=== Local Mode Summary ===\n")
  cat("Phase 1 (anchors):", round(step1_time, 1), "s\n")
  cat("Phase 2 (combined sort+overlap):", round(step2_time, 1), "s\n") 
  cat("Phase 3 (final loops):", round(step3_time, 1), "s\n")
  cat("Total time:", round(total_time, 1), "s\n")
  cat("Output file:", basename(output_file), "\n")
  
  return(output_file)
}





