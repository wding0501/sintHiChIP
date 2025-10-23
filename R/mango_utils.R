#' @importFrom data.table setnames as.data.table setkey rbindlist setDT .SD
#' @importFrom dplyr %>% select arrange full_join
#' @importFrom readr read_delim col_character col_integer
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps GRanges coverage binnedAverage start end
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom methods setClass setGeneric setMethod slot "slot<-" new
#' @importFrom GenomeInfoDb renameSeqlevels sortSeqlevels
#' @importFrom matrixStats colMedians
#' @importFrom stats median complete.cases
#' @importFrom utils head tail
#' @import data.table 
NULL

#' Union class for numeric, integer, and matrix
#'
#' This function is directly from R package diffloop.
#' @name nim-class
#' @keywords internal
setClassUnion("nim", c("numeric", "integer", "matrix"))

#' A class to represent HiChIP interaction data and annotations
#'
#' @slot anchors A GRanges object describing loop anchor locations
#' @slot interactions A matrix. Each row is an interaction between two anchors
#' @slot counts A matrix with the number paired-end reads per loop per sample
#' @slot colData A data.frame with features (columns) for each sample (rows)
#' @slot rowData A data.frame with features (columns) for each loop (rows)
#' 
#' @export
loops <- setClass("loops", slots = c(
  anchors = "GRanges", 
  interactions = "nim", 
  counts = "nim",
  colData = "data.frame", 
  rowData = "data.frame"
))


#' Create Equal Bins
#'
#' This function creates bins for a given vector. This function is directly developed from R package mango.
#'
#' @param vectortobin Vector to bin
#' @param binmethod Method for binning ("equalsize" or "equalocc")
#' @param numberbins Number of bins
#' @return Vector of bin borders
#' @export
#' @examples
#' \dontrun{
#' vector_to_bin <- 1:100
#' bin_borders <- make_equal_bins(vector_to_bin, binmethod = "equalsize", numberbins = 10)
#' }
make_equal_bins <- function(vectortobin, binmethod = "equalsize", numberbins = 30) {
  if (binmethod == "equalsize") {
    sortedvec <- sort(vectortobin)
    borders <- seq(min(sortedvec), max(sortedvec), length.out = numberbins + 1)[2:numberbins]
    return(borders)
  }
    
  if (binmethod == "equalocc") {
    n <- length(vectortobin)
    itemsperbin <- n / numberbins
      
    dt <- data.table::data.table(values = vectortobin)
    data.table::setorder(dt, values)
      
    dt[, bin_id := ceiling(seq_len(.N) / itemsperbin)]
      
    bin_stats <- dt[, .(min_val = min(values), max_val = max(values)), by = bin_id]
    data.table::setorder(bin_stats, bin_id)
      
    borders <- (bin_stats$min_val[-1] + bin_stats$max_val[-nrow(bin_stats)]) / 2
    return(borders)
  }
}



#' Create Loops Object Fast
#'
#' This function creates loop objects from BEDPE files quickly.
#' It is an optimized version of the loopsMake function from R package diffloop.
#'
#' @param beddir A string. The preprocessed data directory
#' @param snames A character vector. Optional list of samples to read in
#' @param type Specifies 'intra', 'inter', or 'all' looping. Default 'all'
#' 
#' @return A loops object
#' @export
#'
#' @examples
#' \dontrun{
#' beddir <- "/home/user/hichip/bedpe_files"
#' sname <- "sample1"
#' snames = paste0(sname, ".filt.intra")
#' loops <- create_loops_object_fast(beddir, snames = snames, type = "all")
#' }
 
create_loops_object_fast<- function(beddir, snames = NA, type = "all") {
  
  # Set data.table parameters for optimal performance
  data.table::setDTthreads(0)
  options(datatable.optimize = Inf)
  
  # Start timing
  total_start <- Sys.time()
  
  # Filter PET interactions based on type
  restrictPets <- function(bt, type) {
    if (type == "intra") {
      return(bt[bt[[1]] == bt[[4]], ])
    } else if (type == "inter") {
      return(bt[bt[[1]] != bt[[4]], ])
    } else {
      return(bt)
    }
  }
  
  # Fast file reading with automatic delimiter detection
  read_bedpe_file <- function(file) {
    # Auto-detect delimiter from first line
    first_line <- readLines(file, n = 1, warn = FALSE)
    delim <- if (grepl("\t", first_line)) "\t" else " "
    
    # Use fread for file reading
    dt <- data.table::fread(
      file, 
      sep = delim,
      header = FALSE,
      colClasses = c("character", "integer", "integer", "character", 
                     "integer", "integer", "character", "integer"),
      showProgress = FALSE,
      verbose = FALSE
    )
    
    # Set column names
    data.table::setnames(dt, paste0("X", 1:8))
    return(dt)
  }
  
  # Handle file paths and sample names
  if (is.null(snames) || all(is.na(snames))) {
    files <- dir(beddir, pattern = "\\.loop_counts\\.bedpe$", full.names = TRUE)
    snames <- tools::file_path_sans_ext(basename(files))
  } else {
    files <- file.path(beddir, paste0(snames, ".loop_counts.bedpe"))
  }
  
  # PHASE 1: File reading
  # cat("PHASE 1: Reading files...\n")
  PHASE1_start <- Sys.time()
  
  file_data <- lapply(seq_along(files), function(i) {
    file_start <- Sys.time()
    dt <- read_bedpe_file(files[i])
    file_time <- as.numeric(difftime(Sys.time(), file_start, units = "secs"))
    cat(sprintf("  File %d: %.2fs (%d rows, %.1f MB)\n", 
                i, file_time, nrow(dt), object.size(dt)/1024^2))
    dt
  })
  
  PHASE1_time <- as.numeric(difftime(Sys.time(), PHASE1_start, units = "secs"))
  # cat(sprintf("PHASE 1 TOTAL: %.2fs\n\n", PHASE1_time))
  
  # PHASE 2: Extract unique anchor regions from all files
  # cat("PHASE 2: Extracting anchors...\n")
  PHASE2_start <- Sys.time()
  
  anchorsData_list <- lapply(seq_along(file_data), function(i) {
    extract_start <- Sys.time()
    dt <- file_data[[i]]
    # Apply intra-chromosomal filtering for anchor extraction
    dt <- restrictPets(dt, type="intra")
    
    # Extract anchor coordinates from both interaction ends
    tmp <- dt[, c("X4", "X5", "X6")]
    names(tmp) <- names(dt[, c("X1", "X2", "X3")])
    result <- rbind(dt[, c("X1", "X2", "X3")], tmp)
    
    extract_time <- as.numeric(difftime(Sys.time(), extract_start, units = "secs"))
    # cat(sprintf("  Extract %d: %.2fs (%d -> %d anchors)\n", 
    #             i, extract_time, nrow(dt), nrow(result)))
    result
  })
  
  # Combine and deduplicate all anchors
  combine_start <- Sys.time()
  anchorsData <- unique(data.table::rbindlist(anchorsData_list))
  combine_time <- as.numeric(difftime(Sys.time(), combine_start, units = "secs"))
  cat(sprintf("  Combine unique: %.2fs (%d unique anchors)\n", 
              combine_time, nrow(anchorsData)))
  
  PHASE2_time <- as.numeric(difftime(Sys.time(), PHASE2_start, units = "secs"))
  # cat(sprintf("PHASE 2 TOTAL: %.2fs\n\n", PHASE2_time))
  
  # PHASE 3: Create GRanges object and precompute lookup structures
  # cat("PHASE 3: Creating anchor GRanges and optimization structures...\n")
  PHASE3_start <- Sys.time()
  
  anchors <- GenomicRanges::makeGRangesFromDataFrame(
    anchorsData, ignore.strand = TRUE,
    seqnames.field = "X1", start.field = "X2", end.field = "X3"
  )
  anchors <- sort(anchors)
  
  # OPTIMIZATION 1: Precompute anchor lookup table for all samples
  lookup_start <- Sys.time()
  anchor_lookup <- data.table::data.table(
    seqnames = as.character(GenomicRanges::seqnames(anchors)),
    start = GenomicRanges::start(anchors),
    end = GenomicRanges::end(anchors),
    anchor_id = seq_along(anchors)
  )
  data.table::setkey(anchor_lookup, seqnames, start, end)
  lookup_time <- as.numeric(difftime(Sys.time(), lookup_start, units = "secs"))
  # cat(sprintf("  Anchor lookup table: %.2fs\n", lookup_time))
  
  # OPTIMIZATION 2: Precompute anchor midpoints for distance calculation
  midpoint_start <- Sys.time()
  anchor_midpoints <- (GenomicRanges::start(anchors) + GenomicRanges::end(anchors)) / 2
  midpoint_time <- as.numeric(difftime(Sys.time(), midpoint_start, units = "secs"))
  cat(sprintf("  Anchor midpoints: %.2fs\n", midpoint_time))
  
  PHASE3_time <- as.numeric(difftime(Sys.time(), PHASE3_start, units = "secs"))
  # cat(sprintf("PHASE 3 TOTAL: %.2fs\n\n", PHASE3_time))
  
  # Optimized getpets function with precomputed lookup
  getpets_optimized <- function(bt_data, sample_name, anchor_lookup) {
    if (nrow(bt_data) == 0) {
      result <- data.table::data.table(left = integer(0), right = integer(0), counts = integer(0))
      data.table::setnames(result, c("left", "right", sample_name))
      return(result)
    }
    
    # Create query data.tables (no need to recreate anchor_dt each time)
    query_start <- Sys.time()
    left_query_dt <- data.table::data.table(
      seqnames = bt_data$X1,
      start = bt_data$X2,
      end = bt_data$X3,
      query_id = seq_len(nrow(bt_data))
    )
    right_query_dt <- data.table::data.table(
      seqnames = bt_data$X4,
      start = bt_data$X5,
      end = bt_data$X6,
      query_id = seq_len(nrow(bt_data))
    )
    
    data.table::setkey(left_query_dt, seqnames, start, end)
    data.table::setkey(right_query_dt, seqnames, start, end)
    query_time <- as.numeric(difftime(Sys.time(), query_start, units = "secs"))
    cat(sprintf("    Query data.tables: %.2fs\n", query_time))
    
    # Use precomputed anchor lookup for overlap detection
    overlap_start <- Sys.time()
    left_overlaps <- left_query_dt[anchor_lookup, nomatch = 0, on = .(seqnames, start, end)]
    right_overlaps <- right_query_dt[anchor_lookup, nomatch = 0, on = .(seqnames, start, end)]
    
    # If no exact matches found, use interval overlap
    if (nrow(left_overlaps) == 0 || nrow(right_overlaps) == 0) {
      # Use foverlaps for interval overlap
      left_overlaps <- data.table::foverlaps(left_query_dt, anchor_lookup, nomatch = 0)
      right_overlaps <- data.table::foverlaps(right_query_dt, anchor_lookup, nomatch = 0)
      
      # Apply minoverlap = 100 filter
      if (nrow(left_overlaps) > 0) {
        left_width <- pmin(left_overlaps$end, left_overlaps$i.end) - 
          pmax(left_overlaps$start, left_overlaps$i.start) + 1
        left_overlaps <- left_overlaps[left_width >= 100]
      }
      if (nrow(right_overlaps) > 0) {
        right_width <- pmin(right_overlaps$end, right_overlaps$i.end) - 
          pmax(right_overlaps$start, right_overlaps$i.start) + 1
        right_overlaps <- right_overlaps[right_width >= 100]
      }
    }
    overlap_time <- as.numeric(difftime(Sys.time(), overlap_start, units = "secs"))
    # cat(sprintf("    Overlap detection: %.2fs\n", overlap_time))
    
    # Create mapping vectors and aggregate
    mapping_start <- Sys.time()
    leftanchor <- rightanchor <- rep(NA_integer_, nrow(bt_data))
    
    if (nrow(left_overlaps) > 0) {
      leftanchor[left_overlaps$query_id] <- left_overlaps$anchor_id
    }
    if (nrow(right_overlaps) > 0) {
      rightanchor[right_overlaps$query_id] <- right_overlaps$anchor_id
    }
    
    # Create interaction data table with counts
    dt <- data.table::data.table(left = leftanchor, right = rightanchor, counts = bt_data[[8]])
    dt <- dt[, .(counts = sum(counts)), by = .(left, right)]
    data.table::setnames(dt, c("left", "right", sample_name))
    mapping_time <- as.numeric(difftime(Sys.time(), mapping_start, units = "secs"))
    cat(sprintf("    Mapping & aggregation: %.2fs\n", mapping_time))
    
    return(dt)
  }
  
  # PHASE 4: Process each sample file with optimized function
  # cat("PHASE 4: Processing samples (optimized)...\n")
  PHASE4_start <- Sys.time()
  
  petlist <- lapply(seq_along(snames), function(i) {
    sample_start <- Sys.time()
    sample_name <- snames[i]
    
    # Use already read data
    bt <- file_data[[i]]
    bt <- restrictPets(bt, type)  # Apply specified interaction type filter
    
    cat(sprintf("  Sample %d (%s): %d interactions\n", i, sample_name, nrow(bt)))
    
    # Use optimized function with precomputed lookup
    result <- getpets_optimized(bt, sample_name, anchor_lookup)
    
    sample_time <- as.numeric(difftime(Sys.time(), sample_start, units = "secs"))
    cat(sprintf("  Sample %d TOTAL: %.2fs\n\n", i, sample_time))
    
    result
  })
  
  PHASE4_time <- as.numeric(difftime(Sys.time(), PHASE4_start, units = "secs"))
  # cat(sprintf("PHASE 4 TOTAL: %.2fs\n\n", PHASE4_time))
  
  # PHASE 5: Merge all sample data
  # cat("PHASE 5: Merging data...\n")
  PHASE5_start <- Sys.time()
  
  pets <- data.table::rbindlist(petlist, fill = TRUE)
  data.table::setkey(pets, left, right)
  
  PHASE5_time <- as.numeric(difftime(Sys.time(), PHASE5_start, units = "secs"))
  # cat(sprintf("PHASE 5 TOTAL: %.2fs\n\n", PHASE5_time))
  
  # PHASE 6: Create standardized interaction matrix
  # cat("PHASE 6: Creating matrices...\n")
  PHASE6_start <- Sys.time()
  
  iraw <- data.table::copy(pets[, .(left, right)])
  # Ensure consistent ordering (left anchor index <= right anchor index)
  iraw[get("left") > get("right"), c("left", "right") := .(get("right"), get("left"))]
  data.table::setkey(iraw, left, right)
  interactions <- as.matrix(iraw)
  
  # Generate count matrix
  counts <- as.matrix(pets[, !c("left", "right"), with = FALSE])
  counts[is.na(counts)] <- 0
  colnames(counts) <- snames
  
  PHASE6_time <- as.numeric(difftime(Sys.time(), PHASE6_start, units = "secs"))
  # cat(sprintf("PHASE 6 TOTAL: %.2fs\n\n", PHASE6_time))
  
  # PHASE 7: Optimized distance calculation using precomputed midpoints
  # cat("PHASE 7: Computing distances (optimized)...\n")
  PHASE7_start <- Sys.time()
  
  # Use precomputed midpoints for highly efficient distance calculation
  w <- abs(anchor_midpoints[interactions[, 2]] - anchor_midpoints[interactions[, 1]])
  w[w < 0] <- 0  # Safety check, though shouldn't be needed
  rowData <- data.frame(loopWidth = as.integer(w))
  
  PHASE7_time <- as.numeric(difftime(Sys.time(), PHASE7_start, units = "secs"))
  # cat(sprintf("PHASE 7 TOTAL: %.2fs\n\n", PHASE7_time))
  
  # PHASE 8: Compute normalization size factors
  # cat("PHASE 8: Computing size factors...\n")
  PHASE8_start <- Sys.time()
  
  groups <- rep("group1", length(snames))
  if(length(snames) == 1) {
    sizeFactor <- 1
  } else {
    # Calculate median-based size factors for normalization
    lc <- log2(counts)
    keep <- rowSums(counts > 0) == ncol(lc)
    lc <- lc[keep, ]
    target <- 2^rowMeans(lc)
    sizeFactor <- matrixStats::colMedians(sweep(2^lc, 1, target, FUN = "/"), na.rm = TRUE)
  }
  dfcd <- data.frame(sizeFactor, groups)
  rownames(dfcd) <- snames
  
  dlo <- methods::new("loops")
  methods::slot(dlo, "anchors", check = TRUE) <- anchors
  methods::slot(dlo, "interactions", check = TRUE) <- interactions
  methods::slot(dlo, "counts", check = TRUE) <- counts
  methods::slot(dlo, "colData", check = TRUE) <- dfcd
  methods::slot(dlo, "rowData", check = TRUE) <- rowData
  
  # Clean up precomputed objects to free memory
  rm(anchor_lookup, anchor_midpoints)
  invisible(gc())
  
  # Final timing summary
  total_time <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  
  cat(paste(rep("=", 50), collapse = ""), "\n")
  cat("TIMING SUMMARY \n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  cat(sprintf("1. File reading:      %8.2fs (%5.1f%%)\n", PHASE1_time, 100*PHASE1_time/total_time))
  cat(sprintf("2. Anchor extraction: %8.2fs (%5.1f%%)\n", PHASE2_time, 100*PHASE2_time/total_time))
  cat(sprintf("3. Anchor GRanges:    %8.2fs (%5.1f%%)\n", PHASE3_time, 100*PHASE3_time/total_time))
  cat(sprintf("4. Sample processing: %8.2fs (%5.1f%%)\n", PHASE4_time, 100*PHASE4_time/total_time))
  cat(sprintf("5. Data merging:      %8.2fs (%5.1f%%)\n", PHASE5_time, 100*PHASE5_time/total_time))
  cat(sprintf("6. Matrix creation:   %8.2fs (%5.1f%%)\n", PHASE6_time, 100*PHASE6_time/total_time))
  cat(sprintf("7. Distance calc:     %8.2fs (%5.1f%%)\n", PHASE7_time, 100*PHASE7_time/total_time))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("TOTAL EXECUTION TIME: %8.2fs\n", total_time))
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  
  return(dlo)
}