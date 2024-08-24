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
#' bin_borders <- create_equal_bins(vector_to_bin, binmethod = "equalsize", numberbins = 10)
#' }
create_equal_bins <- function(vectortobin, binmethod="equalsize", numberbins=30) {
  sortedvec <- sort(vectortobin)
  if (binmethod == "equalsize") {
    borders <- seq(min(sortedvec), max(sortedvec), length.out = numberbins + 1)[2:numberbins]
    return(borders)
  }
  if (binmethod == "equalocc") {
    itemsperbin <- length(vectortobin)/numberbins
    ranks <- rank(sortedvec)
    bins <- ceiling(ranks/itemsperbin)
    mins <- tapply(sortedvec, bins, min)
    maxes <- tapply(sortedvec, bins, max)
    borders <- (mins[-1] + maxes[-length(maxes)])/2
    return(borders)
  }
}

#' Calculate Cut Site Density
#'
#' This function calculates restriction enzyme cut site density information for given bins.
#'
#' @param bins Bins GRanges object
#' @param GRfile Genome GRanges object
#' @return Vector of average scores
#' @export
#' @examples
#' \dontrun{
#' genome_file <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1:100, width = 1)
#' )
#' bins <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = seq(1, 100, by = 10), width = 10)
#' )
#' site_info <- calculate_cut_site_density(bins, genome_file)
#' }
calculate_cut_site_density <- function(bins, GRfile) {
  signal <- sort(GenomeInfoDb::sortSeqlevels(GRfile))
  score <- GenomicRanges::coverage(signal, weight="MeanSite")
  bins <- GenomeInfoDb::renameSeqlevels(bins, names(score))
  binned_data <- GenomicRanges::binnedAverage(bins, score, "signal")
  bin_avg <- S4Vectors::mcols(binned_data)$signal
  return(bin_avg)
}

#' Create Loops Object Fast
#'
#' This function creates binned and peak-merged loops object from BEDPE files quickly.
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
#' snames <- c("sample1", "sample2")
#' loops <- create_loops_object_fast(beddir, snames = snames, type = "all")
#' }
create_loops_object_fast <- function(beddir, snames = NA, type = "all") {
  ct <- list(readr::col_character(), readr::col_integer(), readr::col_integer(),
             readr::col_character(), readr::col_integer(), readr::col_integer(), 
             readr::col_character(), readr::col_integer())
  
  restrictPets <- function(bt, type) {
    if (type == "intra") {
      return(bt[bt[[1]] == bt[[4]], ])
    } else if (type == "inter") {
      return(bt[bt[[1]] != bt[[4]], ])
    } else {
      return(bt)
    }
  }
  
  if (is.null(snames) || all(is.na(snames))) {
    files <- dir(beddir, pattern = "\\.loop_counts\\.bedpe$", full.names = TRUE)
    snames <- tools::file_path_sans_ext(basename(files))
  } else {
    files <- file.path(beddir, paste0(snames, ".loop_counts.bedpe"))
  }
  
  anchorsData <- unique(data.table::rbindlist(lapply(files, function(file) {
    dt <- readr::read_delim(file, delim = " ", col_names = paste0("X", 1:8), col_types = ct)
    dt <- restrictPets(dt, type="intra")
    tmp <- dt[, c("X4", "X5", "X6")]
    names(tmp) <- names(dt[, c("X1", "X2", "X3")])
    rbind(dt[, c("X1", "X2", "X3")], tmp)
  })))
  
  anchors <- GenomicRanges::makeGRangesFromDataFrame(unique(anchorsData), ignore.strand = TRUE,
                                                     seqnames.field = "X1", start.field = "X2", end.field = "X3")
  
  anchors <- sort(anchors)
  
  getpets <- function(left, right, counts, sample_name, anchors) {
    ovl_left <- GenomicRanges::findOverlaps(left, anchors, minoverlap = 100)
    ovl_right <- GenomicRanges::findOverlaps(right, anchors, minoverlap = 100)
    
    leftanchor <- rightanchor <- rep(NA_integer_, length(left))
    
    leftanchor[S4Vectors::queryHits(ovl_left)] <- S4Vectors::subjectHits(ovl_left)
    rightanchor[S4Vectors::queryHits(ovl_right)] <- S4Vectors::subjectHits(ovl_right)
    
    dt <- data.table::data.table(left = leftanchor, right = rightanchor, counts = counts[[8]])
    dt <- dt[, list(counts = sum(counts)), by = list(left, right)]
    data.table::setnames(dt, c("left", "right", sample_name))
    return(dt)
  }
  
  sample_name <- left <- right <- NULL
  petlist <- lapply(snames, function(sample_name) {
    fullfile <- file.path(beddir, paste0(sample_name, ".loop_counts.bedpe"))
    bt <- readr::read_delim(fullfile, delim = " ", col_names = FALSE, col_types = ct)
    bt <- restrictPets(bt, type="intra")
    getpets(GenomicRanges::makeGRangesFromDataFrame(bt[, 1:3], ignore.strand = TRUE,
                                                    seqnames.field = "X1", start.field = "X2", end.field = "X3"),
            GenomicRanges::makeGRangesFromDataFrame(bt[, 4:6], ignore.strand = TRUE,
                                                    seqnames.field = "X4", start.field = "X5", end.field = "X6"),
            bt, sample_name, anchors)
  })
  
  pets <- data.table::rbindlist(petlist, fill = TRUE)
  data.table::setkey(pets, left, right)
  
  iraw <- data.table::copy(pets[, list(left, right)])
  iraw[get("left") > get("right"), c("left", "right") := list(get("right"), get("left"))]
  data.table::setkey(iraw, left, right)
  
  interactions <- as.matrix(iraw)
  
  counts <- as.matrix(pets[, !c("left", "right"), with = FALSE])
  counts[is.na(counts)] <- 0
  colnames(counts) <- snames
  
  w <- abs((GenomicRanges::start(anchors[interactions[, 2]]) + GenomicRanges::end(anchors[interactions[, 2]]))/2 -
             (GenomicRanges::start(anchors[interactions[, 1]]) + GenomicRanges::end(anchors[interactions[, 1]]))/2)
  w[w < 0] <- 0
  rowData <- data.frame(loopWidth = as.integer(w))
  
  groups <- rep("group1", length(snames))
  if(length(snames) == 1){
    sizeFactor <- 1
  } else {
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
  
  return(dlo)
}