#' @useDynLib sintHiChIP, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom readr read_delim
#' @importFrom GenomeInfoDb seqlevels "seqlevels<-"
#' @importFrom stats complete.cases median p.adjust pbinom predict smooth.spline aggregate 
#' @importFrom utils write.table
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame findOverlaps seqnames start end mcols<-
#' @importFrom IRanges IRanges
#' @importFrom data.table fread setDT data.table := setorder setnames as.data.table
#' @importFrom graphics plot
#' @importFrom methods new slot "slot<-"
#' @importFrom S4Vectors mcols
NULL

#' Identify Significant HiChIP Loops
#'
#' This function generates significant loops from HiChIP data.
#'
#' @param sname Sample name. No default value.
#' @param cwd Current working directory. No default value.
#' @param outdir Output directory. No default value.
#' @param normSiteFile Normalization site file. Should be a BED format file containing normalization cut site density information. No default value.
#' @param local Use local method. Default is TRUE.
#' @param FDR False Discovery Rate. Default is 0.01.
#' @param nbins Number of bins. Default is 10.
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' sintHiChIP_loops(sname = "sample1",
#'                  cwd = "/home/user/hichip/working/directory",
#'                  outdir = "/home/user/hichip/output",
#'                  normSiteFile = "/home/user/hichip/normalization/file.bed",
#'                  local = TRUE,
#'                  FDR = 0.01,
#'                  nbins = 10)
#' }
sintHiChIP_loops <- function(sname, cwd, outdir, normSiteFile, local = TRUE, FDR = 0.01, nbins = 10) {
  # Debug print statements
  cat("Parameters:\n")
  cat("Sample:", sname, "\n")
  cat("CWD:", cwd, "\n")
  cat("Output Directory:", outdir, "\n")
  cat("Normalization Site File:", normSiteFile, "\n")
  cat("Local:", local, "\n")
  cat("FDR:", FDR, "\n")
  cat("Number of Bins:", nbins, "\n")
  
  # Suppress R CMD check notes for unbound global variables
  totalAnchorSites_1 <- NULL
  totalAnchorSites_2 <- NULL
  quality <- NULL
  ratio <- NULL
  site <- NULL
  score1 <- NULL
  score2 <- NULL
  loopWidth <- NULL
  P_IAB_distance <- NULL
  P_combos_distance <- NULL
  P_IAB_site <- NULL
  P_combos_site <- NULL
  PETS <- NULL
  p_binom_site <- NULL
  P_site <- NULL
  Q_site <- NULL
  width <- NULL
  strand <- NULL
  # Create loop bins and calculate total PET counts
  lo <- create_loops_object_fast(cwd, snames = paste0(sname, ".filt.intra"), type = "all")
  totalPetCounts <- rowSums(lo@counts)
  
  # Load and process normalization site file
  GenomeBin <- read_delim(normSiteFile, col_names = TRUE)
  GenomeBin <- convert_df_to_granges(GenomeBin, seqnames.field = "seqnames", start.field = "start", end.field = "end")
  genomeBin_seqname <- seqlevels(GenomeBin)
  chrname <- unique(intersect(seqlevels(lo@anchors), seqnames(GenomeBin)))
  seqlevels(GenomeBin, pruning.mode = 'coarse') <- chrname
  
  lAnchor <- lo@anchors[lo@interactions[, 1]]
  rAnchor <- lo@anchors[lo@interactions[, 2]]
  site_l <- calculate_cut_site_density(lAnchor, GenomeBin)
  site_r <- calculate_cut_site_density(rAnchor, GenomeBin)
  lAnchor_site <- cbind(as.data.frame(lo@anchors[lo@interactions[, 1]]), totalAnchorSites = site_l)
  rAnchor_site <- cbind(as.data.frame(lo@anchors[lo@interactions[, 2]]), totalAnchorSites = site_r)
  sAnchors_site <- subset(rbind(lAnchor_site, rAnchor_site), select = -c(width, strand))
  tc.site <- GRanges(stats::aggregate(totalAnchorSites ~ seqnames + start + end, data = sAnchors_site, sum, na.rm = TRUE))
  tc.site <- sort(tc.site)
  
  # Prepare chrpeaks for pairing calculations
  chrpeaks <- as.data.frame(tc.site)
  chrpeaks <- chrpeaks[, c(1, 2, 3, 4, 6, 5)]
  names(chrpeaks) <- c("chr", "start", "end", "name", "score", "strand")
  
  mcols(lo@anchors) <- mcols(tc.site)
  df <- cbind(as.data.frame(lo@anchors[lo@interactions[, 1]])[,-c(4:5)],
              as.data.frame(lo@anchors[lo@interactions[, 2]])[,-c(4:5)])
  setDT(df)
  
  setnames(df, c("chr_1", "start_1", "end_1", "totalAnchorSites_1", "chr_2", "start_2", "end_2", "totalAnchorSites_2"))
  df[, loopWidth := as.integer(lo@rowData$loopWidth)]
  df[, PETS := totalPetCounts]
  df[, site := totalAnchorSites_1 * totalAnchorSites_2]
  df[, ratio := calculate_density_ratio(df)]
  
  # Assign quality tags based on ratio
  df[, quality := cut(ratio,
                      breaks = c(0, 0.1, 0.25, 0.5, 1),
                      labels = c("Low", "Medium", "High", "Very High"),
                      include.lowest = TRUE)]
  
  # Filter data based on quality
  if (!local) {
    df <- df[quality %in% c("High", "Very High")]
  }
  
  totalcombos <- 0
  chromosomes <- sort(unique(as.character(chrpeaks$chr)))
  
  # Create distance and site bins
  distanceborders <- create_equal_bins(df$loopWidth, binmethod = "equalocc", numberbins = nbins)
  siteborders <- create_equal_bins(df$site, binmethod = "equalsize", numberbins = nbins)
  
  # Fit distance and site models
  distance_IAB_model <- model_hichip(df$loopWidth, df$PETS, borders = distanceborders, yvals = TRUE)
  distance_IAB_spline <- smooth.spline(log10(distance_IAB_model[, 1]), distance_IAB_model[, 3], spar = 0.75)
  plot(distance_IAB_spline$x, distance_IAB_spline$y)
  site_IAB_model <- model_hichip(df$site, df$PETS, borders = siteborders, yvals = TRUE)
  site_IAB_model.complete <- site_IAB_model[complete.cases(site_IAB_model), ]
  site_IAB_spline <- smooth.spline(log10(site_IAB_model.complete[, 1]), site_IAB_model.complete[, 3], spar = 0.75)
  plot(site_IAB_spline$x, site_IAB_spline$y)
  
  # Initialize variables
  sumofy_dist <- rep(0, nbins)
  pvals_dist <- rep(0, nbins)
  sumofx_dist <- rep(0, nbins)
  countofx_dist <- rep(0, nbins)
  
  meanofx_site <- rep(0, nbins)
  sumofy_site <- rep(0, nbins)
  pvals_site <- rep(0, nbins)
  sumofx_site <- rep(0, nbins)
  countofx_site <- rep(0, nbins)
  
  for (chrom in chromosomes) {
    cat(sprintf("\n ****** Processing combinations for %s\n", chrom))
    chrom_peaks <- chrpeaks[chrpeaks$chr == chrom, ]
    chrom_peaks <- chrom_peaks[order(chrom_peaks$start), ]
    combos <- makecombos(chrom, chrom_peaks, mindist = 0, maxdist = 2000000)
    combos <- data.table(combos)
    combos[, site := score1 * score2]
    
    # Calculate distance model for chromosome
    distance_combo_model_chrom <- model_hichip(x = combos$dist, y = NA, borders = distanceborders, yvals = FALSE)
    sumofy_dist <- sumofy_dist + distance_combo_model_chrom[, 2]
    sumofx_dist <- sumofx_dist + distance_combo_model_chrom[, 4]
    countofx_dist <- countofx_dist + distance_combo_model_chrom[, 5]
    
    # Calculate site model for chromosome
    site_combo_model_chrom <- model_hichip(x = combos$site, y = NA, borders = siteborders, yvals = FALSE)
    sumofy_site <- sumofy_site + site_combo_model_chrom[, 2]
    sumofx_site <- sumofx_site + site_combo_model_chrom[, 4]
    countofx_site <- countofx_site + site_combo_model_chrom[, 5]
  }
  
  # Combine data from all chromosomes
  site_combo_model <- cbind(sumofx_site / countofx_site, sumofy_site, sumofy_site / sum(sumofy_site))
  site_combo_model.complete <- site_combo_model[complete.cases(site_combo_model),]
  
  distance_combo_model <- cbind(sumofx_dist / countofx_dist, sumofy_dist, sumofy_dist / sum(sumofy_dist))
  distance_combo_model.complete <- distance_combo_model[complete.cases(distance_combo_model),]
  
  # Generate distance and site splines
  distance_combo_spline <- smooth.spline(log10(distance_combo_model.complete[, 1]), distance_combo_model.complete[, 3], spar = 0.75)
  plot(distance_combo_spline$x, distance_combo_spline$y)
  site_combo_spline <- smooth.spline(log10(site_combo_model.complete[, 1]), site_combo_model.complete[, 3], spar = 0.75)
  plot(site_combo_spline$x, site_combo_spline$y)
  
  # Predict and adjust values
  df[, P_IAB_distance := predict(distance_IAB_spline, log10(df$loopWidth))$y]
  df[, P_combos_distance := predict(distance_combo_spline, log10(df$loopWidth))$y]
  df[, P_IAB_site := predict(site_IAB_spline, log10(df$site))$y]
  df[, P_combos_site := predict(site_combo_spline, log10(df$site))$y]
  
  # Limit values to minimum and maximum
  df[P_IAB_distance <= min(distance_IAB_model[, 3]), P_IAB_distance := min(distance_IAB_model[, 3])]
  df[P_IAB_distance >= max(distance_IAB_model[, 3]), P_IAB_distance := max(distance_IAB_model[, 3])]
  df[P_combos_distance <= min(distance_combo_model[, 3]), P_combos_distance := min(distance_combo_model[, 3])]
  df[P_combos_distance >= max(distance_combo_model[, 3]), P_combos_distance := max(distance_combo_model[, 3])]
  df[P_IAB_site <= min(site_IAB_model[, 3]), P_IAB_site := min(site_IAB_model[, 3])]
  df[P_IAB_site >= max(site_IAB_model[, 3]), P_IAB_site := max(site_IAB_model[, 3])]
  df[P_combos_site <= min(site_combo_model[, 3]), P_combos_site := min(site_combo_model[, 3])]
  df[P_combos_site >= max(site_combo_model[, 3]), P_combos_site := max(site_combo_model[, 3])]
  
  # Calculate binomial probability
  totalcombos <- sum(sumofy_dist)
  df[, p_binom_site := (P_IAB_distance * P_IAB_site) / (P_combos_distance * P_combos_site * totalcombos)]
  df[, p_binom_site := pmin(1, pmax(0, p_binom_site))]
  
  # Calculate final interaction P value
  totalIAB <- sum(distance_IAB_model[, 2])
  df[, P_site := 1 - pbinom(PETS - 1, size = totalIAB, prob = p_binom_site)]
  df[, Q_site := p.adjust(P_site, method = "BH")]
  
  df <- df[stats::complete.cases(df), ]
  data.table::setorder(df, Q_site)
  df_filt <- df[, .SD, .SDcols = c("chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "PETS", "Q_site")]
  df_filt <- df_filt[Q_site < FDR]
  data.table::fwrite(df_filt, file = file.path(outdir, paste0(sname, ".interaction.Q", FDR, ".txt")),
                     sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cat("Processing complete. Output saved to: ", outdir, "\n")
  
  invisible(NULL)
}

#' Calculate Density Ratio
#'
#' This function is to calculate the ratio of two anchor cut site density values. 
#' @param interactions Data frame of interactions
#' @return Vector of ratio densities
#' @keywords internal
#' @examples
#' \dontrun{
#' interactions <- data.frame(
#'   totalAnchorSites_1 = c(10, 20, 30, 40),
#'   totalAnchorSites_2 = c(15, 25, 35, 20)
#' )
#' ratios <- calculate_density_ratio(interactions)
#' print(ratios)
#' }
calculate_density_ratio <- function(interactions) {
  pmin(interactions$totalAnchorSites_1, interactions$totalAnchorSites_2) /
    pmax(interactions$totalAnchorSites_1, interactions$totalAnchorSites_2)
}

#' Convert DataFrame to GRanges Object
#'
#' @param df Data frame to convert
#' @param seqnames.field Column name for sequence names
#' @param start.field Column name for start positions
#' @param end.field Column name for end positions
#' @param ignore.strand Whether to ignore strand information
#' @return GRanges object
#' @keywords internal
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' df <- data.frame(
#'   chr = c("chr1", "chr1", "chr2"),
#'   start = c(100, 200, 300),
#'   end = c(150, 250, 350),
#'   score = c(1, 2, 3)
#' )
#' gr <- convert_df_to_granges(df, "chr", "start", "end")
#' print(gr)
#' }
convert_df_to_granges <- function(df, seqnames.field, start.field, end.field, ignore.strand = TRUE) {
  dt <- data.table::as.data.table(df)
  gr <- GenomicRanges::GRanges(seqnames = dt[[seqnames.field]],
                               ranges = IRanges::IRanges(start = dt[[start.field]], end = dt[[end.field]]))
  if (!ignore.strand) {
    gr$strand <- dt[["strand"]]
  }
  other_columns <- setdiff(names(df), c(seqnames.field, start.field, end.field))
  for (col in other_columns) {
    S4Vectors::mcols(gr)[[col]] <- dt[[col]]
  }
  return(gr)
}

