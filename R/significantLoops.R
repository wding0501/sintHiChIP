#' Identify Significant HiChIP Loops
#'
#' This function generates significant loops from HiChIP data with statistical models
#' including overdispersion detection and correction for massive datasets.
#'
#' @param sname Character string. Sample name. No default value.
#' @param cwd Character string. Current working directory. No default value.
#' @param outdir Character string. Output directory. No default value.
#' @param normSiteFile Character string. Normalization site file. Should be a BED format file containing normalization cut site density information. No default value.
#' @param local Logical. Use local method. Default is TRUE.
#' @param FDR Numeric. False Discovery Rate. Default is 0.01.
#' @param nbins Numeric. Number of bins. Default is 10.
#' 
#' @return Invisible NULL
#' 
#' @details 
#' This version includes:
#' \itemize{
#'   \item Anchor site density calculation with unique anchor pre-computation
#'   \item Statistical modeling (Binomial/Negative Binomial/Poisson) based on overdispersion detection
#'   \item Numerical stability with log-space calculations and proper handling of edge cases
#'   \item Memory-efficient batch processing for massive datasets (40M+ interactions)
#'   \item data.table operations throughout
#'   \item Handling of zero and negative site values
#' }
#'
#' @importFrom data.table fread setDT setnames setkey setorder as.data.table foverlaps fwrite
#' @importFrom GenomicRanges GRanges seqnames findOverlaps mcols width start end
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols subjectHits queryHits
#' @importFrom stats smooth.spline predict approx var pbinom pnbinom ppois pnorm p.adjust
#' @importFrom MASS glm.nb
#' 
#' @examples
#' \dontrun{
#' # Basic usage
#' sintHiChIP_sigloops(
#'   sname = "sample1",
#'   cwd = "/path/to/data",
#'   outdir = "/path/to/output",
#'   normSiteFile = "/path/to/norm.bed"
#' )
#' 
#' # With custom parameters
#' sintHiChIP_sigloops(
#'   sname = "sample1",
#'   cwd = "/path/to/data", 
#'   outdir = "/path/to/output",
#'   normSiteFile = "/path/to/norm.bed",
#'   local = FALSE,
#'   FDR = 0.05,
#'   nbins = 15
#' )
#' }
#'
#' @export
sintHiChIP_sigloops <- function(sname, cwd, outdir, normSiteFile, local = TRUE, FDR = 0.01, nbins = 10) {
  
  # ================== Parameter Validation ==================
  if (missing(sname) || is.null(sname) || nchar(sname) == 0) {
    stop("Sample name (sname) is required and cannot be empty")
  }
  if (missing(cwd) || !dir.exists(cwd)) {
    stop("Working directory (cwd) must exist: ", cwd)
  }
  if (missing(outdir) || is.null(outdir) || nchar(outdir) == 0) {
    outdir <- cwd
  }
  if (!file.exists(normSiteFile)) {
    stop("Normalization site file not found: ", normSiteFile)
  }
  if (!is.numeric(FDR) || FDR <= 0 || FDR >= 1) {
    stop("FDR must be a numeric value between 0 and 1")
  }
  if (!is.numeric(nbins) || nbins < 2) {
    stop("nbins must be a numeric value >= 2")
  }
  bedpe_file <- file.path(cwd, paste0(sname, ".filt.intra.loop_counts.bedpe"))
  if (!file.exists(bedpe_file)) {
    stop("Required BEDPE file not found: ", bedpe_file)
  }
  if (missing(local) || is.function(local)) {
    local <- TRUE  
    cat("Fixed local parameter conflict, set to TRUE\n")
  }
  # Suppress R CMD check notes
  totalAnchorSites_1 <- NULL
  totalAnchorSites_2 <- NULL
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
  mu <- NULL
  mu_safe <- NULL
  lambda <- NULL
  
  mode <- if(local) "local" else "global"
  
  cat("=== HiChIP loop detection ===\n")
  cat("Sample:", sname, "\n")
  cat("Mode:", mode, "\n")
  cat("FDR:", FDR, "\n")
  cat("nbins:", nbins, "\n\n")
  
  total_start_time <- Sys.time()
  
  # ================== Phase 1: Load and Process Interaction Data ==================
  step_start <- Sys.time()
  
  lo <- create_loops_object_fast(cwd, snames = paste0(sname, ".filt.intra"), type = "all")
  totalPetCounts <- rowSums(lo@counts)
  
  GenomeBin <- data.table::fread(normSiteFile, header = TRUE, showProgress = FALSE)
  GenomeBin <- sintHiChIP_create_granges_fast(GenomeBin, 
                                              names(GenomeBin)[1],
                                              names(GenomeBin)[2],
                                              names(GenomeBin)[3])
  
  chrname <- unique(intersect(seqlevels(lo@anchors), seqnames(GenomeBin)))
  seqlevels(GenomeBin, pruning.mode = 'coarse') <- chrname
  
  # Anchor site processing
  cat("Computing anchor site densities...\n")
  anchor_start <- Sys.time()
  
  # Get all anchor information
  n_anchors <- length(lo@anchors)
  n_interactions <- nrow(lo@interactions)
  n_genomebins <- length(GenomeBin)
  
  cat(sprintf("Processing %d anchors, %d interactions, %d genome bins\n", 
              n_anchors, n_interactions, n_genomebins))
  
  # Extract anchor coordinates efficiently
  anchors_df <- data.frame(
    seqnames = as.character(seqnames(lo@anchors)),
    start = start(lo@anchors),
    end = end(lo@anchors),
    anchor_id = seq_len(n_anchors),
    stringsAsFactors = FALSE
  )
  
  # Get unique anchors to reduce computation
  unique_anchors_df <- unique(anchors_df[, c("seqnames", "start", "end")])
  unique_anchors_df$unique_id <- seq_len(nrow(unique_anchors_df))
  
  cat(sprintf("Reduced to %d unique anchors (%.1f%% reduction)\n", 
              nrow(unique_anchors_df), 
              100 * (1 - nrow(unique_anchors_df) / n_anchors)))
  
  # Density calculation for unique anchors only
  if (nrow(unique_anchors_df) > 10000 || n_genomebins > 100000) {
    
    # Convert to data.table for maximum speed on large datasets
    anchors_dt <- data.table::as.data.table(unique_anchors_df)
    genomebins_dt <- data.table::as.data.table(as.data.frame(GenomeBin))
    
    # Add unique anchor ID for efficient matching
    anchors_dt[, anchor_id := .I]
    
    # Extract MeanSite weights from genomic bins
    genomebins_dt[, weight := MeanSite]  # Direct assignment since MeanSite exists
    
    # Set keys for efficient genomic overlap
    data.table::setkey(anchors_dt, seqnames, start, end)
    data.table::setkey(genomebins_dt, seqnames, start, end)
    
    # Perform genomic overlap using foverlaps
    # genomebins_dt intervals that fall within anchors_dt intervals
    overlaps <- data.table::foverlaps(genomebins_dt, anchors_dt, 
                                      type = "any", 
                                      nomatch = 0L)
    
    if (nrow(overlaps) > 0) {
      # Efficient aggregation by anchor regions
      # Group by anchor coordinates and sum weights within each anchor
      density_stats <- overlaps[, .(
        total_weight = sum(weight, na.rm = TRUE),
        bin_count = .N
      ), by = .(seqnames, start, end, anchor_id)]
      
      # Calculate density per base pair for each anchor
      density_stats[, width := end - start + 1]
      density_stats[, density := total_weight / width]
      
      # Initialize result vector with zeros
      unique_densities <- numeric(nrow(anchors_dt))
      
      # Efficient assignment using anchor_id as direct index
      unique_densities[density_stats$anchor_id] <- density_stats$density
      
    } else {
      # No overlaps found - all densities are zero
      unique_densities <- numeric(nrow(unique_anchors_df))
    }
    
  } else {
    
    # Method for smaller datasets using GenomicRanges
    unique_gr <- sintHiChIP_create_granges_fast(unique_anchors_df)
    
    # Overlap detection
    overlaps_hits <- GenomicRanges::findOverlaps(unique_gr, GenomeBin, type = "any")
    
    if (length(overlaps_hits) > 0) {
      # Extract MeanSite weights for overlapping bins
      weights <- mcols(GenomeBin)$MeanSite[subjectHits(overlaps_hits)]
      query_indices <- queryHits(overlaps_hits)
      
      # Efficient aggregation using data.table
      overlap_dt <- data.table::data.table(
        anchor_idx = query_indices,
        weight = weights
      )
      
      # Sum weights for each anchor
      weight_sums <- overlap_dt[, .(total_weight = sum(weight, na.rm = TRUE)), 
                                by = anchor_idx]
      
      # Calculate density = total_weight / anchor_width
      anchor_widths <- width(unique_gr)
      unique_densities <- numeric(length(unique_gr))
      unique_densities[weight_sums$anchor_idx] <- weight_sums$total_weight / anchor_widths[weight_sums$anchor_idx]
      
    } else {
      # No overlaps - all densities are zero
      unique_densities <- numeric(length(unique_gr))
    }
  }
  
  # Efficient lookup to map back to all anchors
  anchor_key <- paste(anchors_df$seqnames, anchors_df$start, anchors_df$end, sep = "_")
  unique_key <- paste(unique_anchors_df$seqnames, unique_anchors_df$start, unique_anchors_df$end, sep = "_")
  
  # Use match for O(n) lookup
  density_lookup <- unique_densities[match(anchor_key, unique_key)]
  
  # Extract interaction-specific densities
  site_l <- density_lookup[lo@interactions[, 1]]
  site_r <- density_lookup[lo@interactions[, 2]]
  
  # Efficient aggregation using data.table
  interaction_anchors_dt <- data.table::data.table(
    seqnames = c(anchors_df$seqnames[lo@interactions[, 1]], 
                 anchors_df$seqnames[lo@interactions[, 2]]),
    start = c(anchors_df$start[lo@interactions[, 1]], 
              anchors_df$start[lo@interactions[, 2]]),
    end = c(anchors_df$end[lo@interactions[, 1]], 
            anchors_df$end[lo@interactions[, 2]]),
    totalAnchorSites = c(site_l, site_r)
  )
  
  # Aggregate using data.table (fastest aggregation method)
  tc_dt <- interaction_anchors_dt[, .(totalAnchorSites = sum(totalAnchorSites, na.rm = TRUE)),
                                  by = .(seqnames, start, end)]
  
  # Convert to GRanges using method
  tc.site <- sort(sintHiChIP_create_granges_fast(tc_dt))
  
  anchor_time <- as.numeric(difftime(Sys.time(), anchor_start, units = "secs"))
  # cat(sprintf("Anchor processing completed in %.2f seconds\n", anchor_time))
  
  # Prepare chrpeaks for combination generation
  chrpeaks <- as.data.frame(tc.site)
  if (ncol(chrpeaks) >= 6) {
    chrpeaks <- chrpeaks[, c(1, 2, 3, 4, 6, 5)]
  } else {
    names(chrpeaks)[1:3] <- c("chr", "start", "end")
    chrpeaks$name <- seq_len(nrow(chrpeaks))
    chrpeaks$score <- tc.site$totalAnchorSites
    chrpeaks$strand <- "*"
    chrpeaks <- chrpeaks[, c("chr", "start", "end", "name", "score", "strand")]
  }
  names(chrpeaks) <- c("chr", "start", "end", "name", "score", "strand")
  
  mcols(lo@anchors) <- mcols(tc.site)
  df <- cbind(as.data.frame(lo@anchors[lo@interactions[, 1]])[, -c(4:5)],
              as.data.frame(lo@anchors[lo@interactions[, 2]])[, -c(4:5)])
  df <- data.table::setDT(df)
  data.table::setnames(df, c("chr_1", "start_1", "end_1", "totalAnchorSites_1",
                             "chr_2", "start_2", "end_2", "totalAnchorSites_2"))
  df[, loopWidth := as.integer(lo@rowData$loopWidth)]
  df[, PETS := totalPetCounts]
  df[, site := totalAnchorSites_1 * totalAnchorSites_2]
  
  step1_time <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  # cat(sprintf("Phase 1 completed in %.1f seconds\n", step1_time))
  
  # ================== Phase 2: Fit IAB Models ==================
  step_start <- Sys.time()
  
  distanceborders <- make_equal_bins(df$loopWidth, binmethod = "equalocc", numberbins = nbins)
  siteborders <- make_equal_bins(df$site, binmethod = "equalocc", numberbins = nbins)
  
  distance_IAB_model <- model_hichip(df$loopWidth, df$PETS, borders = distanceborders, yvals = TRUE)
  if (nrow(distance_IAB_model) < 2) stop("Distance IAB model returned insufficient points for spline")
  distance_IAB_spline <- stats::smooth.spline(log10(distance_IAB_model[, 1]), distance_IAB_model[, 3], spar = 0.35)
  
  site_IAB_model <- model_hichip(df$site, df$PETS, borders = siteborders, yvals = TRUE)
  site_IAB_model.complete <- site_IAB_model[complete.cases(site_IAB_model), , drop = FALSE]
  if (nrow(site_IAB_model.complete) < 2) stop("Site IAB model returned insufficient points for spline")
  site_IAB_spline <- stats::smooth.spline(log10(site_IAB_model.complete[, 1]), site_IAB_model.complete[, 3], spar = 0.35)
  
  step2_time <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  # cat(sprintf("Phase 2 completed in %.1f seconds\n", step2_time))
  
  # ================== Phase 3: Combination Processing ==================
  step_start <- Sys.time()
  
  # Initialize result vectors
  sumofy_dist <- rep(0, nbins)
  sumofx_dist <- rep(0, nbins)
  countofx_dist <- rep(0, nbins)
  sumofy_site <- rep(0, nbins)
  sumofx_site <- rep(0, nbins)
  countofx_site <- rep(0, nbins)
  
  chromosomes <- sort(unique(as.character(chrpeaks$chr)))
  # cat(sprintf("Processing %d chromosomes for combinations\n", length(chromosomes)))
  
  # Memory-efficient collection using lists
  collection_start <- Sys.time()
  
  # Use list collection to avoid memory pre-allocation issues
  distance_list <- vector("list", length(chromosomes))
  site_list <- vector("list", length(chromosomes))
  
  # Data structure for chromosome processing
  chrpeaks_dt <- data.table::as.data.table(chrpeaks)
  data.table::setkey(chrpeaks_dt, chr)
  
  # cat("Collecting combinations...\n")
  
  # Process each chromosome efficiently
  for (i in seq_along(chromosomes)) {
    chrom <- chromosomes[i]
    
    # Efficient chromosome data extraction using keyed lookup
    chrom_peaks <- chrpeaks_dt[chrom][order(start)]
    
    if (nrow(chrom_peaks) >= 2) {
      # Use makecombos function
      combos <- makecombos(chrom, chrom_peaks, mindist = 0, maxdist = 2000000)
      
      if (nrow(combos) > 0) {
        # Store directly in list (much more efficient than vector concatenation)
        distance_list[[i]] <- combos$dist
        site_list[[i]] <- combos$score1 * combos$score2
        
        # Report progress for large chromosomes
        if (nrow(combos) > 10000000) {
          cat(sprintf("  %s: %d combos\n", chrom, nrow(combos)))
        }
      }
    }
    
    # Progress reporting
    if (i %% 8 == 0) {
      cat(sprintf("  Progress: %d/%d chromosomes\n", i, length(chromosomes)))
    }
  }
  
  # Efficient merging using unlist
  # cat("Merging results...\n")
  merge_start <- Sys.time()
  
  # Remove empty elements
  distance_list <- distance_list[!sapply(distance_list, is.null)]
  site_list <- site_list[!sapply(site_list, is.null)]
  
  # Efficient merge operation
  all_distances <- unlist(distance_list, use.names = FALSE)
  all_sites <- unlist(site_list, use.names = FALSE)
  
  # Immediate cleanup
  rm(distance_list, site_list)
  invisible(gc())
  
  merge_time <- as.numeric(difftime(Sys.time(), merge_start, units = "secs"))
  collection_time <- as.numeric(difftime(Sys.time(), collection_start, units = "secs"))
  
  actual_combos <- length(all_distances)
  # cat(sprintf("Collection: %.2f s (%.2f s merge), %d combinations\n", collection_time, merge_time, actual_combos))
  
  # Batch model fitting (key: 2 calls instead of 46+ calls)
  modeling_start <- Sys.time()
  
  if (actual_combos > 0) {
    # Distance model (single call replacing multiple individual calls)
    distance_combo_model <- model_hichip(x = all_distances, y = NA, 
                                         borders = distanceborders, yvals = FALSE)
    
    # Site model (single call replacing multiple individual calls)
    site_combo_model <- model_hichip(x = all_sites, y = NA, 
                                     borders = siteborders, yvals = FALSE)
    
    # Extract results
    sumofy_dist <- distance_combo_model[, 2]
    sumofx_dist <- distance_combo_model[, 4]
    countofx_dist <- distance_combo_model[, 5]
    
    sumofy_site <- site_combo_model[, 2]
    sumofx_site <- site_combo_model[, 4]
    countofx_site <- site_combo_model[, 5]
  }
  
  modeling_time <- as.numeric(difftime(Sys.time(), modeling_start, units = "secs"))
  
  # Quick cleanup
  rm(all_distances, all_sites)
  if (exists("distance_combo_model")) rm(distance_combo_model)
  if (exists("site_combo_model")) rm(site_combo_model)
  invisible(gc())
  
  step3_time <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  # cat(sprintf("Phase 3 completed in %.1f seconds\n", step3_time))
  
  # ================== Phase 4: Build Combined Models ==================
  step_start <- Sys.time()
  
  site_combo_model <- cbind(sumofx_site / countofx_site, sumofy_site, sumofy_site / sum(sumofy_site))
  site_combo_model.complete <- site_combo_model[complete.cases(site_combo_model),]
  
  distance_combo_model <- cbind(sumofx_dist / countofx_dist, sumofy_dist, sumofy_dist / sum(sumofy_dist))
  distance_combo_model.complete <- distance_combo_model[complete.cases(distance_combo_model),]
  
  distance_combo_spline <- stats::smooth.spline(log10(distance_combo_model.complete[, 1]), distance_combo_model.complete[, 3], spar = 0.35)
  site_combo_spline <- stats::smooth.spline(log10(site_combo_model.complete[, 1]), site_combo_model.complete[, 3], spar = 0.35)
  
  step4_time <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  # cat(sprintf("Phase 4 completed in %.1f seconds\n", step4_time))
  
  # ================== Phase 5: Statistical Testing ==================
  step_start <- Sys.time()
  
  n_rows <- nrow(df)
  cat(sprintf("Processing %d interactions\n", n_rows))
  
  # Data Preprocessing
  min_positive_site <- 1e-10
  df[site <= 0, site := min_positive_site]
  
  # Spline Predictions
  fast_spline_predict <- function(spline_obj, query_values, grid_size = 3000) {
    x_range <- range(spline_obj$x)
    grid_x <- seq(x_range[1], x_range[2], length.out = grid_size)
    grid_y <- predict(spline_obj, grid_x)$y
    approx(grid_x, grid_y, query_values, rule = 2)$y
  }
  
  # Pre-compute log transformations
  log_site <- log10(df$site)
  log_loopWidth <- log10(df$loopWidth)
  
  # Range constraints
  distance_range <- range(distance_IAB_spline$x)
  site_range <- range(site_IAB_spline$x)
  
  log_site_clamped <- pmax(site_range[1], pmin(site_range[2], log_site))
  log_loopWidth_clamped <- pmax(distance_range[1], pmin(distance_range[2], log_loopWidth))
  
  # All spline predictions
  P_IAB_distance <- fast_spline_predict(distance_IAB_spline, log_loopWidth_clamped)
  P_combos_distance <- fast_spline_predict(distance_combo_spline, log_loopWidth_clamped)
  P_IAB_site <- fast_spline_predict(site_IAB_spline, log_site_clamped, 2500)
  P_combos_site <- fast_spline_predict(site_combo_spline, log_site_clamped, 2500)
  
  # Cleanup
  rm(log_site, log_loopWidth, log_site_clamped, log_loopWidth_clamped)
  
  # Boundary Constraints
  dist_IAB_range <- range(distance_IAB_model[, 3])
  dist_combo_range <- range(distance_combo_model[, 3])
  site_IAB_range <- range(site_IAB_model[, 3])
  site_combo_range <- range(site_combo_model[, 3])
  
  P_IAB_distance <- pmax(dist_IAB_range[1], pmin(dist_IAB_range[2], P_IAB_distance))
  P_combos_distance <- pmax(dist_combo_range[1], pmin(dist_combo_range[2], P_combos_distance))
  P_IAB_site <- pmax(site_IAB_range[1], pmin(site_IAB_range[2], P_IAB_site))
  P_combos_site <- pmax(site_combo_range[1], pmin(site_combo_range[2], P_combos_site))
  
  # Parameter Calculation
  totalcombos <- sum(sumofy_dist)
  totalIAB <- sum(df$PETS)
  
  p_binom_site <- (P_IAB_distance * P_IAB_site) / (P_combos_distance * P_combos_site * totalcombos)
  p_binom_site <- pmax(1e-15, pmin(1, p_binom_site))
  mu <- totalIAB * p_binom_site
  mu_safe <- pmax(mu, 1e-10)
  
  # Dispersion assessment
  observed_var <- var(df$PETS)
  observed_mean <- mean(df$PETS)
  dispersion_ratio <- observed_var / observed_mean
  
  # P-value Computation
  q25_threshold <- quantile(mu_safe, 0.25)
  pets_threshold <- pmax(q25_threshold, mu_safe + 3 * sqrt(mu_safe))
  high_priority_mask <- df$PETS > pets_threshold
  n_high_priority <- sum(high_priority_mask)
  
  P_site_result <- numeric(n_rows)
  
  if (dispersion_ratio > 1.5) {
    # Negative Binomial approach
    set.seed(123)
    sample_size <- min(800, max(200, ceiling(n_rows * 0.00002)))
    sample_indices <- sample(n_rows, sample_size)
    df_sample <- df[sample_indices]
    
    # Create a minimal data frame with matching mu_safe values
    sample_data <- data.frame(
      PETS = df_sample$PETS,
      mu_offset = log(mu_safe[sample_indices])
    )
    
    nb_fit <- MASS::glm.nb(PETS ~ offset(mu_offset), data = sample_data)
    theta_adj <- ifelse(is.null(nb_fit$theta) || !is.finite(nb_fit$theta), 25, max(20, nb_fit$theta))
    
    if (n_high_priority > 0) {
      hp_indices <- which(high_priority_mask)
      P_site_result[hp_indices] <- pnbinom(df$PETS[hp_indices] - 1, 
                                           mu = mu_safe[hp_indices], 
                                           size = theta_adj,
                                           lower.tail = FALSE)
    }
    
    if (n_high_priority < n_rows) {
      bg_indices <- which(!high_priority_mask)
      P_site_result[bg_indices] <- ppois(df$PETS[bg_indices] - 1, 
                                         lambda = mu_safe[bg_indices], 
                                         lower.tail = FALSE)
    }
    
  } else {
    # Binomial approach
    if (n_high_priority > 0) {
      hp_indices <- which(high_priority_mask)
      P_site_result[hp_indices] <- pbinom(df$PETS[hp_indices] - 1, 
                                          size = totalIAB, 
                                          prob = p_binom_site[hp_indices],
                                          lower.tail = FALSE)
    }
    
    if (n_high_priority < n_rows) {
      bg_indices <- which(!high_priority_mask)
      mu_vals <- totalIAB * p_binom_site[bg_indices]
      var_vals <- mu_vals * (1 - p_binom_site[bg_indices])
      z_scores <- (df$PETS[bg_indices] - mu_vals) / sqrt(pmax(var_vals, 1))
      P_site_result[bg_indices] <- pnorm(z_scores, lower.tail = FALSE)
    }
  }
  
  # Data Assignment and Correction
  df[, `:=`(
    P_site = P_site_result,
    p_binom_site = p_binom_site,
    mu_safe = mu_safe
  )]
  
  rm(P_site_result, p_binom_site, mu_safe)
  
  # Fallback check
  if (n_rows > 10000) {
    p_sample_indices <- sample(n_rows, 500)
    p_unique_ratio <- length(unique(df$P_site[p_sample_indices])) / 500
    
    if (p_unique_ratio < 0.2) {
      df[, lambda := totalIAB * p_binom_site]
      poisson_p_values <- 1 - ppois(df$PETS - 1, lambda = df$lambda)
      
      zero_mask <- poisson_p_values == 0
      if (any(zero_mask)) {
        poisson_p_values[zero_mask] <- exp(ppois(df$PETS[zero_mask] - 1,
                                                 lambda = df$lambda[zero_mask],
                                                 log.p = TRUE, lower.tail = FALSE))
      }
      
      df[, P_site := poisson_p_values]
      rm(poisson_p_values)
    }
  }
  
  # BH correction
  df[, Q_site := p.adjust(P_site, method = "BH")]
  
  # Add remaining columns
  df[, `:=`(
    P_IAB_distance = P_IAB_distance,
    P_combos_distance = P_combos_distance,
    P_IAB_site = P_IAB_site,
    P_combos_site = P_combos_site,
    mu = mu
  )]
  
  # Cleanup
  rm(P_IAB_distance, P_combos_distance, P_IAB_site, P_combos_site, mu)
  invisible(gc())
  
  step5_time <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  # cat(sprintf("Phase 5 completed in %.2f seconds\n", step5_time))
  
  # ================== Phase 6: Output Generation ==================
  step_start <- Sys.time()
  
  # Generate final results
  df_final <- df[complete.cases(df)]
  data.table::setorder(df_final, Q_site)
  df_significant <- df_final[Q_site < FDR, .(chr_1, start_1, end_1, chr_2, start_2, end_2, PETS, Q_site)]
  
  # Create output directory and write results
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  output_file <- file.path(outdir, paste0(sname, ".interaction.", mode, ".Q", FDR, ".txt"))
  data.table::fwrite(df_significant, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  step6_time <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  # cat(sprintf("Phase 6 completed in %.1f seconds\n", step6_time))
  
  # ================== Performance Summary ==================
  total_time <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat("=== PROCESSING COMPLETE ===\n")
  cat(sprintf("Total wall time: %.1f seconds (%.1f minutes)\n", total_time, total_time / 60))
  cat(sprintf("Phase 1 (load/process): %.1f s\n", step1_time))
  cat(sprintf("Phase 2 (IAB models): %.1f s\n", step2_time))
  cat(sprintf("Phase 3 (combination processing): %.1f s\n", step3_time))
  cat(sprintf("Phase 4 (combo model building): %.1f s\n", step4_time))
  cat(sprintf("Phase 5 (statistical testing): %.1f s\n", step5_time))
  cat(sprintf("Phase 6 (output generation): %.1f s\n", step6_time))
 
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat(sprintf("Final processed interactions: %d\n", nrow(df_final)))
  cat(sprintf("Processing completed. Found %d significant interactions (FDR < %g)\n", nrow(df_significant), FDR))
  cat(sprintf("Output written to: %s\n", file.path(outdir, paste0(sname, ".interaction.", mode, ".Q", FDR, ".txt"))))
  
  invisible(NULL)
}

