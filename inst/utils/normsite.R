# Optimized R Implementation of normSite
# Using data.table, vectorization, and parallel processing

suppressMessages({
    library(data.table)
    library(GenomicRanges)
    library(IRanges)
    library(readr)
    library(parallel)
    library(foreach)
    library(doParallel)
})

options(scipen = 20)

#' Impute restriction sites using normal distribution
#' 
#' @param variance Variance parameter for normal distribution
#' @param chrseq GRanges object containing chromosome sites
#' @param PosMax Maximum position (chromosome length)
#' @return Numeric vector of imputed density values
imputeSite <- function(variance, chrseq, PosMax) {
    if (length(chrseq) == 0) {
        return(rep(0, PosMax))
    }
    
    # Calculate parameters
    x <- start(chrseq)
    varia <- variance
    delta <- round(sqrt(varia))
    normchr <- rep(0, PosMax)
    threshold <- 5 * delta
    
    # Classify sites into three categories
    X1 <- x[x < threshold]
    X1 <- sort(X1)
    X2 <- x[x >= threshold & x <= PosMax - threshold]
    X2 <- sort(X2)
    X3 <- x[x > PosMax - threshold]
    X3 <- sort(X3)
    
    X1_len <- length(X1)
    X2_len <- length(X2)
    X3_len <- length(X3)
    
    # Precompute normal distribution for X2 sites
    y <- 100000
    start_z <- y - threshold
    end_z <- y + threshold
    Z <- dnorm(start_z:end_z, y, delta)
    
    # Process X1 sites (near chromosome start)
    for (i in 1:X1_len) {
        Pos <- X1[i]
        t <- 1:(Pos + threshold)
        norm <- dnorm(t, Pos, threshold)
        normchr[1:(Pos + threshold)] <- norm + normchr[1:(Pos + threshold)]
    }
    
    # Process X2 sites (middle region)
    for (i in 1:X2_len) {
        Pos <- X2[i]
        normchr[(Pos - threshold):(Pos + threshold)] <- Z + normchr[(Pos - threshold):(Pos + threshold)]
    }
    
    # Process X3 sites (near chromosome end)
    for (i in 1:X3_len) {
        if (X3_len == 0) {
            break
        } else {
            Pos <- X3[i]
            t <- (Pos - threshold):PosMax
            norm <- dnorm(t, Pos, delta)
            normchr[(Pos - threshold):PosMax] <- norm + normchr[(Pos - threshold):PosMax]
        }
    }
    
    return(normchr)
}

#' Calculate window averages
#' 
#' @param data Input numeric vector
#' @param step Window size
#' @return Numeric vector of window averages
WindowFunct <- function(data, step) {
    total <- length(data)
    spots <- seq(from = 1, to = total, by = step)
    tlen <- length(spots)
    result <- vector(length = tlen)
    
    for (i in 1:(tlen - 1)) {
        result[i] <- mean(data[spots[i]:(spots[i] + step - 1)])
    }
    result[tlen] <- mean(data[spots[tlen]:total])
    
    return(result)
}

#' Parallel normalization with intelligent load balancing
#' 
#' @param variance Variance parameter
#' @param resSites GRanges object with restriction sites
#' @param GenomeBin GRanges object with genomic bins
#' @param SeqName Character vector of chromosome names
#' @param binsize Size of genomic bins
#' @param use_parallel Whether to use parallel processing
#' @param ncores Number of cores to use
#' @param chunk_size Number of chromosomes per parallel chunk
#' @return Numeric vector of normalized site densities
normSite <- function(variance = 100000, resSites, GenomeBin, SeqName, binsize, 
                          use_parallel = TRUE, ncores = NULL, chunk_size = NULL) {
    
    # Auto-detect optimal number of cores
    if (is.null(ncores)) {
        total_cores <- parallel::detectCores()
        ncores <- min(4, total_cores - 1)
    }
    
    Chrlen <- length(SeqName)
    cat("Processing", Chrlen, "chromosomes with bin size", binsize, "\n")
    
    # Precompute total result length to avoid dynamic growth
    total_bins <- 0
    chr_bin_counts <- integer(Chrlen)
    chr_site_counts <- integer(Chrlen)
    
    # Calculate workload for intelligent load balancing
    for (i in seq_len(Chrlen)) {
        chr_length <- seqlengths(GenomeBin)[i]
        chr_bins <- ceiling(chr_length / binsize)
        chr_bin_counts[i] <- chr_bins
        total_bins <- total_bins + chr_bins
        
        # Count sites per chromosome for load balancing
        chr_sites <- sum(seqnames(resSites) == SeqName[i])
        chr_site_counts[i] <- chr_sites
    }
    
    cat("Total bins to generate:", total_bins, "\n")
    cat("Total restriction sites:", length(resSites), "\n")
    
    # Pre-allocate final result vector
    GrangeSite <- numeric(total_bins)
    current_offset <- 1
    
    # Setup parallel processing with intelligent strategy
    if (use_parallel && Chrlen > 1 && ncores > 1) {
        cat("Using parallel processing with", ncores, "cores\n")
        
        # Intelligent chunking strategy
        if (is.null(chunk_size)) {
            total_cores <- parallel::detectCores()
            if (total_cores >= 32 && Chrlen <= total_cores) {
                # High-core system: one chromosome per core
                chunk_size <- 1
                ncores <- min(Chrlen, total_cores - 2)
                cat("High-core system detected: using one chromosome per core\n")
            } else {
                # Regular system: balance workload
                avg_sites_per_chr <- mean(chr_site_counts[chr_site_counts > 0])
                if (avg_sites_per_chr > 10000) {
                    chunk_size <- 1
                } else {
                    chunk_size <- max(1, ceiling(Chrlen / (ncores * 2)))
                }
            }
        }
        
        cat("Using chunk size:", chunk_size, "with", ncores, "cores\n")
        
        # Prepare data for parallel processing
        chr_data_list <- list()
        for (i in seq_len(Chrlen)) {
            chr_name <- SeqName[i]
            chrSite <- resSites[seqnames(resSites) == chr_name]
            seqlevels(chrSite, pruning.mode = "coarse") <- chr_name
            
            chr_data_list[[i]] <- list(
                chr_index = i,
                chr_name = chr_name,
                chrSite = chrSite,
                PosMax = seqlengths(GenomeBin)[i],
                expected_bins = chr_bin_counts[i]
            )
        }
        
        # Setup parallel cluster
        cl <- parallel::makeCluster(ncores)
        
        # Export necessary objects to cluster
        parallel::clusterEvalQ(cl, {
            library(data.table)
            library(GenomicRanges)
        })
        
        parallel::clusterExport(cl, c("imputeSite", "WindowFunct", 
                                    "variance", "binsize"), 
                              envir = environment())
        
        # Parallel execution
        cat("Processing chromosomes in parallel...\n")
        
        chr_results <- parallel::parLapply(cl, chr_data_list, function(chr_info) {
            if (length(chr_info$chrSite) == 0) {
                return(list(
                    chr_index = chr_info$chr_index,
                    result = rep(0, chr_info$expected_bins)
                ))
            }
            
            # Execute imputeSite
            normchr <- imputeSite(variance, chr_info$chrSite, chr_info$PosMax)
            
            # Execute WindowFunct
            site <- WindowFunct(normchr, binsize)
            
            return(list(
                chr_index = chr_info$chr_index,
                result = site
            ))
        })
        
        # Stop parallel cluster
        parallel::stopCluster(cl)
        
        # Sort results by chromosome index to maintain order
        chr_results <- chr_results[order(sapply(chr_results, function(x) x$chr_index))]
        
        # Merge results into final vector
        for (i in seq_len(Chrlen)) {
            chr_result <- chr_results[[i]]$result
            
            # Validate result length
            expected_length <- chr_bin_counts[i]
            actual_length <- length(chr_result)
            
            if (actual_length != expected_length) {
                cat("Warning: Length mismatch for chromosome", SeqName[i], 
                    "- expected:", expected_length, "actual:", actual_length, "\n")
                # Pad or trim as needed
                if (actual_length < expected_length) {
                    chr_result <- c(chr_result, rep(0, expected_length - actual_length))
                } else {
                    chr_result <- chr_result[1:expected_length]
                }
            }
            
            result_length <- length(chr_result)
            end_offset <- current_offset + result_length - 1
            
            if (end_offset <= length(GrangeSite)) {
                GrangeSite[current_offset:end_offset] <- chr_result
                current_offset <- end_offset + 1
                
                # Progress reporting
                zero_pct <- sum(chr_result == 0) / length(chr_result) * 100
                max_val <- max(chr_result)
                cat("Completed chromosome", SeqName[i], "- bins:", result_length, 
                    ", zero%:", sprintf("%.1f", zero_pct), 
                    ", max:", sprintf("%.4f", max_val), "\n")
            } else {
                cat("Error: Result vector overflow for chromosome", SeqName[i], "\n")
            }
        }
        
    } else {
        # Serial processing with optimizations
        cat("Using serial processing with optimizations\n")
        
        for (i in seq_len(Chrlen)) {
            cat("Processing chromosome", SeqName[i], "(", i, "/", Chrlen, ")\n")
            
            # Extract chromosome-specific data
            chrSite <- resSites[seqnames(resSites) == SeqName[i]]
            seqlevels(chrSite, pruning.mode = "coarse") <- SeqName[i]
            
            PosMax <- seqlengths(GenomeBin)[i]
            
            if (length(chrSite) == 0) {
                site <- rep(0, chr_bin_counts[i])
            } else {
                # Execute imputeSite with timing
                start_time <- Sys.time()
                normchr <- imputeSite(variance, chrSite, PosMax)
                impute_time <- Sys.time() - start_time
                
                # Execute WindowFunct with timing
                start_time <- Sys.time()
                site <- WindowFunct(normchr, binsize)
                window_time <- Sys.time() - start_time
                
                cat("    Impute time:", format(impute_time), ", Window time:", format(window_time), "\n")
                
                # Memory cleanup
                rm(normchr)
                gc(verbose = FALSE)
            }
            
            # Copy results to final array
            result_length <- length(site)
            end_offset <- current_offset + result_length - 1
            GrangeSite[current_offset:end_offset] <- site
            current_offset <- end_offset + 1
            
            # Report statistics
            zero_count <- sum(site == 0)
            zero_pct <- zero_count / length(site) * 100
            max_val <- max(site)
            cat("    Bins:", length(site), ", Zero%:", sprintf("%.1f", zero_pct), 
                ", Max:", sprintf("%.6f", max_val), "\n")
        }
    }
    
    return(GrangeSite)
}

#' Intelligent BED file preprocessing with filtering and validation
#' 
#' @param bed_file Path to BED file
#' @param target_chroms Vector of target chromosome names
#' @param min_sites_per_chr Minimum sites required per chromosome
#' @return List containing processed data and statistics
preprocess_bed_file <- function(bed_file, target_chroms, min_sites_per_chr = 100) {
    cat("Preprocessing BED file:", bed_file, "\n")
    
    # Fast reading with data.table - handle variable number of columns
    dt <- fread(bed_file, header = FALSE)
    
    # Check number of columns and assign appropriate names
    ncols <- ncol(dt)
    cat("Detected", ncols, "columns in BED file\n")
    
    if (ncols < 3) {
        stop("BED file must have at least 3 columns (chr, start, end)")
    }
    
    # Use only first 3 columns and rename them
    dt <- dt[, 1:3]
    setnames(dt, c("chr", "start", "end"))
    
    cat("Original file:", nrow(dt), "sites across", length(unique(dt$chr)), "chromosomes\n")
    
    # Filter to target chromosomes
    dt_filtered <- dt[chr %in% target_chroms]
    
    # Calculate chromosome statistics
    chr_stats <- dt_filtered[, .(count = .N, max_pos = max(end)), by = chr]
    setorder(chr_stats, -count)
    
    cat("Chromosome statistics:\n")
    print(chr_stats)
    
    # Filter chromosomes with insufficient sites
    valid_chrs <- chr_stats[count >= min_sites_per_chr, chr]
    dt_final <- dt_filtered[chr %in% valid_chrs]
    
    cat("Final data:", nrow(dt_final), "sites across", length(valid_chrs), "chromosomes\n")
    
    # Additional data quality checks
    duplicates <- dt_final[, .N, by = .(chr, start)][N > 1]
    if (nrow(duplicates) > 0) {
        cat("Warning:", nrow(duplicates), "duplicate sites found\n")
        dt_final <- unique(dt_final, by = c("chr", "start"))
        cat("Removed duplicates, final count:", nrow(dt_final), "\n")
    }
    
    return(list(
        data = dt_final,
        valid_chromosomes = valid_chrs,
        chromosome_stats = chr_stats
    ))
}

#' Main turbo normSite function with comprehensive optimization
#' 
#' @param bed_file Input BED file path
#' @param species Species name ("mouse", "human", "rat")
#' @param variance Variance parameter for normal distribution
#' @param binsize Genomic bin size
#' @param use_parallel Enable parallel processing
#' @param ncores Number of cores (NULL for auto-detection)
#' @param output_dir Output directory
#' @return data.table with normalized results
main_normSite_turbo <- function(bed_file, species = "mouse", variance = 100000, 
                               binsize = 5000, use_parallel = TRUE, ncores = NULL, 
                               output_dir = NULL) {
    
    if (is.null(output_dir)) {
        output_dir <- getwd()
    }
    
    # Define chromosome sets for different species
    chr_sets <- list(
        mouse = paste0("chr", c(1:19, "X", "Y", "M")),
        human = paste0("chr", c(1:22, "X", "Y", "M")),
        rat = paste0("chr", c(1:20, "X", "Y", "M"))
    )
    
    target_chroms <- chr_sets[[species]]
    if (is.null(target_chroms)) {
        stop("Unsupported species: ", species, ". Use 'mouse', 'human', or 'rat'")
    }
    
    cat("=== Starting Turbo normSite Analysis ===\n")
    cat("Species:", species, "\n")
    cat("Input file:", bed_file, "\n")
    cat("Variance:", variance, "\n")
    cat("Bin size:", binsize, "\n")
    cat("Parallel processing:", use_parallel, "\n")
    if (use_parallel && !is.null(ncores)) {
        cat("Number of cores:", ncores, "\n")
    }
    cat("\n")
    
    # Data preprocessing with validation
    start_time <- Sys.time()
    processed_data <- preprocess_bed_file(bed_file, target_chroms)
    preprocess_time <- Sys.time() - start_time
    cat("Preprocessing completed in:", format(preprocess_time), "\n\n")
    
    # Convert to GenomicRanges objects
    dt <- processed_data$data
    resSites <- makeGRangesFromDataFrame(
        data.frame(seqnames = dt$chr, start = dt$start, end = dt$start),
        keep.extra.columns = FALSE
    )
    
    # Setup chromosome information
    valid_chroms <- processed_data$valid_chromosomes
    seqlevels(resSites, pruning.mode = 'coarse') <- valid_chroms
    
    # Calculate chromosome lengths
    chr_lengths <- dt[, max(end), by = chr]
    setnames(chr_lengths, c("chr", "length"))
    chr_lengths <- chr_lengths[match(valid_chroms, chr_lengths$chr)]
    
    # Create GenomeBin object
    seqlengths(resSites) <- chr_lengths$length
    GenomeBin <- tileGenome(seqinfo(resSites), tilewidth = binsize, cut.last.tile.in.chrom = TRUE)
    start(GenomeBin) <- start(GenomeBin) - 1
    
    cat("=== Starting normalization ===\n")
    norm_start_time <- Sys.time()
    
    # Execute super-optimized normalization
    SiteCoverage <- normSite(variance, resSites, GenomeBin, valid_chroms, 
                                  binsize, use_parallel = use_parallel, ncores = ncores)
    
    norm_end_time <- Sys.time()
    norm_time <- norm_end_time - norm_start_time
    cat("Normalization completed in:", format(norm_time), "\n\n")
    
    # Generate final results
    mcols(GenomeBin)$MeanSite <- SiteCoverage
    
    NormSiteFile <- data.table(
        seqnames = as.character(seqnames(GenomeBin)),
        start = start(GenomeBin),
        end = end(GenomeBin),
        MeanSite = mcols(GenomeBin)$MeanSite
    )
    
    # Set bins starting at position 0 to have density 0 (as in original code)
    NormSiteFile[start == 0, MeanSite := 0]
    
    # Generate output filename
    base_name <- tools::file_path_sans_ext(basename(bed_file))
    output_file <- file.path(output_dir, sprintf("normsite_%s_%s_%d_turbo.tmp", 
                                                base_name, species, binsize))
    
    # Write results using fast data.table method
    fwrite(NormSiteFile, output_file, sep = "\t", col.names = TRUE)
    
    # Generate comprehensive statistics report
    cat("=== Results Summary ===\n")
    cat("Total processing time:", format(norm_end_time - start_time), "\n")
    cat("Total genomic bins:", nrow(NormSiteFile), "\n")
    cat("Chromosomes processed:", paste(valid_chroms, collapse = ", "), "\n")
    
    zero_count <- sum(NormSiteFile$MeanSite == 0)
    zero_pct <- zero_count / nrow(NormSiteFile) * 100
    cat("Zero density bins:", zero_count, sprintf("(%.1f%%)", zero_pct), "\n")
    
    non_zero_values <- NormSiteFile[MeanSite > 0, MeanSite]
    if (length(non_zero_values) > 0) {
        cat("Non-zero statistics:\n")
        cat("  Min:", sprintf("%.6f", min(non_zero_values)), "\n")
        cat("  Mean:", sprintf("%.6f", mean(non_zero_values)), "\n")
        cat("  Median:", sprintf("%.6f", median(non_zero_values)), "\n")
        cat("  Max:", sprintf("%.6f", max(non_zero_values)), "\n")
        cat("  Standard deviation:", sprintf("%.6f", sd(non_zero_values)), "\n")
    }
    
    # Performance metrics
    total_sites <- nrow(dt)
    sites_per_second <- total_sites / as.numeric(norm_time, units = "secs")
    bins_per_second <- nrow(NormSiteFile) / as.numeric(norm_time, units = "secs")
    
    cat("Performance metrics:\n")
    cat("  Sites processed per second:", sprintf("%.0f", sites_per_second), "\n")
    cat("  Bins generated per second:", sprintf("%.0f", bins_per_second), "\n")
    
    cat("Results written to:", output_file, "\n")
    
    return(NormSiteFile)
}

# ============================================================================
# MAIN EXECUTION SCRIPT
# ============================================================================

# Configuration parameters
setwd("/home/wding/HiChIP/diffloop/normSite")  # Set working directory
BED_FILE <- "mm10_mboi.bed"  # Change to your file
GENOME_BUILD <- "mm10"       # Options: "mm09", "mm10", "hg19", "hg38"
VARIANCE <- 100000
BIN_SIZE <- 5000
USE_PARALLEL <- TRUE         # Set to FALSE for serial processing
NCORES <- NULL              # NULL for auto-detection, or specify number

# Execute analysis using original function structure
cat("Starting Optimized normSite Analysis...\n")

# Read and preprocess BED file
cat("Preprocessing BED file:", BED_FILE, "\n")

# Define chromosome sets for different genome builds
chr_sets <- list(
    mm09 = paste0("chr", c(1:19, "X", "Y", "M")),
    mm10 = paste0("chr", c(1:19, "X", "Y", "M")),
    hg19 = paste0("chr", c(1:22, "X", "Y", "M")),
    hg38 = paste0("chr", c(1:22, "X", "Y", "M"))
)

target_chroms <- chr_sets[[GENOME_BUILD]]
if (is.null(target_chroms)) {
    stop("Unsupported genome build: ", GENOME_BUILD, ". Use 'mm09', 'mm10', 'hg19', or 'hg38'")
}

# Preprocess BED file
processed_data <- preprocess_bed_file(BED_FILE, target_chroms)
dt <- processed_data$data

# Create resSites (GRanges object)
resSites <- makeGRangesFromDataFrame(
    data.frame(seqnames = dt$chr, start = dt$start, end = dt$start),
    keep.extra.columns = FALSE
)

# Setup chromosome information
valid_chroms <- processed_data$valid_chromosomes
seqlevels(resSites, pruning.mode = 'coarse') <- valid_chroms

# Calculate chromosome lengths
chr_lengths <- dt[, max(end), by = chr]
setnames(chr_lengths, c("chr", "length"))
chr_lengths <- chr_lengths[match(valid_chroms, chr_lengths$chr)]

# Create GenomeBin object
seqlengths(resSites) <- chr_lengths$length
GenomeBin <- tileGenome(seqinfo(resSites), tilewidth = BIN_SIZE, cut.last.tile.in.chrom = TRUE)
start(GenomeBin) <- start(GenomeBin) - 1

# Execute normalization using original function pattern
cat("=== Starting normalization ===\n")
start_time <- Sys.time()

SiteCoverage <- normSite(VARIANCE, resSites, GenomeBin, valid_chroms, BIN_SIZE, 
                                  use_parallel = USE_PARALLEL, ncores = NCORES)

end_time <- Sys.time()
processing_time <- end_time - start_time

# Generate final results
mcols(GenomeBin)$MeanSite <- SiteCoverage

NormSiteFile <- data.frame(
    seqnames = seqnames(GenomeBin), 
    start = start(GenomeBin), 
    end = end(GenomeBin),
    MeanSite = mcols(GenomeBin)[,1]
)

# Set bins starting at position 0 to have density 0 (as in original code)
NormSiteFile$MeanSite[NormSiteFile$start == 0] <- 0

# Generate output filename
output_file <- sprintf("normsite_%s_mboi_%d.tmp", GENOME_BUILD, BIN_SIZE)

# Write results
write.table(NormSiteFile, output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Generate statistics report
cat("=== Results Summary ===\n")
cat("Total processing time:", format(processing_time), "\n")
cat("Total genomic bins:", nrow(NormSiteFile), "\n")
cat("Chromosomes processed:", paste(valid_chroms, collapse = ", "), "\n")

zero_count <- sum(NormSiteFile$MeanSite == 0)
zero_pct <- zero_count / nrow(NormSiteFile) * 100
cat("Zero density bins:", zero_count, sprintf("(%.1f%%)", zero_pct), "\n")

non_zero_sites <- NormSiteFile$MeanSite[NormSiteFile$MeanSite > 0]
if (length(non_zero_sites) > 0) {
    cat("Non-zero statistics:\n")
    cat("  Min:", sprintf("%.6f", min(non_zero_sites)), "\n")
    cat("  Mean:", sprintf("%.6f", mean(non_zero_sites)), "\n")
    cat("  Median:", sprintf("%.6f", median(non_zero_sites)), "\n")
    cat("  Max:", sprintf("%.6f", max(non_zero_sites)), "\n")
}

cat("=== Analysis completed successfully ===\n")