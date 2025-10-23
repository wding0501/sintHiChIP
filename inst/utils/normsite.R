# ============================================================================
# normSite.R - Complete Version with generate_normSite_file() Function
# Optimized R Implementation for HiChIP Normalization
# Using data.table, vectorization, and parallel processing
# ============================================================================

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
            
            if (end_offset > total_bins) {
                cat("Warning: Result overflow for chromosome", SeqName[i], "\n")
                end_offset <- total_bins
                result_length <- total_bins - current_offset + 1
                chr_result <- chr_result[1:result_length]
            }
            
            GrangeSite[current_offset:end_offset] <- chr_result
            current_offset <- end_offset + 1
        }
        
    } else {
        # Serial processing
        cat("Using serial processing\n")
        
        for (i in seq_len(Chrlen)) {
            chr_name <- SeqName[i]
            cat("Processing chromosome:", chr_name, "\n")
            
            chrSite <- resSites[seqnames(resSites) == chr_name]
            seqlevels(chrSite, pruning.mode = "coarse") <- chr_name
            PosMax <- seqlengths(GenomeBin)[i]
            
            if (length(chrSite) == 0) {
                site <- rep(0, chr_bin_counts[i])
            } else {
                normchr <- imputeSite(variance, chrSite, PosMax)
                site <- WindowFunct(normchr, binsize)
            }
            
            result_length <- length(site)
            end_offset <- current_offset + result_length - 1
            
            if (end_offset > total_bins) {
                end_offset <- total_bins
                result_length <- total_bins - current_offset + 1
                site <- site[1:result_length]
            }
            
            GrangeSite[current_offset:end_offset] <- site
            current_offset <- end_offset + 1
        }
    }
    
    return(GrangeSite)
}

#' Preprocess BED file with optimized data.table operations
#' 
#' @param bed_file Path to input BED file
#' @param target_chroms Character vector of target chromosome names
#' @return List containing processed data.table and valid chromosomes
preprocess_bed_file <- function(bed_file, target_chroms) {
    cat("Reading BED file:", bed_file, "\n")
    
    # Use data.table's fread for fast reading
    dt <- fread(bed_file, header = FALSE, col.names = c("chr", "start", "end"))
    
    initial_rows <- nrow(dt)
    cat("Initial rows:", initial_rows, "\n")
    
    # Filter for target chromosomes (vectorized operation)
    dt <- dt[chr %in% target_chroms]
    
    filtered_rows <- nrow(dt)
    cat("Rows after chromosome filtering:", filtered_rows, "\n")
    
    if (filtered_rows == 0) {
        stop("No data remaining after filtering for target chromosomes")
    }
    
    # Get valid chromosomes in order
    valid_chroms <- unique(dt$chr)
    valid_chroms <- valid_chroms[order(match(valid_chroms, target_chroms))]
    
    cat("Valid chromosomes found:", length(valid_chroms), "\n")
    cat("Chromosomes:", paste(valid_chroms, collapse = ", "), "\n")
    
    # Convert chr to factor for efficient sorting
    dt[, chr := factor(chr, levels = valid_chroms)]
    
    # Sort by chromosome and position
    setorder(dt, chr, start)
    
    return(list(
        data = dt,
        valid_chromosomes = valid_chroms
    ))
}

# ============================================================================
# MAIN EXPORTED FUNCTION: generate_normSite_file()
# ============================================================================

#' Generate Normalization Site File for sintHiChIP
#' 
#' This function generates a normalization file that accounts for restriction
#' enzyme cut site density across the genome. The normalization file is used
#' by sintHiChIP to correct for biases in chromatin interaction detection.
#' 
#' @param bed_file Path to restriction enzyme cut site BED file (required).
#'   This should be a tab-delimited file with chromosome, start, and end positions.
#' @param species Species name - "mouse", "human", or "rat" (required).
#'   This determines which chromosomes to include in the analysis.
#' @param genome_build Genome build version. If NULL, auto-detected based on species.
#'   Options: "mm09", "mm10" (mouse); "hg19", "hg38" (human); "rn6" (rat).
#' @param variance Variance parameter for Gaussian smoothing (default: 100000).
#'   Higher values result in smoother density profiles.
#' @param binsize Genomic bin size in base pairs (default: 5000).
#'   Should match the resolution used in HiChIP analysis.
#' @param output_dir Output directory path (default: current directory).
#'   The normalization file will be saved here.
#' @param use_parallel Enable parallel processing (default: TRUE).
#'   Recommended for faster processing on multi-core systems.
#' @param ncores Number of cores for parallel processing (default: NULL for auto-detection).
#'   If NULL, uses detectCores() - 1, maximum of 4.
#' 
#' @return Character string containing the path to the generated normalization file.
#' 
#' @details
#' The function performs the following steps:
#' 1. Reads restriction enzyme cut sites from the BED file
#' 2. Filters for chromosomes appropriate to the specified species
#' 3. Applies Gaussian smoothing to estimate site density
#' 4. Bins the genome and calculates mean site density per bin
#' 5. Outputs a tab-delimited file with columns: chr, start, end, MeanSite
#' 
#' The output file naming convention is:
#' normsite_{bed_filename}_{species}_{binsize}_turbo.tmp
#' 
#' @examples
#' \dontrun{
#' # Generate normalization file for mouse genome
#' normsite_file <- generate_normSite_file(
#'   bed_file = "mm10_mboi.bed",
#'   species = "mouse",
#'   binsize = 5000,
#'   output_dir = "./normalization"
#' )
#' 
#' # Generate for human genome with specific parameters
#' normsite_file <- generate_normSite_file(
#'   bed_file = "hg38_mboi.bed",
#'   species = "human",
#'   genome_build = "hg38",
#'   variance = 100000,
#'   binsize = 5000,
#'   output_dir = "./norm_files",
#'   use_parallel = TRUE,
#'   ncores = 8
#' )
#' 
#' # Use the generated file in sintHiChIP
#' run_sintHiChIP(
#'   mode = "local",
#'   normSiteFile = normsite_file,
#'   # ... other parameters
#' )
#' }
#' 
#' @export
generate_normSite_file <- function(bed_file, 
                                   species,
                                   genome_build = NULL,
                                   variance = 100000,
                                   binsize = 5000,
                                   output_dir = ".",
                                   use_parallel = TRUE,
                                   ncores = NULL) {
    
    # Validate required parameters
    if (missing(bed_file) || is.null(bed_file)) {
        stop("bed_file is required")
    }
    if (missing(species) || is.null(species)) {
        stop("species is required. Options: 'mouse', 'human', 'rat'")
    }
    if (!file.exists(bed_file)) {
        stop("BED file does not exist: ", bed_file)
    }
    
    # Validate species
    valid_species <- c("mouse", "human", "rat")
    if (!tolower(species) %in% valid_species) {
        stop("Invalid species. Options: ", paste(valid_species, collapse = ", "))
    }
    species <- tolower(species)
    
    # Auto-detect genome build if not specified
    if (is.null(genome_build)) {
        genome_build <- switch(species,
            "mouse" = "mm10",
            "human" = "hg38",
            "rat" = "rn6",
            stop("Could not auto-detect genome build for species: ", species)
        )
        cat("Auto-detected genome build:", genome_build, "\n")
    }
    
    # Define chromosome sets for different genome builds
    chr_sets <- list(
        # Mouse genomes
        mm09 = paste0("chr", c(1:19, "X", "Y", "M")),
        mm10 = paste0("chr", c(1:19, "X", "Y", "M")),
        # Human genomes
        hg19 = paste0("chr", c(1:22, "X", "Y", "M")),
        hg38 = paste0("chr", c(1:22, "X", "Y", "M")),
        # Rat genome
        rn6 = paste0("chr", c(1:20, "X", "Y", "M"))
    )
    
    target_chroms <- chr_sets[[genome_build]]
    if (is.null(target_chroms)) {
        stop("Unsupported genome build: ", genome_build, 
             ". Options: mm09, mm10 (mouse); hg19, hg38 (human); rn6 (rat)")
    }
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("Created output directory:", output_dir, "\n")
    }
    
    # Print analysis parameters
    cat("\n=== Starting normSite Analysis ===\n")
    cat("BED file:", bed_file, "\n")
    cat("Species:", species, "\n")
    cat("Genome build:", genome_build, "\n")
    cat("Variance:", variance, "\n")
    cat("Bin size:", binsize, "bp\n")
    cat("Output directory:", output_dir, "\n")
    cat("Parallel processing:", use_parallel, "\n")
    if (!is.null(ncores)) {
        cat("Number of cores:", ncores, "\n")
    }
    cat("\n")
    
    # Start timing
    start_time <- Sys.time()
    
    # Preprocess BED file
    cat("=== Step 1: Preprocessing BED file ===\n")
    processed_data <- preprocess_bed_file(bed_file, target_chroms)
    dt <- processed_data$data
    valid_chroms <- processed_data$valid_chromosomes
    
    # Create resSites (GRanges object)
    cat("\n=== Step 2: Creating GRanges object ===\n")
    resSites <- makeGRangesFromDataFrame(
        data.frame(seqnames = dt$chr, start = dt$start, end = dt$start),
        keep.extra.columns = FALSE
    )
    
    # Setup chromosome information
    seqlevels(resSites, pruning.mode = 'coarse') <- valid_chroms
    
    # Calculate chromosome lengths
    chr_lengths <- dt[, max(end), by = chr]
    setnames(chr_lengths, c("chr", "length"))
    chr_lengths <- chr_lengths[match(valid_chroms, chr_lengths$chr)]
    
    # Create GenomeBin object
    cat("Creating genomic bins...\n")
    seqlengths(resSites) <- chr_lengths$length
    GenomeBin <- tileGenome(seqinfo(resSites), tilewidth = binsize, cut.last.tile.in.chrom = TRUE)
    start(GenomeBin) <- start(GenomeBin) - 1
    
    cat("Total genomic bins created:", length(GenomeBin), "\n")
    
    # Execute normalization
    cat("\n=== Step 3: Calculating site density normalization ===\n")
    norm_start_time <- Sys.time()
    
    SiteCoverage <- normSite(
        variance = variance,
        resSites = resSites,
        GenomeBin = GenomeBin,
        SeqName = valid_chroms,
        binsize = binsize,
        use_parallel = use_parallel,
        ncores = ncores
    )
    
    norm_end_time <- Sys.time()
    norm_duration <- norm_end_time - norm_start_time
    cat("Normalization completed in:", format(norm_duration), "\n")
    
    # Generate final results
    cat("\n=== Step 4: Generating output file ===\n")
    mcols(GenomeBin)$MeanSite <- SiteCoverage
    
    NormSiteFile <- data.frame(
        seqnames = seqnames(GenomeBin),
        start = start(GenomeBin),
        end = end(GenomeBin),
        MeanSite = mcols(GenomeBin)[, 1]
    )
    
    # Set bins starting at position 0 to have density 0
    NormSiteFile$MeanSite[NormSiteFile$start == 0] <- 0
    
    # Generate output filename
    base_name <- tools::file_path_sans_ext(basename(bed_file))
    output_file <- file.path(
        output_dir,
        sprintf("normsite_%s_%s_%d_turbo.tmp", base_name, species, binsize)
    )
    
    # Write results
    write.table(NormSiteFile, output_file, sep = "\t",
                row.names = FALSE, quote = FALSE, col.names = TRUE)
    
    # Calculate total processing time
    end_time <- Sys.time()
    total_duration <- end_time - start_time
    
    # Generate summary statistics
    cat("\n=== Results Summary ===\n")
    cat("Total processing time:", format(total_duration), "\n")
    cat("Total genomic bins:", nrow(NormSiteFile), "\n")
    cat("Chromosomes processed:", paste(valid_chroms, collapse = ", "), "\n")
    
    zero_count <- sum(NormSiteFile$MeanSite == 0)
    zero_pct <- zero_count / nrow(NormSiteFile) * 100
    cat("Zero density bins:", zero_count, sprintf("(%.1f%%)", zero_pct), "\n")
    
    non_zero_sites <- NormSiteFile$MeanSite[NormSiteFile$MeanSite > 0]
    if (length(non_zero_sites) > 0) {
        cat("\nNon-zero density statistics:\n")
        cat("  Min:    ", sprintf("%.6f", min(non_zero_sites)), "\n")
        cat("  Mean:   ", sprintf("%.6f", mean(non_zero_sites)), "\n")
        cat("  Median: ", sprintf("%.6f", median(non_zero_sites)), "\n")
        cat("  Max:    ", sprintf("%.6f", max(non_zero_sites)), "\n")
    }
    
    # Performance metrics
    total_sites <- length(resSites)
    total_bins <- nrow(NormSiteFile)
    duration_seconds <- as.numeric(total_duration, units = "secs")
    
    if (duration_seconds > 0) {
        sites_per_second <- total_sites / duration_seconds
        bins_per_second <- total_bins / duration_seconds
        
        cat("\nPerformance metrics:\n")
        cat("  Sites processed per second:", sprintf("%.0f", sites_per_second), "\n")
        cat("  Bins generated per second: ", sprintf("%.0f", bins_per_second), "\n")
    }
    
    cat("\n=== Output File ===\n")
    cat("File:", output_file, "\n")
    cat("Size:", sprintf("%.2f MB", file.size(output_file) / 1024^2), "\n")
    
    cat("\n=== Analysis completed successfully ===\n\n")
    
    # Return the output file path
    invisible(output_file)
}

# ============================================================================
# End of normSite.R
# ============================================================================
