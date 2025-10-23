---
title: "sintHiChIP: detecting significant HiChIP interactions with cut site density correction"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{sintHiChIP: detecting significant HiChIP interactions with cut site density correction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

**Authors:** Weiyue Ding (wyding0501@hotmail.com)  
**Affiliation:** Harbin Institute of Technology

sintHiChIP is a comprehensive R package designed for the processing of HiChIP data, offering both local and global modes of operation. By employing cut site density correction, sintHiChIP facilitates the identification and characterization of significant chromatin interactions from HiChIP experiments.

## Platform

sintHiChIP is primarily designed to run on Linux and Unix-like operating systems. This includes various Linux distributions (such as Ubuntu, CentOS, Fedora) as well as Unix-based systems like macOS. While it may work on other platforms, we strongly recommend using a Linux or Unix environment for performance and compatibility.

## Prerequisites

Before installing and using sintHiChIP, ensure you have the following software and R packages installed.

### Required Software

1. R (version 4.0.0 or higher)
2. bedtools (version 2.29.0 or higher)
   - Used for genomic interval manipulations
   - Installation: `sudo apt-get install bedtools` (Ubuntu/Debian) 

3. tabix (usually comes with samtools, version 1.10 or higher)
   - Used for indexing and querying TAB-delimited genome position files
   - Installation: `sudo apt-get install tabix` (Ubuntu/Debian) 

4. bgzip (usually comes with samtools, version 1.10 or higher)
   - Used for blocking compression of genomic data files
   - Installation: Comes with tabix

5. build_matrix utility
   - Required for matrix generation in global mode
   - Can be obtained from HiC-Pro utilities (version 3.1.0 or higher)
   - Alternative implementations may also be compatible
   - Must be executable and accessible in system PATH or specified directly

### Input Data Requirements

sintHiChIP requires pre-processed HiChIP data from HiC-Pro pipeline. The package expects data in HiC-Pro output format:

1. **Common Requirements (Both Modes):**
   - HiC-Pro output directory containing allValidPairs files
   - Peak files in BED format (from MACS2 or similar peak callers)
   - Normalization restriction enzyme cut site density file (generated using `generate_normSite_file()` function, see [Generating Normalization Files](#generating-normalization-files) section)

2. **Local Mode Specific:**
   - Restriction fragment information in BED format
   - Peak-to-peak interaction within defined genomic regions

3. **Global Mode Specific:**
   - Chromosome size file
   - build_matrix utility for matrix generation
   - Genome-wide matrix-based interaction

### Required R Packages

You can install Bioconductor and the required R packages with:

```r
# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install Bioconductor packages if not already installed
bioc_pkgs <- c("GenomicRanges", "IRanges", "S4Vectors", "GenomeInfoDb")
BiocManager::install(bioc_pkgs[!bioc_pkgs %in% installed.packages()[,"Package"]])

# Install CRAN packages if not already installed
cran_pkgs <- c("data.table", "dplyr", "readr", "Rcpp", "matrixStats", "MASS", 
               "methods", "tools", "graphics", "stats", "utils", "parallel")
install.packages(cran_pkgs[!cran_pkgs %in% installed.packages()[,"Package"]])
```

## Installation

You can install sintHiChIP directly from GitHub using the devtools package:
```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("wding0501/sintHiChIP")
```

## Generating Normalization Files

Before running sintHiChIP analysis, you need to generate a normalization file that accounts for restriction enzyme cut site density across the genome.

### generate_normSite_file()

This function creates the normalization file required for cut site density correction:

```r
library(sintHiChIP)

# Generate normalization file for mouse genome
normsite_file <- generate_normSite_file(
  bed_file = "mm10_mboi.bed",
  species = "mouse",
  variance = 100000,
  binsize = 5000,
  output_dir = "./normalization",
  use_parallel = TRUE,
  ncores = NULL
)
```

#### Parameters:
- **bed_file**: Path to restriction enzyme cut site BED file (required)
- **species**: Species name - "mouse", "human", or "rat" (required)
- **variance**: Variance parameter for Gaussian smoothing (default: 100000)
- **binsize**: Genomic bin size in base pairs (default: 5000)
- **output_dir**: Output directory path (default: current directory)
- **use_parallel**: Enable parallel processing (default: TRUE)
- **ncores**: Number of cores for parallel processing (default: NULL for auto-detection)

#### Examples for different species:

**For Mouse (mm10):**
```r
setwd("/path/to/data")
BED_FILE <- "mm10_mboi.bed"
GENOME_BUILD <- "mm10"
BIN_SIZE <- 5000
USE_PARALLEL <- TRUE
NCORES <- 8
```

**For Human (hg38):**
```r
setwd("/path/to/data")
BED_FILE <- "hg38_mboi.bed"
GENOME_BUILD <- "hg38"
BIN_SIZE <- 5000
USE_PARALLEL <- TRUE
NCORES <- 8
```

**Run:**
```bash
Rscript normSite.R
```

**Output:**
- Mouse: `normsite_mm10_mboi_5000.tmp`
- Human: `normsite_hg38_mboi_5000.tmp`

#### Output:
The function generates a file named `normsite_{filename}_{species}_{binsize}_turbo.tmp` containing:
- Column 1: Chromosome name
- Column 2: Bin start position
- Column 3: Bin end position
- Column 4: Mean restriction site density

**Important:** Generate this file once per genome/enzyme combination. The same normalization file can be reused for all samples with the same genome assembly and restriction enzyme.

## Main Functions

sintHiChIP provides five main functions for HiChIP data processing:

### 1. sintHiChIP_sigloops - Statistical Significance Testing

**`sintHiChIP_sigloops`** is the core statistical engine of the package, responsible for identifying statistically significant chromatin interactions from HiChIP data. This function implements sophisticated statistical modeling with cut site density correction to distinguish genuine interactions from background noise.

#### Key Features:
- **Adaptive Statistical Modeling**: Automatically selects appropriate statistical models (Binomial, Negative Binomial, or Poisson) based on data overdispersion characteristics
- **Cut Site Density Correction**: Accounts for restriction enzyme cut site density to reduce bias in interaction detection
- **Memory-Efficient Processing**: Optimized for large datasets (40M+ interactions) with batch processing capabilities
- **Numerical Stability**: Uses log-space calculations and proper handling of edge cases
- **Comprehensive Error Handling**: Robust processing with detailed progress reporting

#### Statistical Methods:
The function employs a dual-modeling approach:
1. **IAB (Interaction Anchor Bias) Models**: Models actual observed interactions based on distance and site density
2. **Combination Models**: Models all possible interactions to establish background expectations
3. **P-value Calculation**: Uses appropriate statistical tests with FDR correction for multiple testing

#### Usage:
```r
sintHiChIP_sigloops(
  sname = "sample1",
  cwd = "/path/to/data",
  outdir = "/path/to/output", 
  normSiteFile = "/path/to/norm_sites.bed",
  local = TRUE,  # TRUE for local mode, FALSE for global mode
  FDR = 0.01,
  nbins = 10
)
```

#### Parameters:
- **sname**: Sample name for identification
- **cwd**: Working directory containing processed interaction data
- **outdir**: Output directory for results
- **normSiteFile**: BED format file with restriction enzyme cut site density
- **local**: Processing mode (TRUE for peak-to-peak, FALSE for peak-to-all)
- **FDR**: False Discovery Rate threshold for significance testing
- **nbins**: Number of bins for distance and site density modeling

#### Outputs:
- Significant interactions file: `{sample}.interaction.{mode}.Q{FDR}.txt`
- Contains genomic coordinates and statistical measures for significant loops

### 2. run_sintHiChIP

The primary interface for sintHiChIP processing, supporting both local and global modes:

```r
library(sintHiChIP)

# Step 1: Generate normalization file (do this once)
normsite_file <- generate_normSite_file(
  bed_file = "mm10_mboi.bed",
  species = "mouse",
  binsize = 5000,
  output_dir = "./normalization"
)

# Step 2: Local mode processing
run_sintHiChIP(
  mode = "local",
  outdir = "/path/to/output",
  hicpro_output = "/path/to/hicpro_results",
  sample_name = "sample1",
  peaks = "/path/to/peaks.bed",
  resfrags = "/path/to/restriction_fragments.bed",
  normSiteFile = normsite_file,  # Use generated file
  FDR = 0.01,
  min_dist = 20000,
  max_dist = 2000000
)

# Step 3: Global mode processing
run_sintHiChIP(
  mode = "global",
  outdir = "/path/to/output",
  hicpro_output = "/path/to/hicpro_results",
  sample_name = "sample1",
  peaks = "/path/to/peaks.bed",
  chr_size = "/path/to/chrom_sizes.txt",
  build_matrix = "/path/to/build_matrix",
  normSiteFile = normsite_file,  # Use the same generated file
  bin_size = 5000,
  FDR = 0.01,
  min_dist = 20000,
  max_dist = 2000000
)
```

### 3. sintHiChIP_local_mode

Specialized function for local mode processing with peak-to-peak interactions:

```r
sintHiChIP_local(
  outdir = "/path/to/local_output",
  hicpro_output = "/path/to/hicpro_results",
  sample_name = "sample1",
  peaks = "/path/to/peaks.bed",
  resfrags = "/path/to/restriction_fragments.bed",
  normSiteFile = "/path/to/norm_sites.bed",
  FDR = 0.01,
  min_dist = 20000,
  max_dist = 2000000,
  half_length = 73,
  no_merge = FALSE,
  max_anchor_width = 50000,
  keep_temp = FALSE,
  nbins = 10,
  peak_pad = 500,
  merge_gap = 500,
  make_washu = TRUE
)
```

### 4. sintHiChIP_global_mode

Specialized function for global mode processing with genome-wide matrix-based interactions:

```r
sintHiChIP_global(
  outdir = "/path/to/global_output",
  hicpro_output = "/path/to/hicpro_results",
  sample_name = "sample1",
  peaks = "/path/to/peaks.bed",
  chr_size = "/path/to/chrom_sizes.txt",
  build_matrix = "/path/to/build_matrix",
  normSiteFile = "/path/to/norm_sites.bed",
  bin_size = 5000,
  FDR = 0.01,
  min_dist = 20000,
  max_dist = 2000000,
  keep_temp = FALSE
)
```

### 5. sintHiChIP_local_single

Direct local mode processing using allValidPairs file:

```r
sintHiChIP_local_single(
  valid_pairs = "/path/to/sample1.allValidPairs",
  peaks = "/path/to/peaks.bed",
  resfrags = "/path/to/restriction_fragments.bed",
  normSiteFile = "/path/to/norm_sites.bed",
  outdir = "/path/to/local_output",
  FDR = 0.01,
  min_dist = 20000,
  max_dist = 2000000,
  half_length = 73,
  no_merge = FALSE,
  max_anchor_width = 50000,
  nbins = 10
)
```

### 6. sintHiChIP_global_single

Direct global mode processing using allValidPairs file:

```r
sintHiChIP_global_single(
  valid_pairs = "/path/to/sample1.allValidPairs",
  peaks = "/path/to/peaks.bed",
  chr_size = "/path/to/chrom_sizes.txt",
  build_matrix = "/path/to/build_matrix",
  normSiteFile = "/path/to/norm_sites.bed",
  outdir = "/path/to/global_output",
  bin_size = 5000,
  FDR = 0.01,
  min_dist = 20000,
  max_dist = 2000000,
  nbins = 10,
  keep_temp = FALSE
)
```

## Function Selection Guide

### When to use each function:

- **run_sintHiChIP()**: Main interface function, recommended for most users
- **sintHiChIP_local_mode()** or **sintHiChIP_global_mode()**: When you have standard HiC-Pro output structure and want mode-specific functionality
- **sintHiChIP_local_single()** or **sintHiChIP_global_single()**: When you want to specify exact allValidPairs file paths or have non-standard directory structures

### Local vs Global Mode:

- **Local Mode**: Peak-to-peak interactions within specific genomic regions, suitable for targeted interaction detection
- **Global Mode**: Genome-wide matrix-based processing, suitable for comprehensive interaction mapping

## Common Parameters

### Core Parameters:
- **outdir**: Output directory path
- **peaks**: Peak file in BED format (required for both modes)
- **normSiteFile**: Normalization file generated by `generate_normSite_file()` (required for both modes)
- **FDR**: False Discovery Rate threshold (default: 0.01)
- **min_dist**: Minimum interaction distance in bp (default: 20000)
- **max_dist**: Maximum interaction distance in bp (default: 2000000)

### Local Mode Specific:
- **resfrags**: Restriction fragments file in BED format
- **half_length**: Read extension length (default: 73)
- **no_merge**: Skip anchor merging (default: FALSE)
- **max_anchor_width**: Maximum anchor width (default: 50000)
- **peak_pad**: Peak padding in bp (default: 500)
- **merge_gap**: Merge gap for bedtools merge (default: 500)

### Global Mode Specific:
- **chr_size**: Chromosome sizes file
- **build_matrix**: Path to build_matrix utility
- **bin_size**: Genomic bin size for matrix generation (default: 5000)

### Optional Parameters:
- **nbins**: Number of bins for statistical modeling (default: 10)
- **keep_temp**: Keep temporary files (default: FALSE)
- **make_washu**: Create WashU/UCSC compatible files (default: TRUE, local mode only)

## Output Files

sintHiChIP generates several output files:

### 1. Filtered Interaction Files
**File**: `*.filt.intra.loop_counts.bedpe`

Contains filtered intra-chromosomal interactions with raw PET counts.

### 2. Significant Interaction Files
**File**: `*.interaction.[local|global].Q[FDR].txt`

Contains statistically significant interactions after FDR correction.

### 3. WashU Genome Browser Tracks
**File**: `*.interaction.[local|global].Q[FDR].washu.txt.gz`

Browser-compatible format for visualizing significant interactions.

## Statistical Methods

sintHiChIP employs statistical modeling with cut site density correction to identify significant chromatin interactions. The package automatically selects the most appropriate statistical model based on data characteristics.

## Documentation

For detailed function documentation:

```r
?run_sintHiChIP
?sintHiChIP_local
?sintHiChIP_global
?sintHiChIP_local_single
?sintHiChIP_global_single
?generate_normSite_file
```

## Complete Workflow Example

Here is a complete example demonstrating the entire sintHiChIP workflow from normalization file generation to result interpretation:

```r
library(sintHiChIP)

# ============================================================================
# Step 1: Generate Normalization File (one-time setup per genome/enzyme)
# ============================================================================

cat("Generating normalization file...\n")
normsite_file <- generate_normSite_file(
  bed_file = "data/mm10_mboi.bed",
  species = "mouse",
  variance = 100000,
  binsize = 5000,
  output_dir = "normalization",
  use_parallel = TRUE,
  ncores = 8
)

cat("Normalization file created:", normsite_file, "\n\n")

# ============================================================================
# Step 2: Run sintHiChIP Analysis
# ============================================================================

cat("Running sintHiChIP analysis for sample1...\n")

run_sintHiChIP(
  mode = "local",
  outdir = "results/sample1",
  hicpro_output = "hicpro_output",
  sample_name = "sample1",
  peaks = "data/sample1_peaks.bed",
  resfrags = "data/restriction_fragments.bed",
  normSiteFile = normsite_file,
  FDR = 0.01,
  min_dist = 20000,
  max_dist = 2000000,
  make_washu = TRUE
)

cat("Analysis completed!\n\n")

# ============================================================================
# Step 3: Load and Examine Results
# ============================================================================

# Load significant interactions
sig_file <- "results/sample1/sample1.interaction.local.Q0.01.txt"
if (file.exists(sig_file)) {
  sig_interactions <- read.table(sig_file, header = TRUE)
  cat("Found", nrow(sig_interactions), "significant interactions\n")
  cat("\nFirst few interactions:\n")
  print(head(sig_interactions))
}

# ============================================================================
# Step 4: Process Multiple Samples (reuse normalization file)
# ============================================================================

samples <- c("sample2", "sample3")

for (sample in samples) {
  cat("\nProcessing", sample, "...\n")
  
  run_sintHiChIP(
    mode = "local",
    outdir = file.path("results", sample),
    hicpro_output = "hicpro_output",
    sample_name = sample,
    peaks = file.path("data", paste0(sample, "_peaks.bed")),
    resfrags = "data/restriction_fragments.bed",
    normSiteFile = normsite_file,  # Reuse the same normalization file
    FDR = 0.01
  )
}

cat("\nAll samples processed successfully!\n")
```

## Workflow Summary

1. **Normalization File Generation**: Use `generate_normSite_file()` to create the cut site density normalization file (one-time setup per genome/enzyme combination)
2. **Data Preparation**: Ensure HiC-Pro processed data with allValidPairs files
3. **Function Selection**: Choose appropriate function based on your needs
4. **Parameter Configuration**: Set appropriate thresholds and file paths, including the generated normSiteFile
5. **Execute Processing**: Run sintHiChIP with selected function
6. **Results Interpretation**: Process significant interactions and visualization tracks

## Conclusion

sintHiChIP provides a comprehensive framework for HiChIP data processing through five main functions. The unified interface facilitates both targeted (local) and genome-wide (global) interaction detection while ensuring robust statistical significance testing with cut site density correction.
