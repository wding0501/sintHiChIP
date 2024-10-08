# sintHiChIP

sintHiChIP is a comprehensive R package designed for the analysis of HiChIP data, offering both local and global modes of operation. By correcting the bias of restriction enzyme cut site density, sintHiChIP facilitates the identification and characterization of significant chromatin interactions from HiChIP experiments.

## Installation

You can install sintHiChIP directly from GitHub using the devtools package:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("wding0501/sintHiChIP")
```
## Platform

sintHiChIP is primarily designed to run on Linux and Unix-like operating systems. This includes various Linux distributions (such as Ubuntu, CentOS, Fedora) as well as Unix-based systems like macOS. While it may work on other platforms, we strongly recommend using a Linux or Unix environment for optimal performance and compatibility.

## Prerequisites

Before installing and using sintHiChIP, ensure you have the following:

- R (version 4.0.0 or higher)
- Python (version 3.6 or higher)

### For Local Mode:

- hichipper (version 0.7.7 or higher)
- bedtools (version 2.29.0 or higher)
- tabix (version 1.10 or higher)
- bgzip (version 1.10 or higher)

### Required R Packages

```r
# Install Bioconductor packages
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors", "GenomeInfoDb"))

# Install CRAN packages
install.packages(c("data.table", "dplyr", "readr", "yaml", "Rcpp", "matrixStats", "testthat", "knitr", "rmarkdown"))
```

## Usage

sintHiChIP can be run in two modes: local and global.

### Local Mode

```r
library(sintHiChIP)

run_sintHiChIP(
  mode = "local",
  outdir = "P2P",
  yaml = "config.yaml",
  peaks = "EACH,ALL",
  resfrags = "path/to/restriction_fragments.bed",
  hicpro_output = "path/to/hicpro_results",
  normSiteFile = "path/to/norm_sites.bed",
  FDR = 0.01
)
```

### Global Mode

```r
run_sintHiChIP(
  mode = "global",
  outdir = "/path/to/output",
  valid_pairs = "/path/to/valid_pairs.txt",
  chr_size = "/path/to/chrom_sizes.txt",
  build_matrix = "/path/to/build_matrix",
  peaks = "/path/to/peaks.bed",
  normSiteFile = "/path/to/norm_sites.bed",
  bin_size = 5000,
  FDR = 0.01
)
```

## Output Files

sintHiChIP generates the following significant output files:

1. Distance filtered intra interaction files (*.filt.intra.loop_counts.bedpe)
   Columns:
   - chr1: Chromosome of the first anchor
   - start1: Start position of the first anchor
   - end1: End position of the first anchor
   - chr2: Chromosome of the second anchor
   - start2: Start position of the second anchor
   - end2: End position of the second anchor
   - dot:  the dot .
   - couunt: raw loop counts


2. Significant interaction files (*.interaction.Q0.01.txt)
   Columns:
   - chr1: Chromosome of the first anchor
   - start1: Start position of the first anchor
   - end1: End position of the first anchor
   - chr2: Chromosome of the second anchor
   - start2: Start position of the second anchor
   - end2: End position of the second anchor
   - count: raw loop counts
   - q_value: Adjuested P-value of the interaction

3. WashU Genome Browser compatible tracks (*.interaction.Q0.01.washu.txt.gz)


## Documentation

For more detailed information, refer to the package documentation:

```r
?sintHiChIP
?run_sintHiChIP
```

## Citation

If you use sintHiChIP in your research, please cite:

[Your citation information here]

## License

sintHiChIP is licensed under the GNU General Public License (GPL) v3.0. See the [LICENSE](LICENSE) file for details.

## Contact

For questions and feedback, please contact:

Weiyue Ding  
Harbin Institute of Technology  
Email: wyding0501@hotmail.com

