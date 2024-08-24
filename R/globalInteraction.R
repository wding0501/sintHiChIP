#' Process Global (Peak-to-All) Interactions in sintHiChIP Global Mode 
#'
#' This function processes the final intra interactions for the downstream significant analysis.
#'
#' @param peaks Path to peak file
#' @param interaction_file Path to interaction file
#' @param output_directory Output directory
#' @param prefix Prefix for output files
#' @return Invisible NULL
#' @export
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom data.table as.data.table
#' @importFrom dplyr mutate filter select
#' @importFrom readr read_delim write_delim
#' @examples
#' \dontrun{
#' overlap_peaks_global(peaks = "/home/user/hichip/peaks.bed",
#'                      interaction_file = "/home/user/hichip/interactions.txt",
#'                      output_directory = "/home/user/hichip/output",
#'                      prefix = "sample1")
#' }
overlap_peaks_global <- function(peaks, interaction_file, output_directory, prefix) {
  chr1 <- NULL
  start1 <- NULL
  end1 <- NULL
  chr2 <- NULL
  start2 <- NULL
  end2 <- NULL
  dot <- NULL
  cc <- NULL
  peaks <- read_delim(peaks, delim = "\t", col_names = c("chr", "start", "end"))
  interactions <- read_delim(interaction_file, delim = "\t", 
                             col_names = c("chr1", "start1", "end1", "chr2", "start2", "end2", "cc"))
  
  peaks_gr <- GRanges(seqnames = peaks$chr, ranges = IRanges(start = peaks$start, end = peaks$end))
  
  results_all <- list()
  results_filtered <- list()
  
  chrs <- unique(intersect(interactions$chr1, peaks$chr))
  
  for (chr in chrs) {
    chr_interactions <- interactions[interactions$chr1 == chr | interactions$chr2 == chr, ]
    
    left_gr <- GRanges(seqnames = chr_interactions$chr1, 
                       ranges = IRanges(start = chr_interactions$start1, end = chr_interactions$end1))
    right_gr <- GRanges(seqnames = chr_interactions$chr2, 
                        ranges = IRanges(start = chr_interactions$start2, end = chr_interactions$end2))
    
    left_overlaps <- findOverlaps(left_gr, peaks_gr)
    right_overlaps <- findOverlaps(right_gr, peaks_gr)
    
    chr_results <- data.frame(
      left_overlap = as.integer(seq_len(nrow(chr_interactions)) %in% queryHits(left_overlaps)),
      right_overlap = as.integer(seq_len(nrow(chr_interactions)) %in% queryHits(right_overlaps)),
      left_index = match(seq_len(nrow(chr_interactions)), queryHits(left_overlaps)),
      right_index = match(seq_len(nrow(chr_interactions)), queryHits(right_overlaps))
    )
    
    chr_results_all <- cbind(chr_interactions, chr_results)
    chr_results_filtered <- chr_results_all[chr_results_all$left_overlap == 1 | chr_results_all$right_overlap == 1, ]
    
    results_all[[chr]] <- chr_results_all
    results_filtered[[chr]] <- chr_results_filtered
  }
  
  all_results <- do.call(rbind, results_all)
  filtered_results <- do.call(rbind, results_filtered)
  
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  
  
  target_data <- filtered_results %>%
    dplyr::mutate(dot = ".") %>%
    dplyr::select(chr1, start1, end1, chr2, start2, end2, dot, cc)
  
  write_delim(filtered_results,
              file.path(output_directory,  paste0(prefix, ".intra.Idx.bedpe")),
              delim = " ",
              col_names = TRUE)
  message("Filtered intra interactions with Index saved to: ",  file.path(output_directory, paste0(prefix, ".intra.Idx.bedpe")))
  write_delim(target_data,
              file.path(output_directory, paste0(prefix, ".filt.intra.loop_counts.bedpe")),
              delim = " ",
              col_names = FALSE)
 
  message("Filtered intra interactions saved to: ",  file.path(output_directory, paste0(prefix, ".filt.intra.loop_counts.bedpe")))
  invisible(NULL)
}