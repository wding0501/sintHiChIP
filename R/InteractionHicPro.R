#' Create Interaction Matrix from HiC-pro Output
#'
#' This function creates an interaction matrix (chromosome intervals + contact count)
#' from a given HiC-pro output matrix.
#'
#' @param intervals Path to the file storing the indices of individual bins of chromosomes.
#' @param matrix Path to the file storing the interactions between bins (matrix).
#' @param output Path to the output file where the intervals + contact count will be written.
#' @return Invisible NULL
#' @export
#' @importFrom data.table fread setkey fwrite
#' @examples
#' \dontrun{
#' create_interaction_matrix(intervals = "/path/to/intervals.txt",
#'                           matrix = "/path/to/matrix.txt",
#'                           output = "/path/to/output_interactions.txt")
#' }
create_interaction_matrix <- function(intervals, matrix, output) {
  # Declare variables to avoid R CMD check warnings
  idx <- NULL
  idx1 <- NULL
  idx2 <- NULL
  cc <- NULL
  chr1 <- NULL
  s1 <- NULL
  e1 <- NULL
  i.chr1 <- NULL
  i.s1 <- NULL
  i.e1 <- NULL
  IntervalMat <- data.table::fread(intervals, header = FALSE, 
                                   col.names = c("chr1", "s1", "e1", "idx"))
  InpInteraction <- data.table::fread(matrix, header = FALSE, 
                                      col.names = c("idx1", "idx2", "cc"))
  
  data.table::setkey(IntervalMat, idx)
  data.table::setkey(InpInteraction, idx1)
  df1 <- InpInteraction[IntervalMat, nomatch = 0]
  df1 <- df1[, list(idx2, cc, chr1, s1, e1)]
  
  data.table::setkey(df1, idx2)
  df2 <- df1[IntervalMat, nomatch = 0]
  df2 <- df2[, list(chr1, s1, e1, chr2 = i.chr1, s2 = i.s1, e2 = i.e1, cc)]
  
  data.table::fwrite(df2, output, sep = "\t", quote = FALSE)
  
  invisible(NULL)
}