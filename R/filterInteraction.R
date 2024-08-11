#' Filter Distance for Global Analysis
#'
#' This function filters interactions based on distance thresholds for global analysis in sintHiChIP.
#'
#' @param InpBinIntervalFile Character string. Path to input bin interval file.
#' @param InpMatrixFile Character string. Path to input matrix file.
#' @param Interaction_Initial_File Character string. Path to initial interactions output file.
#' @param Interaction_File Character string. Path to filtered interactions output file.
#' @param LowDistThres Numeric. Minimum distance threshold for filtering.
#' @param UppDistThres Numeric. Maximum distance threshold for filtering.
#' @return Invisible NULL
#' @export
#' @importFrom data.table fread
#' @importFrom dplyr filter
#' @examples
#' \dontrun{
#' filter_distance_global(
#'   InpBinIntervalFile = "/home/user/hichip/bin_intervals.bed",
#'   InpMatrixFile = "/home/user/hichip/matrix.txt",
#'   Interaction_Initial_File = "/home/user/hichip/initial_interactions.txt",
#'   Interaction_File = "/home/user/hichip/filtered_interactions.txt",
#'   LowDistThres = 20000,
#'   UppDistThres = 2000000
#' )
#' }
filter_distance_global <- function(InpBinIntervalFile, InpMatrixFile, Interaction_Initial_File, 
                                Interaction_File, LowDistThres, UppDistThres) {
  # Declare variables to avoid R CMD check warnings
  chr1 <- NULL
  chr2 <- NULL
  s1 <- NULL
  s2 <- NULL
  create_interaction_matrix(InpBinIntervalFile, InpMatrixFile, Interaction_Initial_File)
  
  interaction_data <- data.table::fread(Interaction_Initial_File, header = TRUE)
  
  interaction_data$s1 <- as.numeric(as.character(interaction_data$s1))
  interaction_data$s2 <- as.numeric(as.character(interaction_data$s2))
  
  if (any(is.na(interaction_data$s1)) || any(is.na(interaction_data$s2))) {
    warning("Some values in columns 's1' or 's2' could not be converted to numeric and have been turned to NA.")
  }
  
  filtered_data <- dplyr::filter(
    interaction_data,
    chr1 == chr2 & abs(s1 - s2) >= LowDistThres & abs(s1 - s2) <= UppDistThres
  )
  
  write.table(filtered_data, file = Interaction_File, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  message("Filtering completed. Filtered interactions saved to: ", Interaction_File)
  
  invisible(NULL)
}