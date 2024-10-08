% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_global.R
\name{filter_intra_global}
\alias{filter_intra_global}
\title{Filter Interactions for Global Analysis}
\usage{
filter_intra_global(outdir, bin_size, prefix)
}
\arguments{
\item{outdir}{Character string. Output directory.}

\item{bin_size}{Numeric. Bin size used for matrix generation.}

\item{prefix}{Character string. Prefix for output files.}
}
\value{
Invisible NULL
}
\description{
This function filters interactions based on distance thresholds for global analysis in sintHiChIP.
}
\examples{
\dontrun{
filter_intra_global(
  outdir = "/home/user/hichip/global_output",
  bin_size = 5000,
  prefix = "sample1"
)
}
}
\keyword{internal}
