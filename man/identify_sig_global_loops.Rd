% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_global.R
\name{identify_sig_global_loops}
\alias{identify_sig_global_loops}
\title{Identify Significant Global Loops}
\usage{
identify_sig_global_loops(outdir, normSiteFile, FDR, prefix)
}
\arguments{
\item{outdir}{Character string. Output directory.}

\item{normSiteFile}{Character string. Path to normalization site file.}

\item{FDR}{Numeric. False Discovery Rate threshold.}

\item{prefix}{Character string. Prefix for output files.}
}
\value{
Invisible NULL
}
\description{
This function identifies significant loops in the global analysis of sintHiChIP.
}
\examples{
\dontrun{
identify_sig_global_loops(
  outdir = "/home/user/hichip/global_output",
  normSiteFile = "/home/user/hichip/norm_sites.bed",
  FDR = 0.01,
  prefix = "sample1"
)
}
}
\keyword{internal}
