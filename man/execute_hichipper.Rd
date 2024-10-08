% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localInteraction.R
\name{execute_hichipper}
\alias{execute_hichipper}
\title{Execute hichipper}
\usage{
execute_hichipper(yaml_file, outdir, keep_temp = FALSE)
}
\arguments{
\item{yaml_file}{Character string. Path to YAML configuration file.}

\item{outdir}{Character string. Output directory.}

\item{keep_temp}{Logical. Keep temporary files.}
}
\value{
Invisible NULL
}
\description{
This function runs the hichipper tool.
}
\examples{
\dontrun{
execute_hichipper(
  yaml_file = "/home/user/hichip/config.yaml",
  outdir = "hichipper_output",
  keep_temp = FALSE
)
}
}
