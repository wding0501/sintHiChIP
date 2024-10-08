% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/significantLoops.R
\name{calculate_density_ratio}
\alias{calculate_density_ratio}
\title{Calculate Density Ratio}
\usage{
calculate_density_ratio(interactions)
}
\arguments{
\item{interactions}{Data frame of interactions}
}
\value{
Vector of ratio densities
}
\description{
This function is to calculate the ratio of two anchor cut site density values.
}
\examples{
\dontrun{
interactions <- data.frame(
  totalAnchorSites_1 = c(10, 20, 30, 40),
  totalAnchorSites_2 = c(15, 25, 35, 20)
)
ratios <- calculate_density_ratio(interactions)
print(ratios)
}
}
\keyword{internal}
