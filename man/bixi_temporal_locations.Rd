\name{bixi_temporal_locations}
\alias{bixi_temporal_locations}
\docType{data}
\title{
Temporal indices for the 2019 BIXI season
}
\description{
These data represent 196 temporal indices (rows) related to each day
of the 2019 season of Montreal based bike sharing company named BIXI.
The time column represent the label associated to each day and the time_index
column represent the location in time space of each day when compared to each
other. Since no days are missing and they are all spaced by exactly one day,
the time_index is simply a range from 0 to 195.

}
\usage{data("bixi_temporal_locations")}
\format{
  A data frame with 196 observations on the following 2 variables.
  \describe{
    \item{\code{time}}{a IDate}
    \item{\code{time_index}}{a numeric vector}
  }
}
\source{
Lei, M., Labbe, A., & Sun, L. (2021). Scalable Spatiotemporally Varying
Coefficient Modelling with Bayesian Kernelized Tensor Regression.
arXiv preprint arXiv:2109.00046.
}
\keyword{datasets}
