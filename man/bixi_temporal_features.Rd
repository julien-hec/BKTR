\name{bixi_temporal_features}
\alias{bixi_temporal_features}
\docType{data}
\title{
Temporal Features in Montreal applicable to BIXI for 2019
}
\description{
These data represent the temporal features in Montreal applicable to a Montreal
based bike sharing company named BIXI. The data include six features (columns)
for 196 days (rows). The time column represent the label associated to each
captured time for the 2019 season of BIXI. The other columns contain information
about Montral weather and applicable holidays for each day.
}
\usage{data("bixi_temporal_features")}
\format{
  A data frame with 196 observations on the following 6 variables.
  \describe{
    \item{\code{time}}{a IDate}
    \item{\code{humidity}}{a numeric vector}
    \item{\code{max_temp_f}}{a numeric vector}
    \item{\code{mean_temp_c}}{a numeric vector}
    \item{\code{total_precip_mm}}{a numeric vector}
    \item{\code{holiday}}{a numeric vector}
  }
}
\source{
Lei, M., Labbe, A., & Sun, L. (2021). Scalable Spatiotemporally Varying
Coefficient Modelling with Bayesian Kernelized Tensor Regression.
arXiv preprint arXiv:2109.00046.
}
\references{
The weather data is sourced from the Environment and Climate Change
Canada Historical Climate Data website.

The holiday column is specifying if a date is a holiday or not,
according to the Quebec government.
}
\keyword{datasets}
