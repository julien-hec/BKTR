\name{bixi_spatial_locations}
\alias{bixi_spatial_locations}
\docType{data}
\title{
Spatial Locations of Montreal BIXI Stations in 2019
}
\description{
Data points representing the spatial locations of 587 bike sharing stations for
the Montreal based bike sharing company named BIXI. The dataframe contains
a label to identify each station and its associated longitude and latitude
coordinates.
}
\usage{data("bixi_spatial_locations")}
\format{
  A data frame with 587 observations on the following 3 variables.
  \describe{
    \item{\code{location}}{a character vector}
    \item{\code{latitude}}{a numeric vector}
    \item{\code{longitude}}{a numeric vector}
  }
}
\source{
Wang, X., Cheng, Z., Trépanier, M., & Sun, L. (2021).
Modeling bike-sharing demand using a regression model with spatially
varying coefficients. Journal of Transport Geography, 93, 103059.
}
\references{
BIXI Montréal (2023). “Open Data.” Accessed: 2023-07-11, URL
\url{https://bixi.com/en/} open-data.
}
\keyword{datasets}
