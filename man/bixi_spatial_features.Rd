\name{bixi_spatial_features}
\alias{bixi_spatial_features}
\docType{data}
\title{
Spatial Features of Montreal BIXI Stations in 2019
}
\description{
These data represent 14 spatial features (columns) for 587 bike sharing
stations (rows) located at different geographical coordinates (longitude,
latitude) in Montreal. The Montreal based bike sharing company is named BIXI.
The first column contains the descriptive label affected to each station and the
other columns contain information about the infrastructure, points of interests,
walkscore and population surrounding each station for 2019.
}
\usage{data("bixi_spatial_features")}
\format{
  A data frame with 587 observations on the following 14 variables.
  \describe{
    \item{\code{location}}{a character vector}
    \item{\code{area_park}}{a numeric vector}
    \item{\code{len_cycle_path}}{a numeric vector}
    \item{\code{len_major_road}}{a numeric vector}
    \item{\code{len_minor_road}}{a numeric vector}
    \item{\code{num_metro_stations}}{a numeric vector}
    \item{\code{num_other_commercial}}{a numeric vector}
    \item{\code{num_restaurants}}{a numeric vector}
    \item{\code{num_university}}{a numeric vector}
    \item{\code{num_pop}}{a numeric vector}
    \item{\code{num_bus_stations}}{a numeric vector}
    \item{\code{num_bus_routes}}{a numeric vector}
    \item{\code{walkscore}}{a numeric vector}
    \item{\code{capacity}}{a numeric vector}
  }
}
\source{
Wang, X., Cheng, Z., Trépanier, M., & Sun, L. (2021).
Modeling bike-sharing demand using a regression model with spatially
varying coefficients. Journal of Transport Geography, 93, 103059.
}
\references{
Reference for the BIXI station informations:
BIXI Montréal (2023). “Open Data.” Accessed: 2023-07-11, URL
\url{https://bixi.com/en/} open-data.

Reference for point of interests and infrastructure informations:
DMTI Spatial Inc (2019). “Enhanced Point of Interest (DMTI).”
URL \url{https://www.dmtispatial.com}.

Reference for Walkscore:
Walk Score (2023). “Walk Score Methodology.” Accessed: 2023-07-11,
URL \url{https://www.walkscore.com/methodology.shtml}.

The population information comes from the 2016 Canada census data at a
dissemination block level.
}
\keyword{datasets}
