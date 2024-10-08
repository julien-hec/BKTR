% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_spatial_betas}
\alias{plot_spatial_betas}
\title{Plot Spatial Beta Coefficients}
\usage{
plot_spatial_betas(
  bktr_reg,
  plot_feature_labels,
  temporal_point_label,
  nb_cols = 1,
  use_dark_mode = TRUE,
  show_figure = TRUE,
  zoom = 11,
  google_token = NULL,
  stadia_token = NULL,
  fig_width = 8.5,
  fig_height = 5.5,
  fig_resolution = 200
)
}
\arguments{
\item{bktr_reg}{BKTRRegressor: BKTRRegressor object.}

\item{plot_feature_labels}{Array: Array of feature labels to plot.}

\item{temporal_point_label}{String: Temporal point label to plot.}

\item{nb_cols}{Integer: The number of columns to use in the facet grid.}

\item{use_dark_mode}{Boolean: Whether to use a dark mode for the geographic map or not. Defaults to TRUE.}

\item{show_figure}{Boolean: Whether to show the figure. Defaults to True.}

\item{zoom}{Integer: Zoom level for the geographic map. Defaults to 11.}

\item{google_token}{String or NULL: Google API token to use for the geographic map. Defaults to NULL.}

\item{stadia_token}{String or NULL: Stadia API token to use for the geographic map. Defaults to NULL.}

\item{fig_width}{Numeric: Figure width when figure is shown. Defaults to 8.5.}

\item{fig_height}{Numeric: Figure height when figure is shown. Defaults to 5.5.}

\item{fig_resolution}{Numeric: Figure resolution PPI. Defaults to 200.}
}
\value{
ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
}
\description{
Create a plot of beta values through space for a given
temporal point and a set of feature labels. We use ggmap under the hood,
so you need to provide a Google or Stadia API token to plot on a map.
See: https://cran.r-project.org/web/packages/ggmap/readme/README.html for more details
on how to get an API token.
}
\examples{
\dontshow{if (FALSE) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Launch MCMC sampling on a light version of the BIXI dataset
bixi_data <- BixiData$new(is_light = TRUE)
bktr_regressor <- BKTRRegressor$new(
  data_df <- bixi_data$data_df,
  spatial_positions_df = bixi_data$spatial_positions_df,
  temporal_positions_df = bixi_data$temporal_positions_df,
  burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
bktr_regressor$mcmc_sampling()

# Plot spatial beta coefficients for the first time point and the two features
# Use your own Google API token instead of 'GOOGLE_API_TOKEN'
plot_spatial_betas(
  bktr_regressor,
  plot_feature_labels = c('mean_temp_c', 'area_park'),
  temporal_point_label = bixi_data$temporal_positions_df$time[1],
  google_token = 'GOOGLE_API_TOKEN')

# We can also use light mode and plot the maps side by side
# Use your own Stadia API token instead of 'STADIA_API_TOKEN'
plot_spatial_betas(
  bktr_regressor,
  plot_feature_labels = c('mean_temp_c', 'area_park', 'total_precip_mm'),
  temporal_point_label = bixi_data$temporal_positions_df$time[10],
  use_dark_mode = FALSE, nb_cols = 3, stadia_token = 'STADIA_API_TOKEN')
\dontshow{\}) # examplesIf}
}
