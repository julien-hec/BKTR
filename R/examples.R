#' @importFrom R6 R6Class
#' @import data.table

# Private function to normalize a vector to [0, 1]
normalize_0_1 <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
}

#' @title BIXI Data Class
#'
#' @description R6 class encapsulating all BIXI dataframes. It is also
#' possible to use a light version of the dataset by using the \code{is_light}
#' parameter. In this case, the dataset is reduced to its first 25 stations
#' and first 50 days. The light version is only used for testing and short examples.
#'
#' @examples
#' # Create a light BIXI data collection instance containing multiple dataframes
#' # This only uses the first 25 stations and 50 days of the full dataset
#' bixi_data <- BixiData$new(is_light = TRUE)
#' # Dataframe containing the position (latitude and longitude) of M stations
#' bixi_data$spatial_positions_df
#' # Dataframe containing the time position of N days (O to N-1)
#' bixi_data$temporal_positions_df
#' # Dataframe with spatial and temporal features for each day and station (M x N rows)
#' bixi_data$data_df
#'
#' @export
BixiData <- R6::R6Class(
    'BixiData',
    public = list(
        #' @field departure_df The departure dataframe
        departure_df = NULL,
        #' @field spatial_features_df The spatial features dataframe
        spatial_features_df = NULL,
        #' @field temporal_features_df The temporal features dataframe
        temporal_features_df = NULL,
        #' @field spatial_positions_df The spatial positions dataframe
        spatial_positions_df = NULL,
        #' @field temporal_positions_df The temporal positions dataframe
        temporal_positions_df = NULL,
        #' @field data_df The data dataframe
        data_df = NULL,
        #' @field is_light Whether the light version of the dataset is used
        is_light = FALSE,

        #' @description Initialize the BIXI data class
        #'
        #' @param is_light Whether the light version of the dataset is used,
        #' defaults to FALSE.
        #'
        #' @return A new BIXI data instance
        initialize = function(is_light = FALSE) {
            self$is_light <- is_light

            # Normalize departure counts to [0, 1]
            self$departure_df <- BKTR::bixi_station_departures
            max_val <- max(self$departure_df[, !c('location')], na.rm = TRUE)
            cols <- colnames(self$departure_df[, !c('location')])
            self$departure_df[, (cols) := .SD / max_val, .SDcols = cols]

            # Normalize spatial features column wise to [0, 1]
            self$spatial_features_df <- BKTR::bixi_spatial_features
            cols <- colnames(self$spatial_features_df[, !c('location')])
            self$spatial_features_df[, (cols) := lapply(.SD, normalize_0_1), .SDcols = cols]

            # Normalize temporal features column wise to [0, 1]
            self$temporal_features_df <- BKTR::bixi_temporal_features
            cols <- colnames(self$temporal_features_df[, !c('time')])
            self$temporal_features_df[, (cols) := lapply(.SD, normalize_0_1), .SDcols = cols]

            self$spatial_positions_df <- BKTR::bixi_spatial_locations
            self$temporal_positions_df <- BKTR::bixi_temporal_locations

            # Reduce the dataset to its first 25 stations and first 50 days when used for examples
            if (is_light) {
                self$spatial_positions_df <- self$spatial_positions_df[1:25, ]
                self$temporal_positions_df <- self$temporal_positions_df[1:50, ]
                kept_locs <- self$spatial_positions_df$location
                kept_times <- self$temporal_positions_df$time
                self$spatial_features_df <- self$spatial_features_df[self$spatial_features_df$location %in% kept_locs, ]
                self$temporal_features_df <- self$temporal_features_df[self$temporal_features_df$time %in% kept_times, ]
                # Filter rows and columns of departure_df
                use_dep_rows <- self$departure_df$location %in% kept_locs
                use_dep_cols <- colnames(self$departure_df) %in% c('location', as.character(kept_times))
                self$departure_df <- self$departure_df[use_dep_rows, use_dep_cols, with = FALSE]
            }

            self$data_df <- reshape_covariate_dfs(
                spatial_df = self$spatial_features_df,
                temporal_df = self$temporal_features_df,
                y_df = self$departure_df,
                y_column_name = 'nb_departure'
            )
        }
    )
)
