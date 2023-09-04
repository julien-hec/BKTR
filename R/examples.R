#' @importFrom R6 R6Class
#' @import data.table

# Private function to normalize a vector to [0, 1]
normalize_0_1 <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
}

#' @title BIXI Data Class
#'
#' @description R6 class encapsulating all BIXI dataframes
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

        #' @description Initialize the BIXI data class
        #' @return A new BIXI data instance
        initialize = function() {
            self$departure_df <- BKTR::bixi_station_departures
            # Normalize departure counts to [0, 1]
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

            self$data_df <- reshape_covariate_dfs(
                spatial_df = self$spatial_features_df,
                temporal_df = self$temporal_features_df,
                y_df = self$departure_df,
                y_column_name = 'nb_departure'
            )
        }
    )
)
