#' @importFrom R6 R6Class
#' @import data.table

normalize_0_1 <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
}

#' @title BIXI Data Class
#'
#' @description R6 class encapsulating the BIXI data
#'
#' @export
BixiData <- R6::R6Class(
    'BixiData',
    public = list(
        departure_df = NULL,
        spatial_features_df = NULL,
        temporal_features_df = NULL,
        spatial_positions_df = NULL,
        temporal_positions_df = NULL,
        data_df = NULL,

        initialize = function() {
            self$departure_df <- BKTR::bixi_station_departures
            # Normalize departure counts
            max_val <- max(self$departure_df[, !c('location')], na.rm = TRUE)
            self$departure_df[, .SD / max_val, .SDcols = !c('location')]

            # Normalize spatial features
            self$spatial_features_df <- BKTR::bixi_spatial_features
            self$spatial_features_df[, lapply(.SD, normalize_0_1), .SDcols = !c('location')]

            # Normalize temporal features
            self$temporal_features_df <- BKTR::bixi_temporal_features
            self$temporal_features_df[, lapply(.SD, normalize_0_1), .SDcols = !c('time')]

            self$spatial_positions_df <- BKTR::bixi_spatial_locations
            self$temporal_positions_df <- BKTR::bixi_temporal_locations

            # TODO - reshape data
            self$data_df <- NULL
        }
    )
)
