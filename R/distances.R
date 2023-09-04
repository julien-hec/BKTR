#' @title Distance tensor checks
#'
#' @description Check that two tensors are valid for distance computation
#'
#' @noRd
check_dist_tensor_dimensions <- function(x1, x2, expected_nb_dim = 2, expected_last_dim_shape = NULL) {
    if (!(TSR$is_tensor(x1) && TSR$is_tensor(x2))) {
        stop('Distance params must be tensors')
    }
    if (!(x1$ndim == x2$ndim && x2$ndim == expected_nb_dim)) {
        stop(sprintf('Distance params should have %s dimension(s)', expected_nb_dim))
    }
    if (
        !is.null(expected_last_dim_shape) && !(
            expected_last_dim_shape == x1$shape[-1] && expected_last_dim_shape == x2$shape[-1]
        )
    ) {
        stop(
            sprintf('Distance params last dimension should contain %s elements', expected_last_dim_shape)
        )
    }
}

#' @title Function to compute a tensor's euclidean distance
#'
#' @description Function to compute the euclidean distance between a tensor and its transpose.
#'
#' @noRd
get_euclidean_dist_tsr <- function(x) {
    check_dist_tensor_dimensions(x, x)
    x1 <- x$unsqueeze(1)
    x2 <- x$unsqueeze(1)$transpose(1, 2)
    return((x1 - x2)$pow(2)$sum(3)$sqrt())
}

#' @title Class to project coordinates with mercator projection on a 2D plane
#'
#' @description Project coordinates with mercator projection on a 2D plane for
#' a given scale. Keep track of the scale and the center of the projection to
#' be able to project new coordinates which is useful during interpolation.
#'
#' @noRd
GeoMercatorProjector <- R6::R6Class(
    'GeoMercatorProjector',
    public = list(
        ini_df = NULL,
        x_mid_point = NULL,
        y_mid_point = NULL,
        coords_scale = NULL,
        scale = NULL,
        scaled_ini_df = NULL,
        EARTH_RADIUM_KM = 6371,
        initialize = function(df, scale = 10.0) {
            self$ini_df = df
            km_df = private$km_from_coords_df(df)
            lon_x <- km_df$lon_x
            lat_y <- km_df$lat_y
            x_max <- max(lon_x)
            x_min <- min(lon_x)
            y_max <- max(lat_y)
            y_min <- min(lat_y)
            self$x_mid_point <- (x_min + x_max) / 2
            self$y_mid_point <- (y_min + y_max) / 2
            self$coords_scale = max(x_max - x_min, y_max - y_min)
            self$scale = scale
            self$scaled_ini_df = private$scale_and_center_df(km_df)
        },

        project_new_coords = function(df) {
            km_df <- private$km_from_coords_df(df)
            return(private$scale_and_center_df(km_df))
        }
    ),
    private = list(
        scale_and_center_df = function(df) {
                new_df <- df
                scaling_factor <- self$scale / self$coords_scale
                new_df$lon_x <- (new_df$lon_x - self$x_mid_point) * scaling_factor
                new_df$lat_y <- (new_df$lat_y - self$y_mid_point) * scaling_factor
                return(new_df)
        },

        km_from_coords_df = function(df) {
            if (!('latitude' %in% colnames(df) && 'longitude' %in% colnames(df))) {
                stop('Dataframe must have columns "latitude" and "longitude"')
            }
            new_df <- df
            lons <- TSR$tensor(df$longitude)
            lats <- TSR$tensor(df$latitude)
            x <- (self$EARTH_RADIUM_KM / (2 * pi)) * torch_deg2rad(lons)
            merc_n_y <- torch_log(
                torch_tan(pi / 4 + torch_deg2rad(lats) / 2)
            )
            y <- (self$EARTH_RADIUM_KM / (2 * pi)) * merc_n_y
            new_df$lon_x <- as.numeric(x$cpu())
            new_df$lat_y <- as.numeric(y$cpu())
            new_df <- new_df[, -c('latitude', 'longitude')]
            return(new_df)
        }
    )
)
