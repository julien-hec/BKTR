#' @title Function to normalize an R matrix with a min-max method
#'
#' @description Apply a min-max normalization on a R matrix for a given
#' new min and new max
#' @param mat Matrix[numeric]: The matrix on which to apply the normalization
#' @param new_min Numeric: The new mimimal value for the normalized matrix
#' @param new_max Numeric: The new maximal value for the normalized matrix
#' @return A new matrix normalized via the min-max method
#'
#' @export
#' @keywords internal
min_max_normalize <- function(mat, new_min, new_max) {
    normalize_col <- function(x) {
        curr_min <- min(x)
        curr_max <- max(x)
        curr_spread <- curr_max - curr_min
        new_spread <- new_max - new_min
        return((x - curr_min) * new_spread / curr_spread + new_min)
    }
    return(apply(mat, 2, normalize_col))
}

#' @title A helper function to load CSV data as numeric matrix
#'
#' @description Helper utility function to load and transform CSV data into numeric matrix
#' that can be loaded as tensor into the BKTR algorithm
#' @param file_location String: The file path of the source CSV
#' @param has_header Boolean: The new mimimal value for the normalized matrix
#' @param columns_to_drop Vector[String]: A vector of column name that can be removed
#' @param columns_to_keep Vector[String]: A vector of column name that can remain in the matrix
#' @param rows_to_keep Vector[Integer]: A vector of row index that can remain in the matrix
#' @param fill_na_with_zeros Boolean: Replace missing values with zeros
#' @param min_max_normalization Boolean: Apply a (0-1) min-max normalization on the matrix
#' @param transpose_matrix Boolean: Transpose the matrix after all its transformations
#' @return A new matrix on which transformations have been applied
#'
#' @export
load_matrix_from_csv <- function(
    file_location,
    has_header = TRUE,
    columns_to_drop = c(),
    columns_to_keep = NULL,
    rows_to_keep = NULL,
    fill_na_with_zeros = FALSE,
    min_max_normalization = FALSE,
    transpose_matrix = FALSE
) {
    df <- read.csv(file_location, header = has_header)
    df <- df[, !(names(df) %in% columns_to_drop)]
    if (!is.null(columns_to_keep)) {
        df <- df[, columns_to_keep]
    }
    if (!is.null(rows_to_keep)) {
        df <- df[rows_to_keep, ]
    }
    if (fill_na_with_zeros) {
        df[is.na(df)] <- 0
    }
    mat <- as.matrix(df)
    if (transpose_matrix) {
        mat <- t(mat)
    }
    if (min_max_normalization) {
        mat <- min_max_normalize(mat, 0, 1)
    }
    return(mat)
}

#' @title Function to create a distance matrix from a spatial data matrix
#'
#' @description Function that calculate the distance between each spatial entity used in the BKTR
#' algorithm using the haversine distance method of calculation.
#' @param spatial_matrix Matrix[numeric]: The matrix of the spatial data containing at least columns
#' named 'latitude' and 'longitude'
#' @return A new matrix containing the haversine distance between each spatial entity
#'
#' @export
#' @keywords internal
get_haversine_dist_matrix <- function(spatial_matrix) {
    spatial_tensor <- tsr$new_tensor(spatial_matrix)
    # Taking into account the spatial matrix contains a latitude and longitude column
    approx_earth_radius <- 6371.2
    # a spatial_matrix <- tsr$new_tensor(as.matrix(station_loader$df[, c('latitude', 'longitude')]))
    nb_loc <- spatial_tensor$shape[1]
    loc_values_1 <- torch::torch_deg2rad(spatial_tensor$unsqueeze(2)$expand(c(nb_loc, nb_loc, 2)))
    loc_values_2 <- loc_values_1$permute(c(2, 1, 3))
    loc_distances <- torch::torch_abs(loc_values_1 - loc_values_2)

    a <- (
        (loc_distances[, , 1] / 2)$sin() ** 2
        + loc_values_1[, , 1]$cos() * loc_values_2[, , 1]$cos() * (loc_distances[, , 2] / 2)$sin() ** 2
    )
    return(approx_earth_radius * 2 * torch::torch_atan2(a$sqrt(), (1 - a)$sqrt()))
}