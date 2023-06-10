#' @import torch
#' @importFrom R6 R6Class

EARTH_RADIUS_KM <- 6371.2


#' @title Distance Types
#'
#' @description Distance Type Enum. Possibilities of distance functions between
#' two spatial or temporal vectors.
#'
#' @export
DIST_TYPE <- list(
    'EUCLIDEAN' = 'EUCLIDEAN',
    'HAVERSINE' = 'HAVERSINE',
    'MANHATTAN' = 'MANHATTAN',
    'CHEBYSHEV' = 'CHEBYSHEV',
    'NONE' = 'none'
)

#' @title R6 class for Composed Kernels
#'
#' @description R6 class for Composed Kernels
#'
#' @export
DistanceCalculator <- R6::R6Class(
    public = list(
        get_matrix = function(x, distance_type) {
            switch(
                distance_type,
                EUCLIDEAN = {
                    return(self$calc_euclidean_dist(x, x))},
                MANHATTAN = {
                    return(self$calc_manhattan_dist(x, x))},
                CHEBYSHEV = {
                    return(self$calc_chebyshev_dist(x, x))},
                HAVERSINE = {
                    return(self$calc_haversine_dist(x, x))},
                NONE = {
                    return(NULL)
                }
            )
        },

        check_tensor_dimensions = function(x1, x2, expected_nb_dim, expected_last_dim_shape) {
            if (!(tsr$is_tensor(x1) && tsr$is_tensor(x2))) {
                torch:::value_error('Distance params must be tensors')
            }
            if (!(x1$ndim == x2$ndim && x2$ndim == expected_nb_dim)) {
                torch:::value_error(sprintf('Distance params should have %s dimension(s)', expected_nb_dim))
            }
            if (
                expected_last_dim_shape != NULL && !(
                    expected_last_dim_shape == x1$shape[-1] && expected_last_dim_shape == x2$shape[-1]
                )
            ) {
                torch:::value_error(
                    sprintf('Distance params last dimension should contain %s elements', expected_last_dim_shape)
                )
            }
        },

        calc_euclidean_dist = function(x1, x2) {
            self$check_tensor_dimensions(x1, x2, expected_nb_dim = 2)
            xu1 <- x1$unsqueeze(1)
            xu2 <- x2$unsqueeze(1)
            return((xu1 - xu2$transpose(1, 2))$pow(2)$sum(3)$sqrt())
        },

        calc_haversine_dist = function(x1, x2, earth_radius = EARTH_RADIUS_KM) {
            self$check_tensor_dimensions(x1, x2, expected_nb_dim=2, expected_last_dim_shape=2)

            xu1 <- x1$unsqueeze(1)
            xu2 <- x2$unsqueeze(1)
            xu1 <- torch::torch_deg2rad(xu1)
            xu2 <- torch::torch_deg2rad(xu2)
            xu2 <- xu2$transpose(1, 2)

            dist <- (xu1 - xu2)$abs()
            a <- (dist[, , 1] / 2)$sin()$pow(2) + xu1[, , 1]$cos() * xu2[, , 1]$cos() * (dist[, , 2] / 2)$sin()$pow(2)
            return(earth_radius * 2 * torch::torch_atan2(a$sqrt(), (1 - a)$sqrt()))
        },

        calc_manhattan_dist = function(x1, x2) {
            self$check_tensor_dimensions(x1, x2, expected_nb_dim = 2)
            xu1 <- x1$unsqueeze(1)
            xu2 <- x2$unsqueeze(1)
            return((xu1 - xu2$transpose(1, 2))$abs()$sum(3))
        },

        calc_chebyshev_dist = function(x1, x2) {
            self$check_tensor_dimensions(x1, x2, expected_nb_dim = 2)
            xu1 <- x1$unsqueeze(1)
            xu2 <- x2$unsqueeze(1)
            return((xu1 - xu2$transpose(1, 2))$abs()$max(3)[0])
        }
    )
)
