#' @title Distance Types
#'
#' @description Distance Type Enum. Possibilities of distance functions between
#' two spatial or temporal vectors.
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

#' @title Function to compute the euclidean distance between two tensors
#'
#' @description Distance Type Enum. Possibilities of distance functions between
#' two spatial or temporal vectors.
get_euclidean_dist_tsr <- function(x) {
    check_dist_tensor_dimensions(x, x)
    x1 <- x$unsqueeze(1)
    x2 <- x$unsqueeze(1)$transpose(1, 2)
    return((x1 - x2)$pow(2)$sum(3)$sqrt())
}
