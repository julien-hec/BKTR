#' @import torch
#' @importFrom R6 R6Class
#' @importFrom R6P Singleton

#' @title R6 singleton that contains the configuration for the tensor backend
#'
#' @description Tensor backend configuration and methods for all the tensor operations
#' in BKTR
#'
#' @export
TensorOperator <- R6::R6Class(
    'TensorOperator',
    inherit = R6P::Singleton,
    private = list(
        dtype = NULL,
        device = NULL
    ),
    public = list(
        initialize = function(
            dtype = torch::torch_float64,
            device = 'cpu'
        ) {
            private$dtype <- dtype
            private$device <- device
        },

        set_params = function(
            dtype = NULL,
            device = NULL,
            seed = NULL
        ) {
            if (!is.null(dtype)) {
                private$dtype <- dtype
            }
            if (!is.null(device)) {
                private$device <- device
            }
            if (!is.null(seed)) {
                torch::torch_manual_seed(seed)
            }
        },

        tensor = function(tensor_data) {
            return(
                torch::torch_tensor(
                    tensor_data,
                    dtype = private$dtype(),
                    device = private$device
                )
            )
        },

        is_tensor = function(tensor) {
            return(is(tensor, 'torch_tensor'))
        },

        eye = function(eye_dim) {
            return(
                torch::torch_eye(
                    eye_dim,
                    dtype = private$dtype(),
                    device = private$device
                )
            )
        },

        ones = function(tsr_dim) {
            return(
                torch::torch_ones(
                    tsr_dim,
                    dtype = private$dtype(),
                    device = private$device
                )
            )
        },


        zeros = function(tsr_dim) {
            return(
                torch::torch_zeros(
                    tsr_dim,
                    dtype = private$dtype(),
                    device = private$device
                )
            )
        },

        rand = function(tsr_dim) {
            return(torch::torch_rand(tsr_dim, dtype = private$dtype(), device = private$device))
        },

        randn = function(tsr_dim) {
            return(torch::torch_randn(tsr_dim, dtype = private$dtype(), device = private$device))
        },

        randn_like = function(input_tensor) {
            return(torch::torch_randn_like(input_tensor, dtype = private$dtype(), device = private$device))
        },

        arange = function(start, end) {
            return(self$tensor(torch::torch_arange(start, end)))
        },

        kronecker_prod = function(a, b) {
            return(torch::torch_kron(a, b))
        },

        khatri_rao_prod = function(a, b) {
            if (a$shape[2] != b$shape[2]) {
                stop(
                    sprintf(
                        'Matrices must have the same number of columns to perform khatri rao product, got %i and %i',
                        a$shape[2], b$shape[2]
                    )
                )
            }
            return(torch::torch_reshape(
                torch::torch_einsum("ac,bc->abc", c(a, b)),
                c(-1, a$shape[2])
            ))
        }
    ),

    active = list(
        default_jitter = function() {
            if (private$dtype() == torch::torch_float64()) {
                return(1e-8)
            } else if (private$dtype() == torch::torch_float32()) {
                return(1e-4)
            }
            stop('The dtype used by TSR has no default mapped jitter value')
        }
    )
)

#' @title Tensor Operator Singleton
#'
#' @description Singleton instance of the \code{TensorOperator} class that contains
#' all informations related the tensor API; tensor methods, used data type and used device.
#'
#' @export
TSR <- TensorOperator$new()
