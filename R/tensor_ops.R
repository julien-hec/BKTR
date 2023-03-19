#' @import torch
#' @importFrom R6 R6Class
#' @importFrom R6P Singleton

#' @title R6 singleton that contains the configuration for the tensor backend
#'
#' @description Tensor backend configuration and methods for all the tensor operations
#' in BKTR
#'
#' @export
#' @keywords internal
TensorOperator <- R6::R6Class(
    'TensorOperator',
    inherit = R6P::Singleton,
    private = list(
        tensor_dtype = NULL,
        tensor_device = NULL
    ),
    public = list(
        initialize = function(
            tensor_dtype = torch::torch_float64,
            tensor_device = 'cpu'
        ) {
            private$tensor_dtype <- tensor_dtype
            private$tensor_device <- tensor_device
        },

        set_tensor_type = function(tensor_dtype) {
            private$tensor_dtype <- tensor_dtype
            return(self)
        },

        set_device_type = function(tensor_device) {
            private$tensor_device <- tensor_device
            return(self)
        },

        new_tensor = function(tensor_dim) {
            return(
                torch::torch_tensor(
                    tensor_dim,
                    dtype = private$tensor_dtype(),
                    device = private$tensor_device
                )
            )
        },

        is_tensor = function(tensor) {
            return(is(tensor, 'torch_tensor'))
        },

        zeros = function(tensor_dim) {
            return(
                torch::torch_zeros(
                    tensor_dim,
                    dtype = private$tensor_dtype(),
                    device = private$tensor_device
                )
            )
        },
        ones = function(tensor_dim) {
            return(
                torch::torch_ones(
                    tensor_dim,
                    dtype = private$tensor_dtype(),
                    device = private$tensor_device
                )
            )
        },

        eye = function(eye_dim) {
            return(
                torch::torch_eye(
                    eye_dim,
                    dtype = private$tensor_dtype(),
                    device = private$tensor_device
                )
            )
        },

        rand = function(size) {
            return(torch::torch_rand(size, dtype = private$tensor_dtype(), device = private$tensor_device))
        },

        randn = function(size) {
            return(torch::torch_randn(size, dtype = private$tensor_dtype(), device = private$tensor_device))
        },

        randn_like = function(input_tensor) {
            return(torch::torch_randn_like(input_tensor, dtype = private$tensor_dtype(), device = private$tensor_device))
        },

        arange = function(start, end) {
            return(self$new_tensor(torch::torch_arange(start, end)))
        },

        kronecker_prod = function(a, b) {
            kron_prod <- torch::torch_einsum("ab,cd->acbd", c(a, b))
            kron_shape <- c(a$shape[1] * b$shape[1], a$shape[2] * b$shape[2])
            return(torch::torch_reshape(kron_prod, kron_shape))
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
        },

        get_df_tensor_or_null = function(df) {
            if (is.null(df)) {
                return(NULL)
            }
            return(self$new_tensor(as.matrix(df)))
        }
    ),

    active = list(
        default_jitter = function() {
            if (private$tensor_dtype() == torch::torch_float64()) {
                return(1e-8)
            } else if (private$tensor_dtype() == torch::torch_float32()) {
                return(1e-4)
            }
            stop('The dtype used by TSR has no default mapped jitter value')
        }
    )
)

# Singleton containing all the information for the used API, Tensor default type and device
tsr <- TensorOperator$new()
