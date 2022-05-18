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
        tensor_dtype = torch::torch_float64(),
        tensor_device = 'cpu'
    ),
    public = list(
        initialize = function(
            tensor_dtype = torch::torch_float64(),
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
                    dtype = private$tensor_dtype,
                    device = private$tensor_device
                )
            )
        },
        zeros = function(tensor_dim) {
            return(
                torch::torch_zeros(
                    tensor_dim,
                    dtype = private$tensor_dtype,
                    device = private$tensor_device
                )
            )
        },
        ones = function(tensor_dim) {
            return(
                torch::torch_ones(
                    tensor_dim,
                    dtype = private$tensor_dtype,
                    device = private$tensor_device
                )
            )
        },

        eye = function(eye_dim) {
            return(
                torch::torch_eye(
                    eye_dim,
                    dtype = private$tensor_dtype,
                    device = private$tensor_device
                )
            )
        },

        arange = function(start, end) {
            # TODO create an issue on github for arange in torch_float64
            return(self$new_tensor(
                torch::torch_arange(start, end, dtype = torch::torch_float32())
            ))
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
        }
    )
)

# Singleton containing all the information for the used API, Tensor default type and device
tsr <- TensorOperator$new()
