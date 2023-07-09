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
    # Private values used to store the config of the underlying tensor library
    private = list(
        dtype = NULL,
        device = NULL
    ),
    public = list(
        # Public facing config values used streamlined for TensorOperator
        fp_type = NULL,
        fp_device = NULL,

        initialize = function(
            fp_type = 'float64',
            fp_device = 'cpu'
        ) {
            self$set_params(fp_type, fp_device)
        },

        set_params = function(
            fp_type = NULL,
            fp_device = NULL,
            seed = NULL
        ) {
            if (!is.null(fp_type)) {
                if (fp_type == 'float64') {
                    private$dtype <- torch::torch_float64
                } else if (fp_type == 'float32') {
                    private$dtype <- torch::torch_float32
                } else {
                    stop('`fp_type` must be either "float64" or "float32"')
                }
                self$fp_type <- fp_type
            }
            if (!is.null(fp_device)) {
                private$device <- fp_device
                self$fp_device <- fp_device
            }
            if (!is.null(seed)) {
                torch::torch_manual_seed(seed)
                # This is for rWishart until it is implemented in R Torch
                set.seed(seed)
            }
        },

        # I wanted this to be an active binding but it is compiled at build time
        # See: https://github.com/r-lib/R6/issues/152
        get_default_jitter = function() {
            if (private$dtype() == torch::torch_float64()) {
                return(1e-8)
            } else if (private$dtype() == torch::torch_float32()) {
                return(1e-4)
            }
            stop('The dtype used by TSR has no default mapped jitter value')
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
)

#' @title Tensor Operator Singleton
#'
#' @description Singleton instance of the \code{TensorOperator} class that contains
#' all informations related the tensor API; tensor methods, used data type and used device.
#'
#' @export
TSR <- TensorOperator$new()
