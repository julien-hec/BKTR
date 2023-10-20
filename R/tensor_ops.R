#' @import torch
#' @importFrom R6 R6Class
#' @importFrom R6P Singleton

#' @title R6 singleton that contains the configuration for the tensor backend
#'
#' @description Tensor backend configuration and methods for all the tensor operations
#' in BKTR
#'
#' @examplesIf torch::torch_is_installed()
#' # Set the seed, setup the tensor floating point type and device
#' TSR$set_params(fp_type='float64', fp_device='cpu', seed=42)
#' # Create a tensor from a vector
#' TSR$tensor(c(1, 2, 3))
#' # Create a tensor from a matrix
#' TSR$tensor(matrix(c(1, 2, 3, 4), nrow=2))
#' # Create a 3x3 tensor with a diagonal of ones and zeros elsewhere
#' TSR$eye(3)
#' # Create a tensor of ones (with 6 elements, 2 rows and 3 columns)
#' TSR$ones(c(2, 3))
#' # Create a tensor of zeros (with 12 elements, 3 rows and 4 columns)
#' TSR$zeros(c(3, 4))
#' # Create a tensor of random uniform values (with 6 elements)
#' TSR$rand(c(2, 3))
#' # Create a tensor of random normal values (with 6 elements)
#' TSR$randn(c(2, 3))
#' # Create a tensor of random normal values with the same shape as a given tensor
#' tsr_a <- TSR$randn(c(2, 3))
#' TSR$randn_like(tsr_a)
#' # Create a tensor of a range of values (1, 2, 3, 4)
#' TSR$arange(1, 4)
#' # Choose two random values from a given tensor without replacement
#' tsr_b <- TSR$rand(6)
#' TSR$rand_choice(tsr_b, 2)
#' # Use the tensor operator to compute the kronecker product of two 2x2 matrices
#' tsr_c <- TSR$tensor(matrix(c(1, 2, 3, 4), nrow=2))
#' tsr_d <- TSR$tensor(matrix(c(5, 6, 7, 8), nrow=2))
#' TSR$kronecker_prod(tsr_c, tsr_d) # Returns a 4x4 tensor
#' # Use the tensor operator to compute the khatri rao product of two 2x2 matrices
#' TSR$khatri_rao_prod(tsr_c, tsr_d) # Returns a 4x2 tensor
#' # Check if a given object is a tensor
#' TSR$is_tensor(tsr_d) # Returns TRUE
#' TSR$is_tensor(TSR$eye(2)) # Returns TRUE
#' TSR$is_tensor(1) # Returns FALSE
#'
#' @export
TensorOperator <- R6::R6Class(
    'TensorOperator',
    inherit = R6P::Singleton,
    private = list(
        # Private values used to store the config of the underlying tensor library
        dtype = NULL,
        device = NULL
    ),
    public = list(
        # Public facing config values used streamlined for TensorOperator
        #' @field fp_type The floating point type to use for the tensor operations
        fp_type = NULL,
        #' @field fp_device The device to use for the tensor operations
        fp_device = NULL,

        #' @description Initialize the tensor operator with the given floating point type
        #' and device
        #' @param fp_type The floating point type to use for the tensor operations (either
        #' "float64" or "float32")
        #' @param fp_device The device to use for the tensor operations (either "cpu" or
        #' "cuda")
        #' @return A new tensor operator instance
        initialize = function(
            fp_type = 'float64',
            fp_device = 'cpu'
        ) {
            self$set_params(fp_type, fp_device)
        },

        #' @description Set the tensor operator parameters
        #' @param fp_type The floating point type to use for the tensor operations (either
        #' "float64" or "float32")
        #' @param fp_device The device to use for the tensor operations (either "cpu" or
        #' "cuda")
        #' @param seed The seed to use for the random number generator
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

        #' @description Get the default jitter value for the floating point type used by the tensor operator
        #' @return The default jitter value for the floating point type used by the tensor operator
        get_default_jitter = function() {
            # I wanted this to be an active binding but it is compiled at build time
            # See: https://github.com/r-lib/R6/issues/152
            if (private$dtype() == torch::torch_float64()) {
                return(1e-8)
            } else if (private$dtype() == torch::torch_float32()) {
                return(1e-4)
            }
            stop('The dtype used by TSR has no default mapped jitter value')
        },

        #' @description Create a tensor from a vector or matrix of data with the tensor operator dtype and device
        #' @param tensor_data The vector or matrix of data to create the tensor from
        #' @return A new tensor with the tensor operator dtype and device
        tensor = function(tensor_data) {
            return(
                torch::torch_tensor(
                    tensor_data,
                    dtype = private$dtype(),
                    device = private$device
                )
            )
        },

        #' @description Check if a provided object is a tensor
        #' @param tensor The object to check
        #' @return A boolean indicating if the object is a tensor
        is_tensor = function(tensor) {
            return(is(tensor, 'torch_tensor'))
        },

        #' @description Create a tensor with a diagonal of ones and zeros with the tensor operator dtype and device
        #' for a given dimension
        #' @param eye_dim The dimension of the tensor to create
        #' @return A new tensor with a diagonal of ones and zeros with the tensor operator dtype and device
        eye = function(eye_dim) {
            return(
                torch::torch_eye(
                    eye_dim,
                    dtype = private$dtype(),
                    device = private$device
                )
            )
        },

        #' @description Create a tensor of ones with the tensor operator dtype and device for a given dimension
        #' @param tsr_dim The dimension of the tensor to create
        #' @return A new tensor of ones with the tensor operator dtype and device
        ones = function(tsr_dim) {
            return(
                torch::torch_ones(
                    tsr_dim,
                    dtype = private$dtype(),
                    device = private$device
                )
            )
        },

        #' @description Create a tensor of zeros with the tensor operator dtype and device for a given dimension
        #' @param tsr_dim The dimension of the tensor to create
        #' @return A new tensor of zeros with the tensor operator dtype and device
        zeros = function(tsr_dim) {
            return(
                torch::torch_zeros(
                    tsr_dim,
                    dtype = private$dtype(),
                    device = private$device
                )
            )
        },

        #' @description Create a tensor of random uniform values with the tensor operator dtype and
        #' device for a given dimension
        #' @param tsr_dim The dimension of the tensor to create
        #' @return A new tensor of random values with the tensor operator dtype and device
        rand = function(tsr_dim) {
            return(torch::torch_rand(tsr_dim, dtype = private$dtype(), device = private$device))
        },

        #' @description Create a tensor of random normal values with the tensor operator dtype and device
        #' for a given dimension
        #' @param tsr_dim The dimension of the tensor to create
        #' @return A new tensor of random normal values with the tensor operator dtype and device
        randn = function(tsr_dim) {
            return(torch::torch_randn(tsr_dim, dtype = private$dtype(), device = private$device))
        },

        #' @description Create a tensor of random uniform values with the same shape as a given tensor
        #' with the tensor operator dtype and device
        #' @param input_tensor The tensor to use as a shape reference
        #' @return A new tensor of random uniform values with the same shape as a given tensor
        randn_like = function(input_tensor) {
            return(torch::torch_randn_like(input_tensor, dtype = private$dtype(), device = private$device))
        },

        #' @description Create a tensor of a range of values with the tensor operator dtype and device
        #' for a given start and end
        #' @param start The start of the range
        #' @param end The end of the range
        #' @return A new tensor of a range of values with the tensor operator dtype and device
        arange = function(start, end) {
            return(self$tensor(torch::torch_arange(start, end)))
        },

        #' @description Choose random values from a tensor for a given number of samples
        #' @param choices_tsr The tensor to choose values from
        #' @param nb_sample The number of samples to choose
        #' @param use_replace A boolean indicating if the sampling should be done with replacement.
        #' Defaults to FALSE
        #' @param weights_tsr The weights to use for the sampling. If NULL, the sampling is uniform.
        #' Defaults to NULL
        #' @return A new tensor of randomly chosen values from a tensor
        rand_choice = function(choices_tsr, nb_sample, use_replace = FALSE, weights_tsr = NULL) {
            if (is.null(weights_tsr)) {
                weights_tsr <- self$ones(choices_tsr$shape)
            }
            if (choices_tsr$shape != weights_tsr$shape) {
                stop('Choices and weights tensors must have the same shape')
            }
            choices_indx <- torch::torch_multinomial(weights_tsr, nb_sample, use_replace)
            return(choices_tsr[choices_indx])
        },

        #' @description Efficiently compute the kronecker product of two matrices in tensor format
        #' @param a The first tensor
        #' @param b The second tensor
        #' @return The kronecker product of the two matrices
        kronecker_prod = function(a, b) {
            return(torch::torch_kron(a, b))
        },

        #' @description Efficiently compute the khatri rao product of two matrices in tensor format
        #' having the same number of columns
        #' @param a The first tensor
        #' @param b The second tensor
        #' @return The khatri rao product of the two matrices
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
