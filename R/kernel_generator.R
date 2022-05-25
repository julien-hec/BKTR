#' @importFrom R6 "R6Class"
#' @import torch

#' @title R6 class that generates Temporal Kernels
#'
#' @description A TemporalKernelGenerator can create temporal kernels according
#' to provided parameters. A temporal kernel can use a periodic kernel function, a
#' squared exponential (SE) kernel function or a mixture of both.
#'
#' @export
TemporalKernelGenerator <- R6::R6Class(
    'TemporalKernelGenerator',

    public = list(
        time_distances = NULL,
        period_length = NULL,
        kernel = NULL,
        kernel_variance = NULL,
        core_kernel_fn = NULL,
        # Following params are set externally via HyperParamSampler
        periodic_length_scale = NULL,
        decay_time_scale = NULL,

        #' @description Create a new \code{TemporalKernelGenerator} object.
        #' @param kernel_fn_name String ('periodic', 'se', 'periodic_se'): The kernel function's name
        #' @param nb_time_segments Integer: The number of time segments for the kernel (e,g. number of days)
        #' @param period_length Integer: The period's length for the kernel (e,g. 7 weekly patterns)
        #' @param kernel_variance Numeric: The variance used in the kernel
        #' @param ini_period_length_scale Numeric: Initial value for the periodic length scale
        #' @param ini_decay_time_scale Numeric: Initial value for the decay time scale
        #' @return A new \code{TemporalKernelGenerator} object.
        initialize = function(
            kernel_fn_name,
            nb_time_segments,
            period_length,
            kernel_variance,
            ini_period_length_scale = NULL,
            ini_decay_time_scale = NULL
        ) {
            self$core_kernel_fn <- private$get_kernel_fn(kernel_fn_name)
            private$set_time_distance_matrix(nb_time_segments)
            self$period_length <- period_length
            self$kernel_variance <- kernel_variance
            self$periodic_length_scale <- ini_period_length_scale
            self$decay_time_scale <- ini_decay_time_scale
        },

        #' @description Generate and set temporal kernel values using current parameters
        #' @return A new numeric matrix containing the temporal kernel
        kernel_gen = function() {
            self$kernel <- self$kernel_variance * self$core_kernel_fn()
            return(self$kernel)
        }
    ),

    private = list(
        #~ @description Get the kernel generating function according to its name
        #~ @param kernel_fn_name String ('periodic', 'se', 'periodic_se'): The kernel function's name
        #~ @return A function that can generate a kernel
        get_kernel_fn = function(kernel_fn_name) {
            if (kernel_fn_name == 'periodic') {
                return(private$periodic_kernel_gen)
            } else if (kernel_fn_name == 'se') {
                return(private$se_kernel_gen)
            } else if (kernel_fn_name == 'periodic_se') {
                return(private$periodic_se_kernel_gen)
            } else {
                torch:::value_error('Please select a valid temporal kernel fonction name')
            }
        },

        #~ @description Calculate and set the time distance matrix used in the temporal kernel
        #~ @param nb_time_segments Integer: The number of time segments for the kernel
        set_time_distance_matrix = function(nb_time_segments) {
            time_segments <- tsr$arange(1, nb_time_segments)$unsqueeze(2)
            self$time_distances <- time_segments - time_segments$t()
        },

        #~ @description Generate a periodic temporal kernel matrix
        #~ @return A matrix containing the generated periodic kernel
        periodic_kernel_gen = function() {
            return(torch::torch_exp(
                -2 * torch::torch_sin(pi * self$time_distances / self$period_length) ** 2
                / self$periodic_length_scale ** 2
            ))
        },

        #~ @description Generate a squared exponential (SE) temporal kernel matrix
        #~ @return A matrix containing the generated SE kernel
        se_kernel_gen = function() {
            return(torch::torch_exp(
                - self$time_distances ** 2 / (2 * self$decay_time_scale ** 2)
            ))
        },

        #~ @description Generate a squared exponential (SE) temporal kernel matrix
        #~ @return A matrix containing the generated SE kernel
        periodic_se_kernel_gen = function() {
            return(private$periodic_kernel_gen() * private$se_kernel_gen())
        }
    )
)

#' @title R6 class that generates Spatial Kernels
#'
#' @description A SpatialKernelGenerator can create spatial kernels according
#' to provided parameters.
#'
#' @export
SpatialKernelGenerator <- R6::R6Class(
    'SpatiallKernelGenerator',

    public = list(
        distance_matrix = NULL,
        smoothness_factor = NULL, # 2*nu (SKLearn) == ddd (Matlab code)
        kernel_variance = NULL,
        kernel = NULL,
        core_kernel_fn = NULL,
        # Following param is set externally via HyperParamSampler
        spatial_length_scale = NULL,

        #' @description Create a new \code{SpatialKernelGenerator} object.
        #' @param distance_matrix Matrix[numeric]: The distance matrix between each spatial entity
        #' @param smoothness_factor Integer (1,3,5): Factor deciding the matern kernel's type (3 for 3/2)
        #' @param kernel_variance Numeric: The variance used in the kernel
        #' @param ini_spatial_length_scale Numeric: Initial value for the spatial length scale
        #' @return A new \code{SpatialKernelGenerator} object.
        initialize = function(
            distance_matrix,
            smoothness_factor,
            kernel_variance,
            ini_spatial_length_scale = NULL
        ) {
            self$distance_matrix <- distance_matrix
            self$smoothness_factor <- smoothness_factor
            self$kernel_variance <- kernel_variance
            self$spatial_length_scale <- ini_spatial_length_scale
            self$core_kernel_fn <- private$get_core_kernel_fn(smoothness_factor)
        },

        #' @description Generate and set spatial kernel values using current parameters
        #' @return A new numeric matrix containing the spatial kernel
        kernel_gen = function() {
            kernel_results <- self$distance_matrix * sqrt(self$smoothness_factor) / self$spatial_length_scale
            self$kernel <- self$kernel_variance * self$core_kernel_fn(kernel_results) * (-kernel_results)$exp()
            return(self$kernel)
        }
    ),

    private = list(
        #~ @description Get a valid matern kernel function according to a smoothness factor
        #~ @param smoothness_factor Integer (1,3,5): Factor deciding the matern kernel's type (3 for 3/2)
        #~ @return A core kernel matern function
        # TODO missing smoothness factor 7
        get_core_kernel_fn = function(smoothness_factor) {
            if (smoothness_factor == 1)  {
                return(function(t) return(1))
            } else if (smoothness_factor == 3) {
                return(function(t) return(1 + t))
            } else if (smoothness_factor == 5) {
                return(function(t) return(1 + t * (1 + t / 3)))
            } else {
                torch:::value_error('Kernel function for this smoothness factor is not implemented')
            }
        }
    )
)
