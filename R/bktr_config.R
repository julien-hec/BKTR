#' @importFrom R6 R6Class
#' @importFrom torch torch_float64
#' @include sampler_config.R

#' @title R6 class that contains the configuration for a BKTR process
#'
#' @description A BKTRConfig object contains various information that are central to the
#' BKTR sampling process. It is notably used to configure the tensor backend, initialize values
#' for the algorithm (e,g. hyperparameters) or to define key algorithm parameters (e,g. number
#' of iterations or decomposition rank).
#'
#' @export
BKTRConfig <- R6::R6Class(
    public = list(
        # BKTR Params
        #' @field rank_decomp Rank of the CP decomposition (Paper: \eqn{R})
        rank_decomp = NULL,
        #' @field max_iter Maximum number of iteration that the MCMC sampling will use
        #' (Paper: \eqn{K_1 + K_2})
        max_iter = NULL,
        #' @field burn_in_iter Number of iteration before the sampling starts (Paper: \eqn{K_1})
        burn_in_iter = NULL,

        # Kernel Params
        #' @field temporal_period_length Period used in the periodic kernel (Paper: \eqn{T})
        temporal_period_length = NULL,
        #' @field temporal_kernel_fn_name Temporal kernel function used (periodic_se, se or periodic)
        temporal_kernel_fn_name = NULL,
        #' @field kernel_time_segment_duration Duration of one time segment in the kernel
        kernel_time_segment_duration = NULL,
        #' @field has_stabilizing_diag Indicate if we add a stabilizing diagonal in the temporal kernel
        has_stabilizing_diag = NULL,
        #' @field spatial_smoothness_factor Smoothness factor used in the Matern kernel
        spatial_smoothness_factor = NULL,
        #' @field kernel_variance Variance used for all kernels (Paper: \eqn{\sigma^2_s, \sigma^2_t})
        kernel_variance = NULL,
        #' @field sigma_r Variance of the white noise process TODO (Paper: \eqn{\tau^{-1}})
        sigma_r = NULL,
        #' @field a_0 Initial value for the shape (\eqn{\alpha}) in the gamma function generating tau
        a_0 = NULL,
        #' @field b_0 Initial value for the rate (\eqn{\beta}) in the gamma function generating tau
        b_0 = NULL,

        # Periodic length scale Params
        #' @field periodic_scale_config Periodic length scale's KernelSamplerConfig (Paper: \eqn{\gamma_1})
        periodic_scale_config = NULL,
        #' @field decay_scale_config Decay time scale's KernelSamplerConfig (Paper: \eqn{\gamma_2})
        decay_scale_config = NULL,
        #' @field spatial_length_config spatial length-scale's KernelSamplerConfig (Paper: \eqn{\phi})
        spatial_length_config = NULL,

        # Ouput Params
        #' @field sampled_beta_indexes Indexes of beta estimates that need to be sampled through iterations
        sampled_beta_indexes = NULL,
        #' @field sampled_y_indexes Indexes of y estimates that need to be sampled through iterations
        sampled_y_indexes = NULL,
        #' @field results_export_dir Path of the folder where the csv file will be exported
        results_export_dir = NULL,

        # Torch Params
        #' @field torch_seed Seed used by torch to be able to reproduce results for a given seed
        torch_seed = NULL,
        #' @field torch_dtype Type used for floating points in the tensor backend
        torch_dtype = NULL,
        #' @field torch_device Device used by the tensor backend for calculation (\code{'cuda'} or \code{'cpu'})
        torch_device = NULL,

        #' @description Create a new \code{BKTRConfig} object.
        #' @param rank_decomp Integer: CP decomposition's Rank
        #' @param burn_in_iter Integer: Iterations used for burn in during MCMC
        #' @param max_iter Integer: Total number of iterations for MCMC
        #' @param temporal_period_length Integer: Period for periodic kernel (\code{7} represent week's length)
        #' @param temporal_kernel_fn_name String: Temporal kernel function used (periodic_se, se or periodic)
        #' @param kernel_time_segment_duration Numeric: Duration of one time segment in the kernel
        #' @param has_stabilizing_diag Boolean: Indicate if we add a stabilizing diagonal in the temporal kernel
        #' @param spatial_smoothness_factor Integer (1,3,5): Matern kernel's smoothness (3 for Matern 3/2)
        #' @param kernel_variance Numeric: Variance used for all kernels
        #' @param sigma_r Numeric: Initial value of white noise's variance (\eqn{\tau^{-1}})
        #' @param a_0 Numeric: Initial value for the shape (\eqn{\alpha}) in the gamma function generating tau
        #' @param b_0 Numeric: Initial value for the rate (\eqn{\beta}) in the gamma function generating tau
        #' @param period_slice_sampling_scale Numeric: Slice sampling scale of the periodic length scale
        #' @param period_min_hparam_val Numeric: Minimum value for periodic length scale hyperparameter
        #' @param period_max_hparam_val Numeric: Maximum value for periodic length scale hyperparameter
        #' @param period_hparam_mu_prior Numeric: Initial value for periodic length scale mean
        #' @param period_hparam_precision_prior Numeric: Initial value for periodic length scale precision
        #' @param decay_slice_sampling_scale Numeric: Slice sampling scale for decay length scale
        #' @param decay_min_hparam_val Numeric: Minimum value for decay length scale hyperparameter
        #' @param decay_max_hparam_val Numeric: Maximum value for decay length scale hyperparameter
        #' @param decay_hparam_mu_prior Numeric: Initial value for decay length scale mean
        #' @param decay_hparam_precision_prior Numeric: Initial value for decay length scale precision
        #' @param spatial_slice_sampling_scale Numeric: Slice sampling scale for spatial length scale
        #' @param spatial_min_hparam_val Numeric: Minimum value for spatial length scale hyperparameter
        #' @param spatial_max_hparam_val Numeric: Maximum value for spatial length scale hyperparameter
        #' @param spatial_hparam_mu_prior Numeric: Initial value for spatial length scale mean
        #' @param spatial_hparam_precision_prior Numeric: Initial value for spatial length scale precision
        #' @param sampled_beta_indexes Indexes of beta estimates that need to be sampled through iterations
        #' @param sampled_y_indexes Indexes of y estimates that need to be sampled through iterations
        #' @param results_export_dir Path of the folder where the csv file will be exported
        #' @param torch_seed Seed used by torch to be able to reproduce results for a given seed
        #' @param torch_device String ('cpu', 'cuda'): the device used by Torch for computation
        #' @param torch_dtype Torch floating point type (\code{torch_float16()}, \code{torch_float32()},
        #' \code{torch_float64()}): Default \code{torch_float64()} helps for large kernel cholesky decomposition
        #' @return A new \code{BKTRConfig} object.
        initialize = function(
            # BKTR Params
            rank_decomp = NULL,
            burn_in_iter = NULL,
            max_iter = NULL,
            # Kernel Params
            temporal_period_length = 7,
            temporal_kernel_fn_name = 'periodic_se',
            kernel_time_segment_duration = 1,
            has_stabilizing_diag = FALSE,
            spatial_smoothness_factor = 3,
            kernel_variance = 1,
            # Tau params
            sigma_r = 1e-2,
            a_0 = 1e-6,
            b_0 = 1e-6,
            # Periodic length scale Params
            period_slice_sampling_scale = log(10),
            period_min_hparam_val = log(1E-3),
            period_max_hparam_val = log(1E3),
            period_hparam_mu_prior = 0,
            period_hparam_precision_prior = 1,
            # Decay time scale Params
            decay_slice_sampling_scale = log(10),
            decay_min_hparam_val = log(1E-3),
            decay_max_hparam_val = log(1E3),
            decay_hparam_mu_prior = 0,
            decay_hparam_precision_prior = 1,
            # Spatial length scale Params
            spatial_slice_sampling_scale = log(10),
            spatial_min_hparam_val = log(1E-3),
            spatial_max_hparam_val = log(1E3),
            spatial_hparam_mu_prior = 0,
            spatial_hparam_precision_prior = 1,
            # Output Params
            sampled_beta_indexes = c(),
            sampled_y_indexes = c(),
            results_export_dir = '.',
            # Torch Params
            torch_seed = NULL,
            torch_device = 'cpu',
            torch_dtype = torch::torch_float64()
        ) {
            if (burn_in_iter >= max_iter) {
                torch:::value_error('burn_in iterations cannot be greater than the total iterations')
            }

            self$rank_decomp <- rank_decomp
            self$max_iter <- max_iter
            self$burn_in_iter <- burn_in_iter
            self$temporal_period_length <- temporal_period_length
            self$spatial_smoothness_factor <- spatial_smoothness_factor
            self$kernel_variance <- kernel_variance
            self$temporal_kernel_fn_name <- temporal_kernel_fn_name
            self$kernel_time_segment_duration <- kernel_time_segment_duration
            self$has_stabilizing_diag <- has_stabilizing_diag
            self$sigma_r <- sigma_r
            self$a_0 <- a_0
            self$b_0 <- b_0
            self$sampled_beta_indexes <- sampled_beta_indexes
            self$sampled_y_indexes <- sampled_y_indexes
            self$torch_seed <- torch_seed
            self$torch_dtype <- torch_dtype
            self$torch_device <- torch_device
            self$results_export_dir <- results_export_dir

            self$periodic_scale_config <- KernelSamplerConfig$new(
                slice_sampling_scale = period_slice_sampling_scale,
                min_hyper_value = period_min_hparam_val,
                max_hyper_value = period_max_hparam_val,
                hyper_mu_prior = period_hparam_mu_prior,
                hyper_precision_prior = period_hparam_precision_prior
            )

            self$decay_scale_config <- KernelSamplerConfig$new(
                slice_sampling_scale = decay_slice_sampling_scale,
                min_hyper_value = decay_min_hparam_val,
                max_hyper_value = decay_max_hparam_val,
                hyper_mu_prior = decay_hparam_mu_prior,
                hyper_precision_prior = decay_hparam_precision_prior
            )

            self$spatial_length_config <- KernelSamplerConfig$new(
                slice_sampling_scale = spatial_slice_sampling_scale,
                min_hyper_value = spatial_min_hparam_val,
                max_hyper_value = spatial_max_hparam_val,
                hyper_mu_prior = spatial_hparam_mu_prior,
                hyper_precision_prior = spatial_hparam_precision_prior
            )
        }
    )
)
