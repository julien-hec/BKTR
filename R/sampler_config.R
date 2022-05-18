#' @importFrom R6 R6Class

#' @title R6 class representing the configuration used for the sampling of a hyperparameter
#'
#' @description A KernelSamplerConfig object contains all the information necessary to the process
#' of sampling kernel hyperparameters. (Used in \strong{Algorithm 2} of the paper)
#'
#' @export
#' @keywords internal
KernelSamplerConfig <- R6::R6Class(
    public = list(
        #' @field slice_sampling_scale The amplitude of the range used for sampling (Paper: \eqn{\rho})
        slice_sampling_scale = NULL,
        #' @field min_hyper_value The minimal value admissible for the hyperparameter (Paper: \eqn{\phi_{min}})
        min_hyper_value = NULL,
        #' @field max_hyper_value The maximal value admissible for the hyperparameter (Paper: \eqn{\phi_{max}})
        max_hyper_value = NULL,
        #' @field hyper_mu_prior The prior value for the mean of the sampled hyperparameter (Paper: \eqn{\phi})
        hyper_mu_prior = NULL,
        #' @field hyper_precision_prior The prior value attributed to the precision of the sampled hyperparameter
        hyper_precision_prior = NULL,

        #' @description Create a new \code{KernelSamplerConfig} object.
        #' @param slice_sampling_scale Numeric: Slice sampling scale for the hyperparameter
        #' @param min_hyper_value Numeric: Minimum value of the hyperparameter
        #' @param max_hyper_value Numeric: Maximum value of the hyperparameter
        #' @param hyper_mu_prior Numeric: Prior mean of the hyperparameter
        #' @param hyper_precision_prior Numeric: Prior precision of the hyperparameter
        #' @return A new \code{KernelSamplerConfig} object.
        initialize = function(
            slice_sampling_scale = log(10),
            min_hyper_value = log(1E-3),
            max_hyper_value = log(1E3),
            hyper_mu_prior = 0,
            hyper_precision_prior = 1
        ) {
            self$slice_sampling_scale <- slice_sampling_scale
            self$min_hyper_value <- min_hyper_value
            self$max_hyper_value <- max_hyper_value
            self$hyper_mu_prior <- hyper_mu_prior
            self$hyper_precision_prior <- hyper_precision_prior
        }
    )
)
