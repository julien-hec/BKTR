#' @importFrom R6 "R6Class"
#' @import torch


#' @title R6 class for kernel's hyperparameter sampling
#'
#' @description The HyperParamSampler encapsulate all the behavior related to
#' the sampling of the kernel hyperparameters
#'
#' @export
#' @keywords internal
HyperParamSampler <- R6::R6Class(
    'HyperParamSampler',
    public = list(
        # Coming from config
        slice_sampling_scale = NULL,
        min_hyper_value = NULL,
        max_hyper_value = NULL,
        hyper_mu_prior = NULL,
        hyper_precision_prior = NULL,
        # Passed params
        kernel_generator = NULL,
        marginal_ll_eval_fn = NULL,
        kernel_hparam_name = NULL, # Used for reference to the param's name in the kernel generator

        theta_value = NULL,
        theta_min = NULL,
        theta_max = NULL,

        initialize = function(
            config,
            kernel_generator,
            marginal_ll_eval_fn,
            kernel_hparam_name
        ) {
            self$hyper_mu_prior <- config$hyper_mu_prior
            self$hyper_precision_prior <- config$hyper_precision_prior
            self$slice_sampling_scale <- config$slice_sampling_scale
            self$min_hyper_value <- config$min_hyper_value
            self$max_hyper_value <- config$max_hyper_value

            self$kernel_generator <- kernel_generator
            self$marginal_ll_eval_fn <- marginal_ll_eval_fn
            self$kernel_hparam_name <- kernel_hparam_name
            self$set_theta_value(self$hyper_mu_prior)
        },

        set_theta_value = function(theta) {
            self$theta_value <- theta
            self$kernel_generator[[self$kernel_hparam_name]] <- exp(theta)
        },

        initialize_theta_bounds = function() {
            theta_range <- self$slice_sampling_scale * runif(1)
            self$theta_min <- max(self$theta_value - theta_range, self$min_hyper_value)
            self$theta_max <- min(self$theta_min + self$slice_sampling_scale, self$max_hyper_value)
        },

        prior_fn = function(theta) {
            return(-0.5 * self$hyper_precision_prior * (theta - self$hyper_mu_prior) ** 2)
        },

        sample_rand_theta_value = function() {
            return(self$theta_min + (self$theta_max - self$theta_min) * runif(1))
        },

        sample = function() {
            initial_theta <- self$theta_value
            self$initialize_theta_bounds()

            self$kernel_generator$kernel_gen()
            initial_marginal_likelihood <- self$marginal_ll_eval_fn() + self$prior_fn(self$theta_value)

            density_threshold <- runif(1)

            while (TRUE) {
                new_theta <- self$sample_rand_theta_value()
                self$set_theta_value(new_theta)
                self$kernel_generator$kernel_gen()
                # Rarely the generated kernel is not a positive definite matrix (so we resample theta)
                # TODO find how we can circumvent those non PD matrix generation
                new_marginal_likelihood <- self$marginal_ll_eval_fn() + self$prior_fn(new_theta)

                if (exp(new_marginal_likelihood - initial_marginal_likelihood) > density_threshold) {
                    return(exp(new_theta))
                }
                if (new_theta < initial_theta) {
                    self$theta_min <- new_theta;
                }
                else {
                    self$theta_max <- new_theta;
                }
            }
        }
    )
)

#' @title Sample a tensor of random values from a normal multivariate distribution
#'
#' @description The sampling use a tensor of 
#'
#' @export
#' @keywords internal
sample_norm_multivariate <- function(mean_vec, precision_lower_tri) {
    # TODO Open PR & Issue for https://github.com/mlverse/torch/blob/main/R/distributions-multivariate_normal.R L:86
    # Not Able to use the precision matrix because of priority of ops (!is.null(NULL) + !is.null(1) + !is.null(1)) == F
    # ERROR comes from torch::distr_multivariate_normal(torch::torch_zeros(2), precision_matrix = torch::torch_eye(2))
    return(
        torch::linalg_solve(
            precision_lower_tri, tsr$new_tensor(torch::torch_randn_like(mean_vec))
        ) + mean_vec
    )
}

#' @export
#' @keywords internal
get_cov_decomp_chol <- function(
    spatial_decomp, time_decomp, covs, rank_cp, omega, tau, y, wish_precision_tensor
) {
    y_masked <- omega * y
    # TODO Merge some parts with marginal ll of spatial and temporal
    # get corresponding norm multivariate mean
    b <- tsr$khatri_rao_prod(spatial_decomp, time_decomp)$reshape(
        c(spatial_decomp$shape[1], time_decomp$shape[1], rank_cp)
    )
    psi_c <- torch::torch_einsum('ijk,ijl->ijlk', c(covs, b))
    psi_c_mask <- psi_c * omega$unsqueeze(3)$unsqueeze(4)$expand_as(psi_c)
    psi_c_mask <- psi_c_mask$permute(c(2, 1, 3, 4))$reshape(
        c(psi_c$shape[1] * psi_c$shape[2], psi_c$shape[3] * psi_c$shape[4])
    )
    inv_s <- tsr$kronecker_prod(tsr$eye(rank_cp), wish_precision_tensor)
    lambda_c <- tau * psi_c_mask$t()$matmul(psi_c_mask) + inv_s
    chol_lc <- torch::linalg_cholesky(lambda_c)
    cc <- torch::linalg_solve(chol_lc, psi_c_mask$t()$matmul(y_masked$t()$flatten()))
    return(list(chol_lc = chol_lc, cc = cc))
}

#' @title R6 class for the Tau precision hyperparameter sampling
#'
#' @description Encapsulate all the behavior that allows to generate new tau values
#'
#' @export
#' @keywords internal
TauSampler <- R6::R6Class(
    'TauSampler',
    public = list(
        b_0 = NULL,
        a_tau = NULL,

        initialize = function(a_0, b_0, nb_observations) {
            self$b_0 <- b_0
            self$a_tau <- tsr$new_tensor(a_0 + 0.5 * nb_observations)
        },

        sample = function(total_sq_error) {
            b_tau <- self$b_0 + 0.5 * total_sq_error
            return(tsr$new_tensor(
                torch::distr_gamma(self$a_tau$cpu(), b_tau$cpu())$sample()
            ))
        }
    )
)

#' @title R6 class to sample new precision matrices
#'
#' @description Encapsulate all the behavior that allows to sample new precision matrices from
#' a Wishart distribution
#'
#' @export
#' @keywords internal
# TODO create a PR to add rand wishart in R Torch
PrecisionMatrixSampler <- R6::R6Class(
    'PrecisionMatrixSampler',
    public = list(
        nb_covariates = NULL,
        wish_df = NULL,
        wish_precision_tensor = NULL,

        initialize = function(nb_covariates, rank_cp) {
            self$nb_covariates <- nb_covariates
            self$wish_df <- nb_covariates + rank_cp
        },

        sample = function(covs_decomp) {
            w <- covs_decomp$matmul(covs_decomp$t()) + tsr$eye(self$nb_covariates)
            wish_sigma <- as.matrix(((w + w$t()) * 0.5)$inverse()$cpu())
            wish_precision_matrix <- rWishart(1, self$wish_df, wish_sigma) # Lambda C
            self$wish_precision_tensor <- tsr$new_tensor(wish_precision_matrix[, , 1])
            return(self$wish_precision_tensor)
        }
    )
)
