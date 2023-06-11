#' @importFrom R6 "R6Class"
#' @import torch


#' @title R6 class for kernel's hyperparameter sampling
#'
#' @description The KernelParamSampler encapsulate all the behavior related to
#' the sampling of the kernel hyperparameters
#'
#' @export
#' @keywords internal
KernelParamSampler <- R6::R6Class(
    'KernelParamSampler',
    public = list(
        kernel = NULL,
        marginal_ll_eval_fn = NULL,

        initialize = function(
            kernel,
            marginal_ll_eval_fn
        ) {
            self$kernel <- kernel
            self$marginal_ll_eval_fn <- marginal_ll_eval_fn
        },

        set_theta_value = function(theta) {
            self$theta_value <- theta
            self$kernel_generator[[self$kernel_hparam_name]] <- exp(theta)
        },

        initialize_theta_bounds = function(param) {
            theta_range <- param$slice_sampling_scale * as.numeric(tsr$rand(1))
            theta_min <- max(log(param$value) - theta_range, log(param$lower_bound))
            theta_max <- min(theta_min + param$slice_sampling_scale, log(param$upper_bound))
            return(list(min = theta_min, max = theta_max))
        },

        prior_fn = function(param) {
            return(-0.5 * param$hparam_precision * log(param$value) ** 2)
        },

        sample_rand_theta_value = function(theta_min, theta_max) {
            return(theta_min + (theta_max - theta_min) * as.numeric(tsr$rand(1)))
        },

        sample_param = function(param) {
            theta_bounds <- self$initialize_theta_bounds(param)
            theta_min <- theta_bounds$min
            theta_max <- theta_bounds$max
            initial_theta <- log(param$value)
            self$kernel$kernel_gen()
            initial_marginal_likelihood <- self$marginal_ll_eval_fn() + self$prior_fn(param)
            density_threshold <- as.numeric(tsr$rand(1))

            while (TRUE) {
                new_theta <- self$sample_rand_theta_value(theta_min, theta_max)
                param$value <- exp(new_theta)
                self$kernel$kernel_gen()
                new_marginal_likelihood <- self$marginal_ll_eval_fn() + self$prior_fn(param)

                marg_ll_diff <- new_marginal_likelihood - initial_marginal_likelihood
                if (exp(marg_ll_diff) > density_threshold) {
                    return(param$value)
                }
                if (new_theta < initial_theta) {
                    theta_min <- new_theta
                } else {
                    theta_max <- new_theta
                }
            }
        },

        sample = function() {
            for (param in self$kernel$parameters) {
                if (!param$is_fixed) {
                    self$sample_param(param)
                }
            }
        }
    )
)

#' @title Sample a tensor of random values from a normal multivariate distribution
#'
#' @description The sampling use a tensor of mean and the upper triangular portion of the precision matrix
#'
#' @export
#' @keywords internal
sample_norm_multivariate <- function(mean_vec, precision_upper_tri) {
    # TODO Open PR & Issue for https://github.com/mlverse/torch/blob/main/R/distributions-multivariate_normal.R L:86
    # Not Able to use the precision matrix because of priority of ops (!is.null(NULL) + !is.null(1) + !is.null(1)) == F
    # ERROR comes from torch::distr_multivariate_normal(torch::torch_zeros(2), precision_matrix = torch::torch_eye(2))
    return(
        torch::torch_triangular_solve(
            tsr$new_tensor(torch::torch_randn_like(mean_vec))$unsqueeze(2),
            precision_upper_tri,
            upper = TRUE
        )[[1]]$squeeze() + mean_vec
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
                torch::distr_gamma(self$a_tau$cpu(), b_tau)$sample()
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
            wish_precision_matrix <- rWishart(1, self$wish_df, wish_sigma)[, , 1]
            self$wish_precision_tensor <- tsr$new_tensor(wish_precision_matrix)
            return(self$wish_precision_tensor)
        }
    )
)
