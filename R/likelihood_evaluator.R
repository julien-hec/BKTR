#' @import torch
#' @importFrom R6 "R6Class"
#' @include tensor_ops.R

#' @title R6 class to evaluate the marginal likelihood of the hyperparameter
#'
#' @description MarginalLikelihoodEvaluator enable the calculation of the marginal
#' likelihood of the kernel hyperparameters. This likelihood is used during the sampling
#' process.
#'
#' @noRd
MarginalLikelihoodEvaluator <- R6::R6Class(
    'MarginalLikelihoodEvaluator',
    public = list(
        axis_permutation = c(),
        rank_decomp = NULL,
        nb_covariates = NULL,
        covariates = NULL,
        omega = NULL,
        y_masked = NULL,

        inv_k = NULL,
        chol_k = NULL,
        chol_lu = NULL,
        uu = NULL,
        likelihood = NULL,


        initialize = function(rank_decomp, nb_covariates, covariates, omega, y, is_transposed) {
            self$rank_decomp <- rank_decomp
            self$nb_covariates <- nb_covariates
            self$covariates <- covariates
            self$omega <- omega
            self$axis_permutation <- if (is_transposed) c(2, 1) else c(1, 2)
            self$y_masked <- y * omega
        },

        calc_likelihood = function(kernel_values, decomp_values, covs_decomp, tau) {
            rank_decomp <- self$rank_decomp
            kernel_size <- kernel_values$shape[1]
            lambda_size <- kernel_size * self$rank_decomp

            psi_u <- torch::torch_einsum("ijk,jkl->ilj", c(
                self$covariates$permute(c(self$axis_permutation, 3)),
                TSR$khatri_rao_prod(decomp_values, covs_decomp)$reshape(c(-1, self$nb_covariates, rank_decomp))
            ))
            psi_u_mask <- psi_u * self$omega$permute(c(self$axis_permutation))$unsqueeze(2)$expand_as(psi_u)

            self$chol_k <- torch::linalg_cholesky(kernel_values)
            kernel_inverse <- torch::linalg_solve(
                self$chol_k$t(), torch::linalg_solve(self$chol_k, TSR$eye(kernel_size))
            )
            stabilized_kernel_inv <- (kernel_inverse$t() + kernel_inverse) / 2
            self$inv_k <- TSR$kronecker_prod(
                TSR$eye(rank_decomp),
                stabilized_kernel_inv
            ) # I_R Kron inv(Ks)

            lambda_u <- tau * torch::torch_einsum('ijk,ilk->ijl', c(psi_u_mask, psi_u_mask)) # tau * H_T * H_T'
            lambda_u <- (
                lambda_u$transpose(1, -1)$unsqueeze(-1) * TSR$eye(kernel_size)
            )$transpose(2, 3)$reshape(c(lambda_size, lambda_size))
            lambda_u <- lambda_u + self$inv_k
            self$chol_lu <- torch::linalg_cholesky(lambda_u)
            uu <- torch:::torch_linalg_solve_triangular(
                self$chol_lu,
                torch::torch_einsum(
                    'ijk,ik->ji', c(psi_u_mask, self$y_masked$permute(c(self$axis_permutation)))
                )$flatten()$unsqueeze(2),
                upper = FALSE
            )$squeeze()
            self$likelihood <- as.numeric((
                TSR$tensor(0.5 * tau ** 2) * uu$t()$matmul(uu)
                - self$chol_lu$diag()$log()$sum()
                - TSR$tensor(rank_decomp) * self$chol_k$diag()$log()$sum()
            )$cpu())
            self$uu <- uu
            return(self$likelihood)
        }
    )
)
