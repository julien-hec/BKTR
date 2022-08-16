#' @import torch
#' @importFrom R6 "R6Class"
#' @include tensor_ops.R

#' @title R6 class to evaluate the marginal likelihood of the hyperparameter
#'
#' @description MarginalLikelihoodEvaluator enable the calculation of the marginal
#' likelihood of the kernel hyperparameters. This likelihood is used during the sampling
#' process.
#'
#' @export
#' @keywords internal
MarginalLikelihoodEvaluator <- R6::R6Class(
    'MarginalLikelihoodEvaluator',
    public = list(
        permuted_axis = NULL,
        rank_decomp = NULL,
        nb_covariates = NULL,
        covariates = NULL,
        omega = NULL,
        y_masked = NULL,

        inv_k = NULL,
        chol_k = NULL,
        chol_lu = NULL,
        # TODO What is uu
        uu = NULL,
        likelihood = NULL,


        initialize = function(rank_decomp, nb_covariates, covariates, omega, y, is_transposed) {
            self$rank_decomp <- rank_decomp
            self$nb_covariates <- nb_covariates
            self$covariates <- covariates
            self$omega <- omega
            self$permuted_axis <- if (is_transposed) c(2, 1) else c(1, 2)
            self$y_masked <- y * omega
        },

        calc_likelihood = function(kernel_values, decomp_values, covs_decomp, tau) {
            rank_decomp <- self$rank_decomp
            kernel_size <- kernel_values$shape[1]
            lambda_size <- kernel_size * self$rank_decomp

            psi_u <- torch::torch_einsum("ijk,jkl->ilj", c(
                self$covariates$permute(c(self$permuted_axis, 3)),
                tsr$khatri_rao_prod(decomp_values, covs_decomp)$reshape(c(-1, self$nb_covariates, rank_decomp))
            ))
            psi_u_mask <- psi_u * self$omega$permute(c(self$permuted_axis))$unsqueeze(2)$expand_as(psi_u)
            self$chol_k <- torch::linalg_cholesky(kernel_values)
            self$inv_k <- tsr$kronecker_prod(
                tsr$eye(rank_decomp),
                (kernel_values)$inverse()
            ) # I_R Kron inv(Ks)
            lambda_u <- tau * torch::torch_einsum('ijk,ilk->ijl', c(psi_u_mask, psi_u_mask)) # tau * H_T * H_T'
            # TODO Check the broadcasting here, for sure we can simplify it (We just want to diagonalize lambda_u)
            lambda_u <- lambda_u$permute(c(2, 1, 3))$flatten(start_dim = 1, end_dim = 2)$unsqueeze(2)$expand(
                c(-1, kernel_size, -1))$permute(c(1, 3, 2))$reshape(c(lambda_size, lambda_size))
            lambda_u <- tsr$kronecker_prod(
                tsr$ones(c(rank_decomp, rank_decomp)),
                tsr$eye(kernel_size)
            ) * lambda_u
            lambda_u <- lambda_u + self$inv_k
            self$chol_lu <- torch::linalg_cholesky(lambda_u)
            uu <- torch::torch_triangular_solve(
                torch::torch_einsum(
                    'ijk,ik->ji', c(psi_u_mask, self$y_masked$permute(c(self$permuted_axis)))
                )$flatten()$unsqueeze(2),
                self$chol_lu,
                upper = FALSE
            )[[1]]$squeeze()
            self$likelihood <- as.numeric((
                tsr$new_tensor(0.5 * tau ** 2) * uu$t()$matmul(uu)
                - self$chol_lu$diag()$log()$sum()
                - tsr$new_tensor(rank_decomp) * self$chol_k$diag()$log()$sum()
            )$cpu())
            self$uu <- uu
            return(self$likelihood)
        }
    )
)
