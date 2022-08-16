#' @import torch
#' @include kernel_generator.R
#' @include marginal_likelihoods.R
#' @include result_logger.R
#' @include sampler.R
#' @include tensor_ops.R
#' @importFrom R6 R6Class

#' @title R6 class encapsulating the BKTR regression steps
#'
#' @description A BKTRRegressor holds all the key elements to accomplish the MCMC sampling
#' algorithm (\strong{Algorithm 1} of the paper).
#'
#' @export
BKTRRegressor <- R6::R6Class(
    public = list(
        config = NULL,
        spatial_distance_tensor = NULL,
        y = NULL,
        omega = NULL,
        covariates = NULL,
        covariates_dim = NULL,
        logged_params_tensor = NULL,
        tau = NULL,
        # Covariate decompositions (change during iter)
        spatial_decomp = NULL, # U
        temporal_decomp = NULL, # V
        covs_decomp = NULL, # C or W
        # Result Logger
        result_logger = NULL,
        # Kernel factories
        temporal_kernel_factory = NULL,
        spatial_kernel_factory = NULL,
        # Samplers
        spatial_length_sampler = NULL,
        decay_scale_sampler = NULL,
        periodic_length_sampler = NULL,
        tau_sampler = NULL,
        precision_matrix_sampler = NULL,
        # Likelihood evaluators
        spatial_ll_evaluator = NULL,
        temporal_ll_evaluator = NULL,

        #' @description Create a new \code{BKTRRegressor} object.
        #' @param bktr_config BKTRConfig: Configuration used for the BKTR regression
        #' @param temporal_covariate_matrix matrix[numeric]: Temporal Covariates
        #' @param spatial_covariate_matrix matrix[numeric]: Spatial Covariates
        #' @param spatial_distance_matrix matrix[numeric]: Distance between spatial entities
        #' @param y matrix[numeric]: Response variable that we are trying to predict
        #' @param omega matrix[boolean]: Mask showing if a y observation is missing or not
        #' @return A new \code{BKTRRegressor} object.
        initialize = function(
            bktr_config,
            temporal_covariate_matrix,
            spatial_covariate_matrix,
            spatial_distance_matrix,
            y,
            omega
        ) {
            self$config <- bktr_config
            # Set tensor backend according to config
            tsr$set_device_type(self$config$torch_device)
            tsr$set_tensor_type(self$config$torch_dtype)
            # Assignation
            self$spatial_distance_tensor <- tsr$new_tensor(spatial_distance_matrix)
            self$y <- tsr$new_tensor(y)
            self$omega <- tsr$new_tensor(omega)
            self$tau <- 1 / self$config$sigma_r
            private$reshape_covariates(
                tsr$new_tensor(spatial_covariate_matrix),
                tsr$new_tensor(temporal_covariate_matrix)
            )
            private$initialize_params()
        },

        #' @description Launch the MCMC sampling process. \cr
        #' For a predefined number of iterations:
        #' \enumerate{
        #' \item{Sample the spatial length scale hyperparameter}
        #' \item{Sample the decay time scale hyperparameter}
        #' \item{Sample the periodic length scale hyperparameter}
        #' \item{Sample the precision matrix from a wishart distribution}
        #' \item{Sample a new spatial covariate decomposition}
        #' \item{Sample a new covariate decomposition}
        #' \item{Sample a new temporal covariate decomposition}
        #' \item{Calculate respective errors for the iterations}
        #' \item{Sample a new tau value}
        #' \item{Collect all the important data for the iteration}
        #' }
        #' @return A list containing the results of the MCMC sampling.
        mcmc_sampling = function() {
            for (i in 1:self$config$max_iter) {
                print(sprintf('*** Running iter %i ***', i))
                private$sample_kernel_hparam()
                private$sample_precision_wish()
                private$sample_spatial_decomp()
                private$sample_covariate_decomp()
                private$sample_temporal_decomp()
                private$set_errors_and_sample_precision_tau()
                private$collect_iter_values(i)

                # Needed to force garbage collection in Cuda
                gc()
            }
            avg_estimates <- private$calculate_avg_estimates()
            private$log_iter_results()
            return(avg_estimates)
        }
    ),

    private = list(
        #~ @description Reshape the covariate tensors into one single tensor and set this
        #~ tensor into the \code{covariates} property of the BKTRRegressor object
        #~ @param temporal_covariate_tensor tensor: Temporal Covariates
        #~ @param spatial_covariate_tensor tensor: Spatial Covariates
        reshape_covariates = function(spatial_covariate_tensor, temporal_covariate_tensor) {
            time_cov_shape <- temporal_covariate_tensor$shape
            space_cov_shape <- spatial_covariate_tensor$shape
            covs_dim <- list(
                nb_spaces = space_cov_shape[1], # S
                nb_times = time_cov_shape[1], # T
                nb_spatial_covariates = space_cov_shape[2],
                nb_temporal_covariates = time_cov_shape[2],
                nb_covariates = 1 + space_cov_shape[2] + time_cov_shape[2] # P
            )

            intersect_covs <- tsr$ones(c(covs_dim$nb_spaces, covs_dim$nb_times, 1))
            spatial_covs <- spatial_covariate_tensor$unsqueeze(2)$expand(
                c(covs_dim$nb_spaces, covs_dim$nb_times, covs_dim$nb_spatial_covariates)
            )
            time_covs <- temporal_covariate_tensor$unsqueeze(1)$expand(
                c(covs_dim$nb_spaces, covs_dim$nb_times, covs_dim$nb_temporal_covariates)
            )

            self$covariates_dim <- covs_dim
            self$covariates <- torch::torch_dstack(c(intersect_covs, spatial_covs, time_covs))
        },

        #~ @description Initialize the CP decomposed covariate tensors
        #~ using normally distributed random values
        init_covariate_decomp = function() {
            rank_decomp <- self$config$rank_decomp
            covs_dim <- self$covariates_dim

            self$spatial_decomp <- tsr$new_tensor(
                torch::torch_randn(c(covs_dim$nb_spaces, rank_decomp))
            )
            self$temporal_decomp <- tsr$new_tensor(
                torch::torch_randn(c(covs_dim$nb_times, rank_decomp))
            )
            self$covs_decomp <- tsr$new_tensor(
                torch::torch_randn(c(covs_dim$nb_covariates, rank_decomp))
            )
        },

        create_result_logger = function() {
            self$result_logger <- ResultLogger$new(
                y = self$y,
                omega = self$omega,
                covariates = self$covariates,
                nb_iter = self$config$max_iter,
                nb_burn_in_iter = self$config$burn_in_iter,
                export_dir = self$config$results_export_dir,
                seed = self$config$torch_seed,
                sampled_beta_indexes = self$config$sampled_beta_indexes,
                sampled_y_indexes = self$config$sampled_y_indexes
            )
        },

        #~ @description Create and set the kernel factories for the spatial and
        #~ temporal kernels
        create_kernel_factories = function() {
            self$spatial_kernel_factory <- SpatialKernelGenerator$new(
                self$spatial_distance_tensor,
                self$config$spatial_smoothness_factor,
                self$config$kernel_variance
            )
            self$temporal_kernel_factory <- TemporalKernelGenerator$new(
                'periodic_se',
                self$covariates_dim$nb_times,
                self$config$temporal_period_length,
                self$config$kernel_variance
            )
        },

        #~ @description Create and set the evaluators for the spatial and the temporal likelihoods
        create_likelihood_evaluators = function() {
            rank_decomp <- self$config$rank_decomp

            self$spatial_ll_evaluator <- MarginalLikelihoodEvaluator$new(
                rank_decomp, self$covariates_dim$nb_covariates, self$covariates,
                self$omega, self$y, is_transposed = FALSE
            )
            self$temporal_ll_evaluator <- MarginalLikelihoodEvaluator$new(
                rank_decomp, self$covariates_dim$nb_covariates, self$covariates,
                self$omega, self$y, is_transposed = TRUE
            )
        },

        #~ @description Create and set the hyperparameter samplers for the
        #~ spatial length scale, decay time scale, periodic length scale,
        #~ tau and the precision matrix
        create_hparam_samplers = function() {
            self$spatial_length_sampler <- HyperParamSampler$new(
                config = self$config$spatial_length_config,
                kernel_generator = self$spatial_kernel_factory,
                marginal_ll_eval_fn = private$calc_spatial_marginal_ll,
                kernel_hparam_name = 'spatial_length_scale'
            )
            self$decay_scale_sampler <- HyperParamSampler$new(
                config = self$config$decay_scale_config,
                kernel_generator = self$temporal_kernel_factory,
                marginal_ll_eval_fn = private$calc_temporal_marginal_ll,
                kernel_hparam_name = 'decay_time_scale'
            )
            self$periodic_length_sampler <- HyperParamSampler$new(
                config = self$config$periodic_scale_config,
                kernel_generator = self$temporal_kernel_factory,
                marginal_ll_eval_fn = private$calc_temporal_marginal_ll,
                kernel_hparam_name = 'periodic_length_scale'
            )

            self$tau_sampler <- TauSampler$new(
                self$config$a_0, self$config$b_0, self$omega$sum()
            )

            self$precision_matrix_sampler <- PrecisionMatrixSampler$new(
                self$covariates_dim$nb_covariates, self$config$rank_decomp
            )
        },

        #~ @description Calculate the spatial marginal likelihood
        calc_spatial_marginal_ll = function() {
            return(self$spatial_ll_evaluator$calc_likelihood(
                    self$spatial_kernel_factory$kernel, self$temporal_decomp, self$covs_decomp, self$tau
            ))
        },

        #~ @description Calculate the temporal marginal likelihood
        calc_temporal_marginal_ll = function() {
            return(self$temporal_ll_evaluator$calc_likelihood(
                self$temporal_kernel_factory$kernel, self$spatial_decomp, self$covs_decomp, self$tau
            ))
        },

        #~ @description Sample new kernel hyperparameters
        sample_kernel_hparam = function() {
            self$spatial_length_sampler$sample()
            self$decay_scale_sampler$sample()
            self$periodic_length_sampler$sample()
        },

        #~ @description Sample the precision matrix from a Wishart distribution
        sample_precision_wish = function() {
            self$precision_matrix_sampler$sample(self$covs_decomp)
        },

        #~ @description Sample a new covariate decomposition from a mulivariate normal distribution
        #~ @param initial_decomp tensor: Decomposition of the previous iteration
        #~ @param chol_l tensor: The cholesky decomposition of the #TODO
        #~ @param uu tensor: #TODO
        #~ @return A tensor containing the newly sampled covariate decomposition
        sample_decomp_norm = function(initial_decomp, chol_l, uu) {
            precision_mat <- chol_l$t()
            mean_vec <- self$tau * torch::linalg_solve(precision_mat, uu)
            return(
                sample_norm_multivariate(mean_vec, precision_mat)$reshape_as(initial_decomp$t())$t()
            )
        },

        #~ @description Sample a new spatial covariate decomposition
        sample_spatial_decomp = function() {
            ll_eval <- self$spatial_ll_evaluator
            self$spatial_decomp <- private$sample_decomp_norm(self$spatial_decomp, ll_eval$chol_lu, ll_eval$uu)
        },

        #~ @description Sample a new covariate decomposition
        sample_covariate_decomp = function() {
            chol_res <- get_cov_decomp_chol(
                self$spatial_decomp, self$temporal_decomp, self$covariates,
                self$config$rank_decomp, self$omega, self$tau, self$y,
                self$precision_matrix_sampler$wish_precision_tensor
            )
            self$covs_decomp <- private$sample_decomp_norm(self$covs_decomp, chol_res$chol_lc, chol_res$cc)
        },

        #~ @description Sample a new temporal covariate decomposition
        sample_temporal_decomp = function() {
            # Need to recalculate uu and chol_u since covariate decomp changed
            private$calc_temporal_marginal_ll()
            ll_eval <- self$temporal_ll_evaluator
            self$temporal_decomp <- private$sample_decomp_norm(self$temporal_decomp, ll_eval$chol_lu, ll_eval$uu)
        },

        #~ @description Set BKTR error values (MAE, RMSE, Total Sq. Error) and sample a new tau
        set_errors_and_sample_precision_tau = function() {
            self$result_logger$set_y_and_beta_estimates(private$decomposition_tensors())
            error_metrics <- self$result_logger$set_error_metrics()
            self$tau <- self$tau_sampler$sample(error_metrics[['total_sq_error']])
        },

        #~ @description Get a list of current iteration scalar parameters that need to be
        #~ gathered as historical data
        logged_scalar_params = function() {
            return(list(
                tau = as.double(self$tau),
                spatial_length = self$spatial_length_sampler$theta_value,
                decay_scale = self$decay_scale_sampler$theta_value,
                periodic_length = self$periodic_length_sampler$theta_value
            ))
        },

        #~ @description List of all used decomposition tensors
        decomposition_tensors = function() {
            return(
                list(
                    spatial_decomp = self$spatial_decomp,
                    temporal_decomp = self$temporal_decomp,
                    covs_decomp = self$covs_decomp
                )
            )
        },

        #~ @description Collect all necessary iteration values
        collect_iter_values = function(iter) {
            self$result_logger$collect_iter_samples(iter, private$logged_scalar_params)
        },

        #~ @description Calculate the final list of values returned by the MCMC sampling including
        #~ the y estimation, the average estimated betas and the errors
        calculate_avg_estimates = function() {
            self$result_logger$get_avg_estimates()
        },

        #~ @description Log iteration results in via the result logger
        log_iter_results = function() {
            self$result_logger$log_final_results()
        },

        #~ @description Initialize all parameters that are needed before we start the MCMC sampling
        initialize_params = function() {
            private$init_covariate_decomp()
            private$create_result_logger()
            private$create_kernel_factories()
            private$create_likelihood_evaluators()
            private$create_hparam_samplers()
        }
    )
)
