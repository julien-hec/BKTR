#' @import torch
#' @include tensor_ops.R
#' @include kernels.R
#' @include marginal_likelihoods.R
#' @include result_logger.R
#' @include samplers.R
#' @importFrom R6 R6Class

#' @title R6 class encapsulating the BKTR regression steps
#'
#' @description A BKTRRegressor holds all the key elements to accomplish the MCMC sampling
#' algorithm (\strong{Algorithm 1} of the paper).
#'
#' @export
BKTRRegressor <- R6::R6Class(
    public = list(
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
        # Kernels
        temporal_kernel = NULL,
        spatial_kernel = NULL,
        # Samplers
        spatial_params_sampler = NULL,
        temporal_params_sampler = NULL,
        tau_sampler = NULL,
        precision_matrix_sampler = NULL,
        # Likelihood evaluators
        spatial_ll_evaluator = NULL,
        temporal_ll_evaluator = NULL,
        # Params
        rank_decomp = NULL,
        burn_in_iter = NULL,
        sampling_iter = NULL,
        max_iter = NULL,
        a_0 = NULL,
        b_0 = NULL,
        # Export Params
        sampled_beta_indexes = NULL,
        sampled_y_indexes = NULL,
        results_export_dir = NULL,
        results_export_suffix = NULL,
        # Labels
        spatial_labels = NULL,
        temporal_labels = NULL,
        spatial_feature_labels = NULL,
        temporal_feature_labels = NULL,

        #' @description Create a new \code{BKTRRegressor} object.
        #' @param spatial_covariates_df df: Spatial Covariates, 2D dataframe (spatial points x spatial features)
        #' @param temporal_covariates_df df: Temporal Covariates, 2D dataframe (temporal points x temporal features)
        #' @param y_df df: Variable that we try to predict, 2D dataframe (spatial points x temporal points)
        #' @param rank_decomp Integer: Rank of the CP decomposition (Paper -- \eqn{R})
        #' @param burn_in_iter Integer: Number of iteration before sampling (Paper -- :math:`K_1`)
        #' @param sampling_iter Integer: Number of sampling iterations (Paper -- :math:`K_2`)
        #' @param spatial_kernel: Kernel: Spatial kernel Used. Defaults to KernelMatern(smoothness_factor=3).
        #' @param spatial_x_df df: Spatial kernel input tensor used to calculate covariate distance.
        #' Vector of length equal to nb spatial points. Defaults to NULL.
        #' @param spatial_dist_df df: Spatial kernel covariate distance. A two dimensions df
        #' (nb spatial points x nb spatial points).  Should be used instead of `spatial_kernel_x`
        #' if distance was already calculated. Defaults to NULL.
        #' @param temporal_kernel: Kernel: Temporal kernel Used. Defaults to KernelSE().
        #' @param temporal_x_df: df: Temporal kernel input tensor used to calculate covariate distance.
        #' Vector of length equal to nb temporal points. Defaults to Null.
        #' @param temporal_dist_df: df: Spatial kernel covariate distance. A two dimensions df
        #' (nb spatial points x nb spatial points).  Should be used instead of `spatial_kernel_x`
        #' if distance was already calculated. Defaults to NULL.
        #' @param sigma_r Numeric: Initial value of white noise's variance (\eqn{\tau^{-1}})
        #' @param a_0 Numeric: Initial value for the shape (\eqn{\alpha}) in the gamma function generating tau
        #' @param b_0 Numeric: Initial value for the rate (\eqn{\beta}) in the gamma function generating tau
        #' @param sampled_beta_indexes vector(Integer): Indexes of beta estimates that need to be sampled
        #' through iterations. Defaults to c().
        #' @param sampled_y_indexes vector(Integer): Indexes of y estimates that need to be sampled through
        #' iterations. Defaults to c().
        #' @param results_export_dir (str): Path of the folder where the csv file will be exported
        #' (if NULL it is only printed). Defaults to NULL.
        #' @param results_export_suffix (str): Suffix added at the end of the csv file name
        #' (if NULL, no suffix is added). Defaults to NULL.
        #' @return A new \code{BKTRRegressor} object.
        initialize = function(
            spatial_covariates_df,
            temporal_covariates_df,
            y_df,
            rank_decomp,
            burn_in_iter,
            sampling_iter,
            spatial_kernel = KernelMatern$new(smoothness_factor=3),
            spatial_x_df = NULL,
            spatial_dist_df = NULL,
            temporal_kernel = KernelSE$new(),
            temporal_x_df = NULL,
            temporal_dist_df = NULL,
            sigma_r = 1E-2,
            a_0 = 1E-6,
            b_0 = 1E-6,
            sampled_beta_indexes = c(),
            sampled_y_indexes = c(),
            results_export_dir = NULL,
            results_export_suffix = NULL
        ) {
            private$verify_input_labels(
                y_df,
                spatial_covariates_df,
                temporal_covariates_df,
                spatial_x_df,
                spatial_dist_df,
                temporal_x_df,
                temporal_dist_df
            )

            # Set labels
            self$spatial_labels <- rownames(spatial_covariates_df)
            self$temporal_labels <- rownames(temporal_covariates_df)
            self$spatial_feature_labels <- colnames(spatial_covariates_df)
            self$temporal_feature_labels <- colnames(temporal_covariates_df)

            # Tensor Assignation
            y_matrix <- as.matrix(y_df)
            # Omega is 1 if y is not NA, 0 otherwise
            self$omega <- tsr$new_tensor(ifelse(is.na(y_matrix), 0.0, 1.0))
            # Y should replace all NA values by 0
            y_matrix[is.na(y_matrix)] <- 0.0
            self$y <- tsr$new_tensor(y_matrix)
            spatial_covariates <- tsr$new_tensor(as.matrix(spatial_covariates_df))
            temporal_covariates <- tsr$new_tensor(as.matrix(temporal_covariates_df))
            self$tau <- 1 / tsr$new_tensor(sigma_r)
            temporal_x_tsr <- tsr$get_df_tensor_or_null(temporal_x_df)
            temporal_dist_tsr <- tsr$get_df_tensor_or_null(temporal_dist_df)
            spatial_x_tsr <- tsr$get_df_tensor_or_null(spatial_x_df)
            spatial_dist_tsr <- tsr$get_df_tensor_or_null(spatial_dist_df)

            # Params Assignation
            self$rank_decomp <- rank_decomp
            self$burn_in_iter <- burn_in_iter
            self$sampling_iter <- sampling_iter
            self$max_iter <- burn_in_iter + sampling_iter
            self$a_0 <- a_0
            self$b_0 <- b_0
            self$sampled_beta_indexes <- sampled_beta_indexes
            self$sampled_y_indexes <- sampled_y_indexes
            self$results_export_dir <- results_export_dir
            self$results_export_suffix <- results_export_suffix

            # Reshape Covariates
            private$reshape_covariates(spatial_covariates, temporal_covariates)

            #Kernel Assignation
            self$spatial_kernel <- spatial_kernel
            self$temporal_kernel <- temporal_kernel
            self$spatial_kernel$set_distance_matrix(spatial_x_tsr, spatial_dist_tsr)
            self$temporal_kernel$set_distance_matrix(temporal_x_tsr, temporal_dist_tsr)

            # Create First Kernels
            self$spatial_kernel$kernel_gen()
            self$temporal_kernel$kernel_gen()
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
            private$initialize_params()
            for (i in 1:self$max_iter) {
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
            return(private$log_iter_results())
        }
    ),

    private = list(
        #~ @description Verify if kernel inputs are valid and align with covariates labels.
        verify_kernel_labels = function(
            kernel_x,
            kernel_dist,
            expected_labels,
            kernel_type
        ) {
            if (is.null(kernel_x) == is.null(kernel_dist)) {
                torch:::value_error(
                    sprintf(
                        'The %s kernel requires either a kernel_x or a kernel_dist input.',
                        kernel_type
                    )
                )
            }
            if (!is.null(kernel_x) && !identical(rownames(kernel_x), expected_labels)) {
                torch:::value_error(
                    sprintf(
                        'The %s kernel_x input should have the same index as the %s covariates dataframe.',
                        kernel_type, kernel_type
                    )
                )
            }
            if (
                !is.null(kernel_dist) && !(
                    identical(rownames(kernel_dist), colnames(kernel_dist))
                    && identical(colnames(kernel_dist), expected_labels)
                )
            ) {
                torch:::value_error(
                    sprintf(
                        'The %s_dist input should have the same index and columns as the %s covariates dataframe.',
                        kernel_type, kernel_type
                    )
                )
            }
        },

        #~ @description Verify validity of BKTR dataframe input labels
        verify_input_labels = function(
            y,
            spatial_covariates,
            temporal_covariates,
            spatial_kernel_x,
            spatial_kernel_dist,
            temporal_kernel_x,
            temporal_kernel_dist
        ) {
            spatial_labels <- rownames(spatial_covariates)
            temporal_labels <- rownames(temporal_covariates)
            if (!identical(spatial_labels, rownames(y))) {
                torch:::value_error('The spatial_covariates and y dataframes must have the same index labels.')
            }
            if (!identical(temporal_labels, colnames(y))) {
                torch:::value_error(
                    'The temporal_covariates index should hold the same values as the y dataframe column names.'
                )
            }
            private$verify_kernel_labels(spatial_kernel_x, spatial_kernel_dist, spatial_labels, 'spatial')
            private$verify_kernel_labels(temporal_kernel_x, temporal_kernel_dist, temporal_labels, 'temporal')
        },

        #~ @description Reshape the covariate tensors into one single tensor and set this
        #~ tensor into the \code{covariates} property of the BKTRRegressor object
        #~ @param temporal_covariate_tensor tensor: Temporal Covariates
        #~ @param spatial_covariate_tensor tensor: Spatial Covariates
        reshape_covariates = function(spatial_covariate_tensor, temporal_covariate_tensor) {
            space_cov_shape <- spatial_covariate_tensor$shape
            time_cov_shape <- temporal_covariate_tensor$shape
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
            rank_decomp <- self$rank_decomp
            covs_dim <- self$covariates_dim

            self$spatial_decomp <- tsr$new_tensor(
                tsr$randn(c(covs_dim$nb_spaces, rank_decomp))
            )
            self$temporal_decomp <- tsr$new_tensor(
                tsr$randn(c(covs_dim$nb_times, rank_decomp))
            )
            self$covs_decomp <- tsr$new_tensor(
                tsr$randn(c(covs_dim$nb_covariates, rank_decomp))
            )
        },

        create_result_logger = function() {
            self$result_logger <- ResultLogger$new(
                y = self$y,
                omega = self$omega,
                covariates = self$covariates,
                nb_iter = self$max_iter,
                nb_burn_in_iter = self$burn_in_iter,
                sampled_beta_indexes = self$sampled_beta_indexes,
                sampled_y_indexes = self$sampled_y_indexes,
                results_export_dir = self$results_export_dir,
                results_export_suffix = self$results_export_suffix
            )
        },

        #~ @description Create and set the evaluators for the spatial and the temporal likelihoods
        create_likelihood_evaluators = function() {
            rank_decomp <- self$rank_decomp

            self$spatial_ll_evaluator <- MarginalLikelihoodEvaluator$new(
                rank_decomp, self$covariates_dim$nb_covariates, self$covariates,
                self$omega, self$y, is_transposed = FALSE
            )
            self$temporal_ll_evaluator <- MarginalLikelihoodEvaluator$new(
                rank_decomp, self$covariates_dim$nb_covariates, self$covariates,
                self$omega, self$y, is_transposed = TRUE
            )
        },

        #~ @description Create and set the hyperparameter samplers
        #~ for spatial and temporal kernels, tau and the precision matrix
        create_hparam_samplers = function() {
            self$spatial_params_sampler <- KernelParamSampler$new(
                kernel=self$spatial_kernel,
                marginal_ll_eval_fn=private$calc_spatial_marginal_ll
            )
            self$temporal_params_sampler <- KernelParamSampler$new(
                kernel=self$temporal_kernel,
                marginal_ll_eval_fn=private$calc_temporal_marginal_ll
            )

            self$tau_sampler <- TauSampler$new(
                self$a_0, self$b_0, self$omega$sum()
            )

            self$precision_matrix_sampler <- PrecisionMatrixSampler$new(
                self$covariates_dim$nb_covariates, self$rank_decomp
            )
        },

        #~ @description Calculate the spatial marginal likelihood
        calc_spatial_marginal_ll = function() {
            return(self$spatial_ll_evaluator$calc_likelihood(
                    self$spatial_kernel$kernel, self$temporal_decomp, self$covs_decomp, self$tau
            ))
        },

        #~ @description Calculate the temporal marginal likelihood
        calc_temporal_marginal_ll = function() {
            return(self$temporal_ll_evaluator$calc_likelihood(
                self$temporal_kernel$kernel, self$spatial_decomp, self$covs_decomp, self$tau
            ))
        },

        #~ @description Sample new kernel hyperparameters
        sample_kernel_hparam = function() {
            self$spatial_params_sampler$sample()
            self$temporal_params_sampler$sample()
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
            mean_vec <- self$tau * torch::torch_triangular_solve(
                uu$unsqueeze(2),
                precision_mat,
                upper = TRUE
            )[[1]]$squeeze()
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
                self$rank_decomp, self$omega, self$tau, self$y,
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
            scalar_list <- list(tau = as.double(self$tau))
            for (p in self$spatial_kernel$parameters) {
                scalar_list[p$full_name] <- p$value
            }
            for (p in self$temporal_kernel$parameters) {
                scalar_list[p$full_name] <- p$value
            }
            return(scalar_list)
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

        #~ @description Log iteration results in via the result logger
        log_iter_results = function() {
            self$result_logger$log_iter_results()
        },

        #~ @description Initialize all parameters that are needed before we start the MCMC sampling
        initialize_params = function() {
            private$init_covariate_decomp()
            private$create_result_logger()
            private$create_likelihood_evaluators()
            private$create_hparam_samplers()

            # Calculate first likelihoods
            private$calc_spatial_marginal_ll()
            private$calc_temporal_marginal_ll()
        }
    ),
    active = list(
        beta_estimates = function(value) {
            if (is.null(self$result_logger)) {
                torch:::value_error('Beta estimates can only be accessed after running the MCMC sampling.')
            }
            return(self$result_logger$beta_estimates)
        },
        beta_stdev = function() {
            if (is.null(self$result_logger)) {
                torch:::value_error('Beta standard dev can only be accessed after running the MCMC sampling.')
            }
            return(self$result_logger$beta_stdev)
        },
        y_estimates = function() {
            if (self$result_logger == NULL) {
                torch:::value_error('Y estimates can only be accessed after running the MCMC sampling.')
            }
            return(self$result_logger$y_estimates)
        }
    )
)
