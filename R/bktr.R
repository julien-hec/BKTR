#' @import torch
#' @include tensor_ops.R
#' @include kernels.R
#' @include likelihood_evaluator.R
#' @include result_logger.R
#' @include samplers.R
#' @importFrom R6 R6Class

#' @title R6 class encapsulating the BKTR regression elements
#'
#' @description A BKTRRegressor holds all the key elements to accomplish the MCMC sampling
#' algorithm (\strong{Algorithm 1} of the paper).
#'
#' @examplesIf torch::torch_is_installed()
#' # Create a BIXI data collection instance containing multiple dataframes
#' bixi_data <- BixiData$new(is_light = TRUE) # Use light version for example
#'
#' # Create a BKTRRegressor instance
#' bktr_regressor <- BKTRRegressor$new(
#'   formula = nb_departure ~ 1 + mean_temp_c + area_park,
#'   data_df <- bixi_data$data_df,
#'   spatial_positions_df = bixi_data$spatial_positions_df,
#'   temporal_positions_df = bixi_data$temporal_positions_df,
#'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
#'
#' # Launch the MCMC sampling
#' bktr_regressor$mcmc_sampling()
#'
#' # Get the summary of the bktr regressor
#' summary(bktr_regressor)
#'
#' # Get estimated response variables for missing values
#' bktr_regressor$imputed_y_estimates
#'
#' # Get the list of sampled betas for given spatial, temporal and feature labels
#' bktr_regressor$get_iterations_betas(
#'   spatial_label = bixi_data$spatial_positions_df$location[1],
#'   temporal_label = bixi_data$temporal_positions_df$time[1],
#'   feature_label = 'mean_temp_c')
#'
#' # Get the summary of all betas for the 'mean_temp_c' feature
#' bktr_regressor$get_beta_summary_df(feature_labels = 'mean_temp_c')
#'
#' @export
BKTRRegressor <- R6::R6Class(
    classname = 'BKTRRegressor',
    public = list(
        #' @field data_df The dataframe containing all the covariates through time and space (including
        #' the response variable)
        data_df = NULL,
        #' @field y The response variable tensor
        y = NULL,
        #' @field omega The tensor indicating which response values are not missing
        omega = NULL,
        #' @field covariates The tensor containing all the covariates
        covariates = NULL,
        #' @field covariates_dim The dimensions of the covariates tensor
        covariates_dim = NULL,
        #' @field logged_params_tensor The tensor containing all the sampled hyperparameters
        logged_params_tensor = NULL,
        #' @field tau The precision hyperparameter
        tau = NULL,
        #' @field spatial_decomp The spatial covariate decomposition
        spatial_decomp = NULL, # U
        #' @field temporal_decomp The temporal covariate decomposition
        temporal_decomp = NULL, # V
        #' @field covs_decomp The feature covariate decomposition
        covs_decomp = NULL, # C or W
        # Result Logger
        #' @field result_logger The result logger instance used to store the results of the MCMC sampling
        result_logger = NULL,
        #' @field has_completed_sampling Boolean showing wheter the MCMC sampling has been completed
        has_completed_sampling = FALSE,
        # Kernels
        #' @field spatial_kernel The spatial kernel used
        spatial_kernel = NULL,
        #' @field temporal_kernel The temporal kernel used
        temporal_kernel = NULL,
        #' @field spatial_positions_df The dataframe containing the spatial positions
        spatial_positions_df = NULL,
        #' @field temporal_positions_df The dataframe containing the temporal positions
        temporal_positions_df = NULL,
        # Samplers
        #' @field spatial_params_sampler The spatial kernel hyperparameter sampler
        spatial_params_sampler = NULL,
        #' @field temporal_params_sampler The temporal kernel hyperparameter sampler
        temporal_params_sampler = NULL,
        #' @field tau_sampler The tau hyperparameter sampler
        tau_sampler = NULL,
        #' @field precision_matrix_sampler The precision matrix sampler
        precision_matrix_sampler = NULL,
        # Likelihood evaluators
        #' @field spatial_ll_evaluator The spatial likelihood evaluator
        spatial_ll_evaluator = NULL,
        #' @field temporal_ll_evaluator The temporal likelihood evaluator
        temporal_ll_evaluator = NULL,
        # Params
        #' @field rank_decomp The rank of the CP decomposition
        rank_decomp = NULL,
        #' @field burn_in_iter The number of burn in iterations
        burn_in_iter = NULL,
        #' @field sampling_iter The number of sampling iterations
        sampling_iter = NULL,
        #' @field max_iter The total number of iterations
        max_iter = NULL,
        #' @field a_0 The initial value for the shape in the gamma function generating tau
        a_0 = NULL,
        #' @field b_0 The initial value for the rate in the gamma function generating tau
        b_0 = NULL,
        #' @field formula The formula used to specify the relation between the response variable and the covariates
        formula = NULL,
        #' @field spatial_labels The spatial labels
        spatial_labels = NULL,
        #' @field temporal_labels The temporal labels
        temporal_labels = NULL,
        #' @field feature_labels The feature labels
        feature_labels = NULL,
        #' @field geo_coords_projector The geographic coordinates projector
        geo_coords_projector = NULL,

        #' @description Create a new \code{BKTRRegressor} object.
        #' @param data_df data.table: A dataframe containing all the covariates
        #' through time and space. It is important that the dataframe has a two
        #' indexes named `location` and `time` respectively. The dataframe should
        #' also contain every possible combinations of `location` and `time`
        #' (i.e. even missing rows should be filled present but filled with NaN).
        #' So if the dataframe has 10 locations and 5 time points, it should have
        #' 50 rows (10 x 5). If formula is None, the dataframe should contain
        #' the response variable `Y` as the first column. Note that the covariate
        #' columns cannot contain NaN values, but the response variable can.
        #' @param formula A Wilkinson R formula to specify the relation
        #' between the response variable `Y` and the covariates. If Null, the first
        #' column of the data frame will be used as the response variable and all the
        #' other columns will be used as the covariates.  Defaults to Null.
        #' @param rank_decomp Integer: Rank of the CP decomposition (Paper -- \eqn{R}). Defaults to 10.
        #' @param burn_in_iter Integer: Number of iteration before sampling (Paper -- \eqn{K_1}). Defaults to 500.
        #' @param sampling_iter Integer: Number of sampling iterations (Paper -- \eqn{K_2}). Defaults to 500.
        #' @param spatial_positions_df data.table: Spatial kernel input tensor used
        #' to calculate covariates' distance. Vector of length equal to the number of location points.
        #' @param temporal_positions_df data.table: Temporal kernel input tensor used to
        #' calculate covariate distance. Vector of length equal to the number of time points.
        #' @param spatial_kernel Kernel: Spatial kernel Used. Defaults to
        #' a KernelMatern(smoothness_factor=3).
        #' @param temporal_kernel Kernel: Temporal kernel used. Defaults to KernelSE().
        #' @param sigma_r Numeric:  Variance of the white noise process (\eqn{\tau^{-1}})
        #' defaults to 1E-2.
        #' @param a_0 Numeric: Initial value for the shape (\eqn{\alpha}) in the gamma function
        #' generating tau defaults to 1E-6.
        #' @param b_0 Numeric: Initial value for the rate (\eqn{\beta}) in the gamma function
        #' generating tau defaults to 1E-6.
        #' @param has_geo_coords Boolean: Whether the spatial positions df use geographic coordinates
        #' (latitude, longitude). Defaults to TRUE.
        #' @param geo_coords_scale Numeric: Scale factor to convert geographic coordinates to euclidean
        #' 2D space via Mercator projection using x & y domains of [-scale/2, +scale/2]. Only used if
        #' has_geo_coords is TRUE. Defaults to 10.
        #' @return A new \code{BKTRRegressor} object.
        initialize = function(
            data_df,
            spatial_positions_df,
            temporal_positions_df,
            rank_decomp = 10,
            burn_in_iter = 500,
            sampling_iter = 500,
            formula = NULL,
            spatial_kernel = KernelMatern$new(smoothness_factor = 3),
            temporal_kernel = KernelSE$new(),
            sigma_r = 1E-2,
            a_0 = 1E-6,
            b_0 = 1E-6,
            has_geo_coords = TRUE,
            geo_coords_scale = 10
        ) {
            self$has_completed_sampling <- FALSE
            private$verify_input_labels(data_df, spatial_positions_df, temporal_positions_df)

            # We don't need to sort since keys in data.table already do so
            self$data_df <- data_df
            self$temporal_positions_df <- temporal_positions_df
            if (has_geo_coords) {
                self$geo_coords_projector <- GeoMercatorProjector$new(spatial_positions_df, geo_coords_scale)
                self$spatial_positions_df <- self$geo_coords_projector$scaled_ini_df
            } else {
                self$spatial_positions_df <- spatial_positions_df
            }


            # Set formula and get model's matrix
            xy_df_list <- private$get_x_and_y_dfs_from_formula(self$data_df[, -c('location', 'time')], formula)
            y_df <- xy_df_list$y_df
            x_df <- xy_df_list$x_df

            # Set labels
            self$spatial_labels <- self$spatial_positions_df$location
            self$temporal_labels <- self$temporal_positions_df$time
            self$feature_labels <- colnames(x_df)

            # Tensor Assignation
            y_matrix <- matrix(y_df[[1]], ncol = length(self$temporal_labels), byrow = TRUE)
            # Omega is 1 if y is not NA, 0 otherwise
            self$omega <- TSR$tensor(ifelse(is.na(y_matrix), 0.0, 1.0))
            # Y should replace all NA values by 0
            y_matrix[is.na(y_matrix)] <- 0.0
            self$y <- TSR$tensor(y_matrix)
            covariates <- TSR$tensor(as.matrix(x_df))
            self$tau <- 1 / TSR$tensor(sigma_r)

            # Params Assignation
            self$rank_decomp <- rank_decomp
            self$burn_in_iter <- burn_in_iter
            self$sampling_iter <- sampling_iter
            self$max_iter <- burn_in_iter + sampling_iter
            self$a_0 <- a_0
            self$b_0 <- b_0

            # Reshape Covariates
            private$reshape_covariates(covariates, length(self$spatial_labels), length(self$temporal_labels))

            #Kernel Assignation
            self$spatial_kernel <- spatial_kernel
            self$temporal_kernel <- temporal_kernel
            self$spatial_kernel$set_positions(self$spatial_positions_df)
            self$temporal_kernel$set_positions(self$temporal_positions_df)
            # Create First Kernels
            self$spatial_kernel$kernel_gen()
            self$temporal_kernel$kernel_gen()
        },

        #' @description Launch the MCMC sampling process. \cr
        #' For a predefined number of iterations:
        #' \enumerate{
        #' \item{Sample spatial kernel hyperparameters}
        #' \item{Sample temporal kernel hyperparameters}
        #' \item{Sample the precision matrix from a wishart distribution}
        #' \item{Sample a new spatial covariate decomposition}
        #' \item{Sample a new feature covariate decomposition}
        #' \item{Sample a new temporal covariate decomposition}
        #' \item{Calculate respective errors for the iterations}
        #' \item{Sample a new tau value}
        #' \item{Collect all the important data for the iteration}
        #' }
        #' @return NULL Results are stored and can be accessed via summary()
        mcmc_sampling = function() {
            private$initialize_params()
            for (i in 1:self$max_iter) {
                private$sample_kernel_hparam()
                private$sample_precision_wish()
                private$sample_spatial_decomp()
                private$sample_covariate_decomp()
                private$sample_temporal_decomp()
                private$set_errors_and_sample_precision_tau(i)
                private$collect_iter_values(i)
            }
            private$log_final_iter_results()
        },

        #' @description Use interpolation to predict betas and response values for new data.
        #' @param new_data_df data.table: New covariates. Must have the same columns as
        #'   the covariates used to fit the model. The index should contain the combination
        #'   of all old spatial coordinates with all new temporal coordinates, the combination
        #'   of all new spatial coordinates with all old temporal coordinates, and the
        #'   combination of all new spatial coordinates with all new temporal coordinates.
        #' @param new_spatial_positions_df data.table or NULL: A data frame containing the new
        #'   spatial positions. Defaults to NULL.
        #' @param new_temporal_positions_df data.table or NULL: A data frame containing the new
        #'   temporal positions. Defaults to NULL.
        #' @param jitter Numeric or NULL: A small value to add to the diagonal of the precision matrix.
        #'   Defaults to NULL.
        #'
        #' @examplesIf torch::torch_is_installed()
        #' ## PREDICTION EXAMPLE ##
        #' # Create a light version of the BIXI data collection instance
        #' bixi_data <- BixiData$new(is_light = TRUE)
        #' # Simplify variable names
        #' data_df <- bixi_data$data_df
        #' spa_pos_df <- bixi_data$spatial_positions_df
        #' temp_pos_df <- bixi_data$temporal_positions_df
        #'
        #' # Keep some data aside for prediction
        #' new_spa_pos_df <- spa_pos_df[1:2, ]
        #' new_temp_pos_df <- temp_pos_df[1:5, ]
        #' reg_spa_pos_df <- spa_pos_df[-(1:2), ]
        #' reg_temp_pos_df <- temp_pos_df[-(1:5), ]
        #' reg_data_df_mask <- data_df$location %in% reg_spa_pos_df$location &
        #'   data_df$time %in% reg_temp_pos_df$time
        #' reg_data_df <- data_df[reg_data_df_mask, ]
        #' new_data_df <- data_df[!reg_data_df_mask, ]
        #'
        #' # Launch mcmc sampling on regression data
        #' bktr_regressor <- BKTRRegressor$new(
        #'   formula = nb_departure ~ 1 + mean_temp_c + area_park,
        #'   data_df = reg_data_df,
        #'   spatial_positions_df = reg_spa_pos_df,
        #'   temporal_positions_df = reg_temp_pos_df,
        #'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
        #' bktr_regressor$mcmc_sampling()
        #'
        #' # Predict response values for new data
        #' bktr_regressor$predict(
        #'   new_data_df = new_data_df,
        #'   new_spatial_positions_df = new_spa_pos_df,
        #'   new_temporal_positions_df = new_temp_pos_df)
        #'
        #' @return List: A list of two dataframes. The first represents the beta
        #'   forecasted for all new spatial locations or temporal points.
        #'   The second represents the forecasted response for all new spatial
        #'   locations or temporal points.
        predict = function(
            new_data_df,
            new_spatial_positions_df = NULL,
            new_temporal_positions_df = NULL,
            jitter = 1e-5
        ) {
            # private$pred_valid_and_sort_data(
            #     new_data_df, new_spatial_positions_df, new_temporal_positions_df
            # )
            ini_fp_type <- TSR$fp_type
            TSR$set_params(fp_type = 'float64')

            if (!is.null(new_spatial_positions_df) && !is.null(self$geo_coords_projector)) {
                new_spatial_positions_df <- self$geo_coords_projector$project_new_coords(new_spatial_positions_df)
            }

            spatial_positions_df <- (
                if (!is.null(new_spatial_positions_df)) rbind(self$spatial_positions_df, new_spatial_positions_df)
                else self$spatial_positions_df
            )
            temporal_positions_df <- (
                if (!is.null(new_temporal_positions_df)) rbind(self$temporal_positions_df, new_temporal_positions_df)
                else self$temporal_positions_df
            )
            data_df <- rbind(self$data_df, new_data_df)
            spa_order_df <- data.table(
                location = spatial_positions_df$location, spa_order = seq_len(nrow(spatial_positions_df))
            )
            temp_order_df <- data.table(
                time = temporal_positions_df$time, temp_order = seq_len(nrow(temporal_positions_df))
            )
            data_df <- data_df[spa_order_df, on = 'location'][temp_order_df, on = 'time']
            setorder(data_df, spa_order, temp_order)
            data_df[, c('spa_order', 'temp_order') := list(NULL, NULL)]

            # private$verify_input_labels(
            #     data_df,
            #     spatial_positions_df,
            #     temporal_positions_df
            # )
            all_betas <- TSR$zeros(
                c(
                    nrow(spatial_positions_df),
                    nrow(temporal_positions_df),
                    length(self$feature_labels),
                    self$sampling_iter
                )
            )
            for (i in seq_len(self$sampling_iter)) {
                new_spa_decomp <- private$pred_simu_new_decomp(
                    'spatial', i, spatial_positions_df, new_spatial_positions_df, jitter
                )
                new_temp_decomp <- private$pred_simu_new_decomp(
                    'temporal', i, temporal_positions_df, new_temporal_positions_df, jitter
                )
                covs_decomp <- TSR$tensor(self$result_logger$covs_decomp_per_iter[, , i])
                all_betas[, , , i] <- torch::torch_einsum(
                    'il,jl,kl->ijk', c(new_spa_decomp, new_temp_decomp, covs_decomp)
                )
            }

            new_betas <- all_betas$mean(dim = -1)
            x_df <- private$get_x_and_y_dfs_from_formula(data_df[, -c('location', 'time')], self$formula)$x_df
            covariates <- TSR$tensor(as.matrix(x_df))$reshape(
                c(nrow(spatial_positions_df), nrow(temporal_positions_df), -1)
            )
            new_y_est <- torch::torch_einsum('ijk,ijk->ij', c(new_betas, covariates))
            new_index_df <- data_df[, c('location', 'time')]
            new_beta_df <- cbind(
                new_index_df,
                data.table(as.matrix(new_betas$flatten(start_dim = 1, end_dim = 2)$cpu()))
            )
            setnames(new_beta_df, c('location', 'time', self$feature_labels))
            new_y_df <- cbind(new_index_df, data.table(as.matrix(new_y_est$flatten()$cpu(), byrow = TRUE)))
            setnames(new_y_df, c('location', 'time', 'y_est'))
            new_locs <- unique(new_spatial_positions_df$location)
            new_times <- unique(new_temporal_positions_df$time)
            new_beta_df <- new_beta_df[, new_beta_df[new_beta_df[, .I[location %in% new_locs | time %in% new_times]], ]]
            new_y_df <- new_y_df[, new_y_df[new_y_df[, .I[location %in% new_locs | time %in% new_times]], ]]
            TSR$set_params(fp_type = ini_fp_type)
            return(list(new_y_df = new_y_df, new_beta_df = new_beta_df))
        },

        #' @description Return all sampled betas through sampling iterations for a given
        #' set of spatial, temporal and feature labels. Useful for plotting the
        #' distribution of sampled beta values.
        #' @param spatial_label String: The spatial label for which we want to get the betas
        #' @param temporal_label String: The temporal label for which we want to get the betas
        #' @param feature_label String: The feature label for which we want to get the betas
        #' @return A list containing the sampled betas through iteration for the given labels
        get_iterations_betas = function(spatial_label, temporal_label, feature_label) {
            if (!self$has_completed_sampling) {
                stop('Beta values can only be accessed after MCMC sampling.')
            }
            beta_per_iter_tensor <- self$result_logger$get_iteration_betas_tensor(
                c(spatial_label), c(temporal_label), c(feature_label)
            )[1]
            return(as.array(beta_per_iter_tensor))
        },

        #' @description Get a summary of estimated beta values. If no labels are given,
        #' then the summary is for all the betas. If labels are given, then the summary
        #' is for the given labels.
        #' @param spatial_labels vector: The spatial labels used in summary. If NULL,
        #' then all spatial labels are used. Defaults to NULL.
        #' @param temporal_labels vector: The temporal labels used in summary. If NULL,
        #' then all temporal labels are used. Defaults to NULL.
        #' @param feature_labels vector: The feature labels used in summary. If NULL,
        #' then all feature labels are used. Defaults to NULL.
        #' @return A new data.table with the beta summary for the given labels.
        get_beta_summary_df = function(
            spatial_labels = NULL,
            temporal_labels = NULL,
            feature_labels = NULL
        ) {
            if (!self$has_completed_sampling) {
                stop('Beta values can only be accessed after MCMC sampling.')
            }
            return(self$result_logger$get_beta_summary_df(spatial_labels, temporal_labels, feature_labels))
        }
    ),

    private = list(
        #~ @description Verify if kernel inputs are valid and align with covariates labels.
        verify_kernel_labels = function(
            kernel_positions,
            expected_labels,
            kernel_type
        ) {
            cov_related_indx_name <- ifelse(kernel_type == 'spatial', 'location', 'time')
            if (key(kernel_positions) != cov_related_indx_name) {
                stop(sprintf(
                    '`%s_positions_df` must have a `%s` key.',
                    kernel_type, cov_related_indx_name
                ))
            }
            if (!identical(expected_labels, unique(kernel_positions[[cov_related_indx_name]]))) {
                stop(paste0(
                    '`', kernel_type, '_positions_df` must contain in its ', cov_related_indx_name,
                    ' index the unique values located in `data_df` ', cov_related_indx_name, ' index.'
                ))
            }
        },

        #~ @description Verify validity of BKTR dataframe input labels
        verify_input_labels = function(
            data_df,
            spatial_positions_df,
            temporal_positions_df
        ) {
            if (!identical(key(data_df), c('location', 'time'))) {
                stop(paste(
                    'The data_df dataframe must have a multi index on location and time.',
                    'Set the keys on the table with `setkey(data_df, location, time)``.'
                ))
            }
            loc_set <- unique(data_df$location)
            time_set <- unique(data_df$time)
            product_set <- CJ(location = loc_set, time = time_set)
            data_df_index_set <- unique(data_df[, c('location', 'time')])

            if (!identical(data_df_index_set, product_set)) {
                stop(paste(
                    'The data_df dataframe must have a row for every possible combination of location and time.',
                    'Even if response values are missing (NaN).'
                ))
            }
            private$verify_kernel_labels(spatial_positions_df, loc_set, 'spatial')
            private$verify_kernel_labels(temporal_positions_df, time_set, 'temporal')
        },

        #~ @description Use formula to get x and y dataframes.
        #~ @param data_df data.table: The initial dataframe used to obtain the x and y dataframes.
        #~ @param formula: Formula to give the y and X dataframes matrix. If formula is
        #~     None, use the first column as y and all other columns as covariates.
        #~ @return A list containing the y and x dataframes.
        get_x_and_y_dfs_from_formula = function(data_df, formula = NULL) {
            if (is.null(formula)) {
                formula <- paste0(colnames(data_df)[1], ' ~ .')
            }
            self$formula <- as.formula(formula)
            formula_y_name <- as.character(self$formula[[2]])
            if (length(formula_y_name) != 1) {
                stop(paste(
                    'The formula provided to the regressor is not valid.',
                    'It must contain one and only one response variable.'
                ))
            }
            mf <- model.frame(self$formula, data = data_df, na.action=na.pass)
            x_df <- data.table(model.matrix(self$formula, mf))
            y_df <- data.table(model.response(mf))
            x_colnames <- colnames(x_df)
            if ('(Intercept)' %in% x_colnames) {
                colnames(x_df)[x_colnames == '(Intercept)'] <- 'Intercept'
            }
            setnames(y_df, formula_y_name)
            return(list(y_df=y_df, x_df=x_df))
        },

        #~ @description Reshape the covariate tensors into one single tensor and set this
        #~ tensor into the \code{covariates} property of the BKTRRegressor object
        #~ @param temporal_covariate_tensor tensor: Temporal Covariates
        #~ @param spatial_covariate_tensor tensor: Spatial Covariates
        reshape_covariates = function(covariate_tensor, nb_locations, nb_times) {
            nb_covariates <- covariate_tensor$shape[2]
            self$covariates_dim <- list(
                nb_spaces = nb_locations,  # S
                nb_times = nb_times,  # T
                nb_covariates = nb_covariates  # P
            )
            self$covariates <- covariate_tensor$reshape(c(nb_locations, nb_times, nb_covariates))
        },

        #~ @description Initialize the CP decomposed covariate tensors
        #~ using normally distributed random values
        init_covariate_decomp = function() {
            rank_decomp <- self$rank_decomp
            covs_dim <- self$covariates_dim

            self$spatial_decomp <- TSR$randn(c(covs_dim$nb_spaces, rank_decomp))
            self$temporal_decomp <- TSR$randn(c(covs_dim$nb_times, rank_decomp))
            self$covs_decomp <- TSR$randn(c(covs_dim$nb_covariates, rank_decomp))
        },

        create_result_logger = function() {
            self$result_logger <- ResultLogger$new(
                y = self$y,
                omega = self$omega,
                covariates = self$covariates,
                nb_burn_in_iter = self$burn_in_iter,
                nb_sampling_iter = self$sampling_iter,
                rank_decomp = self$rank_decomp,
                formula = self$formula,
                spatial_labels = self$spatial_labels,
                temporal_labels = self$temporal_labels,
                feature_labels = self$feature_labels,
                spatial_kernel = self$spatial_kernel,
                temporal_kernel = self$temporal_kernel
            )
        },

        #~ @description Create and set the evaluators for the spatial and the temporal likelihoods
        create_likelihood_evaluators = function() {
            self$spatial_ll_evaluator <- MarginalLikelihoodEvaluator$new(
                self$rank_decomp,
                self$covariates_dim$nb_covariates,
                self$covariates,
                self$omega,
                self$y,
                is_transposed = FALSE
            )
            self$temporal_ll_evaluator <- MarginalLikelihoodEvaluator$new(
                self$rank_decomp,
                self$covariates_dim$nb_covariates,
                self$covariates,
                self$omega,
                self$y,
                is_transposed = TRUE
            )
        },

        #~ @description Create and set the hyperparameter samplers
        #~ for spatial and temporal kernels, tau and the precision matrix
        create_hparam_samplers = function() {
            self$spatial_params_sampler <- KernelParamSampler$new(
                kernel = self$spatial_kernel,
                marginal_ll_eval_fn = private$calc_spatial_marginal_ll
            )
            self$temporal_params_sampler <- KernelParamSampler$new(
                kernel = self$temporal_kernel,
                marginal_ll_eval_fn = private$calc_temporal_marginal_ll
            )

            self$tau_sampler <- TauSampler$new(self$a_0, self$b_0, self$omega$sum())

            self$precision_matrix_sampler <- PrecisionMatrixSampler$new(
                self$covariates_dim$nb_covariates, self$rank_decomp
            )
        },

        #~ @description Calculate the spatial marginal likelihood
        calc_spatial_marginal_ll = function() {
            return(self$spatial_ll_evaluator$calc_likelihood(
                    self$spatial_kernel$covariance_matrix, self$temporal_decomp, self$covs_decomp, self$tau
            ))
        },

        #~ @description Calculate the temporal marginal likelihood
        calc_temporal_marginal_ll = function() {
            return(self$temporal_ll_evaluator$calc_likelihood(
                self$temporal_kernel$covariance_matrix, self$spatial_decomp, self$covs_decomp, self$tau
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
        #~ @param chol_l tensor: The cholesky decomposition of the l tensor
        #~ @param uu tensor: uu decomposition
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
            self$spatial_decomp <- private$sample_decomp_norm(
                self$spatial_decomp, ll_eval$chol_lu, ll_eval$uu
            )
        },

        #~ @description Sample a new covariate decomposition
        sample_covariate_decomp = function() {
            chol_res <- get_cov_decomp_chol(
                self$spatial_decomp,
                self$temporal_decomp,
                self$covariates,
                self$rank_decomp,
                self$omega,
                self$tau,
                self$y,
                self$precision_matrix_sampler$wish_precision_tensor
            )
            self$covs_decomp <- private$sample_decomp_norm(
                self$covs_decomp, chol_res$chol_lc, chol_res$cc
            )
        },

        #~ @description Sample a new temporal covariate decomposition
        sample_temporal_decomp = function() {
            # Need to recalculate uu and chol_u since covariate decomp changed
            private$calc_temporal_marginal_ll()
            ll_eval <- self$temporal_ll_evaluator
            self$temporal_decomp <- private$sample_decomp_norm(
                self$temporal_decomp, ll_eval$chol_lu, ll_eval$uu
            )
        },

        #~ @description Set BKTR error values (MAE, RMSE, Total Sq. Error) and sample a new tau
        set_errors_and_sample_precision_tau = function(iter) {
            self$result_logger$set_y_and_beta_estimates(self$decomposition_tensors, iter)
            error_metrics <- self$result_logger$set_error_metrics()
            self$tau <- self$tau_sampler$sample(self$result_logger$total_sq_error)
        },

        #~ @description Collect all necessary iteration values
        collect_iter_values = function(iter) {
            self$result_logger$collect_iter_samples(iter, as.numeric(self$tau$cpu()))
        },

        #~ @description Log final iteration results in via the result logger
        log_final_iter_results = function() {
            self$result_logger$log_final_iter_results()
            self$has_completed_sampling <- TRUE
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
        },

        pred_simu_new_decomp = function(
            pred_type,
            iter_no,
            position_df,
            new_position_df,
            jitter
        ) {
            old_decomp <- (
                if (pred_type == 'spatial') self$result_logger$spatial_decomp_per_iter[, , iter_no]
                else self$result_logger$temporal_decomp_per_iter[, , iter_no]
            )
            old_decomp <- TSR$tensor(old_decomp)
            if (is.null(new_position_df)) {
                return(old_decomp)
            }
            nb_pos <- nrow(position_df)
            nb_new_pos <- nrow(new_position_df)
            nb_old_pos <- nb_pos - nb_new_pos

            old_kernel <- if (pred_type == 'spatial') self$spatial_kernel else self$temporal_kernel
            new_kernel <- old_kernel
            for (param in new_kernel$parameters){
                if (!param$is_fixed) {
                    param_full_repr <- paste(capitalize_str(pred_type), param$full_name, sep = ' - ')
                    param$value <- as.numeric(
                        self$result_logger$hyperparameters_per_iter_df[iter_no, ..param_full_repr]
                    )
                }
            }
            new_kernel$set_positions(position_df)
            cov_mat <- new_kernel$kernel_gen()
            old_cov <- cov_mat[1:nb_old_pos, 1:nb_old_pos]
            new_old_cov <- cov_mat[-nb_new_pos:nb_pos, 1:nb_old_pos]
            old_new_cov <- cov_mat[1:nb_old_pos, -nb_new_pos:nb_pos]
            new_cov <- cov_mat[-nb_new_pos:nb_pos, -nb_new_pos:nb_pos]
            new_decomp_mus <- new_old_cov$matmul(old_cov$inverse())$matmul(old_decomp)
            new_decomp_cov <- new_cov - new_old_cov$matmul(old_cov$inverse())$matmul(old_new_cov)
            new_decomp_cov <- (new_decomp_cov + new_decomp_cov$t()) / 2
            if (!is.null(jitter)) {
                new_decomp_cov <- new_decomp_cov + jitter * TSR$eye(new_decomp_cov$shape[1])
            }
            new_decomp <- (
                torch::distr_multivariate_normal(new_decomp_mus$t(), new_decomp_cov)$sample()$t()
            )
            return(torch::torch_cat(c(old_decomp, new_decomp), dim = 1))
        }
    ),

    active = list(
        #' @field summary A summary of the BKTRRegressor instance
        summary = function() {
            if (!self$has_completed_sampling) {
                stop('Summary can only be accessed after running the MCMC sampling.')
            }
            return(self$result_logger$summary())
        },
        #' @field beta_covariates_summary A dataframe containing the summary of the beta covariates
        beta_covariates_summary = function() {
            if (!self$has_completed_sampling) {
                stop('Beta covariates summary can only be accessed after running the MCMC sampling.')
            }
            return(self$result_logger$beta_covariates_summary_df)
        },
        #' @field y_estimates A dataframe containing the y estimates
        y_estimates = function() {
            if (!self$has_completed_sampling) {
                stop('Y estimates can only be accessed after running the MCMC sampling.')
            }
            y_est <- bktr_regressor$result_logger$y_estimates_df
            y_est[as.array(bktr_regressor$omega$flatten()$cpu()) == 0, 3] <- NaN
            return(y_est)
        },
        #' @field imputed_y_estimates A dataframe containing the imputed y estimates
        imputed_y_estimates = function() {
            if (!self$has_completed_sampling) {
                stop('Imputed Y estimates can only be accessed after running the MCMC sampling.')
            }
            return(self$result_logger$y_estimates_df)
        },
        #' @field beta_estimates A dataframe containing the beta estimates
        beta_estimates = function() {
            if (!self$has_completed_sampling) {
                stop('Beta estimates can only be accessed after running the MCMC sampling.')
            }
            return(self$result_logger$beta_estimates_df)
        },
        #' @field hyperparameters_per_iter_df A dataframe containing the beta estimates per iteration
        hyperparameters_per_iter_df = function() {
            if (!self$has_completed_sampling) {
                stop('Hyperparameters trace can only be accessed after running the MCMC sampling.')
            }
            return(self$result_logger$hyperparameters_per_iter_df)
        },
        #' @field decomposition_tensors List of all used decomposition tensors
        decomposition_tensors = function() {
            return(
                list(
                    spatial_decomp = self$spatial_decomp,
                    temporal_decomp = self$temporal_decomp,
                    covs_decomp = self$covs_decomp
                )
            )
        }

    )
)

#' @title Summarize a BKTRRegressor instance
#' @param object A BKTRRegressor instance
#' @param ... Additional arguments to comply with generic function
#' @export
summary.BKTRRegressor <- function(object, ...) {
    cat(object$summary)
}


#' @title Print the summary of a BKTRRegressor instance
#' @param x A BKTRRegressor instance
#' @param ... Additional arguments to comply with generic function
#' @export
print.BKTRRegressor <- function(x, ...) {
    cat(x$summary)
}
