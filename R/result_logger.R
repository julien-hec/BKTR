#' @importFrom R6 "R6Class"
#' @import torch
#' @include tensor_ops.R


#' @title R6 class for Logging BKTR Results
#'
#' @description The ResultLogger encapsulate all the behavior related to
#' logging iteration information during the BKTR algorithm
#'
#' @export
#' @keywords internal
ResultLogger <- R6::R6Class(
    'ResultLogger',
    public = list(
        y = NULL,
        omega = NULL,
        covariates = NULL,
        nb_burn_in_iter = NULL,
        nb_sampling_iter = NULL,
        logged_params_map = NULL,
        beta_estimates = NULL,
        y_estimates = NULL,
        total_elapsed_time = NULL,
        formula = NULL,
        rank_decomp = NULL,
        spatial_labels = NULL,
        temporal_labels = NULL,
        feature_labels = NULL,
        spatial_kernel = NULL,
        temporal_kernel = NULL,
        hparam_labels = NULL,
        hparam_per_iter = NULL,
        spatial_decomp_per_iter = NULL,
        temporal_decomp_per_iter = NULL,
        covs_decomp_per_iter = NULL,
        sum_beta_est = NULL,
        sum_y_est = NULL,
        beta_estimates_df = NULL,
        y_estimates_df = NULL,
        beta_covariates_summary_df = NULL,
        hyperparameters_per_iter_df = NULL,
        last_time_stamp = NULL,
        error_metrics = NULL,
        total_sq_error = NULL,
            # Metrics used to create beta summaries
        moment_metrics = c('Mean', 'SD'),
        quantile_metrics = c(
            'Min',
            'Q1',
            'Median',
            'Q3',
            'Max',
            'Low2.5p',
            'Up97.5p'
        ),
        quantile_values = c(0, 0.25, 0.5, 0.75, 1, 0.025, 0.975),

        # Summary parameters
        # LINE_NCHAR = 70
        # TAB_STR = '  '
        # LINE_SEPARATOR = LINE_NCHAR * '='
        # COL_WIDTH = 18
        # DF_DISTRIB_STR_PARAMS = {
        #     'float_format': '{:,.3f}'.format,
        #     'col_space': 8,
        #     'line_width': LINE_NCHAR,
        #     'max_colwidth': COL_WIDTH,
        #     'formatters': {'__index__': lambda x, w=COL_WIDTH: f'{x:{w}}'},
        # }
        MAIN_SUMMARY_COLS = c('Mean', 'Median', 'SD'),
        DISTRIB_COLS = c('Mean', 'Median', 'SD', 'Low2.5p', 'Up97.5p'),

        initialize = function(
            y,
            omega,
            covariates,
            nb_burn_in_iter,
            nb_sampling_iter,
            rank_decomp,
            formula,
            spatial_labels,
            temporal_labels,
            feature_labels,
            spatial_kernel,
            temporal_kernel
        ) {
            # Create a tensor dictionary holding scalar data gathered through all iterations
            self$logged_params_map <- list()

            # Create tensors that accumulate values needed for estimates
            self$spatial_labels <- spatial_labels
            self$temporal_labels <- temporal_labels
            self$feature_labels <- feature_labels
            nofix_spa_params <- Filter(function(p) !p$is_fixed, spatial_kernel$parameters)
            nofix_temp_params <- Filter(function(p) !p$is_fixed, temporal_kernel$parameters)
            self$hparam_labels <- c(
                'Tau',
                paste('Spatial', sapply(nofix_spa_params, function(x) x$full_name), sep = ' - '),
                paste('Temporal', sapply(nofix_temp_params, function(x) x$full_name), sep = ' - ')
            )

            self$spatial_decomp_per_iter <- TSR$zeros(
                c(length(spatial_labels), rank_decomp, nb_sampling_iter)
            )
            self$temporal_decomp_per_iter <- TSR$zeros(
                c(length(temporal_labels), rank_decomp, nb_sampling_iter)
            )
            self$covs_decomp_per_iter <- TSR$zeros(
                c(length(feature_labels), rank_decomp, nb_sampling_iter)
            )

            self$hparam_per_iter <- TSR$zeros(c(length(self$hparam_labels), nb_sampling_iter))
            self$sum_beta_est <- TSR$zeros(covariates$shape)
            self$sum_y_est <- TSR$zeros(y$shape)
            self$total_elapsed_time <- 0

            self$y <- y
            self$omega <- omega
            self$covariates <- covariates
            self$formula <- formula
            self$rank_decomp <- rank_decomp
            self$spatial_kernel <- spatial_kernel
            self$temporal_kernel <- temporal_kernel
            self$nb_burn_in_iter <- nb_burn_in_iter
            self$nb_sampling_iter <- nb_sampling_iter

            # Set initial timer value to calculate iterations' processing time
            self$last_time_stamp <- Sys.time()
        },

        #~ @description Collect current iteration values inside the historical data tensor
        #~list. Note that errors have already been calculated before tau sampling.
        collect_iter_samples = function(iter, tau_value) {
            elapsed_time <- private$get_elapsed_time()

            if (iter > self$nb_burn_in_iter) {
                self$sum_beta_est <- self$sum_beta_est + self$beta_estimates
                self$sum_y_est <- self$sum_y_est + self$y_estimates

                # Collect hyperparameters
                s_iter <- iter - self$nb_burn_in_iter
                s_params <- Filter(function(p) !p$is_fixed, self$spatial_kernel$parameters)
                t_params <- Filter(function(p) !p$is_fixed, self$temporal_kernel$parameters)
                self$hparam_per_iter[1, s_iter] <- tau_value
                self$hparam_per_iter[2:(1 + length(s_params)), s_iter] <- sapply(
                    s_params, function(p) p$value
                )
                self$hparam_per_iter[(2 + length(s_params)):length(self$hparam_labels), s_iter] <- sapply(
                    t_params, function(p) p$value
                )
            }

            total_logged_params <- c(
                list(
                    iter = iter,
                    is_burn_in = ifelse(iter <= self$nb_burn_in_iter, 1, 0),
                    elapsed_time = elapsed_time
                ),
                self$error_metrics
            )

            for (p_name in names(total_logged_params)) {
                self$logged_params_map[[p_name]] <- c(
                    self$logged_params_map[[p_name]], total_logged_params[[p_name]]
                )
            }

            private$print_iter_result(iter, elapsed_time)
        },

        set_error_metrics = function() {
            nb_observ <- self$omega$sum()
            err_matrix <- (self$y_estimates - self$y) * self$omega
            total_sq_error <- err_matrix$norm() ** 2
            mae <- err_matrix$abs()$sum() / nb_observ
            rmse <- (total_sq_error / nb_observ)$sqrt()
            self$total_sq_error <- as.numeric(total_sq_error$cpu())
            self$error_metrics <- list(
                MAE = as.double(mae$cpu()),
                RMSE = as.double(rmse$cpu())
            )
        },

        set_y_and_beta_estimates = function(decomp_tensors_map, iter) {
            # Calculate Coefficient Estimation
            if (iter > self$nb_burn_in_iter) {
                iter_indx <- iter - self$nb_burn_in_iter
                self$spatial_decomp_per_iter[, , iter_indx] <- decomp_tensors_map[['spatial_decomp']]
                self$temporal_decomp_per_iter[, , iter_indx] <- decomp_tensors_map[['temporal_decomp']]
                self$covs_decomp_per_iter[, , iter_indx] <- decomp_tensors_map[['covs_decomp']]
            }

            self$beta_estimates <- torch::torch_einsum(
                'im,jm,km->ijk', c(
                    decomp_tensors_map[['spatial_decomp']],
                    decomp_tensors_map[['temporal_decomp']],
                    decomp_tensors_map[['covs_decomp']]
                )
            )
            self$y_estimates <- torch::torch_einsum('ijk,ijk->ij', c(self$covariates, self$beta_estimates))
        },

        log_final_iter_results = function() {
            self$beta_estimates <- self$sum_beta_est / self$nb_sampling_iter
            self$y_estimates <- self$sum_y_est / self$nb_sampling_iter
            beta_covariates_summary <- private$create_distrib_values_summary(
                self$beta_estimates$reshape(c(-1, length(self$feature_labels)))$cpu(), dim = 1
            )
            self$beta_covariates_summary_df <- cbind(
                data.table(self$feature_labels),
                data.table(as.matrix(beta_covariates_summary$t()))
            )
            setnames(self$beta_covariates_summary_df, c('feature', self$moment_metrics, self$quantile_metrics))
            y_beta_index <- CJ(location=self$spatial_labels, time=self$temporal_labels)
            self$y_estimates_df <- cbind(
                y_beta_index,
                data.table(as.matrix(self$y_estimates$cpu()$flatten()))
            )
            setnames(self$y_estimates_df, c('location', 'time', 'y_est'))
            self$beta_estimates_df <- cbind(
                y_beta_index,
                data.table(as.matrix(self$beta_estimates$reshape(c(-1, length(self$feature_labels)))$cpu()))
            )
            setnames(self$beta_estimates_df, c('location', 'time', self$feature_labels))
            self$hyperparameters_per_iter_df <- cbind(
                data.table(1:self$nb_sampling_iter),
                data.table(as.matrix(self$hparam_per_iter$t()$cpu()))
            )
            setnames(self$hyperparameters_per_iter_df, c('iter', self$hparam_labels))
            self$set_error_metrics()
            private$print_iter_result('TOTAL', self$total_elapsed_time)
        }
    ),

    private = list(
        get_elapsed_time = function() {
            iter_elapsed_time <- Sys.time() - self$last_time_stamp
            self$total_elapsed_time <- self$total_elapsed_time + iter_elapsed_time
            self$last_time_stamp <- Sys.time()
            return(iter_elapsed_time)
        },

        print_iter_result = function(iter, elapsed_time) {
            formatted_err_vals <- sprintf('%7.4f', unlist(self$error_metrics))
            formatted_errors <- paste0(names(self$error_metrics), ': ', formatted_err_vals, collapse = ' | ')
            iter_format <- paste('Iter', ifelse(is.character(iter), '%s', '%-5d'))
            result_items <- c(
                sprintf(iter_format, iter),
                sprintf('Elapsed Time: %8.2fs', elapsed_time),
                formatted_errors
            )
            print(paste0(result_items, collapse = ' | '))
        },

        # Create a summary for a given tensor of beta values across a given dimension
        # for the metrics set in the class.
        # Args:
        #    values (Tensor): Values to summarize
        #    dim (Integer): Dimension of the tensor we want to summarize. If None,
        #        we want to summarize the whole tensor and flatten it. Defaults to None.
        # Returns:
        #    Tensor: A tensor with summaries for the given beta values
        create_distrib_values_summary = function(values, dim) {
            all_metrics <- c(self$moment_metrics, self$quantile_metrics)
            summary_shape <- c(length(all_metrics))
            if (!is.null(dim)) {
                beta_val_shape <- values$shape
                summary_shape <- c(
                    summary_shape,
                    beta_val_shape[-dim]
                )
            }
            beta_summaries <- TSR$zeros(summary_shape)
            # In the advent of having no values, we return the empty tensor
            if (values$numel() == 0) {
                return(beta_summaries)
            }
            # Dimension for moment calculations are a bit different than for quantile
            moment_dim <- ifelse(is.null(dim), c(), dim)
            beta_summaries[1] <- values$mean(dim = moment_dim)
            beta_summaries[2] <- values$std(dim = moment_dim)
            beta_summaries[(length(self$moment_metrics) + 1):length(all_metrics)] <- torch::torch_quantile(
                values, TSR$tensor(self$quantile_values), dim = dim
            )
            return(beta_summaries)
        }
    )
)
