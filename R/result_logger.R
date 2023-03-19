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
        nb_iter = NULL,
        nb_burn_in_iter = NULL,
        sampled_beta_indexes = NULL,
        sampled_y_indexes = NULL,
        logged_params_map = NULL,
        beta_estimates = NULL,
        beta_stdev = NULL,
        y_estimates = NULL,
        total_elapsed_time = NULL,
        sum_beta_est = NULL,
        sum_sq_beta_est = NULL,
        sum_y_est = NULL,
        last_time_stamp = NULL,
        export_path = NULL,
        export_suffix = NULL,
        error_metrics = NULL,

        initialize = function(
            y,
            omega,
            covariates,
            nb_iter,
            nb_burn_in_iter,
            results_export_dir = NULL,
            results_export_suffix = NULL,
            sampled_beta_indexes = c(),
            sampled_y_indexes = c()
        ) {
            self$logged_params_map <- list()

            if (is.null(results_export_dir) && !is.null(results_export_suffix)) {
                torch:::value_error(
                    'Cannot set a suffix for the export file if no export directory is provided.'
                )
            }
            if (!is.null(results_export_dir) && !dir.exists(results_export_dir)) {
                torch:::value_error(
                    sprintf('The selected export directory `%s` does not exists.', export_dir)
                )
            }
            self$export_path <- results_export_suffix
            self$export_suffix <- results_export_suffix

            self$sum_beta_est <- tsr$zeros(covariates$shape)
            self$sum_sq_beta_est <- tsr$zeros(covariates$shape)
            self$sum_y_est <- tsr$zeros(y$shape)
            self$total_elapsed_time <- 0

            self$y <- y
            self$omega <- omega
            self$covariates <- covariates
            self$nb_iter <- nb_iter
            self$nb_burn_in_iter <- nb_burn_in_iter

            # Indexes that need to be logged in the result output
            self$sampled_beta_indexes <- sampled_beta_indexes
            self$sampled_y_indexes <- sampled_y_indexes

            # Set initial timer value to calculate iterations' processing time
            self$last_time_stamp <- Sys.time()
        },

        #~ @description Collect current iteration values inside the historical data tensor
        #~list. Note that errors have already been calculated before tau sampling.
        collect_iter_samples = function(iter, iter_logged_params) {
            elapsed_time_dict <- private$get_elapsed_time_dict()
            y_beta_sampled_values <- private$get_y_and_beta_sampled_values(
                self$y_estimates, self$beta_estimates
            )

            if (iter > self$nb_burn_in_iter) {
                self$sum_beta_est <- self$sum_beta_est + self$beta_estimates
                self$sum_sq_beta_est <- self$sum_sq_beta_est + self$beta_estimates ** 2
                self$sum_y_est <- self$sum_y_est + self$y_estimates
            }

            all_par_iter <- c(
                list(
                    iter = iter,
                    is_burn_in = ifelse(iter <= self$nb_burn_in_iter, 1, 0)
                ),
                iter_logged_params,
                self$error_metrics,
                y_beta_sampled_values,
                elapsed_time_dict
            )
            for (p_name in names(all_par_iter)) {
                self$logged_params_map[[p_name]] <- c(
                    self$logged_params_map[[p_name]], all_par_iter[[p_name]]
                )
            }

            private$print_iter_result(iter, append(elapsed_time_dict, self$error_metrics))
        },

        set_error_metrics = function() {
            nb_observ <- self$omega$sum()
            err_matrix <- (self$y_estimates - self$y) * self$omega
            total_sq_error <- err_matrix$norm() ** 2
            mae <- err_matrix$abs()$sum() / nb_observ
            rmse <- (total_sq_error / nb_observ)$sqrt()
            self$error_metrics <- list(
                total_sq_error = as.double(total_sq_error$cpu()),
                mae = as.double(mae$cpu()),
                rmse = as.double(rmse$cpu())
            )
            return(self$error_metrics)
        },

        set_y_and_beta_estimates = function(decomp_tensors_map) {
            # Calculate Coefficient Estimation
            self$beta_estimates <- torch::torch_einsum(
                'im,jm,km->ijk', c(
                    decomp_tensors_map[['spatial_decomp']],
                    decomp_tensors_map[['temporal_decomp']],
                    decomp_tensors_map[['covs_decomp']]
                )
            )
            self$y_estimates <- torch::torch_einsum('ijk,ijk->ij', c(self$covariates, self$beta_estimates))
        },

        log_iter_results = function() {
            iter_results_df <- as.data.frame(do.call(cbind, self$logged_params_map))
            avg_estimates <- private$get_and_print_avg_estimates()
            if (!is.null(self$export_path)) {
                write.csv(iter_results_df, private$get_file_name('iter_results'), row.names = FALSE)
                y_est <- as.array(avg_estimates[['y_est']]$cpu()$flatten())
                beta_est <- as.array(avg_estimates[['beta_est']]$cpu()$flatten())
                write.table(y_est, private$get_file_name('y_estimates'), row.names = FALSE, col.names = FALSE)
                write.table(beta_est, private$get_file_name('beta_estimates'), row.names = FALSE, col.names = FALSE)
            }
            return(avg_estimates)
        }
    ),

    private = list(
        get_and_print_avg_estimates = function() {
            nb_sampling_iter <- self$nb_iter - self$nb_burn_in_iter
            self$beta_estimates <- self$sum_beta_est / nb_sampling_iter
            self$beta_stdev <- (
                (self$sum_sq_beta_est / nb_sampling_iter) - self$beta_estimates ** 2
            )$sqrt()
            self$y_estimates <- self$sum_y_est / nb_sampling_iter
            error_metrics <- self$set_error_metrics()
            iters_summary_dict <- c(
                list(elapsed_time = self$total_elapsed_time),
                error_metrics
            )
            private$print_iter_result('TOTAL', iters_summary_dict)
            avg_estimates <- c(
                error_metrics,
                list(y_est = self$y_estimates, beta_est = self$beta_estimates)
            )
            return(avg_estimates)
        },

        get_elapsed_time_dict = function() {
            iter_elapsed_time <- Sys.time() - self$last_time_stamp
            self$total_elapsed_time <- self$total_elapsed_time + iter_elapsed_time
            self$last_time_stamp <- Sys.time()
            return(list(elapsed_time = iter_elapsed_time))
        },

        get_y_and_beta_sampled_values = function(y_estimates, beta_estimates) {
            y_beta_sampled <- list()
            flat_y_est <- y_estimates$cpu()$flatten()
            flat_beta_est <- beta_estimates$cpu()$flatten()

            for (idx in self$sampled_y_indexes) {
                y_beta_sampled[paste0('y_', idx)] <- as.double(flat_y_est[idx])
            }
            for (idx in self$sampled_beta_indexes) {
                y_beta_sampled[paste0('beta_', idx)] <- as.double(flat_beta_est[idx])
            }

            return(y_beta_sampled)
        },

        print_iter_result = function(iter, result_dict) {
            iter_result_str <- sprintf('Results for iter %5s', toString(iter))
            for (i in seq_len(length(result_dict))) {
                result_item_str <- sprintf('%s is %.4f', names(result_dict)[[i]], result_dict[[i]])
                iter_result_str <- paste(iter_result_str, '||', result_item_str)
            }
            print(iter_result_str)
        },

        get_file_name = function(export_type_name) {
            time_str <- format(Sys.time(), format = "%Y%m%d_%H%M")
            file_name <- paste0(export_type_name, '_', time_str)
            if (!is.null(self$export_suffix)) {
                file_name <- paste0(file_name, '_', self$export_suffix)
            }
            return(file.path(self$export_dir, paste0(file_name, '.csv')))
        }
    )
)
