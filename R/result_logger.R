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
        LINE_NCHAR = 70,
        TAB_STR = '  ',
        MAIN_COL_WIDTH = 18,
        OTHER_COL_WIDTH = 8,
        OTHER_COL_FMT = '%8.3f',
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
                self$beta_estimates$reshape(c(-1, length(self$feature_labels))), dim = 1
            )
            self$beta_covariates_summary_df <- cbind(
                data.table(self$feature_labels),
                data.table(as.matrix(beta_covariates_summary$t()$cpu()))
            )
            setnames(self$beta_covariates_summary_df, c('feature', self$moment_metrics, self$quantile_metrics))
            y_beta_index <- CJ(location = self$spatial_labels, time = self$temporal_labels)
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
        },

        # Print a summary of the BKTR regressor instance after MCMC sampling.
        summary = function() {
            line_sep <- strrep('=', self$LINE_NCHAR)
            summary_str <- c(
                '',
                line_sep,
                format('BKTR Regressor Summary', width = self$LINE_NCHAR, justify = 'centre'),
                line_sep,
                private$get_formula_str(),
                '',
                sprintf('Burn-in iterations: %i', self$nb_burn_in_iter),
                sprintf('Sampling iterations: %i', self$nb_sampling_iter),
                sprintf('Rank decomposition: %i', self$rank_decomp),
                sprintf('Nb Spatial Locations: %i', length(self$spatial_labels)),
                sprintf('Nb Temporal Points: %i', length(self$temporal_labels)),
                sprintf('Nb Covariates: %i', length(self$feature_labels)),
                line_sep,
                'In Sample Errors:',
                sprintf('%sRMSE: %.3f', self$TAB_STR, self$error_metrics['RMSE']),
                sprintf('%sMAE: %.3f', self$TAB_STR, self$error_metrics['MAE']),
                sprintf('Computation time: %.2fs.', self$total_elapsed_time),
                line_sep,
                '-- Spatial Kernel --',
                private$kernel_summary(self$spatial_kernel, 'spatial'),
                '',
                '-- Temporal Kernel --',
                private$kernel_summary(self$temporal_kernel, 'temporal'),
                line_sep,
                private$beta_summary(),
                line_sep,
                ''
            )
            return(paste(summary_str, collapse = '\n'))
        },

        get_beta_summary_df = function(
            spatial_labels = NULL,
            temporal_labels = NULL,
            feature_labels = NULL
        ) {
            spatial_labs <- if (is.null(spatial_labels)) self$spatial_labels else spatial_labels
            temporal_labs <- if (is.null(temporal_labels)) self$temporal_labels else temporal_labels
            feature_labs <- if (is.null(feature_labels)) self$feature_labels else feature_labels
            iteration_betas <- self$get_iteration_betas_tensor(spatial_labs, temporal_labs, feature_labs)
            beta_summary <- private$create_distrib_values_summary(iteration_betas, dim = 2)$t()$cpu()

            index_cols <- CJ(location = spatial_labs, time = temporal_labs, feature = feature_labs)
            df <- cbind(
                index_cols,
                data.table(as.matrix(beta_summary))
            )
            setnames(df, c('location', 'time', 'feature', self$moment_metrics, self$quantile_metrics))
            return(df)
        },

        get_iteration_betas_tensor = function(spatial_labels, temporal_labels, feature_labels) {
            spatial_indexes <- get_label_indexes(spatial_labels, self$spatial_labels, 'spatial')
            temporal_indexes <- get_label_indexes(temporal_labels, self$temporal_labels, 'temporal')
            feature_indexes <- get_label_indexes(feature_labels, self$feature_labels, 'feature')
            betas_per_iterations <- torch::torch_einsum(
                'sri,tri,cri->stci',
                c(
                    self$spatial_decomp_per_iter[spatial_indexes, , , drop = FALSE],
                    self$temporal_decomp_per_iter[temporal_indexes, , , drop = FALSE],
                    self$covs_decomp_per_iter[feature_indexes, , , drop = FALSE]
                )
            )
            return(betas_per_iterations$reshape(c(-1, self$nb_sampling_iter)))
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

        #~ @description Create a summary for a given tensor of beta values across a given dimension
        #~ for the metrics set in the class.
        #~ @param values Tensor: Values to summarize
        #~ @param dim Integer: Dimension of the tensor we want to summarize. If NULL,
        #~     we want to summarize the whole tensor and flatten it. Defaults to NULL.
        #~ @return A tensor with summaries for the given beta values
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
        },

        get_formula_str = function() {
            formula_str <- paste(self$formula[2], self$formula[3], sep = " ~ ")
            formula_str <- paste('Formula:', formula_str)
            f_wrap <- strwrap(formula_str, width = self$LINE_NCHAR)
            return(paste(f_wrap, collapse = paste0('\n', self$TAB_STR)))
        },

        #~ @description Get a string representation of a given kernel. Since the kernel can be
        #~ composed, this function needs to be recursive.
        #~ @param kernel Kernel: The kernel we want summarize.
        #~ @param kernel_type ('spatial' or 'temporal'): The type of kernel.
        #~ @param indent_count Integer: Indentation level (related to the depth of composition).
        #~    Defaults to 0.
        #~ @return A string representation of the kernel. Containing the name of the kernel,
        #~   the estimated parameters distribution and the fixed parameters.
        kernel_summary = function(kernel, kernel_type, indent_count = 0) {
            params <- kernel$parameters
            if (class(kernel)[1] %in% c('KernelAddComposed', 'KernelMulComposed')) {
                new_ind_nb <- indent_count + 1
                op_str <- capitalize_str(kernel$composition_operation)
                kernel_elems <- c(
                    paste0('Composed Kernel (', op_str, ')'),
                    private$kernel_summary(kernel$left_kernel, kernel_type, new_ind_nb),
                    paste0(self$TAB_STR, ifelse(op_str == 'Add', '+', '*')),
                    private$kernel_summary(kernel$right_kernel, kernel_type, new_ind_nb)
                )
            } else {
                fixed_params <- params[sapply(params, function(x) x$is_fixed)]
                sampled_params <- params[sapply(params, function(x) !x$is_fixed)]
                sampled_par_indexes <- sapply(
                    sampled_params,
                    function(x) which(self$hparam_labels == paste0(capitalize_str(kernel_type), ' - ', x$full_name))
                )
                sampled_par_tsr <- self$hparam_per_iter[sampled_par_indexes, drop = FALSE]
                sampled_par_summary <- private$create_distrib_values_summary(sampled_par_tsr, dim = 2)
                sampled_par_df <- cbind(
                    data.table(sapply(sampled_params, function(x) x$name)),
                    data.table(as.matrix(sampled_par_summary$t()$cpu()))
                )
                setnames(sampled_par_df, c('Parameter', self$moment_metrics, self$quantile_metrics))
                out_cols <- c('Parameter', self$DISTRIB_COLS)
                sampled_par_df <- sampled_par_df[, ..out_cols]
                sampled_par_strs <- private$get_formatted_df_rows(sampled_par_df)
                fixed_par_strs <- sapply(
                    fixed_params,
                    function(x) sprintf('%-20s   Fixed Value: %.3f', x$name, x$value)
                )
                kernel_elems <- c(
                    kernel$name,
                    'Parameter(s):',
                    sampled_par_strs,
                    fixed_par_strs
                )
            }
            kernel_elems <- paste0(rep(self$TAB_STR, indent_count), kernel_elems)
            return(paste(kernel_elems, collapse = '\n'))
        },

        #~ @description Get a string representation of the beta estimates aggregated per
        #~    covariates. (This shows the distribution of the beta hats per covariates)
        #~ @return A string representation of the beta estimates.
        beta_summary = function() {
            beta_est_cols <- c('feature', self$MAIN_SUMMARY_COLS)
            distrib_df <- self$beta_covariates_summary_df[, ..beta_est_cols]
            beta_distrib_str_rows <- private$get_formatted_df_rows(distrib_df)
            beta_summary_strs <- c(
                'Beta Estimates Summary (Aggregated Per Covariates)',
                '',
                beta_distrib_str_rows
            )
            return(paste(beta_summary_strs, collapse = '\n'))
        },

        format_df_row = function(df_row) {
            df_cols <- c(
                format(trunc_str(df_row[, 1], self$MAIN_COL_WIDTH), width = self$MAIN_COL_WIDTH, justify = 'left'),
                sapply(unlist(df_row[, -1]), function(x) sprintf(self$OTHER_COL_FMT, x))
            )
            return(paste(df_cols, collapse = ' '))
        },

        # Format a dataframe to be printed in a table.
        get_formatted_df_rows = function(df) {
            df_header_cols <- c(
                strrep(' ', self$MAIN_COL_WIDTH),
                sapply(colnames(df)[-1], function(x) format(x, width = self$OTHER_COL_WIDTH, justify = 'right'))
            )
            df_header_row <- paste(df_header_cols, collapse = ' ')
            df_str_rows <- by(df, seq_len(nrow(df)), private$format_df_row)
            return(c(df_header_row, df_str_rows))
        }
    )
)
