#' @import ggplot2
#' @import ggmap
#' @importFrom stats reshape

#' @title Print ggplot figure
#' @description Utility function to print a ggplot figure
#' @param fig ggplot: ggplot figure to print
#' @param fig_width Numeric: Figure width.
#' @param fig_height Numeric: Figure height.
#' @param fig_resolution Numeric: Figure resolution PPI.
#' @return NULL
#'
#' @noRd
print_ggplot_fig <- function(fig, fig_width, fig_height, fig_resolution) {
    # The following options are mainly for notebooks rendering (like Colab)
    options(repr.plot.width = fig_width, repr.plot.height = fig_height, repr.plot.res = fig_resolution)
    # The following sizes are for RStudio and other rendering
    print(fig, vp = grid::viewport(width = unit(fig_width, 'inches'), height = unit(fig_height, 'inches')))
}

#' @title Plot Temporal Beta Coefficients
#' @description Create a plot of the beta values through time for a given
#' spatial point and a set of feature labels.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param plot_feature_labels Array: Array of feature labels to plot.
#' @param spatial_point_label String: Spatial point label to plot.
#' @param date_format String: Format of the date to use in bktr dataframes for the time.
#'   Defaults to '\%Y-\%m-\%d'.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Numeric: Figure width when figure is shown. Defaults to 8.5.
#' @param fig_height Numeric: Figure height when figure is shown. Defaults to 5.5.
#' @param fig_resolution Numeric: Figure resolution PPI. Defaults to 200.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @examplesIf torch::torch_is_installed()
#' # Launch MCMC sampling on a light version of the BIXI dataset
#' bixi_data <- BixiData$new(is_light = TRUE)
#' bktr_regressor <- BKTRRegressor$new(
#'   data_df <- bixi_data$data_df,
#'   spatial_positions_df = bixi_data$spatial_positions_df,
#'   temporal_positions_df = bixi_data$temporal_positions_df,
#'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
#' bktr_regressor$mcmc_sampling()
#'
#' # Plot temporal beta coefficients for the first station and the two features
#' plot_temporal_betas(
#'   bktr_regressor,
#'   plot_feature_labels = c('mean_temp_c', 'area_park'),
#'   spatial_point_label = bixi_data$spatial_positions_df$location[1])
#'
#' @export
plot_temporal_betas <- function(
    bktr_reg,
    plot_feature_labels,
    spatial_point_label,
    date_format = '%Y-%m-%d',
    show_figure = TRUE,
    fig_width = 8.5,
    fig_height = 5.5,
    fig_resolution = 200
) {
    if (!bktr_reg$has_completed_sampling) {
        stop('Plots can only be accessed after MCMC sampling.')
    }
    # Verify all labels are valid
    get_label_index_or_raise(spatial_point_label, bktr_reg$spatial_labels, 'spatial')
    sapply(plot_feature_labels, function(feature_label) {
        get_label_index_or_raise(feature_label, bktr_reg$feature_labels, 'feature')
    })
    beta_est_df <- bktr_reg$get_beta_summary_df(
        c(spatial_point_label),
        NULL,
        plot_feature_labels
    )
    plot_title <- paste('Location:', spatial_point_label)

    if (!is.null(date_format)) beta_est_df$time <- as.Date(beta_est_df$time, format = date_format)
    fig <- (
        ggplot(beta_est_df, aes(.data$time, .data$Mean, group = .data$feature, color = .data$feature))
        + geom_line()
        + geom_ribbon(
            aes(ymin = .data$Low2.5p, ymax = .data$Up97.5p, fill = .data$feature),
            alpha = 0.3,
            color = NA
        )
        + ggtitle(plot_title)
        + theme_bw()
        + ylab('Beta Value')
        + xlab('Time')
        + labs(fill = 'Feature', color = 'Feature')

    )
    if (!show_figure) {
        return(fig)
    }
    print_ggplot_fig(fig, fig_width, fig_height, fig_resolution)
}

#' @title Plot Spatial Beta Coefficients
#' @description Create a plot of beta values through space for a given
#' temporal point and a set of feature labels. We use ggmap under the hood,
#' so you need to provide a Google or Stadia API token to plot on a map.
#' See: https://cran.r-project.org/web/packages/ggmap/readme/README.html for more details
#' on how to get an API token.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param plot_feature_labels Array: Array of feature labels to plot.
#' @param temporal_point_label String: Temporal point label to plot.
#' @param google_token String or NULL: Google API token to use for the geographic map. Defaults to NULL.
#' @param stadia_token String or NULL: Stadia API token to use for the geographic map. Defaults to NULL.
#' @param nb_cols Integer: The number of columns to use in the facet grid.
#' @param use_dark_mode Boolean: Whether to use a dark mode for the geographic map or not. Defaults to TRUE.
#' @param zoom Integer: Zoom level for the geographic map. Defaults to 11.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Numeric: Figure width when figure is shown. Defaults to 8.5.
#' @param fig_height Numeric: Figure height when figure is shown. Defaults to 5.5.
#' @param fig_resolution Numeric: Figure resolution PPI. Defaults to 200.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @examplesIf FALSE
#' # Launch MCMC sampling on a light version of the BIXI dataset
#' bixi_data <- BixiData$new(is_light = TRUE)
#' bktr_regressor <- BKTRRegressor$new(
#'   data_df <- bixi_data$data_df,
#'   spatial_positions_df = bixi_data$spatial_positions_df,
#'   temporal_positions_df = bixi_data$temporal_positions_df,
#'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
#' bktr_regressor$mcmc_sampling()
#'
#' # Plot spatial beta coefficients for the first time point and the two features
#' # Use your own Google API token instead of 'GOOGLE_API_TOKEN'
#' plot_spatial_betas(
#'   bktr_regressor,
#'   plot_feature_labels = c('mean_temp_c', 'area_park'),
#'   temporal_point_label = bixi_data$temporal_positions_df$time[1],
#'   google_token = 'GOOGLE_API_TOKEN')
#'
#' # We can also use light mode and plot the maps side by side
#' # Use your own Stadia API token instead of 'STADIA_API_TOKEN'
#' plot_spatial_betas(
#'   bktr_regressor,
#'   plot_feature_labels = c('mean_temp_c', 'area_park', 'total_precip_mm'),
#'   temporal_point_label = bixi_data$temporal_positions_df$time[10],
#'   use_dark_mode = FALSE, nb_cols = 3, stadia_token = 'STADIA_API_TOKEN')
#'
#' @export
plot_spatial_betas <- function(
    bktr_reg,
    plot_feature_labels,
    temporal_point_label,
    nb_cols = 1,
    use_dark_mode = TRUE,
    show_figure = TRUE,
    zoom = 11,
    google_token = NULL,
    stadia_token = NULL,
    fig_width = 8.5,
    fig_height = 5.5,
    fig_resolution = 200
) {
    if (!bktr_reg$has_completed_sampling) {
        stop('Plots can only be accessed after MCMC sampling.')
    }
    # Verify all labels are valid
    point_label_index <- get_label_index_or_raise(temporal_point_label, bktr_reg$temporal_labels, 'temporal')
    sapply(plot_feature_labels, function(feature_label) {
        get_label_index_or_raise(feature_label, bktr_reg$feature_labels, 'feature')
    })

    # Get only spatial estimates
    beta_df <- bktr_reg$beta_estimates
    out_col_names <- c('location', plot_feature_labels)
    beta_df <- beta_df[beta_df$time == temporal_point_label, out_col_names, with = FALSE]

    coords_projector <- bktr_reg$geo_coords_projector
    is_map <- !is.null(coords_projector)
    coord_df <- if (is_map) coords_projector$ini_df else bktr_reg$spatial_positions_df
    if (ncol(coord_df[, -1]) != 2) {
        stop('Spatial coordinates must be 2 dimensions to be plotted.')
    }

    full_df <- beta_df[coord_df, on = 'location', nomatch = NULL]
    full_df <- reshape(
        full_df,
        direction = 'long',
        idvar = colnames(coord_df),
        v.names = 'value',
        timevar = 'feature',
        varying = plot_feature_labels,
        times = plot_feature_labels
    )

    plot_title <- paste0('Estimated Beta at Time Point : ', temporal_point_label)
    longitude <- latitude <- value <- NULL # Used for CRAN global binding checks
    if (is_map) {
        if (!is.null(google_token)) {
            ggmap::register_google(google_token)
            map_source <- 'google'
            map_type <- 'roadmap'
            map_color <- ifelse(use_dark_mode, 'bw', 'color')
        } else if (!is.null(stadia_token)) {
            ggmap::register_stadiamaps(stadia_token)
            map_source <- 'stadia'
            map_type <- ifelse(use_dark_mode, 'stamen_toner', 'stamen_toner_lite')
            map_color <- 'bw'
        } else {
            stop('You must provide a Google or Stadia API token to plot on a map.')
        }
        fig <- (
            ggmap::qmplot(
              x = longitude,
              y = latitude,
              color = value,
              data = full_df,
              source = map_source,
              maptype = map_type,
              mapcolor = map_color,
              zoom = zoom
            )
            + facet_wrap(~feature, ncol = nb_cols)
            + theme_bw()
            + theme(
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            )
            + labs(x = NULL, y = NULL)
            + scale_color_viridis_c()
            + ggtitle(plot_title)
        )
    } else {
        x_col_name <- colnames(full_df)[2]
        y_col_name <- colnames(full_df)[3]
        fig <- (
            ggplot(full_df, aes(x = .data[[x_col_name]], y = .data[[y_col_name]], color = .data$value))
            + geom_point()
            + facet_wrap(~feature, ncol = nb_cols)
            + theme_bw()
            + scale_color_viridis_c()
            + ggtitle(plot_title)
        )
    }
    if (!show_figure) {
        return(fig)
    }
    print_ggplot_fig(fig, fig_width, fig_height, fig_resolution)
}



#' @title Plot Beta Coefficients Distribution
#' @description Plot the distribution of beta values for a given list of labels.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param labels_list List: List of labels tuple (spatial, temporal, feature) for
#'  which to plot the beta distribution through iterations
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Integer: Figure width. Defaults to 9.
#' @param fig_height Integer: Figure height. Defaults to 6.
#' @param fig_resolution Numeric: Figure resolution PPI. Defaults to 200.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @examplesIf torch::torch_is_installed()
#' # Launch MCMC sampling on a light version of the BIXI dataset
#' bixi_data <- BixiData$new(is_light = TRUE)
#' bktr_regressor <- BKTRRegressor$new(
#'   data_df <- bixi_data$data_df,
#'   spatial_positions_df = bixi_data$spatial_positions_df,
#'   temporal_positions_df = bixi_data$temporal_positions_df,
#'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
#' bktr_regressor$mcmc_sampling()
#'
#' # Plot temporal beta coefficients for the first station and the first feature
#' spa_lab <- bixi_data$spatial_positions_df$location[3]
#' plot_beta_dists(
#'   bktr_regressor,
#'   labels_list = list(
#'      c(spa_lab, '2019-04-15', 'area_park'),
#'      c(spa_lab, '2019-04-16', 'area_park'),
#'      c(spa_lab, '2019-04-16', 'mean_temp_c')
#'   ),
#' )
#'
#' @export
plot_beta_dists <- function(
    bktr_reg,
    labels_list,
    show_figure = TRUE,
    fig_width = 9,
    fig_height = 6,
    fig_resolution = 200
) {
    if (!bktr_reg$has_completed_sampling) {
        stop('Plots can only be accessed after MCMC sampling.')
    }
    plot_title <- 'Posterior distribution of beta values per given spatial point, temporal point and feature'
    df <- data.table(sapply(labels_list, function(x) bktr_reg$get_iterations_betas(x[1], x[2], x[3])))
    col_names <- sapply(labels_list, function(x) paste(x, collapse = '\n'))
    setnames(df, col_names)
    df <- reshape(
        df,
        direction = 'long',
        v.names = 'value',
        timevar = 'labels',
        varying = col_names,
        times = col_names
    )
    fig <- (
        ggplot(df, aes(x = .data$labels, y = .data$value, fill = .data$labels))
        + geom_violin(trim = FALSE)
        + ggtitle(plot_title)
        + ylab('Beta Value')
        + xlab('Labels')
        + labs(fill = 'Labels', color = 'Labels')
    )
    if (!show_figure) {
        return(fig)
    }
    print_ggplot_fig(fig, fig_width, fig_height, fig_resolution)
}


#' @title Plot Beta Coefficients Distribution Regrouped by Covariates
#' @description Plot the distribution of beta estimates regrouped by covariates.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param feature_labels Array or NULL: Array of feature labels for
#'   which to plot the beta estimates distribution. If NULL plot for all features.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Integer: Figure width. Defaults to 9.
#' @param fig_height Integer: Figure height. Defaults to 6.
#' @param fig_resolution Numeric: Figure resolution PPI. Defaults to 200.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @examplesIf torch::torch_is_installed()
#' # Launch MCMC sampling on a light version of the BIXI dataset
#' bixi_data <- BixiData$new(is_light = TRUE)
#' bktr_regressor <- BKTRRegressor$new(
#'   formula = 'nb_departure ~ 1 + area_park + mean_temp_c + total_precip_mm',
#'   data_df <- bixi_data$data_df,
#'   spatial_positions_df = bixi_data$spatial_positions_df,
#'   temporal_positions_df = bixi_data$temporal_positions_df,
#'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
#' bktr_regressor$mcmc_sampling()
#'
#' # Plot beta estimates distribution for all features
#' plot_covariates_beta_dists(bktr_regressor)
#' # Or plot for a subset of features
#' plot_covariates_beta_dists(bktr_regressor, c('area_park', 'mean_temp_c'))
#'
#' @export
plot_covariates_beta_dists <- function(
    bktr_reg,
    feature_labels = NULL,
    show_figure = TRUE,
    fig_width = 9,
    fig_height = 6,
    fig_resolution = 200
) {
    if (!bktr_reg$has_completed_sampling) {
        stop('Plots can only be accessed after MCMC sampling.')
    }

    if (is.null(feature_labels)) {
        feature_labels <- bktr_reg$feature_labels
    } else {
        sapply(feature_labels, function(feature_label) {
            get_label_index_or_raise(feature_label, bktr_reg$feature_labels, 'feature')
        })
    }

    plot_title <- 'Distribution of beta estimates by feature across time and space'
    full_df <- reshape(
        bktr_reg$beta_estimates[, feature_labels, with = FALSE],
        direction = 'long',
        v.names = 'value',
        timevar = 'feature',
        varying = feature_labels,
        times = feature_labels
    )
    fig <- (
        ggplot(full_df, aes(x = .data$feature, y = .data$value, fill = .data$feature))
        + geom_violin()
        + ggtitle(plot_title)
        + ylab('Beta Value')
        + xlab('Feature')
        + labs(fill = 'Feature', color = 'Feature')
    )
    if (!show_figure) {
        return(fig)
    }
    print_ggplot_fig(fig, fig_width, fig_height, fig_resolution)
}

#' @title Plot Hyperparameters Distributions
#' @description Plot the distribution of hyperparameters through iterations
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param hyperparameters Array or NULL: Array of hyperparameters to plot.
#'     If NULL, plot all hyperparameters. Defaults to NULL.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Integer: Figure width. Defaults to 9.
#' @param fig_height Integer: Figure height. Defaults to 6.
#' @param fig_resolution Numeric: Figure resolution PPI. Defaults to 200.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @examplesIf torch::torch_is_installed()
#' # Launch MCMC sampling on a light version of the BIXI dataset
#' bixi_data <- BixiData$new(is_light = TRUE)
#' k_matern <- KernelMatern$new()
#' k_periodic <- KernelPeriodic$new()
#' bktr_regressor <- BKTRRegressor$new(
#'   data_df <- bixi_data$data_df,
#'   spatial_kernel = k_matern,
#'   temporal_kernel = k_periodic,
#'   spatial_positions_df = bixi_data$spatial_positions_df,
#'   temporal_positions_df = bixi_data$temporal_positions_df,
#'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
#' bktr_regressor$mcmc_sampling()
#'
#' # Plot the distribution of all hyperparameters
#' plot_hyperparams_dists(bktr_regressor)
#'
#' # Plot the distribution of the spatial kernel hyperparameters
#' spa_par_name <- paste0('Spatial - ', k_matern$parameters[[1]]$full_name)
#' plot_hyperparams_dists(bktr_regressor, spa_par_name)
#'
#' # Plot the distribution of the temporal kernel hyperparameters
#' temp_par_names <- sapply(k_periodic$parameters, function(x) x$full_name)
#' temp_par_names <- paste0('Temporal - ', temp_par_names)
#' plot_hyperparams_dists(bktr_regressor, temp_par_names)
#'
#' @export
plot_hyperparams_dists <- function(
    bktr_reg,
    hyperparameters = NULL,
    show_figure = TRUE,
    fig_width = 9,
    fig_height = 6,
    fig_resolution = 200
) {
    if (!bktr_reg$has_completed_sampling) {
        stop('Plots can only be accessed after MCMC sampling.')
    }

    df <- bktr_reg$result_logger$hyperparameters_per_iter_df
    all_hparams <- colnames(df[, -1])
    hparams <- if (is.null(hyperparameters)) all_hparams else hyperparameters
    hparam_diff <- setdiff(hparams, all_hparams)
    if (length(hparam_diff) > 0) {
        formatted_available_params <- paste(all_hparams, collapse = ',\n\t')
        formatted_hparam_diff <- paste(hparam_diff, collapse = ', ')
        stop(sprintf(
            'Hyperparameter(s) %s not found. Available hyperparameters are:\n\t%s',
            formatted_hparam_diff, formatted_available_params
        ))
    }
    df <- reshape(
        df[, -1],
        direction = 'long',
        v.names = 'value',
        timevar = 'hyperparameter',
        varying = hparams,
        times = hparams
    )
    fig <- (
        ggplot(df, aes(x = .data$hyperparameter, y = .data$value, fill = .data$hyperparameter))
        + geom_violin(trim = FALSE)
        + ggtitle('Posterior Distribution of BKTR Hyperparameters')
        + ylab('Hyperparameter Value')
        + xlab('Hyperparameter')
        + labs(fill = 'Hyperparameter', color = 'Hyperparameter')
    )
    if (!show_figure) {
        return(fig)
    }
    print_ggplot_fig(fig, fig_width, fig_height, fig_resolution)
}

#' @title Plot Hyperparameters Traceplot
#' @description Plot the evolution of hyperparameters through iterations. (Traceplot)
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param hyperparameters Array or NULL: Array of hyperparameters to plot.
#'     If NULL, plot all hyperparameters. Defaults to NULL.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Integer: Figure width. Defaults to 9.
#' @param fig_height Integer: Figure height. Defaults to 5.5.
#' @param fig_resolution Numeric: Figure resolution PPI. Defaults to 200.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @examplesIf torch::torch_is_installed()
#' # Launch MCMC sampling on a light version of the BIXI dataset
#' bixi_data <- BixiData$new(is_light = TRUE)
#' k_matern <- KernelMatern$new()
#' k_periodic <- KernelPeriodic$new()
#' bktr_regressor <- BKTRRegressor$new(
#'   data_df <- bixi_data$data_df,
#'   spatial_kernel = k_matern,
#'   temporal_kernel = k_periodic,
#'   spatial_positions_df = bixi_data$spatial_positions_df,
#'   temporal_positions_df = bixi_data$temporal_positions_df,
#'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
#' bktr_regressor$mcmc_sampling()
#'
#' # Plot the traceplot of all hyperparameters
#' plot_hyperparams_traceplot(bktr_regressor)
#'
#' # Plot the traceplot of the spatial kernel hyperparameters
#' spa_par_name <- paste0('Spatial - ', k_matern$parameters[[1]]$full_name)
#' plot_hyperparams_traceplot(bktr_regressor, spa_par_name)
#'
#' # Plot the traceplot of the temporal kernel hyperparameters
#' temp_par_names <- sapply(k_periodic$parameters, function(x) x$full_name)
#' temp_par_names <- paste0('Temporal - ', temp_par_names)
#' plot_hyperparams_traceplot(bktr_regressor, temp_par_names)
#'
#' @export
plot_hyperparams_traceplot <- function(
    bktr_reg,
    hyperparameters = NULL,
    show_figure = TRUE,
    fig_width = 9,
    fig_height = 5.5,
    fig_resolution = 200
) {
    if (!bktr_reg$has_completed_sampling) {
        stop('Plots can only be accessed after MCMC sampling.')
    }

    df <- bktr_reg$result_logger$hyperparameters_per_iter_df
    all_hparams <- colnames(df[, -1])
    hparams <- if (is.null(hyperparameters)) all_hparams else hyperparameters
    hparam_diff <- setdiff(hparams, all_hparams)
    if (length(hparam_diff) > 0) {
        formatted_available_params <- paste(all_hparams, collapse = ',\n\t')
        formatted_hparam_diff <- paste(hparam_diff, collapse = ', ')
        stop(sprintf(
            'Hyperparameter(s) %s not found. Available hyperparameters are:\n\t%s',
            formatted_hparam_diff, formatted_available_params
        ))
    }
    df_cols <- c('iter', hparams)
    df <- reshape(
        df[, df_cols, with = FALSE],
        direction = 'long',
        v.names = 'value',
        idvar = 'iter',
        timevar = 'hyperparameter',
        varying = hparams,
        times = hparams
    )
    fig <- (
        ggplot(df, aes(.data$iter, .data$value, group = .data$hyperparameter, color = .data$hyperparameter))
        + geom_line()
        + ggtitle('Hyperparameter values through sampling iterations (Traceplot)')
        + theme_bw()
        + ylab('Hyperparameter Value')
        + xlab('Sampling Iter')
        + labs(fill = 'Hyperparameter', color = 'Hyperparameter')
        + theme(legend.position = 'bottom')
    )
    if (!show_figure) {
        return(fig)
    }
    print_ggplot_fig(fig, fig_width, fig_height, fig_resolution)
}


#' @title Plot Y Estimates
#' @description Plot y estimates vs observed y values.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Numeric: Figure width when figure is shown. Defaults to 5.
#' @param fig_height Numeric: Figure height when figure is shown. Defaults to 5.
#' @param fig_resolution Numeric: Figure resolution PPI when figure is shown. Defaults to 200.
#' @param fig_title String or NULL: Figure title if provided. Defaults to 'y estimates vs observed y values'
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @examplesIf torch::torch_is_installed()
#' # Launch MCMC sampling on a light version of the BIXI dataset
#' bixi_data <- BixiData$new(is_light = TRUE)
#' bktr_regressor <- BKTRRegressor$new(
#'   data_df <- bixi_data$data_df,
#'   spatial_positions_df = bixi_data$spatial_positions_df,
#'   temporal_positions_df = bixi_data$temporal_positions_df,
#'   burn_in_iter = 5, sampling_iter = 10) # For example only (too few iterations)
#' bktr_regressor$mcmc_sampling()
#'
#' # Plot Y estimates vs observed y values
#' plot_y_estimates(bktr_regressor)
#'
#' @export
plot_y_estimates <- function(
    bktr_reg,
    show_figure = TRUE,
    fig_width = 5,
    fig_height = 5,
    fig_resolution = 200,
    fig_title = 'y estimates vs observed y values'
) {
    if (!bktr_reg$has_completed_sampling) {
        stop('Plots can only be accessed after MCMC sampling.')
    }
    # Verify all labels are valid
    omega_list <- as.numeric(bktr_reg$omega$flatten()$cpu()) != 0
    y_est_list <- bktr_reg$y_estimates[omega_list]$y_est
    y_list <- as.numeric(bktr_reg$y$flatten()$cpu()[omega_list])
    min_y <- min(y_list)
    max_y <- max(y_list)
    df <- data.table(y = y_list, y_est = y_est_list)
    fig <- (
        ggplot(df, aes(x = .data$y, y = .data$y_est))
        + geom_point(color = '#39a7d0', alpha = 0.6, shape = 21, fill = '#20a0d0')
        + geom_segment(aes(x = min_y, y = min_y, xend = max_y, yend = max_y), color = 'black',
                       linetype = 'twodash', linewidth = 1)
        + theme_bw()
        + ylab('Estimated y')
        + xlab('Observed y')
    )
    if (!is.null(fig_title)) {
        fig <- fig + ggtitle(fig_title)
    }
    if (!show_figure) {
        return(fig)
    }
    print_ggplot_fig(fig, fig_width, fig_height, fig_resolution)
}
