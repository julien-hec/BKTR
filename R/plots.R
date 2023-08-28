#' @import ggplot2
#' @import ggmap

#' @description Create a plot of the beta values through time for a given
#' spatial point and a set of feature labels.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param plot_feature_labels Array: Array of feature labels to plot.
#' @param spatial_point_label String: Spatial point label to plot.
#' @param date_format String: Format of the date to use in bktr dataframes for the time.
#'   Defaults to '%Y-%m-%d'.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Numeric: Figure width when figure is shown. Defaults to 8.5.
#' @param fig_height Numeric: Figure height when figure is shown. Defaults to 5.5.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @export
plot_temporal_betas <- function(
    bktr_reg,
    plot_feature_labels,
    spatial_point_label,
    date_format = '%Y-%m-%d',
    show_figure = TRUE,
    fig_width = 8.5,
    fig_height = 5.5
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
        ggplot(beta_est_df, aes(time, Mean, group = feature, color = feature))
        + geom_line()
        + geom_ribbon(
            aes(ymin = Low2.5p, ymax = Up97.5p, fill = feature),
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
    print(fig, vp = grid::viewport(width = unit(fig_width, 'inches'), height = unit(fig_height, 'inches')))
}

#' @description Create a plot of beta values through space for a given
#' temporal point and a set of feature labels.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param plot_feature_labels Array: Array of feature labels to plot.
#' @param temporal_point_label String: Temporal point label to plot.
#' @param nb_cols Integer: The number of columns to use in the facet grid.
#' @param use_dark_mode Boolean: Whether to use a dark mode for the geographic map or not.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Numeric: Figure width when figure is shown. Defaults to 8.5.
#' @param fig_height Numeric: Figure height when figure is shown. Defaults to 5.5.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @export
plot_spatial_betas <- function(
    bktr_reg,
    plot_feature_labels,
    temporal_point_label,
    nb_cols = 1,
    use_dark_mode = TRUE,
    show_figure = TRUE,
    fig_width = 8.5,
    fig_height = 5.5
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
    beta_df <- beta_df[beta_df$time == temporal_point_label, ..out_col_names]

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
    if (is_map) {
        map_type <- ifelse(use_dark_mode, 'toner', 'toner-lite')
        fig <- (
            qmplot(x = longitude, y = latitude, color = value, data = full_df, maptype = map_type)
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
            ggplot(full_df, aes(x = .data[[x_col_name]], y = .data[[y_col_name]], color = value))
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
    print(fig, vp = grid::viewport(width = unit(fig_width, 'inches'), height = unit(fig_height, 'inches')))
}



#' @description Plot the distribution of beta values for a given list of labels.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param labels_list List: List of labels tuple (spatial, temporal, feature) for
#'  which to plot the beta distribution through iterations
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Integer: Figure width. Defaults to 9.
#' @param fig_height Integer: Figure height. Defaults to 6.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @export
plot_beta_dists <- function(
    bktr_reg,
    labels_list,
    show_figure = TRUE,
    fig_width = 9,
    fig_height = 6
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
        ggplot(df, aes(x = labels, y = value, fill = labels))
        + geom_violin(trim = FALSE)
        + ggtitle(plot_title)
        + ylab('Beta Value')
        + xlab('Labels')
        + labs(fill = 'Labels', color = 'Labels')
    )
    if (!show_figure) {
        return(fig)
    }
    print(fig, vp = grid::viewport(width = unit(fig_width, 'inches'), height = unit(fig_height, 'inches')))
}


#' @description Plot the distribution of beta estimates regrouped by covariates.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param feature_labels Array or NULL: Array of feature labels for
#'   which to plot the beta estimates distribution. If NULL plot for all features.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Integer: Figure width. Defaults to 9.
#' @param fig_height Integer: Figure height. Defaults to 6.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @export
plot_covariates_beta_dists <- function(
    bktr_reg,
    feature_labels,
    show_figure = TRUE,
    fig_width = 9,
    fig_height = 6
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
        bktr_reg$beta_estimates[, ..feature_labels],
        direction = 'long',
        v.names = 'value',
        timevar = 'feature',
        varying = feature_labels,
        times = feature_labels
    )
    fig <- (
        ggplot(full_df, aes(x=feature, y=value, fill=feature))
        + geom_violin()
        + ggtitle(plot_title)
        + ylab('Beta Value')
        + xlab('Feature')
        + labs(fill = 'Feature', color = 'Feature')
    )
    if (!show_figure) {
        return(fig)
    }
    print(fig, vp = grid::viewport(width = unit(fig_width, 'inches'), height = unit(fig_height, 'inches')))
}

#' @description Plot the distribution of hyperparameters through iterations
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param hyperparameters Array or NULL: Array of hyperparameters to plot.
#'     If NULL, plot all hyperparameters. Defaults to NULL.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Integer: Figure width. Defaults to 9.
#' @param fig_height Integer: Figure height. Defaults to 6.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @export
plot_hyperparams_dists <- function(
    bktr_reg,
    hyperparameters = NULL,
    show_figure = TRUE,
    fig_width = 9,
    fig_height = 6
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
        ggplot(df, aes(x = hyperparameter, y = value, fill = hyperparameter))
        + geom_violin(trim = FALSE)
        + ggtitle('Posterior Distribution of BKTR Hyperparameters')
        + ylab('Hyperparameter Value')
        + xlab('Hyperparameter')
        + labs(fill = 'Hyperparameter', color = 'Hyperparameter')
    )
    if (!show_figure) {
        return(fig)
    }
    print(fig, vp = grid::viewport(width = unit(fig_width, 'inches'), height = unit(fig_height, 'inches')))
}

#' @description Plot the evolution of hyperparameters through iterations. (Traceplot)
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param hyperparameters Array or NULL: Array of hyperparameters to plot.
#'     If NULL, plot all hyperparameters. Defaults to NULL.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Integer: Figure width. Defaults to 9.
#' @param fig_height Integer: Figure height. Defaults to 5.5.
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @export
plot_hyperparams_traceplot <- function(
    bktr_reg,
    hyperparameters = NULL,
    show_figure = TRUE,
    fig_width = 9,
    fig_height = 5.5
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
        df[, ..df_cols],
        direction = 'long',
        v.names = 'value',
        idvar = 'iter',
        timevar = 'hyperparameter',
        varying = hparams,
        times = hparams
    )
    fig <- (
        ggplot(df, aes(iter, value, group = hyperparameter, color = hyperparameter))
        + geom_line()
        + ggtitle('Hyperparameter values through sampling iterations (Traceplot)')
        + theme_bw()
        + ylab('Hyperparameter Value')
        + xlab('Sampling Iter')
        + labs(fill = 'Hyperparameter', color = 'Hyperparameter')
        + theme(legend.position='bottom')
    )
    if (!show_figure) {
        return(fig)
    }
    print(fig, vp = grid::viewport(width = unit(fig_width, 'inches'), height = unit(fig_height, 'inches')))
}


#' @description Plot y estimates vs observed y values.
#' @param bktr_reg BKTRRegressor: BKTRRegressor object.
#' @param show_figure Boolean: Whether to show the figure. Defaults to True.
#' @param fig_width Numeric: Figure width when figure is shown. Defaults to 5.
#' @param fig_height Numeric: Figure height when figure is shown. Defaults to 5.
#' @param fig_title String or NULL: Figure title if provided. Defaults to 'y estimates vs observed y values'
#' @return ggplot or NULL: ggplot object or NULL if show_figure is set to FALSE.
#'
#' @export
plot_y_estimates <- function(
    bktr_reg,
    show_figure = TRUE,
    fig_width = 5,
    fig_height = 5,
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
        ggplot(df, aes(x = y, y = y_est))
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
    print(fig, vp = grid::viewport(width = unit(fig_width, 'inches'), height = unit(fig_height, 'inches')))
}
