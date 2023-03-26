#' @import ggplot2
#' @import ggmap


#' @title R6 class for Creating plots of the BKTR estimates
#'
#' @description R6 class for creating plots of the BKTR estimates. Can be used
#' to create plots of the beta estimates for the spatial and temporal features.
#'
#' @export
#' @keywords internal
BKTRBetaPlotMaker <- R6::R6Class(
    public = list(
        bktr_regressor = NULL,
        spatial_labels = NULL,
        temporal_labels = NULL,
        feature_labels = NULL,
        initialize = function(bktr_regressor) {
            self$bktr_regressor <- bktr_regressor
            self$spatial_labels <- bktr_regressor$spatial_labels
            self$temporal_labels <- bktr_regressor$temporal_labels
            self$feature_labels <- c(
                '_INTERSECT_',
                bktr_regressor$spatial_feature_labels,
                bktr_regressor$temporal_feature_labels
            )
        },
        create_temporal_beta_plot = function(
          plot_feature_labels,
          plot_point_label,
          date_format = '%Y-%m-%d'
        ) {
          beta_est_df <- self$get_long_temporal_est_df(
            self$bktr_regressor$beta_estimates,
            plot_feature_labels,
            plot_point_label
          )
          beta_stdev_df <- self$get_long_temporal_est_df(
            self$bktr_regressor$beta_stdev,
            plot_feature_labels,
            plot_point_label
          )
          plot_df <- cbind(beta_est_df, stdev = beta_stdev_df$value)
          plot_title <- paste0('Estimated Beta Values for ', plot_point_label)

          if (!is.null(date_format)) plot_df$time <- as.Date(plot_df$time, format = date_format)
          return(
            ggplot(plot_df, aes(time, value, group = covariate, color = covariate))
            + geom_line()
            + geom_ribbon(
              aes(ymin = value - stdev, ymax = value + stdev, fill = covariate),
              alpha = 0.3,
              color = NA
            )
            + ggtitle(plot_title)
            + theme_bw()
          )
        },
        get_long_temporal_est_df = function(estimates_tensor, plot_feature_labels, plot_point_label) {
          feature_labels_indexes <- match(plot_feature_labels, self$feature_labels)
          point_label_index <- match(plot_point_label, self$spatial_labels)
          est_values <- estimates_tensor[point_label_index, , feature_labels_indexes]
          est_df <- as.data.frame(as.matrix(est_values), row.names = self$temporal_labels)
          colnames(est_df) <- plot_feature_labels
          est_df <- cbind(time = rownames(est_df), est_df)
          est_df <- reshape(
            est_df,
            direction = 'long',
            idvar = 'time',
            v.names = 'value',
            timevar = 'covariate',
            varying = plot_feature_labels,
            times = plot_feature_labels
          )
          rownames(est_df) <- 1:nrow(est_df)
          return(est_df)
        },
        create_spatial_beta_plots = function(
            plot_feature_labels,
            plot_point_label,
            geo_coordinates,
            nb_cols = 1,
            mapbox_zoom = 9,
            use_dark_mode = TRUE
        ) {
            beta_est <- self$bktr_regressor$beta_estimates
            feature_labels_indexes <- match(plot_feature_labels, self$feature_labels)
            # Get only spatial estimates
            point_label_index <- match(plot_point_label, self$temporal_labels)
            beta_est_values <- beta_est[, point_label_index, feature_labels_indexes]
            beta_est_df <- as.data.frame(as.matrix(beta_est_values), row.names = self$spatial_labels)
            colnames(beta_est_df) <- plot_feature_labels

            map_df <- cbind(geo_coordinates, beta_est_df[c(plot_feature_labels)])
            map_df <- cbind(location = rownames(map_df), map_df)
            map_df <- reshape(
              map_df,
              direction = 'long',
              idvar = c('location', 'latitude', 'longitude'),
              v.names = 'value',
              timevar = 'covariate',
              varying = plot_feature_labels,
              times = plot_feature_labels
            )

            map_type <- ifelse(use_dark_mode, 'toner', 'toner-lite')
            plot_title <- paste0('Estimated Beta Values for ', plot_point_label)
            return(
              qmplot(x = longitude, y = latitude, color = value, data = map_df, maptype = map_type)
              + facet_wrap(~covariate, ncol = nb_cols)
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
        }
    )
)
