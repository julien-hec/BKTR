#' @import torch

# Private utility function to do a cross join between two data.tables.
# Taken from https://github.com/Rdatatable/data.table/issues/1717
cross_join_dt <- function(...) {
    rows <- do.call(CJ, lapply(
        list(...),
        function(x) if (is.data.frame(x)) seq_len(nrow(x)) else seq_along(x)
    ))
    do.call(data.table, Map(function(x, y) x[y], list(...), rows))
}

#' Function used to transform covariates coming from two dataframes one for spatial and
#' one for temporal into a single dataframe with the right shape for the BKTR Regressor.
#' This is useful when the temporal covariates do not vary trough space and the spatial
#' covariates do not vary trough time (Like in the BIXI example). The function also adds
#' a column for the target variable at the beginning of the dataframe.
#' @param spatial_df data.table: Spatial covariates dataframe with an index named
#'       location and a shape of (n_locations, n_spatial_covariates)
#' @param temporal_df data.table: Temporal covariates dataframe with an index named
#'       time and a shape of (n_times, n_temporal_covariates)
#' @param y_df data.table: The dataframe containing the target variable. It should have
#'       a shape of (n_locations, n_times). The columns and index names of this dataframe
#'       should be correspond to the one of the spatial_df and temporal_df.
#' @param y_column_name string: The name of the target variable column in y_df. Default
#'       to 'y'.

#' @return data.table: The reshaped covariates dataframe with a shape of
#'       (n_locations * n_times, 1 + n_spatial_covariates + n_temporal_covariates).
#'       The first two columns are the indexes (location, time), the following column
#'       is the target variable and the other columns are the covariates.
#'
#' @export
reshape_covariate_dfs <- function(
    spatial_df,
    temporal_df,
    y_df,
    y_column_name = 'y'
) {
    spa_index_name <- 'location'
    temp_index_name <- 'time'
    if (indices(spatial_df) != c(spa_index_name) || indices(y_df) != c(spa_index_name)) {
        stop(paste('Index names of spatial_df and y_df must be', spa_index_name))
    }
    if (indices(temporal_df) != c(temp_index_name)) {
        stop(paste('Index name of temporal_df must be', temp_index_name))
    }
    spatial_df_cp <- spatial_df
    temporal_df_cp <- temporal_df
    y_df_cp <- y_df

    # Sort indexes and columns
    for (df in list(spatial_df_cp, temporal_df_cp, y_df_cp)) {
        setorderv(df, indices(df)[1])
    }
    y_df_col_names <- sort(colnames(y_df_cp)[!colnames(y_df_cp) == spa_index_name])
    setcolorder(y_df_cp, c(spa_index_name, y_df_col_names))

    # Compare indexes values
    if (!identical(spatial_df_cp[[indices(spatial_df_cp)[1]]], y_df_cp[[indices(y_df_cp)[1]]])) {
        stop('Index values of spatial_df and y_df must be the same')
    }
    if (!identical(as.character(temporal_df_cp[[indices(temporal_df_cp)[1]]]), y_df_col_names)) {
        stop('temporal_df index values and y_df columns names must be the same')
    }
    data_df <- cross_join_dt(spatial_df_cp, temporal_df_cp)
    setindexv(data_df, c(spa_index_name, temp_index_name))
    setorderv(data_df, c(spa_index_name, temp_index_name))
    y_flat_values <- as.vector(unlist(y_df_cp[, ..y_df_col_names]))
    data_df[, (y_column_name) := y_flat_values]
    setcolorder(data_df, c(spa_index_name, temp_index_name, y_column_name))
    return(data_df)
}

#' Simulate spatiotemporal data using kernel covariances.
#'
#' @param nb_spatial_locations Integer: Number of spatial locations
#' @param nb_time_points Integer: Number of time points
#' @param nb_spatial_dimensions Integer: Number of spatial dimensions
#' @param spatial_scale Numeric: Spatial scale
#' @param time_scale Numeric: Time scale
#' @param spatial_covariates_means Vector: Spatial covariates means
#' @param temporal_covariates_means Vector: Temporal covariates means
#' @param spatial_kernel Kernel: Spatial kernel
#' @param temporal_Kernel Kernel: Temporal kernel
#' @param noise_variance_scale Numeric: Noise variance scale
#
#' @return A list containing 4 dataframes:
#'      - `data_df` contains the response variable and the covariates
#'      - `spatial_positions_df` contains the spatial locations and their coordinates
#'      - `temporal_positions_df` contains the time points and their coordinates
#'      - `beta_df` contains the true beta coefficients
#'
#' @export
simulate_spatiotemporal_data <- function(
    nb_locations,
    nb_time_points,
    nb_spatial_dimensions,
    spatial_scale,
    time_scale,
    spatial_covariates_means,
    temporal_covariates_means,
    spatial_kernel,
    temporal_kernel,
    noise_variance_scale
) {
    spa_pos <- BKTR:::tsr$rand(c(nb_locations, nb_spatial_dimensions)) * spatial_scale
    temp_pos <- BKTR:::tsr$arange(0, nb_time_points - 1) * time_scale / (nb_time_points - 1)
    temp_pos <- temp_pos$reshape(c(nb_time_points, 1))

    # Dimension labels
    s_dims <- get_dim_labels('s_dim', nb_spatial_dimensions)
    s_locs <- get_dim_labels('s_loc', nb_locations)
    t_points <- get_dim_labels('t_point', nb_time_points)
    s_covs <- get_dim_labels('s_cov', length(spatial_covariates_means))
    t_covs <- get_dim_labels('t_cov', length(temporal_covariates_means))

    spa_pos_df <- cbind(data.table(s_locs), data.table(as.matrix(spa_pos)))
    setnames(spa_pos_df, c('location', s_dims))
    setindexv(spa_pos_df, 'location')
    temp_pos_df <- cbind(data.table(t_points), data.table(as.matrix(temp_pos)))
    setnames(temp_pos_df, c('time', 'time_val'))
    setindexv(temp_pos_df, 'time')

    spa_means <- BKTR:::tsr$new_tensor(spatial_covariates_means)
    nb_spa_covariates <- length(spa_means)
    spa_covariates <- BKTR:::tsr$randn(c(nb_locations, nb_spa_covariates))
    spa_covariates <- spa_covariates + spa_means

    temp_means <- BKTR:::tsr$new_tensor(temporal_covariates_means)
    nb_temp_covariates <- length(temp_means)
    temp_covariates <- BKTR:::tsr$randn(c(nb_time_points, nb_temp_covariates))
    temp_covariates <- temp_covariates + temp_means
    intercept_covariates <- BKTR:::tsr$ones(c(nb_locations, nb_time_points, 1))
    covs <- torch::torch_cat(
        c(
            intercept_covariates,
            spa_covariates$unsqueeze(2)$expand(c(nb_locations, nb_time_points, nb_spa_covariates)),
            temp_covariates$unsqueeze(1)$expand(c(nb_locations, nb_time_points, nb_temp_covariates))
        ),
        dim = 3
    )
    nb_covs <- 1 + nb_spa_covariates + nb_temp_covariates

    covs_covariance_mat <- rWishart(1, nb_covs, diag(nb_covs))[,,1]
    covs_covariance <- BKTR:::tsr$new_tensor(covs_covariance_mat)

    spatial_kernel$set_positions(spa_pos_df)
    spatial_covariance <- spatial_kernel$kernel_gen()
    temporal_kernel$set_positions(temp_pos_df)
    temporal_covariance <- temporal_kernel$kernel_gen()

    # Use Matrix Normal distribution to sample beta values (to reduce memory usage)
    # the second covariance matrix is the Kronecker product of temporal and covariates covariances
    chol_spa <- torch::linalg_cholesky(spatial_covariance)
    chol_temp_covs <- torch::linalg_cholesky(
        BKTR:::tsr$kronecker_prod(temporal_covariance, covs_covariance)
    )
    beta_values <- (
        chol_spa$matmul(BKTR:::tsr$randn(c(nb_locations, nb_time_points * nb_covs)))$matmul(chol_temp_covs$t())
    )$reshape(c(nb_locations, nb_time_points, nb_covs))
    y_val <- torch:::torch_einsum('ijk,ijk->ij', c(covs, beta_values))
    err <- BKTR:::tsr$randn(c(nb_locations, nb_time_points)) * (noise_variance_scale ** 0.5)
    y_val <- y_val + err
    y_val <- y_val$reshape(c(nb_locations * nb_time_points, 1))
    # We remove the intercept from the covariates
    covs <- covs$reshape(c(nb_locations * nb_time_points, nb_covs))[, 2:nb_covs]

    index_cols_df <- CJ(spa_pos_df[['location']], temp_pos_df[['time']])
    setnames(index_cols_df, c('location', 'time'))
    data_df <- data.table(cbind(as.matrix(y_val), as.matrix(covs)))
    setnames(data_df, c('y', s_covs, t_covs))
    data_df <- cbind(index_cols_df, data_df)
    setindexv(data_df, c('location', 'time'))
    beta_df <- data.table(as.matrix(beta_values$reshape(c(nb_locations * nb_time_points, nb_covs))))
    setnames(beta_df, c('Intercept', s_covs, t_covs))
    beta_df <- cbind(index_cols_df, beta_df)
    setindexv(beta_df, c('location', 'time'))

    return(list(
        data_df = data_df,
        spatial_positions_df = spa_pos_df,
        temporal_positions_df = temp_pos_df,
        beta_df = beta_df
    ))
}


# Private utility function to get the dimension labels for a
#     given dimension prefix and max value.
# Args:
#     dim_prefix (str): The prefix of the dimension labels
#     max_value (int): The maximum value of the dimension labels
# Returns:
#     str: The dimension labels
get_dim_labels <- function(dim_prefix, max_value) {
    max_digits <- nchar(as.character(max_value - 1))
    formatted_numbers <- formatC(0:(max_value - 1), width = max_digits, flag = "0")
    return(paste(dim_prefix, formatted_numbers, sep = "_"))
}


# Get the index of a label in a list of labels. If the label is not in the list,
# raise an error.
#
# @param label Any: The label for which we want to get the index
# @param label_list Vector[Any]: The list of labels
# label_type String: The label type either 'spatial', 'temporal', 'feature'.
#
# @return Integer: The index of the label in the list
get_label_index_or_raise <- function(label, label_list, label_type) {
    match_indx <- match(label, label_list)
    if (is.na(match_indx)) {
        stop(sprintf('Label `%s` does not exist in %s labels.', label, label_type))
    }
    return(match_indx)
}