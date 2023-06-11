#' @import ggplot2


DEFAULT_LBOUND <- 1e-3
DEFAULT_UBOUND <- 1e3


#' @title R6 class for kernel's hyperparameter
#'
#' @description KernelParameter contains all information and behaviour related to a kernel parameters.
#'
#' @export
KernelParameter <- R6::R6Class(
    public = list(
        value = 0,
        is_fixed = FALSE,
        lower_bound = DEFAULT_LBOUND,
        upper_bound = DEFAULT_UBOUND,
        slice_sampling_scale = log(10),
        hparam_precision = 1.0,
        #' @field The kernel associated with the parameter (it is set at kernel instanciation)
        kernel = NULL,
        #' @field The kernel parameter's name
        name = NULL,
        #' @description Create a new \code{KernelParameter} object.
        #' @param value Numeric: The hyperparameter mean's prior value (Paper - \eqn{\phi}) or its constant value
        #' @param is_fixed Boolean: Says if the kernel parameter is fixed or not (if fixed, there is no sampling)
        #' @param lower_bound Numeric: Hyperparameter's minimal value during sampling (Paper - \eqn{\phi_{min}})
        #' @param upper_bound Numeric: Hyperparameter's maximal value during sampling (Paper - \eqn{\phi_{max}})
        #' @param slice_sampling_scale Numeric: The sampling range's amplitude (Paper - \eqn{\rho})
        #' @param hparam_precision Numeric: WHAT IS THAT? TODO
        #' @return A new \code{KernelParameter} object.
        initialize = function(
            value,
            is_fixed = FALSE,
            lower_bound = DEFAULT_LBOUND,
            upper_bound = DEFAULT_UBOUND,
            slice_sampling_scale = log(10),
            hparam_precision = 1.0
        ){
            self$value <- value
            self$lower_bound <- lower_bound
            self$upper_bound <- upper_bound
            self$is_fixed <- is_fixed
            self$slice_sampling_scale <- slice_sampling_scale
            self$hparam_precision <- hparam_precision
        },

        #' @description Set \code{Kernel} for a given \code{KernelParameter} object.
        #' @param kernel Kernel: The kernel to associate with the parameter
        #' @return NULL, set a new kernel for the parameter
        set_kernel = function(kernel, param_name) {
            self$kernel <- kernel
            self$name <- param_name
            self$kernel$parameters <- c(self$kernel$parameters, self)
        }
    ),
    active = list(
        full_name = function() {
            if (self$kernel == NULL) {
                return(self$name)
            }
            return(sprintf('%s - %s', self$kernel$name, self$name))
        }
    )
)


Kernel <- R6::R6Class(
    'Kernel',
    public = list(
        kernel_variance = 1,
        jitter_value = NULL,
        distance_matrix = NULL,
        name = NULL,
        parameters = c(),
        covariance_matrix = NULL,
        positions_df = NULL,
        distance_type = NULL,

        #' @description Kernel abstract base constructor
        #' @param kernel_variance Numeric: The variance of the kernel
        #' @param distance_type String: The distance type to use for calculating the kernel distance matrix
        #' @param jitter_value Numeric: The jitter value to add to the kernel matrix
        #' @return A new \code{Kernel} object.
        initialize = function(kernel_variance, distance_type, jitter_value) {
            self$parameters <- c()
            self$kernel_variance <- kernel_variance
            self$distance_type <- distance_type
            self$jitter_value <- jitter_value
        },

        core_kernel_fn = function() {
            stop('core_kernel_fn() is not implemented')
        },

        add_jitter_to_kernel = function() {
            has_null_jitter <- is.null(self$jitter_value)
            if (!has_null_jitter && self$jitter_value == 0) {
                return()
            }
            jitter_val <- ifelse(has_null_jitter, BKTR:::tsr$default_jitter, self$jitter_value)
            self$covariance_matrix <- self$covariance_matrix + jitter_val * BKTR:::tsr$eye(nrow(self$covariance_matrix))
        },

        kernel_gen = function() {
            if (is.null(self$positions_df)) {
                stop('Set `positions_df` via `set_positions` before kernel evaluation.')
            }
            self$covariance_matrix <- self$kernel_variance * self$core_kernel_fn()
            self$add_jitter_to_kernel()
            return(self$covariance_matrix)
        },

        set_positions = function(positions_df) {
            pos_df_indx <- indices(positions_df)
            if (is.null(pos_df_indx) || length(pos_df_indx) != 1) {
                stop('`positions_df` must have one and only index set via setindex.')
            }
            self$positions_df <- positions_df
            positions_tensor <- BKTR:::tsr$new_tensor(as.matrix(positions_df[, !..pos_df_indx]))
            # TODO: check to transform that into a function `get_distance_matrix` #13
            self$distance_matrix <- DistanceCalculator$new()$get_matrix(positions_tensor, self$distance_type)
        },

        plot = function(self, show_figure = TRUE) {
            x_name <- indices(self$positions_df)[1]
            y_name <- paste0(x_name, "'")
            df <- data.table(self$covariance_matrix)
            pos_labels <- as.character(self$positions_df[[x_name]])
            colnames(df) <- pos_labels
            df[[x_name]] <- pos_labels
            df <- melt(df, id.vars = c(x_name), variable.name = y_name, value.name = 'covariance')
            fig <- ggplot(df, aes(.data[[x_name]], .data[[y_name]], fill = covariance)) +
                geom_tile() + theme_minimal() + scale_x_discrete(limits = pos_labels) +
                scale_y_discrete(limits = rev(pos_labels))
            if (show_figure) {
                print(fig)
                return(NULL)
            }
            return(fig)
        }
    )
)

#' @title R6 class for White Noise Kernels
#'
#' @description R6 class for White Noise Kernels
#'
#' @export
KernelWhiteNoise <- R6::R6Class(
    inherit = Kernel,
    public = list(
        variance = NULL,
        distance_matrix = NULL,
        name = 'White Noise Kernel',
        initialize = function(
            kernel_variance = 1,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, DIST_TYPE$NONE, jitter_value)
        },
        core_kernel_fn = function() {
            return(BKTR:::tsr$eye(nrow(self$positions_df)))
        }
    )
)


#' @title R6 class for Square Exponential Kernels
#'
#' @description R6 class for Square Exponential Kernels
#'
#' @export
KernelSE <- R6::R6Class(
    inherit = Kernel,
    public = list(
        lengthscale = NULL,
        distance_matrix = NULL,
        name = 'SE Kernel',
        initialize = function(
            lengthscale = KernelParameter$new(log(2)),
            kernel_variance = 1,
            distance_type = DIST_TYPE$EUCLIDEAN,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$lengthscale <- lengthscale
            self$lengthscale$set_kernel(self, 'lengthscale')
        },
        core_kernel_fn = function() {
            return(torch::torch_exp(-self$distance_matrix^2 / (2 * self$lengthscale$value^2)))
        }
    )
)

#' @title R6 class for Square Exponential Kernels
#'
#' @description R6 class for Square Exponential Kernels
#'
#' @export
KernelRQ <- R6::R6Class(
    inherit = Kernel,
    public = list(
        lengthscale = NULL,
        alpha = NULL,
        distance_matrix = NULL,
        name = 'RQ Kernel',
        initialize = function(
            lengthscale = KernelParameter$new(log(2)),
            alpha = KernelParameter$new(log(2)),
            kernel_variance = 1,
            distance_type = DIST_TYPE$EUCLIDEAN,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$lengthscale <- lengthscale
            self$lengthscale$set_kernel(self, 'lengthscale')
            self$alpha <- alpha
            self$alpha$set_kernel(self, 'alpha')
        },
        core_kernel_fn = function() {
            return(
                1 + self$distance_matrix^2 / (2 * self$lengthscale$value^2 * self$alpha$value)
            ) ** -self$alpha$value
        }
    )
)

#' @title R6 class for Periodic Kernels
#'
#' @description R6 class for Periodic Kernels
#'
#' @export
KernelPeriodic <- R6::R6Class(
    inherit = Kernel,
    public = list(
        lengthscale = NULL,
        period_length = NULL,
        distance_matrix = NULL,
        name = 'Periodic Kernel',
        initialize = function(
            lengthscale = KernelParameter$new(log(2)),
            period_length = KernelParameter$new(log(2)),
            kernel_variance = 1,
            distance_type = DIST_TYPE$EUCLIDEAN,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$lengthscale <- lengthscale
            self$lengthscale$set_kernel(self, 'lengthscale')
            self$period_length <- period_length
            self$period_length$set_kernel(self, 'period length')
        },
        core_kernel_fn = function() {
            return(torch::torch_exp(
                -2 * torch::torch_sin(pi * self$distance_matrix / self$period_length$value)^2 /
                    self$lengthscale$value^2
            ))
        }
    )
)

#' @title R6 class for Matern Kernels
#'
#' @description R6 class for Matern Kernels
#'
#' @export
KernelMatern <- R6::R6Class(
    inherit = Kernel,
    public = list(
        lengthscale = NULL,
        smoothness_factor = NULL,
        distance_matrix = NULL,
        initialize = function(
            smoothness_factor = 1,
            lengthscale = KernelParameter$new(log(2)),
            kernel_variance = 1,
            distance_type = DIST_TYPE$HAVERSINE,
            jitter_value = NULL
        ) {
            if (smoothness_factor %in% c(1, 3, 5) == FALSE) {
                stop('Smoothness factor should be one of the following values 1, 3 or 5')
            }
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$name <- paste0('Matern ', smoothness_factor, '/2 Kernel')
            self$smoothness_factor <- smoothness_factor
            self$lengthscale <- lengthscale
            self$lengthscale$set_kernel(self, 'lengthscale')
        },
        get_smoothness_kernel_fn = function() {
            if (self$smoothness_factor == 1)  {
                return(function(t) return(1))
            } else if (self$smoothness_factor == 3) {
                return(function(t) return(1 + t))
            } else if (self$smoothness_factor == 5) {
                return(function(t) return(1 + t * (1 + t / 3)))
            } else {
                torch:::value_error('Kernel function for this smoothness factor is not implemented')
            }
        },
        core_kernel_fn = function() {
            temp_kernel <- (
                self$distance_matrix * torch::torch_sqrt(self$smoothness_factor) / self$lengthscale$value
            )
            return(self$get_smoothness_kernel_fn()(temp_kernel) * torch::torch_exp(-temp_kernel))
        }
    ),
)

#' @title Kernel Composition Operations
#'
#' @description Kernel Composition Operations Enum. Possibilities of operation between
#' two kernels to generate a new composed kernel.
#'
#' @export
CompositionOps = list(
    'MUL' = 'MUL',
    'ADD' = 'ADD'
)

#' @title R6 class for Composed Kernels
#'
#' @description R6 class for Composed Kernels
#'
#' @export
#' @keywords internal
KernelComposed <- R6::R6Class(
    inherit = Kernel,
    public = list(
        name = '',
        parameters = c(),
        left_kernel = NULL,
        right_kernel = NULL,
        composition_operation = NULL,
        initialize = function(
            left_kernel,
            right_kernel,
            new_name,
            composition_operation
        ) {
            composed_variance <- 1
            if (left_kernel$distance_type != right_kernel$distance_type) {
                stop('Composed kernel must have the same distance type')
            }
            new_jitter_val <- max(
                left_kernel$jitter_value,
                right_kernel$jitter_value,
                BKTR:::tsr$default_jitter
            )
            super$initialize(composed_variance, left_kernel$distance_type, new_jitter_val)
            self$left_kernel <- left_kernel
            self$right_kernel <- right_kernel
            self$name <- new_name
            self$parameters <- c(
                left_kernel$parameters,
                right_kernel$parameters
            )
            self$composition_operation <- composition_operation
        },
        core_kernel_fn = function() {
            if (self$composition_operation == CompositionOps$MUL) {
                return(self$left_kernel$core_kernel_fn() * self$right_kernel$core_kernel_fn())
            } else if (self$composition_operation == CompositionOps$ADD) {
                return(self$left_kernel$core_kernel_fn() + self$right_kernel$core_kernel_fn())
            } else {
                torch:::value_error('Composition operation is not implemented')
            }
        },
        set_positions = function(positions_df) {
            super$set_positions(positions_df)
            self$left_kernel$set_positions(positions_df)
            self$right_kernel$set_positions(positions_df)
        }
    )
)

#' @export
`+.Kernel` = function(k1, k2) {
    composed_kernel <- KernelComposed$new(
        k1,
        k2,
        paste0(k1$name, ' + ', k2$name),
        CompositionOps$ADD
    )
    return(composed_kernel)
}

#' @export
`*.Kernel` = function(e1, e2) {
    composed_kernel <- KernelComposed$new(
        e1,
        e2,
        paste0(e1$name, ' * ', e2$name),
        CompositionOps$MUL
    )
    return(composed_kernel)
}
