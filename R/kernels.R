#' @import ggplot2
#' @include tensor_ops.R
#' @include distances.R


DEFAULT_LBOUND <- 1e-3
DEFAULT_UBOUND <- 1e3


#' @title R6 class for kernel's hyperparameter
#'
#' @description KernelParameter contains all information and behaviour related to a kernel parameters.
#'
#' @export
KernelParameter <- R6::R6Class(
    public = list(
        #' @field value The hyperparameter mean's prior value or its constant value
        value = 0,
        #' @field is_fixed Says if the kernel parameter is fixed or not (if fixed, there is no sampling)
        is_fixed = FALSE,
        #' @field lower_bound The hyperparameter's minimal value during sampling
        lower_bound = DEFAULT_LBOUND,
        #' @field upper_bound The hyperparameter's maximal value during sampling
        upper_bound = DEFAULT_UBOUND,
        #' @field slice_sampling_scale The sampling range's amplitude
        slice_sampling_scale = log(10),
        #' @field hparam_precision Precision of the hyperparameter
        hparam_precision = 1.0,
        #' @field kernel The kernel associated with the parameter (it is set at kernel instanciation)
        kernel = NULL,
        #' @field name The kernel parameter's name
        name = NULL,
        #' @description Create a new \code{KernelParameter} object.
        #' @param value Numeric: The hyperparameter mean's prior value (Paper - \eqn{\phi}) or its constant value
        #' @param is_fixed Boolean: Says if the kernel parameter is fixed or not (if fixed, there is no sampling)
        #' @param lower_bound Numeric: Hyperparameter's minimal value during sampling (Paper - \eqn{\phi_{min}})
        #' @param upper_bound Numeric: Hyperparameter's maximal value during sampling (Paper - \eqn{\phi_{max}})
        #' @param slice_sampling_scale Numeric: The sampling range's amplitude (Paper - \eqn{\rho})
        #' @param hparam_precision Numeric: The hyperparameter's precision
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
        #' @param param_name String: The parameter's name
        #' @return NULL, set a new kernel for the parameter
        set_kernel = function(kernel, param_name) {
            self$kernel <- kernel
            self$name <- param_name
            self$kernel$parameters <- c(self$kernel$parameters, self)
        }
    ),
    active = list(
        #' @field full_name The kernel parameter's full name
        full_name = function() {
            if (is.null(self$kernel)) {
                return(self$name)
            }
            return(sprintf('%s - %s', self$kernel$name, self$name))
        }
    )
)


#' @title Base R6 class for Kernels
#' @description Abstract base class for kernels
#' @export
Kernel <- R6::R6Class(
    'Kernel',
    public = list(
        #' @field kernel_variance The variance of the kernel
        kernel_variance = 1,
        #' @field jitter_value The jitter value to add to the kernel matrix
        jitter_value = NULL,
        #' @field distance_matrix The distance matrix between points in a tensor format
        distance_matrix = NULL,
        #' @field name The kernel's name
        name = NULL,
        #' @field parameters The parameters of the kernel (list of \code{KernelParameter})
        parameters = c(),
        #' @field covariance_matrix The covariance matrix of the kernel in a tensor format
        covariance_matrix = NULL,
        #' @field positions_df The positions of the points in a dataframe format
        positions_df = NULL,
        #' @field has_dist_matrix Identify if the kernel has a distance matrix or not
        has_dist_matrix = NULL,

        #' @description Kernel abstract base constructor
        #' @param kernel_variance Numeric: The variance of the kernel
        #' @param jitter_value Numeric: The jitter value to add to the kernel matrix
        #' @return A new \code{Kernel} object.
        initialize = function(kernel_variance, jitter_value) {
            self$parameters <- c()
            self$kernel_variance <- kernel_variance
            self$jitter_value <- jitter_value
        },

        #' @description Abstract method to compute the core kernel's covariance matrix
        core_kernel_fn = function() {
            stop('core_kernel_fn() is not implemented')
        },

        #' @description Method to add jitter to the kernel's covariance matrix
        add_jitter_to_kernel = function() {
            has_null_jitter <- is.null(self$jitter_value)
            if (!has_null_jitter && self$jitter_value == 0) {
                return()
            }
            jitter_val <- ifelse(has_null_jitter, TSR$get_default_jitter(), self$jitter_value)
            self$covariance_matrix <- self$covariance_matrix + jitter_val * TSR$eye(nrow(self$covariance_matrix))
        },

        #' @description Method to compute the kernel's covariance matrix
        kernel_gen = function() {
            if (is.null(self$positions_df)) {
                stop('Set `positions_df` via `set_positions` before kernel evaluation.')
            }
            self$covariance_matrix <- self$kernel_variance * self$core_kernel_fn()
            self$add_jitter_to_kernel()
            return(self$covariance_matrix)
        },

        #' @description Method to set the kernel's positions and compute the distance matrix
        #' @param positions_df Dataframe: The positions of the points in a dataframe format
        set_positions = function(positions_df) {
            if (ncol(positions_df) < 2) {
                stop('`positions_df` must have at least two columns.')
            }
            self$positions_df <- positions_df
            positions_tensor <- TSR$tensor(as.matrix(positions_df[, -1]))
            if (self$has_dist_matrix) {
                self$distance_matrix <- get_euclidean_dist_tsr(positions_tensor)
            }
        },

        #' @description Method to plot the kernel's covariance matrix
        #' @param show_figure Boolean: If TRUE, the figure is shown, otherwise it is returned
        #' @return If \code{show_figure} is TRUE, the figure is shown, otherwise it is returned
        plot = function(show_figure = TRUE) {
            x_name <- colnames(self$positions_df)[1]
            y_name <- paste0(x_name, "'")
            df <- data.table(as.matrix(self$covariance_matrix$cpu()))
            pos_labels <- sapply(self$positions_df[, 1], as.character)
            colnames(df) <- pos_labels
            df[[x_name]] <- pos_labels
            df <- melt(df, id.vars = c(x_name), variable.name = y_name, value.name = 'covariance')
            fig <- ggplot(df, aes(.data[[x_name]], .data[[y_name]], fill = covariance)) +
                geom_tile() + theme_minimal() + scale_x_discrete(limits = pos_labels) +
                scale_y_discrete(limits = rev(pos_labels)) + ggtitle(self$name)
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
    'KernelWhiteNoise',
    inherit = Kernel,
    public = list(
        #' @field has_dist_matrix Identify if the kernel has a distance matrix or not
        has_dist_matrix = FALSE,
        #' @field name The kernel's name
        name = 'White Noise Kernel',
        # @description Create a new \code{KernelWhiteNoise} object.
        #' @param kernel_variance Numeric: The variance of the kernel
        #' @param jitter_value Numeric: The jitter value to add to the kernel matrix
        #' @return A new \code{KernelWhiteNoise} object.
        initialize = function(
            kernel_variance = 1,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, jitter_value)
        },
        #' @description Method to compute the core kernel's covariance matrix
        #' @return The core kernel's covariance matrix
        core_kernel_fn = function() {
            return(TSR$eye(nrow(self$positions_df)))
        }
    )
)


#' @title R6 class for Square Exponential Kernels
#'
#' @description R6 class for Square Exponential Kernels
#'
#' @export
KernelSE <- R6::R6Class(
    'KernelSE',
    inherit = Kernel,
    public = list(
        #' @field lengthscale The lengthscale parameter instance of the kernel
        lengthscale = NULL,
        #' @field has_dist_matrix Identify if the kernel has a distance matrix or not
        has_dist_matrix = TRUE,
        #' @field name The kernel's name
        name = 'SE Kernel',
        #' @description Create a new \code{KernelSE} object.
        #' @param lengthscale KernelParameter: The lengthscale parameter instance of the kernel
        #' @param kernel_variance Numeric: The variance of the kernel
        #' @param jitter_value Numeric: The jitter value to add to the kernel matrix
        #' @return A new \code{KernelSE} object.
        initialize = function(
            lengthscale = KernelParameter$new(2),
            kernel_variance = 1,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, jitter_value)
            self$lengthscale <- lengthscale
            self$lengthscale$set_kernel(self, 'lengthscale')
        },
        #' @description Method to compute the core kernel's covariance matrix
        #' @return The core kernel's covariance matrix
        core_kernel_fn = function() {
            return(torch::torch_exp(-self$distance_matrix^2 / (2 * self$lengthscale$value^2)))
        }
    )
)

#' @title R6 class for Rational Quadratic Kernels
#'
#' @description R6 class for Rational Quadratic Kernels
#'
#' @export
KernelRQ <- R6::R6Class(
    'KernelRQ',
    inherit = Kernel,
    public = list(
        #' @field lengthscale The lengthscale parameter instance of the kernel
        lengthscale = NULL,
        #' @field alpha The alpha parameter instance of the kernel
        alpha = NULL,
        #' @field has_dist_matrix The distance matrix between points in a tensor format
        has_dist_matrix = TRUE,
        #' @field name The kernel's name
        name = 'RQ Kernel',
        #' @description Create a new \code{KernelRQ} object.
        #' @param lengthscale KernelParameter: The lengthscale parameter instance of the kernel
        #' @param alpha KernelParameter: The alpha parameter instance of the kernel
        #' @param kernel_variance Numeric: The variance of the kernel
        #' @param jitter_value Numeric: The jitter value to add to the kernel matrix
        #' @return A new \code{KernelRQ} object.
        initialize = function(
            lengthscale = KernelParameter$new(2),
            alpha = KernelParameter$new(2),
            kernel_variance = 1,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, jitter_value)
            self$lengthscale <- lengthscale
            self$lengthscale$set_kernel(self, 'lengthscale')
            self$alpha <- alpha
            self$alpha$set_kernel(self, 'alpha')
        },
        #' @description Method to compute the core kernel's covariance matrix
        #' @return The core kernel's covariance matrix
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
    'KernelPeriodic',
    inherit = Kernel,
    public = list(
        #' @field lengthscale The lengthscale parameter instance of the kernel
        lengthscale = NULL,
        #' @field period_length The period length parameter instance of the kernel
        period_length = NULL,
        #' @field has_dist_matrix Identify if the kernel has a distance matrix or not
        has_dist_matrix = TRUE,
        #' @field name The kernel's name
        name = 'Periodic Kernel',
        #' @description Create a new \code{KernelPeriodic} object.
        #' @param lengthscale KernelParameter: The lengthscale parameter instance of the kernel
        #' @param period_length KernelParameter: The period length parameter instance of the kernel
        #' @param kernel_variance Numeric: The variance of the kernel
        #' @param jitter_value Numeric: The jitter value to add to the kernel matrix
        #' @return A new \code{KernelPeriodic} object.
        initialize = function(
            lengthscale = KernelParameter$new(2),
            period_length = KernelParameter$new(2),
            kernel_variance = 1,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, jitter_value)
            self$lengthscale <- lengthscale
            self$lengthscale$set_kernel(self, 'lengthscale')
            self$period_length <- period_length
            self$period_length$set_kernel(self, 'period length')
        },
        #' @description Method to compute the core kernel's covariance matrix
        #' @return The core kernel's covariance matrix
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
    'KernelMatern',
    inherit = Kernel,
    public = list(
        #' @field lengthscale The lengthscale parameter instance of the kernel
        lengthscale = NULL,
        #' @field smoothness_factor The smoothness factor of the kernel
        smoothness_factor = NULL,
        #' @field has_dist_matrix Identify if the kernel has a distance matrix or not
        has_dist_matrix = TRUE,
        #' @description Create a new \code{KernelMatern} object.
        #' @param smoothness_factor Numeric: The smoothness factor of the kernel (1, 3 or 5)
        #' @param lengthscale KernelParameter: The lengthscale parameter instance of the kernel
        #' @param kernel_variance Numeric: The variance of the kernel
        #' @param jitter_value Numeric: The jitter value to add to the kernel matrix
        initialize = function(
            smoothness_factor = 5,
            lengthscale = KernelParameter$new(2),
            kernel_variance = 1,
            jitter_value = NULL
        ) {
            if (smoothness_factor %in% c(1, 3, 5) == FALSE) {
                stop('Smoothness factor should be one of the following values 1, 3 or 5')
            }
            super$initialize(kernel_variance, jitter_value)
            self$name <- paste0('Matern ', smoothness_factor, '/2 Kernel')
            self$smoothness_factor <- smoothness_factor
            self$lengthscale <- lengthscale
            self$lengthscale$set_kernel(self, 'lengthscale')
        },
        #' @description Method to the get the smoothness kernel function for a given integer smoothness factor
        #' @return The smoothness kernel function
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
        #' @description Method to compute the core kernel's covariance matrix
        #' @return The core kernel's covariance matrix
        core_kernel_fn = function() {
            temp_kernel <- (
                self$distance_matrix * sqrt(self$smoothness_factor) / self$lengthscale$value
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
KernelComposed <- R6::R6Class(
    'KernelComposed',
    inherit = Kernel,
    public = list(
        #' @field name The kernel's name
        name = '',
        #' @field parameters The parameters of the kernel (list of \code{KernelParameter})
        parameters = c(),
        #' @field left_kernel The left kernel to use for composition
        left_kernel = NULL,
        #' @field right_kernel The right kernel to use for composition
        right_kernel = NULL,
        #' @field composition_operation The operation to use for composition
        composition_operation = NULL,
        #' @field has_dist_matrix Identify if the kernel has a distance matrix or not
        has_dist_matrix = TRUE,
        #' @description Create a new \code{KernelComposed} object.
        #' @param left_kernel Kernel: The left kernel to use for composition
        #' @param right_kernel Kernel: The right kernel to use for composition
        #' @param new_name String: The name of the composed kernel
        #' @param composition_operation CompositionOps: The operation to use for composition
        initialize = function(
            left_kernel,
            right_kernel,
            new_name,
            composition_operation
        ) {
            composed_variance <- 1
            new_jitter_val <- max(
                left_kernel$jitter_value,
                right_kernel$jitter_value,
                TSR$get_default_jitter()
            )
            super$initialize(composed_variance, new_jitter_val)
            self$left_kernel <- left_kernel
            self$right_kernel <- right_kernel
            self$name <- new_name
            self$parameters <- c(
                left_kernel$parameters,
                right_kernel$parameters
            )
            self$composition_operation <- composition_operation
        },
        #' @description Method to compute the core kernel's covariance matrix
        #' @return The core kernel's covariance matrix
        core_kernel_fn = function() {
            if (self$composition_operation == CompositionOps$MUL) {
                return(self$left_kernel$core_kernel_fn() * self$right_kernel$core_kernel_fn())
            } else if (self$composition_operation == CompositionOps$ADD) {
                return(self$left_kernel$core_kernel_fn() + self$right_kernel$core_kernel_fn())
            } else {
                torch:::value_error('Composition operation is not implemented')
            }
        },
        #' @description Method to set the kernel's positions and compute the distance matrix
        #' @param positions_df Dataframe: The positions of the points in a dataframe format
        #' @return NULL, set the kernel's positions and compute the distance matrix
        set_positions = function(positions_df) {
            super$set_positions(positions_df)
            self$left_kernel$set_positions(positions_df)
            self$right_kernel$set_positions(positions_df)
        }
    )
)

#' @title R6 class for Kernels Composed via Addition
#'
#' @description R6 class automatically generated when
#' adding two kernels together.
#'
#' @export
KernelAddComposed <- R6::R6Class(
    'KernelAddComposed',
    inherit = KernelComposed,
    public = list(
        #' @description Create a new \code{KernelAddComposed} object.
        #' @param left_kernel Kernel: The left kernel to use for composition
        #' @param right_kernel Kernel: The right kernel to use for composition
        #' @param new_name String: The name of the composed kernel
        #' @return A new \code{KernelAddComposed} object.
        initialize = function(left_kernel, right_kernel, new_name) {
            super$initialize(left_kernel, right_kernel, new_name, CompositionOps$ADD)
        }
    )
)

#' @title R6 class for Kernels Composed via Multiplication
#'
#' @description R6 class automatically generated when
#' multiplying two kernels together.
#'
#' @export
KernelMulComposed <- R6::R6Class(
    'KernelMulComposed',
    inherit = KernelComposed,
    public = list(
        #' @description Create a new \code{KernelMulComposed} object.
        #' @param left_kernel Kernel: The left kernel to use for composition
        #' @param right_kernel Kernel: The right kernel to use for composition
        #' @param new_name String: The name of the composed kernel
        #' @return A new \code{KernelMulComposed} object.
        initialize = function(left_kernel, right_kernel, new_name) {
            super$initialize(left_kernel, right_kernel, new_name, CompositionOps$MUL)
        }
    )
)

#' @title Operator overloading for kernel addition
#' @description Operator overloading for kernel addition
#' @param k1 Kernel: The left kernel to use for composition
#' @param k2 Kernel: The right kernel to use for composition
#' @return A new \code{KernelAddComposed} object.
#'
#' @export
`+.Kernel` <- function(k1, k2) {
    composed_kernel <- KernelAddComposed$new(k1, k2, paste0(k1$name, ' + ', k2$name))
    return(composed_kernel)
}

#' @title Operator overloading for kernel multiplication
#' @description Operator overloading for kernel multiplication
#' @param k1 Kernel: The left kernel to use for composition
#' @param k2 Kernel: The right kernel to use for composition
#' @return A new \code{KernelMulComposed} object.
#'
#' @export
`*.Kernel` <- function(k1, k2) {
    composed_kernel <- KernelMulComposed$new(k1, k2, paste0(k1$name, ' * ', k2$name))
    return(composed_kernel)
}
