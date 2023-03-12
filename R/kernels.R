DEFAULT_LBOUND = log(1e-3)
DEFAULT_UBOUND = log(1e3)


KernelParameter <- R6::R6Class(
    public = list(
        value = 0,
        name = '',
        is_constant = FALSE,
        lower_bound = DEFAULT_LBOUND,
        upper_bound = DEFAULT_UBOUND,
        slice_sampling_scale = log(10),
        hparam_precision = 1.0,
        #' @field The kernel associated with the parameter (it is set at kernel instanciation
        kernel = NULL,
        #' @description Create a new \code{KernelParameter} object.
        #' @param value Numeric: The hyperparameter mean's prior value (Paper - \eqn{\phi}) or its constant value
        #' @param name String: The name of the paramater (used in logging and in kernel representation)
        #' @param is_constant Boolean: Says if the kernel parameter is constant or not (if it is constant, there is no sampling)
        #' @param lower_bound Numeric: The hyperparameter's minimal admissible value in sampling (Paper - \eqn{\phi_{min}})
        #' @param upper_bound Numeric: The hyperparameter's maximal admissible value in sampling (Paper - \eqn{\phi_{max}})
        #' @param slice_sampling_scale Numeric: The sampling range's amplitude (Paper - \eqn{\rho})
        #' @param hparam_precision Numeric: WHAT IS THAT? TODO
        #' @return A new \code{KernelParameter} object.
        initialize = function(
            value,
            name,
            is_constant=FALSE,
            lower_bound=DEFAULT_LBOUND,
            upper_bound=DEFAULT_UBOUND,
            slice_sampling_scale=log(10),
            hparam_precision=1.0
        ){
            self$value <- value
            self$name <- name
            self$lower_bound <- lower_bound
            self$upper_bound <- upper_bound
            self$is_constant <- is_constant
            self$slice_sampling_scale <- slice_sampling_scale
            self$hparam_precision <- hparam_precision
        },

        #' @description Set \code{Kernel} for a given \code{KernelParameter} object.
        #' @param kernel Kernel: The kernel to associate with the parameter
        #' @return NULL, set a new kernel for the parameter
        set_kernel = function(kernel) {
            self$kernel <- kernel
            if (!self$is_constant) {
                self$kernel$parameters <- c(self$kernel$parameters, self)
            }
        },

        get_full_name = function() {
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
        parameters = c(),
        kernel = NULL,
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

        get_name = function() {
            stop('get_name() is not implemented')
        },

        core_kernel_fn = function() {
            stop('core_kernel_fn() is not implemented')
        },

        add_jitter_to_kernel = function() {
            has_null_jitter <- is.null(self$jitter_value)
            if (!has_null_jitter && self$jitter_value == 0) {
                return()
            }
            jitter_val <- ifelse(has_null_jitter, tsr$default_jitter, self$jitter_value)
            self$kernel <- self$kernel + jitter_val * tsr$eye(nrow(self$kernel))
        },

        kernel_gen = function() {
            if (is.null(self$distance_matrix)) {
                stop('Set kernel distance via `set_distance_matrix` before kernel evaluation.')
            }
            self$kernel <- self$kernel_variance * self$core_kernel_fn()
            self$add_jitter_to_kernel()
            return(self$kernel)
        },

        set_distance_matrix = function(x=NULL, distance_matrix=NULL) {
            if (is.null(x) == is.null(distance_matrix)) {
                stop('Either `x` or `distance_matrix` must be provided.')
            } else if (!is.null(x)) {
                self$distance_matrix <- DistanceCalculator$get_matrix(x, self$distance_type)
            } else {
                self$distance_matrix <- distance_matrix
            }
        }
    )
)

KernelWhiteNoise <- R6::R6Class(
    inherit = Kernel,
    public = list(
        variance = NULL,
        distance_matrix = NULL,
        name = 'White Noise Kernel',
        initialize = function(
            variance = KernelParameter$new(1, 'variance', is_constant=TRUE),
            kernel_variance = 1,
            distance_type = DIST_TYPE$LINEAR,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$variance = variance
            self$variance$set_kernel(self)
        },
        core_kernel_fn = function() {
            return(tsr$eye(nrow(self$distance_matrix)))
        }
    )
)

KernelDotProduct <- R6::R6Class(
    inherit = Kernel,
    public = list(
        variance = NULL,
        distance_matrix = NULL,
        name = 'Linear Kernel',
        initialize = function(
            variance = KernelParameter$new(1, 'variance', lower_bound=0),
            kernel_variance = 1,
            jitter_value = NULL
        ) {
            distance_type = DIST_TYPE$DOT_PRODUCT
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$variance = variance
            self$variance$set_kernel(self)
        },
        core_kernel_fn = function() {
            return(self$distance_matrix + self$variance$value)
        }
    )
)

KernelSE <- R6::R6Class(
    inherit = Kernel,
    public = list(
        lengthscale = NULL,
        distance_matrix = NULL,
        name = 'Squared Exponential Kernel',
        initialize = function(
            lengthscale = KernelParameter$new(log(2), 'lengthscale'),
            kernel_variance = 1,
            distance_type = DIST_TYPE$LINEAR,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$lengthscale = lengthscale
            self$lengthscale$set_kernel(self)
        },
        core_kernel_fn = function() {
            return(torch::torch_exp(-self$distance_matrix^2 / (2 * self$lengthscale$value^2)))
        }
    )
)

KernelRQ <- R6::R6Class(
    inherit = Kernel,
    public = list(
        lengthscale = NULL,
        alpha = NULL,
        distance_matrix = NULL,
        name = 'Rational Quadratic Kernel',
        initialize = function(
            lengthscale = KernelParameter$new(log(2), 'lengthscale'),
            alpha = KernelParameter$new(log(2), 'alpha'),
            kernel_variance = 1,
            distance_type = DIST_TYPE$LINEAR,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$lengthscale = lengthscale
            self$lengthscale$set_kernel(self)
            self$alpha = alpha
            self$alpha$set_kernel(self)
        },
        core_kernel_fn = function() {
            return(
                1 + self$distance_matrix^2 / (2 * self$lengthscale$value^2 * self$alpha$value)
            ) ** -self$alpha$value
        }
    )
)

KernelPeriodic <- R6::R6Class(
    inherit = Kernel,
    public = list(
        lengthscale = NULL,
        period_length = NULL,
        distance_matrix = NULL,
        name = 'Periodic Kernel',
        initialize = function(
            lengthscale = KernelParameter$new(log(2), 'lengthscale'),
            period_length = KernelParameter$new(log(2), 'period length'),
            kernel_variance = 1,
            distance_type = DIST_TYPE$LINEAR,
            jitter_value = NULL
        ) {
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$lengthscale = lengthscale
            self$lengthscale$set_kernel(self)
            self$period_length = period_length
            self$period_length$set_kernel(self)
        },
        core_kernel_fn = function() {
            return(torch::torch_exp(
                -2 * torch::torch_sin(pi * self$distance_matrix / self$period_length$value)^2 /
                    self$lengthscale$value^2
            ))
        }
    )
)

KernelMatern <- R6::R6Class(
    inherit = Kernel,
    public = list(
        lengthscale = NULL,
        smoothness_factor = NULL,
        distance_matrix = NULL,
        name = 'Matern Kernel',
        initialize = function(
            smoothness_factor = 1,
            lengthscale = KernelParameter$new(log(2), 'lengthscale'),
            kernel_variance = 1,
            distance_type = DIST_TYPE$HAVERSINE,
            jitter_value = NULL
        ) {
            if (smoothness_factor %in% c(1, 3, 5) == FALSE) {
                stop('Smoothness factor should be one of the following values 1, 3 or 5')
            }
            super$initialize(kernel_variance, distance_type, jitter_value)
            self$smoothness_factor = smoothness_factor
            self$lengthscale = lengthscale
            self$lengthscale$set_kernel(self)
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

CompositionOps = list(
    'MUL' = 'MUL',
    'ADD' = 'ADD'
)

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
                tsr$default_jitter
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
        set_distance_matrix = function(x=NULL, distance_matrix=NULL) {
            super$set_distance_matrix(x, distance_matrix)
            self$left_kernel$set_distance_matrix(x, distance_matrix)
            self$right_kernel$set_distance_matrix(x, distance_matrix)
        }
    )
)

#' @export
`+.Kernel` = function(k1, k2) {
    composed_kernel = KernelComposed$new(
        k1,
        k2,
        paste0(k1$name, ' + ', k2$name),
        CompositionOps$ADD
    )
    return (composed_kernel)
}

#' @export
`*.Kernel` = function(e1, e2) {
    composed_kernel = KernelComposed$new(
        e1,
        e2,
        paste0(e1$name, ' * ', e2$name),
        CompositionOps$MUL
    )
    return (composed_kernel)
}
