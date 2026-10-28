# Default analysis-prior constructors for Bayesian t tests

#' Create the default two-sided analysis prior
#'
#' Constructs a zero-centred normal prior for the standardized effect with
#' support on the entire real line.
#'
#' @param sd Positive prior standard deviation.
#'
#' @return An object of class `"bf_ttest_prior"`.
#'
#' @examples
#' make_analysis_prior_two_sided()
#'
#' @export
make_analysis_prior_two_sided <- function(sd = 0.5) {
  make_prior(
    family = "normal",
    location = 0,
    scale = sd,
    support = c(-Inf, Inf),
    name = sprintf("N(0, %s^2)", format(sd, trim = TRUE))
  )
}

#' Create the default positive-direction analysis prior
#'
#' Constructs a zero-centred normal prior for the standardized effect,
#' truncated to positive effects. This is the default alternative prior for
#' `hypothesis = "greater.pointnull"` and the H1 prior for
#' `hypothesis = "greater.composite"`.
#'
#' @param sd Positive prior standard deviation.
#'
#' @return An object of class `"bf_ttest_prior"` with support `c(0, Inf)`.
#'
#' @examples
#' make_analysis_prior_greater(sd = 0.5)
#'
#' @export
make_analysis_prior_greater <- function(sd = 0.5) {
  make_prior(
    family = "normal",
    location = 0,
    scale = sd,
    support = c(0, Inf),
    name = sprintf("TN(0, %s^2; [0, Inf))", format(sd, trim = TRUE))
  )
}

#' Create the default non-positive analysis prior
#'
#' Constructs a zero-centred normal prior for the standardized effect,
#' truncated to non-positive effects. This is the default H0 prior for
#' `hypothesis = "greater.composite"`.
#'
#' @param sd Positive prior standard deviation.
#'
#' @return An object of class `"bf_ttest_prior"` with support `c(-Inf, 0)`.
#'
#' @examples
#' make_analysis_prior_less_equal(sd = 0.5)
#'
#' @export
make_analysis_prior_less_equal <- function(sd = 0.5) {
  make_prior(
    family = "normal",
    location = 0,
    scale = sd,
    support = c(-Inf, 0),
    name = sprintf("TN(0, %s^2; (-Inf, 0])", format(sd, trim = TRUE))
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Create a prior for the standardized effect
#'
#' @param family Prior family.
#' @param location Location parameter.
#' @param scale Scale parameter (> 0).
#' @param df Degrees of freedom for Student prior.
#' @param support Numeric vector `c(lower, upper)` defining support.
#' @param log_density Optional custom log-density.
#' @param density Optional custom density.
#' @param name Optional display name.
#'
#' @return Object of class `bf_ttest_prior`.
#' @export
make_prior <- function(
    family = c("normal", "student", "cauchy", "jzs", "custom"),
    location = 0,
    scale = 1,
    df = 1,
    support = c(-Inf, Inf),
    log_density = NULL,
    density = NULL,
    name = NULL
) {
  family <- match.arg(family)
  
  if (!is.numeric(location) || length(location) != 1L || !is.finite(location)) {
    stop("`location` must be a finite scalar.", call. = FALSE)
  }
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("`scale` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) {
    stop("`df` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(support) || length(support) != 2L || support[1L] >= support[2L]) {
    stop("`support` must be a numeric vector `c(lower, upper)` with lower < upper.",
         call. = FALSE)
  }
  if (family == "custom" && is.null(log_density) && is.null(density)) {
    stop("For `family = 'custom'`, provide `log_density` or `density`.",
         call. = FALSE)
  }
  
  structure(
    list(
      family = family,
      location = as.numeric(location),
      scale = as.numeric(scale),
      df = as.numeric(df),
      support = as.numeric(support),
      log_density = log_density,
      density = density,
      name = name %||% family
    ),
    class = "bf_ttest_prior"
  )
}

#' Point-mass prior at a fixed effect
#'
#' @param effect Scalar effect value.
#'
#' @return Object of class `bf_ttest_point_prior`.
#' @export
make_point_prior <- function(effect = 0) {
  if (!is.numeric(effect) || length(effect) != 1L || !is.finite(effect)) {
    stop("`effect` must be one finite number.", call. = FALSE)
  }
  
  structure(
    list(effect = as.numeric(effect)),
    class = "bf_ttest_point_prior"
  )
}

#' Create a zero-centred normal design-prior object
#'
#' Convenience wrapper around [make_prior()] for a normal prior on the
#' standardized effect, with support on the entire real line.
#'
#' @param sd Positive prior standard deviation.
#'
#' @return An object of class `"bf_ttest_prior"`.
#'
#' @examples
#' make_design_prior_normal(sd = 0.5)
#'
#' @export
make_design_prior_normal <- function(sd = 0.5) {
  make_prior(
    family = "normal",
    location = 0,
    scale = sd,
    support = c(-Inf, Inf),
    name = sprintf(
      "N(0, %s^2)",
      format(sd, trim = TRUE)
    )
  )
}

#' Backward-compatible alias for `make_point_prior()`
#'
#' @param effect Point effect value.
#'
#' @return Object of class `bf_ttest_point_prior`.
#' @export
make_design_prior_point <- function(effect = 0) {
  make_point_prior(effect = effect)
}

#' @export
print.bf_ttest_prior <- function(x, ...) {
  cat("Prior:", x$name, "\n")
  cat("  family:   ", x$family, "\n", sep = "")
  cat("  location: ", format(x$location), "\n", sep = "")
  cat("  scale:    ", format(x$scale), "\n", sep = "")
  cat("  support:  [", format(x$support[1L]), ", ", format(x$support[2L]), "]\n",
      sep = "")
  invisible(x)
}

#' @export
print.bf_ttest_point_prior <- function(x, ...) {
  cat("Point prior: delta = ", format(x$effect), "\n", sep = "")
  invisible(x)
}

.validate_prior_bf_ttest <- function(prior, allow_point = TRUE) {
  ok_cont <- inherits(prior, "bf_ttest_prior")
  ok_point <- allow_point && inherits(prior, "bf_ttest_point_prior")
  if (!(ok_cont || ok_point)) {
    stop("Invalid prior object.", call. = FALSE)
  }
  invisible(TRUE)
}

prior_mass <- function(prior, rel.tol = 1e-10, subdivisions = 500L) {
  stopifnot(inherits(prior, "bf_ttest_prior"))
  
  if (identical(prior$support, c(-Inf, Inf))) {
    return(1)
  }
  
  base_pdf <- function(x) {
    switch(
      prior$family,
      normal = stats::dnorm(x, mean = prior$location, sd = prior$scale),
      student = stats::dt((x - prior$location) / prior$scale, df = prior$df) / prior$scale,
      cauchy = stats::dcauchy(x, location = prior$location, scale = prior$scale),
      jzs = stats::dcauchy(x, location = prior$location, scale = prior$scale),
      custom = if (!is.null(prior$density)) prior$density(x) else exp(prior$log_density(x))
    )
  }
  
  stats::integrate(
    base_pdf,
    lower = prior$support[1L],
    upper = prior$support[2L],
    rel.tol = rel.tol,
    subdivisions = subdivisions,
    stop.on.error = FALSE
  )$value
}

prior_logpdf <- function(delta, prior, rel.tol = 1e-10) {
  stopifnot(inherits(prior, "bf_ttest_prior"))
  
  out <- rep(-Inf, length(delta))
  in_support <- is.finite(delta) &
    delta >= prior$support[1L] &
    delta <= prior$support[2L]
  
  if (!any(in_support)) return(out)
  
  x <- delta[in_support]
  
  log_base <- switch(
    prior$family,
    normal = stats::dnorm(x, mean = prior$location, sd = prior$scale, log = TRUE),
    student = stats::dt((x - prior$location) / prior$scale, df = prior$df, log = TRUE) -
      log(prior$scale),
    cauchy = stats::dcauchy(x, location = prior$location, scale = prior$scale, log = TRUE),
    jzs = stats::dcauchy(x, location = prior$location, scale = prior$scale, log = TRUE),
    custom = if (!is.null(prior$log_density)) {
      prior$log_density(x)
    } else {
      log(prior$density(x))
    }
  )
  
  if (is.finite(prior$support[1L]) || is.finite(prior$support[2L])) {
    norm_const <- prior_mass(prior, rel.tol = rel.tol)
    log_base <- log_base - log(norm_const)
  }
  
  out[in_support] <- log_base
  out
}

prior_pdf <- function(delta, prior, rel.tol = 1e-10) {
  exp(prior_logpdf(delta, prior, rel.tol = rel.tol))
}