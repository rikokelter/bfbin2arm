# ============================================================
# Deterministic Bayes factor calculations for t-tests
# ============================================================

pnt_stable <- function(q, df, ncp, lower.tail = TRUE, log.p = FALSE) {
  q <- as.numeric(q)
  ncp <- as.numeric(ncp)
  
  if (length(q) == 1L && length(ncp) > 1L) q <- rep(q, length(ncp))
  if (length(ncp) == 1L && length(q) > 1L) ncp <- rep(ncp, length(q))
  
  if (length(q) != length(ncp)) {
    stop("`q` and `ncp` must have compatible lengths.", call. = FALSE)
  }
  
  out <- suppressWarnings(
    stats::pt(q = q, df = df, ncp = ncp, lower.tail = lower.tail, log.p = log.p)
  )
  
  bad <- !is.finite(out) | is.na(out)
  
  if (!any(bad)) {
    return(out)
  }
  
  fallback_one <- function(qq, lambda, lower.tail = TRUE, log.p = FALSE) {
    integrand <- function(u) {
      stats::pnorm(qq * sqrt(u) - lambda) *
        stats::dchisq(df * u, df = df) * df
    }
    
    value <- tryCatch(
      stats::integrate(
        integrand,
        lower = 0,
        upper = Inf,
        rel.tol = 1e-8,
        subdivisions = 1000L,
        stop.on.error = FALSE
      )$value,
      error = function(e) NA_real_
    )
    
    if (!is.finite(value)) {
      return(NA_real_)
    }
    
    value <- min(max(value, 0), 1)
    
    if (!lower.tail) {
      value <- 1 - value
    }
    
    value <- min(max(value, 0), 1)
    
    if (log.p) {
      if (value <= 0) return(-Inf)
      return(log(value))
    }
    
    value
  }
  
  out[bad] <- vapply(
    which(bad),
    function(i) fallback_one(
      qq = q[i],
      lambda = ncp[i],
      lower.tail = lower.tail,
      log.p = log.p
    ),
    numeric(1)
  )
  
  out
}

make_log_bf01_greater_pointnull_fun <- function(
    df,
    n_eff,
    prior_h1,
    rel.tol = 1e-9,
    subdivisions = 500L,
    cache = TRUE
) {
  .validate_prior_bf_ttest(prior_h1, allow_point = FALSE)
  
  store <- if (isTRUE(cache)) {
    new.env(hash = TRUE, parent = emptyenv())
  } else {
    NULL
  }
  
  one <- function(tt) {
    key <- sprintf("%.17g", tt)
    
    if (!is.null(store) && exists(key, envir = store, inherits = FALSE)) {
      return(get(key, envir = store, inherits = FALSE))
    }
    
    log_m0 <- stats::dt(tt, df = df, log = TRUE)
    
    log_m1 <- log_m1_t_generic(
      t = tt,
      df = df,
      n_eff = n_eff,
      prior = prior_h1,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    )
    
    out <- log_m0 - log_m1
    
    if (!is.null(store)) {
      assign(key, out, envir = store)
    }
    
    out
  }
  
  function(t) {
    t <- as.numeric(t)
    if (any(!is.finite(t))) {
      stop("`t` must contain only finite numeric values.", call. = FALSE)
    }
    vapply(t, one, numeric(1))
  }
}

make_log_bf01_composite_fun <- function(
    df,
    n_eff,
    prior_h0,
    prior_h1,
    rel.tol = 1e-9,
    subdivisions = 500L
) {
  .validate_prior_bf_ttest(prior_h0, allow_point = FALSE)
  .validate_prior_bf_ttest(prior_h1, allow_point = FALSE)
  
  function(t) {
    vapply(
      t,
      function(tt) {
        log_m0 <- log_m1_t_generic(
          t = tt,
          df = df,
          n_eff = n_eff,
          prior = prior_h0,
          rel.tol = rel.tol,
          subdivisions = subdivisions
        )
        
        log_m1 <- log_m1_t_generic(
          t = tt,
          df = df,
          n_eff = n_eff,
          prior = prior_h1,
          rel.tol = rel.tol,
          subdivisions = subdivisions
        )
        
        log_m0 - log_m1
      },
      numeric(1)
    )
  }
}

logspace_add <- function(logx) {
  m <- max(logx)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(logx - m)))
}

n_eff_one_sample <- function(n) {
  n
}

n_eff_two_sample <- function(n1, n2) {
  1 / (1 / n1 + 1 / n2)
}

log_m0_t <- function(t, df) {
  stats::dt(t, df = df, log = TRUE)
}

log_d_nct_stable <- function(
    t,
    df,
    ncp,
    rel.tol = 1e-8,
    subdivisions = 200L,
    direct_cutoff_t = 6,
    direct_cutoff_ncp = 6
) {
  t <- as.numeric(t)
  ncp <- as.numeric(ncp)
  
  if (length(t) != 1L) {
    return(vapply(
      t,
      function(tt) {
        log_d_nct_stable(
          t = tt,
          df = df,
          ncp = ncp,
          rel.tol = rel.tol,
          subdivisions = subdivisions,
          direct_cutoff_t = direct_cutoff_t,
          direct_cutoff_ncp = direct_cutoff_ncp
        )
      },
      numeric(length(ncp))
    ))
  }
  
  direct <- suppressWarnings(
    stats::dt(t, df = df, ncp = ncp, log = TRUE)
  )
  out <- as.numeric(direct)
  
  lim <- .log_u_limits_nct(df)
  
  fallback_one <- function(lambda) {
    integrand <- function(v) {
      v <- as.numeric(v)
      u <- exp(v)
      
      log_kernel <-
        stats::dnorm(t * sqrt(u) - lambda, log = TRUE) +
        0.5 * log(u) +
        stats::dchisq(df * u, df = df, log = TRUE) +
        log(df) + v
      
      ans <- exp(log_kernel)
      ans[!is.finite(ans) | ans < 0] <- 0
      ans
    }
    
    value <- tryCatch(
      stats::integrate(
        f = integrand,
        lower = lim[1L],
        upper = lim[2L],
        rel.tol = rel.tol,
        subdivisions = as.integer(subdivisions),
        stop.on.error = FALSE
      )$value,
      error = function(e) NA_real_
    )
    
    if (!is.finite(value) || value <= 0) {
      return(-Inf)
    }
    
    log(value)
  }
  
  force_fallback <- (!is.finite(out)) |
    (abs(t) >= direct_cutoff_t) |
    (abs(ncp) >= direct_cutoff_ncp)
  
  if (any(force_fallback)) {
    out[force_fallback] <- vapply(ncp[force_fallback], fallback_one, numeric(1))
  }
  
  out
}

.log_u_limits_nct <- function(df, eps = 1e-12) {
  if (length(df) != 1L || !is.finite(df) || df <= 0) {
    stop("`df` must be one positive finite number.", call. = FALSE)
  }
  
  lower <- log(stats::qchisq(eps, df = df) / df)
  upper <- log(stats::qchisq(eps, df = df, lower.tail = FALSE) / df)
  
  c(lower, upper)
}

log_integrand_m1 <- function(delta, t, df, n_eff, prior) {
  log_d_nct_stable(
    t = t,
    df = df,
    ncp = sqrt(n_eff) * delta
  ) + prior_logpdf(delta, prior)
}

log_m1_t_generic <- function(
    t,
    df,
    n_eff,
    prior,
    rel.tol = 1e-9,
    subdivisions = 500L,
    tail_sd = 8
) {
  if (length(t) != 1L || !is.finite(t)) {
    stop("`t` must be one finite numeric value.", call. = FALSE)
  }
  
  .validate_prior_bf_ttest(prior, allow_point = FALSE)
  
  lower <- prior$support[1L]
  upper <- prior$support[2L]
  
  if (identical(prior$family, "normal")) {
    lower <- max(lower, prior$location - tail_sd * prior$scale)
    upper <- min(upper, prior$location + tail_sd * prior$scale)
  }
  
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop(
      "Could not construct finite integration bounds for the analysis prior.",
      call. = FALSE
    )
  }
  
  log_kernel <- function(delta) {
    log_d_nct_stable(
      t = t,
      df = df,
      ncp = sqrt(n_eff) * delta,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    ) + prior_logpdf(delta, prior)
  }
  
  opt <- stats::optimize(
    f = function(delta) -log_kernel(delta),
    interval = c(lower, upper),
    maximum = FALSE,
    tol = .Machine$double.eps^0.25
  )
  
  delta_mode <- opt$minimum
  log_scale <- log_kernel(delta_mode)
  
  if (!is.finite(log_scale)) {
    return(-Inf)
  }
  
  scaled_integrand <- function(delta) {
    z <- log_kernel(delta) - log_scale
    ans <- exp(z)
    ans[!is.finite(ans)] <- 0
    ans
  }
  
  value <- tryCatch(
    stats::integrate(
      f = scaled_integrand,
      lower = lower,
      upper = upper,
      rel.tol = rel.tol,
      subdivisions = as.integer(subdivisions),
      stop.on.error = FALSE
    )$value,
    error = function(e) NA_real_
  )
  
  if (!is.finite(value) || value <= 0) {
    return(-Inf)
  }
  
  log_scale + log(value)
}

log_bf01_t_generic <- function(
    t,
    df,
    n_eff,
    prior,
    rel.tol = 1e-9,
    subdivisions = 500L
) {
  t <- as.numeric(t)
  
  if (any(!is.finite(t))) {
    stop("`t` must contain only finite numeric values.", call. = FALSE)
  }
  
  vapply(
    t,
    function(tt) {
      log_m0 <- stats::dt(tt, df = df, log = TRUE)
      
      log_m1 <- log_m1_t_generic(
        t = tt,
        df = df,
        n_eff = n_eff,
        prior = prior,
        rel.tol = rel.tol,
        subdivisions = subdivisions
      )
      
      log_m0 - log_m1
    },
    numeric(1)
  )
}

bf01_t_generic <- function(...) {
  exp(log_bf01_t_generic(...))
}

make_log_bf01_fun <- function(
    df, n_eff, prior,
    rel.tol = 1e-9,
    subdivisions = 500L
) {
  force(df)
  force(n_eff)
  force(prior)
  force(rel.tol)
  force(subdivisions)
  
  function(t) {
    log_bf01_t_generic(
      t = t,
      df = df,
      n_eff = n_eff,
      prior = prior,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    )
  }
}