#' Bayesian and frequentist operating characteristics for a single-arm two-stage BF design
#'
#' Computes naive fixed-sample and corrected two-stage operating characteristics
#' for a single-arm binomial design with one interim analysis for futility. The
#' Bayes factor is oriented as `BF01`, so efficacy corresponds to small values
#' (`BF01 <= k`) and futility corresponds to large values (`BF01 >= kf`).
#'
#' Bayesian operating characteristics are computed under separate design priors:
#' - for `type = "direction"`, Bayesian power averages over `p > p0` under the
#'   H1 design prior truncated to `(p0, 1]`, and Bayesian type-I error averages
#'   over `p <= p0` under the H0 design prior truncated to `[0, p0]`;
#' - for `type = "point"`, Bayesian power averages under the H1 design prior on
#'   `(0, 1)`, and Bayesian type-I error is evaluated at the point null `p = p0`.
#'
#' If `dp` is supplied, additional frequentist power under `H1` is computed at
#' the fixed point alternative `p = dp`. Frequentist type-I error is computed at
#' `p = p0`.
#'
#' @param n1 Integer scalar. Interim sample size.
#' @param n2 Integer scalar. Final sample size, with `n1 < n2`.
#' @param k Numeric scalar. Efficacy threshold on the `BF01` scale.
#' @param kf Numeric scalar. Futility threshold on the `BF01` scale.
#' @param p0 Numeric scalar in `(0, 1)`. Null response probability.
#' @param a0,b0 Numeric scalars. Beta analysis-prior parameters under H0.
#' @param a1,b1 Numeric scalars. Beta analysis-prior parameters under H1.
#' @param da0,db0 Numeric scalars. Beta design-prior parameters under H0.
#' @param da1,db1 Numeric scalars. Beta design-prior parameters under H1.
#' @param dp Optional numeric scalar in `(0,1)`. If supplied, frequentist power
#'   under `H1` is computed at `p = dp`.
#' @param type Character string. One of `\"point\"` or `\"direction\"`.
#' @param k_ce Optional numeric scalar greater than 1. Threshold for compelling
#'   evidence in favour of `H0` on the `BF01` scale.
#' @param grid_size Integer number of grid points used for numerical averaging.
#'
#' @return A list with Bayesian and frequentist operating characteristics.
#' @export
powerbinbf01seq <- function(
    n1, n2, k, kf, p0,
    a0 = 1, b0 = 1,
    a1 = 1, b1 = 1,
    da0 = 1, db0 = 1,
    da1 = 1, db1 = 1,
    dp = NA_real_,
    type = c("point", "direction"),
    k_ce = NULL,
    grid_size = 801L
) {
  type <- match.arg(type)
  
  if (!is.numeric(n1) || !is.numeric(n2) || length(n1) != 1L || length(n2) != 1L) {
    stop("'n1' and 'n2' must be numeric scalars.")
  }
  if (n1 < 1L || n2 < 1L || n1 != floor(n1) || n2 != floor(n2)) {
    stop("'n1' and 'n2' must be positive integers.")
  }
  if (n1 >= n2) {
    stop("Need n1 < n2.")
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("'k' must be a positive number (BF01 threshold for efficacy).")
  }
  if (!is.numeric(kf) || length(kf) != 1L || !is.finite(kf) || kf <= 1) {
    stop("'kf' must be > 1 (BF01 threshold for futility).")
  }
  if (!is.numeric(p0) || length(p0) != 1L || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
    stop("'p0' must be in (0,1).")
  }
  if (!is.numeric(a0) || !is.numeric(b0) || a0 <= 0 || b0 <= 0) {
    stop("'a0' and 'b0' must be positive.")
  }
  if (!is.numeric(a1) || !is.numeric(b1) || a1 <= 0 || b1 <= 0) {
    stop("'a1' and 'b1' must be positive.")
  }
  if (!is.numeric(da0) || !is.numeric(db0) || da0 <= 0 || db0 <= 0) {
    stop("'da0' and 'db0' must be positive.")
  }
  if (!is.numeric(da1) || !is.numeric(db1) || da1 <= 0 || db1 <= 0) {
    stop("'da1' and 'db1' must be positive.")
  }
  if (!is.na(dp) && (!is.numeric(dp) || length(dp) != 1L || !is.finite(dp) || dp <= 0 || dp >= 1)) {
    stop("'dp' must be NA or a number in (0,1).")
  }
  if (!is.null(k_ce)) {
    if (!is.numeric(k_ce) || length(k_ce) != 1L || !is.finite(k_ce) || k_ce <= 1) {
      stop("'k_ce' must be a number > 1 (BF01 threshold for compelling evidence).")
    }
  }
  if (!is.numeric(grid_size) || length(grid_size) != 1L || grid_size < 101L) {
    stop("'grid_size' must be an integer >= 101.")
  }
  
  grid_size <- as.integer(grid_size)
  
  x1_vals <- 0:n1
  y_vals  <- 0:n2
  z_vals  <- 0:(n2 - n1)
  
  bf_int <- singlearm_bf01(
    x = x1_vals,
    n = n1,
    p0 = p0,
    a0 = a0, b0 = b0,
    a1 = a1, b1 = b1,
    type = type
  )
  
  bf_fin <- singlearm_bf01(
    x = y_vals,
    n = n2,
    p0 = p0,
    a0 = a0, b0 = b0,
    a1 = a1, b1 = b1,
    type = type
  )
  
  fut_int <- bf_int >= kf
  eff_fin <- bf_fin <= k
  fut_x   <- x1_vals[fut_int]
  
  calc_two_stage_at_p <- function(p) {
    pmf_int <- dbinom(x1_vals, size = n1, prob = p)
    pmf_fin <- dbinom(y_vals,  size = n2, prob = p)
    
    pnaive  <- sum(pmf_fin[eff_fin])
    pintfut <- sum(pmf_int[fut_int])
    
    erased <- 0
    for (x in fut_x) {
      for (z in z_vals) {
        y <- x + z
        if (y >= 0 && y <= n2 && bf_fin[y + 1L] <= k) {
          erased <- erased +
            dbinom(x, size = n1, prob = p) *
            dbinom(z, size = n2 - n1, prob = p)
        }
      }
    }
    
    pfineff <- pnaive - erased
    nexp <- n1 * pintfut + n2 * (1 - pintfut)
    
    list(
      pnaive = pnaive,
      pintfut = pintfut,
      erased = erased,
      pfineff = pfineff,
      nexp = nexp
    )
  }
  
  integrate_region <- function(lower, upper, da, db, fun) {
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      return(NA_real_)
    }
    
    mass <- pbeta(upper, shape1 = da, shape2 = db) -
      pbeta(lower, shape1 = da, shape2 = db)
    
    if (!is.finite(mass) || mass <= 0) {
      return(NA_real_)
    }
    
    p_grid <- seq(lower, upper, length.out = grid_size)
    dens <- dbeta(p_grid, shape1 = da, shape2 = db)
    dens[p_grid < lower | p_grid > upper] <- 0
    
    if (!all(is.finite(dens)) || sum(dens) <= 0) {
      return(NA_real_)
    }
    
    w <- dens / sum(dens)
    vals <- vapply(p_grid, fun, numeric(1))
    sum(vals * w)
  }
  
  integrate_region_oc <- function(lower, upper, da, db) {
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      return(list(
        pnaive = NA_real_,
        pintfut = NA_real_,
        erased = NA_real_,
        pfineff = NA_real_,
        nexp = NA_real_
      ))
    }
    
    mass <- pbeta(upper, shape1 = da, shape2 = db) -
      pbeta(lower, shape1 = da, shape2 = db)
    
    if (!is.finite(mass) || mass <= 0) {
      return(list(
        pnaive = NA_real_,
        pintfut = NA_real_,
        erased = NA_real_,
        pfineff = NA_real_,
        nexp = NA_real_
      ))
    }
    
    p_grid <- seq(lower, upper, length.out = grid_size)
    dens <- dbeta(p_grid, shape1 = da, shape2 = db)
    dens[p_grid < lower | p_grid > upper] <- 0
    w <- dens / sum(dens)
    
    vals <- lapply(p_grid, calc_two_stage_at_p)
    
    avg_component <- function(name) {
      v <- vapply(vals, function(z) z[[name]], numeric(1))
      sum(v * w)
    }
    
    list(
      pnaive = avg_component("pnaive"),
      pintfut = avg_component("pintfut"),
      erased = avg_component("erased"),
      pfineff = avg_component("pfineff"),
      nexp = avg_component("nexp")
    )
  }
  
  ## Bayesian H1
  if (type == "direction") {
    lower_h1 <- p0
    upper_h1 <- 1
    
    if (lower_h1 >= upper_h1) {
      bayes_H1 <- list(
        pnaive = NA_real_,
        pintfut = NA_real_,
        erased = NA_real_,
        pfineff = NA_real_,
        nexp = NA_real_
      )
    } else {
      lower_h1_open <- min(upper_h1 - 1e-10, lower_h1 + 1e-10)
      bayes_H1 <- integrate_region_oc(
        lower = lower_h1_open, upper = upper_h1,
        da = da1, db = db1
      )
    }
  } else {
    pmf_int_H1 <- singlearm_priorpred_pmf(
      x = x1_vals, n = n1,
      dp = NA_real_,
      da = da1, db = db1,
      dl = 0, du = 1
    )
    
    pmf_fin_H1 <- singlearm_priorpred_pmf(
      x = y_vals, n = n2,
      dp = NA_real_,
      da = da1, db = db1,
      dl = 0, du = 1
    )
    
    pfin_H1_naive <- sum(pmf_fin_H1[eff_fin])
    pint_fut_H1   <- sum(pmf_int_H1[fut_int])
    
    erased_H1 <- 0
    for (x in fut_x) {
      for (z in z_vals) {
        y <- x + z
        if (y >= 0 && y <= n2 && bf_fin[y + 1L] <= k) {
          erased_H1 <- erased_H1 + dbinbin_truncbeta(
            x = x,
            z = z,
            n1 = n1,
            n2 = n2,
            a = da1,
            b = db1,
            lower = 0,
            upper = 1
          )
        }
      }
    }
    
    bayes_H1 <- list(
      pnaive = pfin_H1_naive,
      pintfut = pint_fut_H1,
      erased = erased_H1,
      pfineff = pfin_H1_naive - erased_H1,
      nexp = n1 * pint_fut_H1 + n2 * (1 - pint_fut_H1)
    )
  }
  
  ## Bayesian H0
  if (type == "direction") {
    lower_h0 <- 0
    upper_h0 <- p0
    
    if (lower_h0 >= upper_h0) {
      bayes_H0 <- list(
        pnaive = NA_real_,
        pintfut = NA_real_,
        erased = NA_real_,
        pfineff = NA_real_,
        nexp = NA_real_
      )
    } else {
      bayes_H0 <- integrate_region_oc(
        lower = lower_h0, upper = upper_h0,
        da = da0, db = db0
      )
    }
  } else {
    bayes_H0 <- calc_two_stage_at_p(p0)
  }
  
  ## Frequentist H1 at dp
  freq_H1 <- NULL
  if (!is.na(dp)) {
    freq_H1 <- calc_two_stage_at_p(dp)
  }
  
  ## Frequentist H0 at p0
  freq_H0 <- calc_two_stage_at_p(p0)
  freq0_argmax <- p0
  
  ## Compelling evidence for H0
  pce0_naive <- NA_real_
  pce0_corr  <- NA_real_
  
  if (!is.null(k_ce)) {
    ce_fin <- bf_fin >= k_ce
    
    if (type == "direction") {
      pce0_naive <- integrate_region(
        lower = 0,
        upper = p0,
        da = da0,
        db = db0,
        fun = function(p) {
          pmf_fin <- dbinom(y_vals, size = n2, prob = p)
          sum(pmf_fin[ce_fin])
        }
      )
      
      pce0_corr <- integrate_region(
        lower = 0,
        upper = p0,
        da = da0,
        db = db0,
        fun = function(p) {
          pmf_fin <- dbinom(y_vals, size = n2, prob = p)
          pce_naive_p <- sum(pmf_fin[ce_fin])
          
          delta_ce0 <- 0
          for (x in 0:n1) {
            bf_x <- bf_int[x + 1L]
            if (bf_x >= k_ce) {
              for (z in z_vals) {
                y <- x + z
                if (y >= 0 && y <= n2 && bf_fin[y + 1L] < k_ce) {
                  delta_ce0 <- delta_ce0 +
                    dbinom(x, size = n1, prob = p) *
                    dbinom(z, size = n2 - n1, prob = p)
                }
              }
            }
          }
          
          pce_naive_p + delta_ce0
        }
      )
    } else {
      pmf_fin_H0 <- dbinom(y_vals, size = n2, prob = p0)
      pce0_naive <- sum(pmf_fin_H0[ce_fin])
      
      delta_ce0 <- 0
      for (x in 0:n1) {
        bf_x <- bf_int[x + 1L]
        if (bf_x >= k_ce) {
          for (z in z_vals) {
            y <- x + z
            if (y >= 0 && y <= n2 && bf_fin[y + 1L] < k_ce) {
              delta_ce0 <- delta_ce0 +
                dbinom(x, size = n1, prob = p0) *
                dbinom(z, size = n2 - n1, prob = p0)
            }
          }
        }
      }
      
      pce0_corr <- pce0_naive + delta_ce0
    }
  }
  
  list(
    n1 = n1,
    n2 = n2,
    k = k,
    kf = kf,
    p0 = p0,
    dp = dp,
    type = type,
    bf_interim = bf_int,
    bf_final = bf_fin,
    
    ## Bayesian H1
    pintfut = bayes_H1$pintfut,
    pfineff = bayes_H1$pfineff,
    pnaive = bayes_H1$pnaive,
    perased = bayes_H1$erased,
    nexp = bayes_H1$nexp,
    
    ## Bayesian H0
    pintfut0 = bayes_H0$pintfut,
    pfineff0 = bayes_H0$pfineff,
    pnaive0 = bayes_H0$pnaive,
    perased0 = bayes_H0$erased,
    nexp0 = bayes_H0$nexp,
    
    ## Frequentist H1 at dp
    pintfut_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pintfut,
    pfineff_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pfineff,
    pnaive_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pnaive,
    perased_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$erased,
    nexp_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$nexp,
    
    ## Frequentist H0 at p0
    pintfut_freq0 = freq_H0$pintfut,
    pfineff_freq0 = freq_H0$pfineff,
    pnaive_freq0 = freq_H0$pnaive,
    perased_freq0 = freq_H0$erased,
    nexp_freq0 = freq_H0$nexp,
    freq0_argmax = freq0_argmax,
    
    k_ce = k_ce,
    pce0_naive = pce0_naive,
    pce0_corr = pce0_corr
  )
}

#' Bayesian and frequentist operating characteristics for a fixed-sample
#' single-arm BF design
#'
#' Computes operating characteristics for a genuine fixed-sample single-arm
#' binomial design with final efficacy decision based on `BF01 <= k`.
#'
#' Bayesian operating characteristics are computed under separate design priors:
#' - for `type = "direction"`, Bayesian power averages over `p > p0` under the
#'   H1 design prior truncated to `(p0, 1]`, Bayesian type-I error averages over
#'   `p <= p0` under the H0 design prior truncated to `[0, p0]`, and CE(H0) is
#'   averaged over the same truncated H0 design prior;
#' - for `type = "point"`, Bayesian power averages under the H1 design prior on
#'   `(0, 1)`, Bayesian type-I error is evaluated at the point null `p = p0`,
#'   and CE(H0) is also evaluated at `p = p0`.
#'
#' @param n Integer scalar. Total sample size.
#' @param k Numeric scalar. Efficacy threshold on the `BF01` scale.
#' @param p0 Numeric scalar in `(0, 1)`. Null response probability.
#' @param a0,b0 Numeric scalars. Beta analysis-prior parameters under H0.
#' @param a1,b1 Numeric scalars. Beta analysis-prior parameters under H1.
#' @param da0,db0 Numeric scalars. Beta design-prior parameters under H0.
#' @param da1,db1 Numeric scalars. Beta design-prior parameters under H1.
#' @param dp Optional numeric scalar in `(0,1)`. If supplied, frequentist power
#'   under `H1` is computed at `p = dp`.
#' @param type Character string. One of `\"point\"` or `\"direction\"`.
#' @param k_ce Optional numeric scalar greater than 1. Threshold for compelling
#'   evidence in favour of `H0` on the `BF01` scale.
#' @param grid_size Integer number of grid points used for numerical averaging.
#'
#' @return A list with Bayesian and frequentist operating characteristics for the
#' fixed-sample design.
#' @export
powerbinbf01_fixed <- function(
    n,
    k,
    p0,
    a0 = 1, b0 = 1,
    a1 = 1, b1 = 1,
    da0 = 1, db0 = 1,
    da1 = 1, db1 = 1,
    dp = NA_real_,
    type = c("point", "direction"),
    k_ce = NULL,
    grid_size = 801L
) {
  type <- match.arg(type)
  
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1L || n != floor(n)) {
    stop("'n' must be a positive integer.")
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("'k' must be a positive number (BF01 threshold for efficacy).")
  }
  if (!is.numeric(p0) || length(p0) != 1L || !is.finite(p0) ||
      p0 <= 0 || p0 >= 1) {
    stop("'p0' must be in (0,1).")
  }
  if (!is.numeric(a0) || !is.numeric(b0) || a0 <= 0 || b0 <= 0) {
    stop("'a0' and 'b0' must be positive.")
  }
  if (!is.numeric(a1) || !is.numeric(b1) || a1 <= 0 || b1 <= 0) {
    stop("'a1' and 'b1' must be positive.")
  }
  if (!is.numeric(da0) || !is.numeric(db0) || da0 <= 0 || db0 <= 0) {
    stop("'da0' and 'db0' must be positive.")
  }
  if (!is.numeric(da1) || !is.numeric(db1) || da1 <= 0 || db1 <= 0) {
    stop("'da1' and 'db1' must be positive.")
  }
  if (!is.na(dp) && (!is.numeric(dp) || length(dp) != 1L || !is.finite(dp) ||
                     dp <= 0 || dp >= 1)) {
    stop("'dp' must be NA or a number in (0,1).")
  }
  if (!is.null(k_ce)) {
    if (!is.numeric(k_ce) || length(k_ce) != 1L || !is.finite(k_ce) || k_ce <= 1) {
      stop("'k_ce' must be a number > 1 (BF01 threshold for compelling evidence).")
    }
  }
  if (!is.numeric(grid_size) || length(grid_size) != 1L || grid_size < 101L) {
    stop("'grid_size' must be an integer >= 101.")
  }
  
  grid_size <- as.integer(grid_size)
  x_vals <- 0:n
  
  bf_fin <- singlearm_bf01(
    x = x_vals,
    n = n,
    p0 = p0,
    a0 = a0, b0 = b0,
    a1 = a1, b1 = b1,
    type = type
  )
  
  eff_fin <- bf_fin <= k
  
  calc_fixed_at_p <- function(p) {
    pmf_fin <- dbinom(x_vals, size = n, prob = p)
    p_eff <- sum(pmf_fin[eff_fin])
    
    list(
      pnaive = p_eff,
      pintfut = 0,
      erased = 0,
      pfineff = p_eff,
      nexp = n
    )
  }
  
  integrate_region <- function(lower, upper, da, db, fun) {
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      return(NA_real_)
    }
    
    mass <- pbeta(upper, shape1 = da, shape2 = db) -
      pbeta(lower, shape1 = da, shape2 = db)
    
    if (!is.finite(mass) || mass <= 0) {
      return(NA_real_)
    }
    
    p_grid <- seq(lower, upper, length.out = grid_size)
    dens <- dbeta(p_grid, shape1 = da, shape2 = db)
    dens[p_grid < lower | p_grid > upper] <- 0
    
    if (!all(is.finite(dens)) || sum(dens) <= 0) {
      return(NA_real_)
    }
    
    w <- dens / sum(dens)
    vals <- vapply(p_grid, fun, numeric(1))
    sum(vals * w)
  }
  
  integrate_region_oc <- function(lower, upper, da, db) {
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      return(list(
        pnaive = NA_real_,
        pintfut = NA_real_,
        erased = NA_real_,
        pfineff = NA_real_,
        nexp = NA_real_
      ))
    }
    
    mass <- pbeta(upper, shape1 = da, shape2 = db) -
      pbeta(lower, shape1 = da, shape2 = db)
    
    if (!is.finite(mass) || mass <= 0) {
      return(list(
        pnaive = NA_real_,
        pintfut = NA_real_,
        erased = NA_real_,
        pfineff = NA_real_,
        nexp = NA_real_
      ))
    }
    
    p_grid <- seq(lower, upper, length.out = grid_size)
    dens <- dbeta(p_grid, shape1 = da, shape2 = db)
    dens[p_grid < lower | p_grid > upper] <- 0
    w <- dens / sum(dens)
    
    vals <- lapply(p_grid, calc_fixed_at_p)
    
    avg_component <- function(name) {
      v <- vapply(vals, function(z) z[[name]], numeric(1))
      sum(v * w)
    }
    
    list(
      pnaive = avg_component("pnaive"),
      pintfut = avg_component("pintfut"),
      erased = avg_component("erased"),
      pfineff = avg_component("pfineff"),
      nexp = avg_component("nexp")
    )
  }
  
  ## Bayesian H1
  if (type == "direction") {
    lower_h1 <- p0
    upper_h1 <- 1
    
    if (lower_h1 >= upper_h1) {
      bayes_H1 <- list(
        pnaive = NA_real_,
        pintfut = NA_real_,
        erased = NA_real_,
        pfineff = NA_real_,
        nexp = NA_real_
      )
    } else {
      lower_h1_open <- min(upper_h1 - 1e-10, lower_h1 + 1e-10)
      bayes_H1 <- integrate_region_oc(
        lower = lower_h1_open, upper = upper_h1,
        da = da1, db = db1
      )
    }
  } else {
    pmf_fin_H1 <- singlearm_priorpred_pmf(
      x = x_vals, n = n,
      dp = NA_real_,
      da = da1, db = db1,
      dl = 0, du = 1
    )
    
    pfin_H1 <- sum(pmf_fin_H1[eff_fin])
    
    bayes_H1 <- list(
      pnaive = pfin_H1,
      pintfut = 0,
      erased = 0,
      pfineff = pfin_H1,
      nexp = n
    )
  }
  
  ## Bayesian H0
  if (type == "direction") {
    lower_h0 <- 0
    upper_h0 <- p0
    
    if (lower_h0 >= upper_h0) {
      bayes_H0 <- list(
        pnaive = NA_real_,
        pintfut = NA_real_,
        erased = NA_real_,
        pfineff = NA_real_,
        nexp = NA_real_
      )
    } else {
      bayes_H0 <- integrate_region_oc(
        lower = lower_h0, upper = upper_h0,
        da = da0, db = db0
      )
    }
  } else {
    bayes_H0 <- calc_fixed_at_p(p0)
  }
  
  ## Frequentist H1 at dp
  freq_H1 <- NULL
  if (!is.na(dp)) {
    freq_H1 <- calc_fixed_at_p(dp)
  }
  
  ## Frequentist H0 at p0
  freq_H0 <- calc_fixed_at_p(p0)
  
  ## Compelling evidence for H0
  pce0_naive <- NA_real_
  pce0_corr  <- NA_real_
  
  if (!is.null(k_ce)) {
    ce_fin <- bf_fin >= k_ce
    
    if (type == "direction") {
      pce0_naive <- integrate_region(
        lower = 0,
        upper = p0,
        da = da0,
        db = db0,
        fun = function(p) {
          pmf_fin <- dbinom(x_vals, size = n, prob = p)
          sum(pmf_fin[ce_fin])
        }
      )
      pce0_corr <- pce0_naive
    } else {
      pmf_fin_H0 <- dbinom(x_vals, size = n, prob = p0)
      pce0_naive <- sum(pmf_fin_H0[ce_fin])
      pce0_corr <- pce0_naive
    }
  }
  
  list(
    n = n,
    k = k,
    p0 = p0,
    dp = dp,
    type = type,
    bf_final = bf_fin,
    
    ## Bayesian H1
    pintfut = bayes_H1$pintfut,
    pfineff = bayes_H1$pfineff,
    pnaive = bayes_H1$pnaive,
    perased = bayes_H1$erased,
    nexp = bayes_H1$nexp,
    
    ## Bayesian H0
    pintfut0 = bayes_H0$pintfut,
    pfineff0 = bayes_H0$pfineff,
    pnaive0 = bayes_H0$pnaive,
    perased0 = bayes_H0$erased,
    nexp0 = bayes_H0$nexp,
    
    ## Frequentist H1 at dp
    pintfut_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pintfut,
    pfineff_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pfineff,
    pnaive_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pnaive,
    perased_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$erased,
    nexp_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$nexp,
    
    ## Frequentist H0 at p0
    pintfut_freq0 = freq_H0$pintfut,
    pfineff_freq0 = freq_H0$pfineff,
    pnaive_freq0 = freq_H0$pnaive,
    perased_freq0 = freq_H0$erased,
    nexp_freq0 = freq_H0$nexp,
    
    freq0_argmax = p0,
    k_ce = k_ce,
    pce0_naive = pce0_naive,
    pce0_corr = pce0_corr
  )
}