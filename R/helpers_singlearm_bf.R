# ============================================================
# Internal helpers for single-arm Bayes factor calculations
# ============================================================

# Canonical single-arm BF01
# type = "point": H0: p = p0 vs H1: p ~ Beta(a1,b1)
# type = "direction": H0: p <= p0 vs H1: p > p0
# with separate Beta analysis priors under H0 and H1,
# truncated accordingly.
singlearm_bf01 <- function(
    x, n, p0,
    a0 = 1, b0 = 1,
    a1 = 1, b1 = 1,
    type = c("point", "direction")
) {
  type <- match.arg(type)
  
  if (length(p0) != 1 || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
    stop("p0 must be a single number in (0,1).")
  }
  if (any(x < 0 | x > n | x != floor(x))) {
    stop("x must contain integers between 0 and n.")
  }
  if (length(n) != 1 || n < 1 || n != floor(n)) {
    stop("n must be a positive integer.")
  }
  if (a0 <= 0 || b0 <= 0 || a1 <= 0 || b1 <= 0) {
    stop("a0, b0, a1, and b1 must be positive.")
  }
  
  x <- as.integer(x)
  
  if (type == "point") {
    # H0: p = p0
    m0 <- dbinom(x = x, size = n, prob = p0)
    # H1: p ~ Beta(a1, b1)
    m1 <- exp(
      lchoose(n, x) + lbeta(a1 + x, b1 + n - x) - lbeta(a1, b1)
    )
    return(m0 / m1)
  }
  
  # type == "direction":
  # H0: p <= p0, Beta(a0,b0) truncated to [0, p0]
  # H1: p >  p0, Beta(a1,b1) truncated to (p0, 1]
  m0 <- dbetabinom_truncbeta(
    x = x, n = n,
    a = a0, b = b0,
    lower = 0, upper = p0
  )
  m1 <- dbetabinom_truncbeta(
    x = x, n = n,
    a = a1, b = b1,
    lower = p0, upper = 1
  )
  m0 / m1
}

# Exact probability for BF event by summing over counts
# lower.tail = TRUE  -> P(BF01 <= k)
# lower.tail = FALSE -> P(BF01 >  k)
pbinbf01_singlearm_exact <- function(
    k, n, p0,
    a0 = 1, b0 = 1,
    a1 = 1, b1 = 1,
    dp = NA_real_,
    da0 = 1, db0 = 1,
    da1 = 1, db1 = 1,
    type = c("point", "direction"),
    lower.tail = TRUE
) {
  type <- match.arg(type)
  
  if (length(k) != 1 || !is.finite(k) || k <= 0) {
    stop("k must be a positive number.")
  }
  
  xvals <- 0:n
  bfvals <- singlearm_bf01(
    x = xvals, n = n, p0 = p0,
    a0 = a0, b0 = b0,
    a1 = a1, b1 = b1,
    type = type
  )
  
  if (type == "point") {
    if (is.na(dp)) {
      pmf <- singlearm_priorpred_pmf(
        x = xvals, n = n,
        dp = NA_real_,
        da = da1, db = db1,
        dl = 0, du = 1
      )
    } else {
      pmf <- dbinom(xvals, size = n, prob = dp)
    }
  } else {
    if (is.na(dp)) {
      pmf <- singlearm_priorpred_pmf(
        x = xvals, n = n,
        dp = NA_real_,
        da = da1, db = db1,
        dl = p0, du = 1
      )
    } else {
      pmf <- dbinom(xvals, size = n, prob = dp)
    }
  }
  
  if (lower.tail) {
    sum(pmf[bfvals <= k])
  } else {
    sum(pmf[bfvals > k])
  }
}

# Fixed-sample operating characteristics for a single-arm BF design
#'
#' @keywords internal
singlearm_fixed_oc <- function(
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
    stop("'k' must be a positive number.")
  }
  
  if (!is.numeric(p0) || length(p0) != 1L || !is.finite(p0) ||
      p0 <= 0 || p0 >= 1) {
    stop("'p0' must be in (0, 1).")
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
    stop("'dp' must be NA or in (0, 1).")
  }
  
  if (!is.null(k_ce)) {
    if (!is.numeric(k_ce) || length(k_ce) != 1L || !is.finite(k_ce) || k_ce <= 1) {
      stop("'k_ce' must be NULL or a numeric value > 1.")
    }
  }
  
  if (!is.numeric(grid_size) || length(grid_size) != 1L ||
      !is.finite(grid_size) || grid_size < 101L) {
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
  
  integrate_region <- function(lower, upper, da, db) {
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
      bayes_H1 <- integrate_region(
        lower = lower_h1_open,
        upper = upper_h1,
        da = da1,
        db = db1
      )
    }
  } else {
    pmf_fin_H1 <- singlearm_priorpred_pmf(
      x = x_vals,
      n = n,
      dp = NA_real_,
      da = da1,
      db = db1,
      dl = 0,
      du = 1
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
      bayes_H0 <- integrate_region(
        lower = lower_h0,
        upper = upper_h0,
        da = da0,
        db = db0
      )
    }
  } else {
    bayes_H0 <- calc_fixed_at_p(p0)
  }
  
  ## Frequentist H1
  freq_H1 <- NULL
  if (!is.na(dp)) {
    freq_H1 <- calc_fixed_at_p(dp)
  }
  
  ## Frequentist H0
  freq_H0 <- calc_fixed_at_p(p0)
  
  ## Compelling evidence in favour of H0 at p0
  pce0_naive <- NA_real_
  pce0_corr <- NA_real_
  if (!is.null(k_ce)) {
    ce_fin <- bf_fin >= k_ce
    pmf_fin_H0 <- dbinom(x_vals, size = n, prob = p0)
    pce0_naive <- sum(pmf_fin_H0[ce_fin])
    pce0_corr <- pce0_naive
  }
  
  list(
    n = n,
    k = k,
    p0 = p0,
    dp = dp,
    type = type,
    bf_final = bf_fin,
    
    ## Bayesian H1
    pnaive = bayes_H1$pnaive,
    pintfut = bayes_H1$pintfut,
    perased = bayes_H1$erased,
    pfineff = bayes_H1$pfineff,
    nexp = bayes_H1$nexp,
    
    ## Bayesian H0
    pnaive0 = bayes_H0$pnaive,
    pintfut0 = bayes_H0$pintfut,
    perased0 = bayes_H0$erased,
    pfineff0 = bayes_H0$pfineff,
    nexp0 = bayes_H0$nexp,
    
    ## Frequentist H1 at dp
    pnaive_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pnaive,
    pintfut_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pintfut,
    perased_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$erased,
    pfineff_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$pfineff,
    nexp_freq = if (is.null(freq_H1)) NA_real_ else freq_H1$nexp,
    
    ## Frequentist H0 at p0
    pnaive_freq0 = freq_H0$pnaive,
    pintfut_freq0 = freq_H0$pintfut,
    perased_freq0 = freq_H0$erased,
    pfineff_freq0 = freq_H0$pfineff,
    nexp_freq0 = freq_H0$nexp,
    
    ## Compelling evidence for H0
    pce0_naive = pce0_naive,
    pce0_corr = pce0_corr,
    freq0_argmax = p0
  )
}

#' Check sustained feasibility over future n
#'
#' Given vectors of operating characteristics over n, check whether
#' power, type-I-error, and CE(H0) satisfy their thresholds at n and
#' for at least sustain_n subsequent sample sizes.
#'
#' @param n_vec Integer vector of sample sizes.
#' @param power_vec Numeric vector of power values (same length as n_vec).
#' @param type1_vec Numeric vector of type-I-error values.
#' @param ce_vec Numeric vector of CE(H0) values (may contain NA).
#' @param target_power Numeric target power.
#' @param target_type1 Numeric target type-I-error.
#' @param target_ce Numeric target CE(H0); if <= 0, CE constraint is ignored.
#' @param sustain_n Integer, number of subsequent sample sizes that must
#'   also satisfy the constraints.
#'
#' @return Logical vector of length length(n_vec): TRUE if n_i and the next
#'   sustain_n sample sizes satisfy all active constraints.
#'
.sustained_singlearm_feasibility <- function(
    n_vec,
    power_vec,
    type1_vec,
    ce_vec,
    target_power,
    target_type1,
    target_ce,
    sustain_n
) {
  stopifnot(length(n_vec) == length(power_vec),
            length(n_vec) == length(type1_vec),
            length(n_vec) == length(ce_vec))
  
  len <- length(n_vec)
  ok <- rep(FALSE, len)
  
  for (i in seq_len(len)) {
    # indices to check: i ... i + sustain_n, clipped at len
    j_max <- min(len, i + sustain_n)
    idx <- i:j_max
    
    # basic Bayes constraints
    cond_power <- power_vec[idx] >= target_power
    cond_type1 <- type1_vec[idx] <= target_type1
    
    cond_ce <- rep(TRUE, length(idx))
    if (isTRUE(target_ce > 0)) {
      # CE(H0) only active if target_ce > 0; NA values fail the condition
      cond_ce <- !is.na(ce_vec[idx]) & ce_vec[idx] >= target_ce
    }
    
    ok[i] <- all(cond_power & cond_type1 & cond_ce)
  }
  
  ok
}