# ============================================================
# Internal helpers for single-arm predictive densities
# ============================================================

# Normalizing mass of Beta(a,b) on [lower, upper]
truncbeta_mass <- function(a, b, lower = 0, upper = 1) {
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("a and b must be positive.")
  }
  if (!is.finite(lower) || !is.finite(upper) || lower < 0 || upper > 1 || lower >= upper) {
    stop("Need 0 <= lower < upper <= 1.")
  }
  pbeta(upper, shape1 = a, shape2 = b) - pbeta(lower, shape1 = a, shape2 = b)
}

# Exact prior-predictive pmf under truncated Beta prior:
# X | p ~ Bin(n,p), p ~ Beta(a,b) truncated to [lower, upper]
dbetabinom_truncbeta <- function(x, n, a, b, lower = 0, upper = 1, log = FALSE) {
  if (any(x < 0 | x > n | x != floor(x))) {
    out <- rep(if (log) -Inf else 0, length(x))
    return(out)
  }
  
  z <- truncbeta_mass(a = a, b = b, lower = lower, upper = upper)
  
  if (z <= 0) {
    stop("Truncated Beta prior has zero mass on [lower, upper].")
  }
  
  # Integral:
  # choose(n,x) / z / B(a,b) * \int_lower^upper p^(x+a-1) (1-p)^(n-x+b-1) dp
  # = choose(n,x) * B(a+x, b+n-x) / B(a,b) *
  #   [I_upper(a+x,b+n-x)-I_lower(a+x,b+n-x)] / z
  log_choose <- lchoose(n, x)
  log_beta_ratio <- lbeta(a + x, b + n - x) - lbeta(a, b)
  
  inc_part <- pbeta(upper, shape1 = a + x, shape2 = b + n - x) -
    pbeta(lower, shape1 = a + x, shape2 = b + n - x)
  
  dens <- exp(log_choose + log_beta_ratio) * inc_part / z
  
  if (log) {
    out <- log(dens)
    out[dens <= 0] <- -Inf
    return(out)
  }
  dens
}

# Unified prior-predictive pmf:
# either point prior p = dp, or truncated Beta design prior
singlearm_priorpred_pmf <- function(x, n, dp = NA_real_,
                                    da = 1, db = 1,
                                    dl = 0, du = 1,
                                    log = FALSE) {
  if (!is.na(dp)) {
    out <- dbinom(x = x, size = n, prob = dp, log = log)
    return(out)
  }
  
  dbetabinom_truncbeta(
    x = x, n = n,
    a = da, b = db,
    lower = dl, upper = du,
    log = log
  )
}



# joint prior-predictive integral for (X,Z) in closed form for the single-arm
# truncated-Beta design prior
dbinbin_truncbeta <- function(x, z, n1, n2, a, b, lower = 0, upper = 1, log = FALSE) {
  zmass <- truncbeta_mass(a, b, lower, upper)
  if (zmass <= 0) stop("Truncated Beta prior has zero mass on [lower, upper].")
  
  logp <- lchoose(n1, x) +
    lchoose(n2 - n1, z) +
    lbeta(a + x + z, b + n2 - x - z) -
    lbeta(a, b)
  
  inc_part <- pbeta(upper, a + x + z, b + n2 - x - z) -
    pbeta(lower, a + x + z, b + n2 - x - z)
  
  dens <- exp(logp) * inc_part / zmass
  if (log) {
    out <- log(dens)
    out[dens <= 0] <- -Inf
    return(out)
  }
  dens
}