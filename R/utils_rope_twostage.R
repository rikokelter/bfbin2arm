# =============================================================================
# Utility helpers for single-arm two-stage ROPE designs
#
# All computations reuse the shared helpers already defined in utils_rope.R:
# beta_binom_pmf_rope(y, n, a, b)
# posterior_rope_prob(y, n, p0, delta, analysis_prior)
# .validate_beta_prior(), .validate_probability(), .validate_count()
# =============================================================================

# -----------------------------------------------------------------------------
# Continuation region for stage 1
# -----------------------------------------------------------------------------
.continuation_region_twostage <- function(n1, p0, delta, gamma_1,
                                          analysis_prior) {
  y1_vals <- 0:n1
  post <- vapply(
    y1_vals,
    posterior_rope_prob,
    numeric(1L),
    n = n1,
    p0 = p0,
    delta = delta,
    analysis_prior = analysis_prior
  )
  
  ## Stop for futility when Pr(H0 supported | y1, n1) >= gamma_1,
  ## i.e. post <= 1 - gamma_1. Continue otherwise.
  y1_vals[post > (1 - gamma_1)]
}

# -----------------------------------------------------------------------------
# Predictive probability of declaring equivalence — two-stage design
# -----------------------------------------------------------------------------
.predictive_equiv_twostage <- function(
    n1, n2, p0, delta, gamma_1, gamma_eq, analysis_prior, design_prior
) {
  cont <- .continuation_region_twostage(n1, p0, delta, gamma_1, analysis_prior)
  if (length(cont) == 0L) return(0)
  
  n <- n1 + n2
  aD <- design_prior[1]
  bD <- design_prior[2]
  
  out <- 0
  for (y1 in cont) {
    p_y1 <- beta_binom_pmf_rope(y1, n1, aD, bD)
    dp_post <- c(aD + y1, bD + n1 - y1)
    
    for (y2 in 0:n2) {
      y <- y1 + y2
      p_y2_y1 <- beta_binom_pmf_rope(y2, n2, dp_post[1], dp_post[2])
      post_h1 <- posterior_rope_prob(y, n, p0, delta, analysis_prior)
      
      if (post_h1 >= gamma_eq) {
        out <- out + p_y1 * p_y2_y1
      }
    }
  }
  
  out
}

# -----------------------------------------------------------------------------
# Predictive probability of compelling evidence against equivalence — two-stage
# -----------------------------------------------------------------------------
.predictive_pce_twostage <- function(
    n1, n2, p0, delta, gamma_1, gamma_diff, analysis_prior, design_prior
) {
  n <- n1 + n2
  aD <- design_prior[1]
  bD <- design_prior[2]
  
  y1_vals <- 0:n1
  post1 <- vapply(
    y1_vals,
    posterior_rope_prob,
    numeric(1L),
    n = n1,
    p0 = p0,
    delta = delta,
    analysis_prior = analysis_prior
  )
  
  out <- 0
  
  for (i in seq_along(y1_vals)) {
    y1 <- y1_vals[i]
    p_y1 <- beta_binom_pmf_rope(y1, n1, aD, bD)
    post1_y1 <- post1[i]
    
    ## Stage 1: compelling evidence for H0
    if ((1 - post1_y1) >= gamma_diff) {
      out <- out + p_y1
      next
    }
    
    ## Stage 1: stop for futility without compelling H0 evidence
    if (post1_y1 <= (1 - gamma_1)) {
      next
    }
    
    ## Stage 2: continue, then check final compelling evidence for H0
    dp_post <- c(aD + y1, bD + n1 - y1)
    for (y2 in 0:n2) {
      y <- y1 + y2
      p_y2_y1 <- beta_binom_pmf_rope(y2, n2, dp_post[1], dp_post[2])
      post_h1 <- posterior_rope_prob(y, n, p0, delta, analysis_prior)
      
      if ((1 - post_h1) >= gamma_diff) {
        out <- out + p_y1 * p_y2_y1
      }
    }
  }
  
  out
}

# -----------------------------------------------------------------------------
# Expected sample size under a design prior — two-stage design
# -----------------------------------------------------------------------------
.expected_n_twostage <- function(n1, n2, p0, delta, gamma_1,
                                 analysis_prior, design_prior) {
  cont <- .continuation_region_twostage(n1, p0, delta, gamma_1, analysis_prior)
  aD <- design_prior[1]
  bD <- design_prior[2]
  
  p_cont <- if (length(cont) == 0L) {
    0
  } else {
    sum(vapply(cont, beta_binom_pmf_rope, numeric(1L), n = n1, a = aD, b = bD))
  }
  
  n1 + n2 * p_cont
}

# -----------------------------------------------------------------------------
# One-stage predictive probability of declaring equivalence
# -----------------------------------------------------------------------------
.predictive_equiv_onestage <- function(n, p0, delta, gamma_eq,
                                       analysis_prior, design_prior) {
  y_vals <- 0:n
  post_h1 <- vapply(
    y_vals, posterior_rope_prob, numeric(1L),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  
  aD <- design_prior[1]
  bD <- design_prior[2]
  
  sum(vapply(
    y_vals[post_h1 >= gamma_eq],
    beta_binom_pmf_rope, numeric(1L),
    n = n, a = aD, b = bD
  ))
}

# -----------------------------------------------------------------------------
# One-stage predictive probability of compelling evidence against equivalence
# -----------------------------------------------------------------------------
.predictive_pce_onestage <- function(n, p0, delta, gamma_diff,
                                     analysis_prior, design_prior) {
  y_vals <- 0:n
  post_h1 <- vapply(
    y_vals, posterior_rope_prob, numeric(1L),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  
  aD <- design_prior[1]
  bD <- design_prior[2]
  
  sum(vapply(
    y_vals[(1 - post_h1) >= gamma_diff],
    beta_binom_pmf_rope, numeric(1L),
    n = n, a = aD, b = bD
  ))
}

# -----------------------------------------------------------------------------
# Frequentist operating characteristics under fixed true p
# -----------------------------------------------------------------------------

# One-stage frequentist probability of declaring H1
.frequentist_equiv_onestage <- function(n, p_true, p0, delta, gamma_eq,
                                        analysis_prior) {
  y_vals <- 0:n
  post_h1 <- vapply(
    y_vals, posterior_rope_prob, numeric(1L),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  sum(dbinom(y_vals[post_h1 >= gamma_eq], size = n, prob = p_true))
}

# One-stage frequentist probability of compelling evidence for H0
.frequentist_pce_onestage <- function(n, p_true, p0, delta, gamma_diff,
                                      analysis_prior) {
  y_vals <- 0:n
  post_h1 <- vapply(
    y_vals, posterior_rope_prob, numeric(1L),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  sum(dbinom(y_vals[(1 - post_h1) >= gamma_diff], size = n, prob = p_true))
}

# Two-stage frequentist probability of declaring H1
.frequentist_equiv_twostage <- function(n1, n2, p_true, p0, delta, gamma_1,
                                        gamma_eq, analysis_prior) {
  cont <- .continuation_region_twostage(n1, p0, delta, gamma_1, analysis_prior)
  if (length(cont) == 0L) return(0)
  
  n <- n1 + n2
  out <- 0
  
  for (y1 in cont) {
    p_y1 <- dbinom(y1, size = n1, prob = p_true)
    for (y2 in 0:n2) {
      y <- y1 + y2
      p_y2_y1 <- dbinom(y2, size = n2, prob = p_true)
      post_h1 <- posterior_rope_prob(y, n, p0, delta, analysis_prior)
      
      if (post_h1 >= gamma_eq) {
        out <- out + p_y1 * p_y2_y1
      }
    }
  }
  
  out
}

# Two-stage frequentist probability of compelling evidence for H0
.frequentist_pce_twostage <- function(n1, n2, p_true, p0, delta, gamma_1,
                                      gamma_diff, analysis_prior) {
  n <- n1 + n2
  y1_vals <- 0:n1
  post1 <- vapply(
    y1_vals,
    posterior_rope_prob,
    numeric(1L),
    n = n1,
    p0 = p0,
    delta = delta,
    analysis_prior = analysis_prior
  )
  
  out <- 0
  
  for (i in seq_along(y1_vals)) {
    y1 <- y1_vals[i]
    p_y1 <- dbinom(y1, size = n1, prob = p_true)
    post1_y1 <- post1[i]
    
    if ((1 - post1_y1) >= gamma_diff) {
      out <- out + p_y1
      next
    }
    
    if (post1_y1 <= (1 - gamma_1)) {
      next
    }
    
    for (y2 in 0:n2) {
      y <- y1 + y2
      p_y2_y1 <- dbinom(y2, size = n2, prob = p_true)
      post_h1 <- posterior_rope_prob(y, n, p0, delta, analysis_prior)
      
      if ((1 - post_h1) >= gamma_diff) {
        out <- out + p_y1 * p_y2_y1
      }
    }
  }
  
  out
}

# Two-stage frequentist expected sample size under fixed true p
.frequentist_en_twostage <- function(n1, n2, p_true, p0, delta, gamma_1,
                                     analysis_prior) {
  cont <- .continuation_region_twostage(n1, p0, delta, gamma_1, analysis_prior)
  p_cont <- if (length(cont) == 0L) 0 else sum(dbinom(cont, size = n1, prob = p_true))
  n1 + n2 * p_cont
}

# -----------------------------------------------------------------------------
# Evaluate all operating characteristics for one candidate two-stage design
# -----------------------------------------------------------------------------
.evaluate_rope_twostage_oc <- function(
    n1, n2, p0, delta, gamma_1, gamma_eq, gamma_diff, analysis_prior,
    design_prior_h0, design_prior_h1, p_t1e = NULL, p_power = NULL
) {
  n <- n1 + n2
  
  out <- list(
    n1 = n1,
    n2 = n2,
    n = n,
    type1_1st = .predictive_equiv_onestage(
      n, p0, delta, gamma_eq, analysis_prior, design_prior_h0
    ),
    power_1st = .predictive_equiv_onestage(
      n, p0, delta, gamma_eq, analysis_prior, design_prior_h1
    ),
    pce_1st = .predictive_pce_onestage(
      n, p0, delta, gamma_diff, analysis_prior, design_prior_h0
    ),
    type1_2st = .predictive_equiv_twostage(
      n1, n2, p0, delta, gamma_1, gamma_eq, analysis_prior, design_prior_h0
    ),
    power_2st = .predictive_equiv_twostage(
      n1, n2, p0, delta, gamma_1, gamma_eq, analysis_prior, design_prior_h1
    ),
    pce_2st = .predictive_pce_twostage(
      n1, n2, p0, delta, gamma_1, gamma_diff, analysis_prior, design_prior_h0
    ),
    EN0 = .expected_n_twostage(
      n1, n2, p0, delta, gamma_1, analysis_prior, design_prior_h0
    ),
    EN1 = .expected_n_twostage(
      n1, n2, p0, delta, gamma_1, analysis_prior, design_prior_h1
    ),
    cont_region = .continuation_region_twostage(
      n1, p0, delta, gamma_1, analysis_prior
    )
  )
  
  ## Frequentist OC at p_t1e (type-I error, PCE, EN)
  if (!is.null(p_t1e)) {
    out$freq_type1_1st <- .frequentist_equiv_onestage(
      n, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_eq = gamma_eq, analysis_prior = analysis_prior
    )
    out$freq_type1_2st <- .frequentist_equiv_twostage(
      n1, n2, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_1 = gamma_1, gamma_eq = gamma_eq, analysis_prior = analysis_prior
    )
    out$freq_pce_1st <- .frequentist_pce_onestage(
      n, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_diff = gamma_diff, analysis_prior = analysis_prior
    )
    out$freq_pce_2st <- .frequentist_pce_twostage(
      n1, n2, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_1 = gamma_1, gamma_diff = gamma_diff, analysis_prior = analysis_prior
    )
    out$freq_EN_t1e <- .frequentist_en_twostage(
      n1, n2, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_1 = gamma_1, analysis_prior = analysis_prior
    )
  }
  
  ## Frequentist OC at p_power (power, EN)
  if (!is.null(p_power)) {
    out$freq_power_1st <- .frequentist_equiv_onestage(
      n, p_true = p_power, p0 = p0, delta = delta,
      gamma_eq = gamma_eq, analysis_prior = analysis_prior
    )
    out$freq_power_2st <- .frequentist_equiv_twostage(
      n1, n2, p_true = p_power, p0 = p0, delta = delta,
      gamma_1 = gamma_1, gamma_eq = gamma_eq, analysis_prior = analysis_prior
    )
    out$freq_EN_power <- .frequentist_en_twostage(
      n1, n2, p_true = p_power, p0 = p0, delta = delta,
      gamma_1 = gamma_1, analysis_prior = analysis_prior
    )
  }
  
  out
}