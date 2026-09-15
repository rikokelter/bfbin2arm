# =============================================================================
# Utility helpers for single-arm two-stage ROPE designs
#
# Prerequisite in utils_rope.R:
# posterior_rope_prob(y, n, p0, delta, analysis_prior,
#                     direction = c("equivalence", "noninferiority", "superiority"))
# must return posterior support for H1 in the selected direction.
# =============================================================================


# -----------------------------------------------------------------------------
# Direction-specific labels for print and summary methods
# -----------------------------------------------------------------------------
.rope_twostage_direction_labels <- function(
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  
  switch(
    direction,
    
    equivalence = list(
      title = "equivalence",
      h1 = "equivalence",
      h0 = "non-equivalence",
      h0_short = "non-equivalence",
      interim = "non-equivalence",
      success = "Declare equivalence",
      null_evidence = "Compelling evidence for non-equivalence"
    ),
    
    noninferiority = list(
      title = "non-inferiority",
      h1 = "non-inferiority",
      h0 = "clinically relevant inferiority",
      h0_short = "inferiority",
      interim = "clinically relevant inferiority",
      success = "Declare non-inferiority",
      null_evidence = "Compelling evidence for inferiority"
    ),
    
    superiority = list(
      title = "superiority",
      h1 = "superiority",
      h0 = "non-superiority",
      h0_short = "non-superiority",
      interim = "non-superiority",
      success = "Declare superiority",
      null_evidence = "Compelling evidence for non-superiority"
    )
  )
}

.posterior_rope_vec <- function(y, n, p0, delta, analysis_prior, direction) {
  direction <- match.arg(
    direction,
    c("equivalence", "noninferiority", "superiority")
  )
  
  vapply(
    y,
    posterior_rope_prob,
    numeric(1L),
    n = n,
    p0 = p0,
    delta = delta,
    analysis_prior = analysis_prior,
    direction = direction
  )
}

# -----------------------------------------------------------------------------
# Continuation region for stage 1
# -----------------------------------------------------------------------------
.continuation_region_twostage <- function(
    n1, p0, delta, gamma_1, analysis_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  y1_vals <- 0:n1
  
  post_h1 <- .posterior_rope_vec(
    y1_vals, n1, p0, delta, analysis_prior, direction
  )
  
  ## Stop for futility when posterior support for H0 is at least gamma_1.
  ## Continue when posterior support for H1 is strictly greater than 1-gamma_1.
  y1_vals[post_h1 > (1 - gamma_1)]
}

# -----------------------------------------------------------------------------
# Predictive probability of declaring H1 -- two-stage design
# -----------------------------------------------------------------------------
.predictive_equiv_twostage <- function(
    n1, n2, p0, delta, gamma_1, gamma_eq, analysis_prior, design_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  
  cont <- .continuation_region_twostage(
    n1, p0, delta, gamma_1, analysis_prior, direction
  )
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
      post_h1 <- posterior_rope_prob(
        y, n, p0, delta, analysis_prior, direction = direction
      )
      
      if (post_h1 >= gamma_eq) {
        out <- out + p_y1 * p_y2_y1
      }
    }
  }
  
  out
}

# -----------------------------------------------------------------------------
# Predictive probability of compelling evidence for H0 -- two-stage design
# -----------------------------------------------------------------------------
.predictive_pce_twostage <- function(
    n1, n2, p0, delta, gamma_1, gamma_diff, analysis_prior, design_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  n <- n1 + n2
  aD <- design_prior[1]
  bD <- design_prior[2]
  y1_vals <- 0:n1
  
  post1 <- .posterior_rope_vec(
    y1_vals, n1, p0, delta, analysis_prior, direction
  )
  
  out <- 0
  
  for (i in seq_along(y1_vals)) {
    y1 <- y1_vals[i]
    p_y1 <- beta_binom_pmf_rope(y1, n1, aD, bD)
    post_h1_interim <- post1[i]
    post_h0_interim <- 1 - post_h1_interim
    
    ## Interim futility stop with compelling evidence for H0.
    if (post_h0_interim >= max(gamma_1, gamma_diff)) {
      out <- out + p_y1
      next
    }
    
    ## Interim futility stop without compelling evidence for H0.
    if (post_h0_interim >= gamma_1) {
      next
    }
    
    ## Continuation to the final analysis.
    dp_post <- c(aD + y1, bD + n1 - y1)
    
    for (y2 in 0:n2) {
      y <- y1 + y2
      p_y2_y1 <- beta_binom_pmf_rope(y2, n2, dp_post[1], dp_post[2])
      post_h1_final <- posterior_rope_prob(
        y, n, p0, delta, analysis_prior, direction = direction
      )
      
      if ((1 - post_h1_final) >= gamma_diff) {
        out <- out + p_y1 * p_y2_y1
      }
    }
  }
  
  out
}

# -----------------------------------------------------------------------------
# Expected sample size under a design prior -- two-stage design
# -----------------------------------------------------------------------------
.expected_n_twostage <- function(
    n1, n2, p0, delta, gamma_1, analysis_prior, design_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  
  cont <- .continuation_region_twostage(
    n1, p0, delta, gamma_1, analysis_prior, direction
  )
  aD <- design_prior[1]
  bD <- design_prior[2]
  
  p_cont <- if (length(cont) == 0L) {
    0
  } else {
    sum(vapply(
      cont, beta_binom_pmf_rope, numeric(1L),
      n = n1, a = aD, b = bD
    ))
  }
  
  n1 + n2 * p_cont
}

# -----------------------------------------------------------------------------
# One-stage predictive probability of declaring H1
# -----------------------------------------------------------------------------
.predictive_equiv_onestage <- function(
    n, p0, delta, gamma_eq, analysis_prior, design_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  y_vals <- 0:n
  
  post_h1 <- .posterior_rope_vec(
    y_vals, n, p0, delta, analysis_prior, direction
  )
  aD <- design_prior[1]
  bD <- design_prior[2]
  
  sum(vapply(
    y_vals[post_h1 >= gamma_eq],
    beta_binom_pmf_rope,
    numeric(1L),
    n = n,
    a = aD,
    b = bD
  ))
}

# -----------------------------------------------------------------------------
# One-stage predictive probability of compelling evidence for H0
# -----------------------------------------------------------------------------
.predictive_pce_onestage <- function(
    n, p0, delta, gamma_diff, analysis_prior, design_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  y_vals <- 0:n
  
  post_h1 <- .posterior_rope_vec(
    y_vals, n, p0, delta, analysis_prior, direction
  )
  aD <- design_prior[1]
  bD <- design_prior[2]
  
  sum(vapply(
    y_vals[(1 - post_h1) >= gamma_diff],
    beta_binom_pmf_rope,
    numeric(1L),
    n = n,
    a = aD,
    b = bD
  ))
}

# -----------------------------------------------------------------------------
# Frequentist operating characteristics under fixed true p
# -----------------------------------------------------------------------------

.frequentist_equiv_onestage <- function(
    n, p_true, p0, delta, gamma_eq, analysis_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  y_vals <- 0:n
  
  post_h1 <- .posterior_rope_vec(
    y_vals, n, p0, delta, analysis_prior, direction
  )
  
  sum(dbinom(y_vals[post_h1 >= gamma_eq], size = n, prob = p_true))
}

.frequentist_pce_onestage <- function(
    n, p_true, p0, delta, gamma_diff, analysis_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  y_vals <- 0:n
  
  post_h1 <- .posterior_rope_vec(
    y_vals, n, p0, delta, analysis_prior, direction
  )
  
  sum(dbinom(
    y_vals[(1 - post_h1) >= gamma_diff],
    size = n,
    prob = p_true
  ))
}

.frequentist_equiv_twostage <- function(
    n1, n2, p_true, p0, delta, gamma_1, gamma_eq, analysis_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  
  cont <- .continuation_region_twostage(
    n1, p0, delta, gamma_1, analysis_prior, direction
  )
  if (length(cont) == 0L) return(0)
  
  n <- n1 + n2
  out <- 0
  
  for (y1 in cont) {
    p_y1 <- dbinom(y1, size = n1, prob = p_true)
    
    for (y2 in 0:n2) {
      y <- y1 + y2
      p_y2_y1 <- dbinom(y2, size = n2, prob = p_true)
      post_h1 <- posterior_rope_prob(
        y, n, p0, delta, analysis_prior, direction = direction
      )
      
      if (post_h1 >= gamma_eq) {
        out <- out + p_y1 * p_y2_y1
      }
    }
  }
  
  out
}

.frequentist_pce_twostage <- function(
    n1, n2, p_true, p0, delta, gamma_1, gamma_diff, analysis_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  n <- n1 + n2
  y1_vals <- 0:n1
  
  post1 <- .posterior_rope_vec(
    y1_vals, n1, p0, delta, analysis_prior, direction
  )
  out <- 0
  
  for (i in seq_along(y1_vals)) {
    y1 <- y1_vals[i]
    p_y1 <- dbinom(y1, size = n1, prob = p_true)
    post_h1_interim <- post1[i]
    post_h0_interim <- 1 - post_h1_interim
    
    ## Interim futility stop with compelling evidence for H0.
    if (post_h0_interim >= max(gamma_1, gamma_diff)) {
      out <- out + p_y1
      next
    }
    
    ## Interim futility stop without compelling evidence for H0.
    if (post_h0_interim >= gamma_1) {
      next
    }
    
    ## Continuation to final analysis.
    for (y2 in 0:n2) {
      y <- y1 + y2
      p_y2_y1 <- dbinom(y2, size = n2, prob = p_true)
      post_h1_final <- posterior_rope_prob(
        y, n, p0, delta, analysis_prior, direction = direction
      )
      
      if ((1 - post_h1_final) >= gamma_diff) {
        out <- out + p_y1 * p_y2_y1
      }
    }
  }
  
  out
}

.frequentist_en_twostage <- function(
    n1, n2, p_true, p0, delta, gamma_1, analysis_prior,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  
  cont <- .continuation_region_twostage(
    n1, p0, delta, gamma_1, analysis_prior, direction
  )
  
  p_cont <- if (length(cont) == 0L) {
    0
  } else {
    sum(dbinom(cont, size = n1, prob = p_true))
  }
  
  n1 + n2 * p_cont
}

# -----------------------------------------------------------------------------
# Evaluate all operating characteristics for one candidate two-stage design
# -----------------------------------------------------------------------------
.evaluate_rope_twostage_oc <- function(
    n1, n2, p0, delta, gamma_1, gamma_eq, gamma_diff, analysis_prior,
    design_prior_h0, design_prior_h1, p_t1e = NULL, p_power = NULL,
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  n <- n1 + n2
  
  out <- list(
    n1 = n1,
    n2 = n2,
    n = n,
    type1_1st = .predictive_equiv_onestage(
      n, p0, delta, gamma_eq, analysis_prior, design_prior_h0, direction
    ),
    power_1st = .predictive_equiv_onestage(
      n, p0, delta, gamma_eq, analysis_prior, design_prior_h1, direction
    ),
    pce_1st = .predictive_pce_onestage(
      n, p0, delta, gamma_diff, analysis_prior, design_prior_h0, direction
    ),
    type1_2st = .predictive_equiv_twostage(
      n1, n2, p0, delta, gamma_1, gamma_eq, analysis_prior,
      design_prior_h0, direction
    ),
    power_2st = .predictive_equiv_twostage(
      n1, n2, p0, delta, gamma_1, gamma_eq, analysis_prior,
      design_prior_h1, direction
    ),
    pce_2st = .predictive_pce_twostage(
      n1, n2, p0, delta, gamma_1, gamma_diff, analysis_prior,
      design_prior_h0, direction
    ),
    EN0 = .expected_n_twostage(
      n1, n2, p0, delta, gamma_1, analysis_prior, design_prior_h0, direction
    ),
    EN1 = .expected_n_twostage(
      n1, n2, p0, delta, gamma_1, analysis_prior, design_prior_h1, direction
    ),
    cont_region = .continuation_region_twostage(
      n1, p0, delta, gamma_1, analysis_prior, direction
    )
  )
  
  if (!is.null(p_t1e)) {
    out$freq_type1_1st <- .frequentist_equiv_onestage(
      n, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_eq = gamma_eq, analysis_prior = analysis_prior,
      direction = direction
    )
    out$freq_type1_2st <- .frequentist_equiv_twostage(
      n1, n2, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_1 = gamma_1, gamma_eq = gamma_eq,
      analysis_prior = analysis_prior, direction = direction
    )
    out$freq_pce_1st <- .frequentist_pce_onestage(
      n, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_diff = gamma_diff, analysis_prior = analysis_prior,
      direction = direction
    )
    out$freq_pce_2st <- .frequentist_pce_twostage(
      n1, n2, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_1 = gamma_1, gamma_diff = gamma_diff,
      analysis_prior = analysis_prior, direction = direction
    )
    out$freq_EN_t1e <- .frequentist_en_twostage(
      n1, n2, p_true = p_t1e, p0 = p0, delta = delta,
      gamma_1 = gamma_1, analysis_prior = analysis_prior,
      direction = direction
    )
  }
  
  if (!is.null(p_power)) {
    out$freq_power_1st <- .frequentist_equiv_onestage(
      n, p_true = p_power, p0 = p0, delta = delta,
      gamma_eq = gamma_eq, analysis_prior = analysis_prior,
      direction = direction
    )
    out$freq_power_2st <- .frequentist_equiv_twostage(
      n1, n2, p_true = p_power, p0 = p0, delta = delta,
      gamma_1 = gamma_1, gamma_eq = gamma_eq,
      analysis_prior = analysis_prior, direction = direction
    )
    out$freq_EN_power <- .frequentist_en_twostage(
      n1, n2, p_true = p_power, p0 = p0, delta = delta,
      gamma_1 = gamma_1, analysis_prior = analysis_prior,
      direction = direction
    )
  }
  
  out
}