test_that(
  "targetpce_h0 is applied to Bayesian feasibility",
  {
    prior_eq_test <- list(
      type = "baseline_difference",
      
      baseline = list(
        family = "beta",
        shape1 = 37,
        shape2 = 13
      ),
      
      difference = list(
        family = "normal",
        mean = 0,
        sd = 0.03,
        truncation = c(-0.10, 0.10)
      )
    )
    
    result <- design_twoarm_onestage_rope(
      nmin = 50L,
      nmax = 60L,
      nstep = 5L,
      
      delta = 0.10,
      gammaeq = 0.80,
      gammadiff = 0.975,
      
      analysis_prior_C = c(1, 1),
      analysis_prior_T = c(1, 1),
      
      design_prior_eq = prior_eq_test,
      design_prior_ne_C = c(37, 13),
      design_prior_ne_T = c(42, 8),
      
      calibration = "Bayesian",
      
      targetpower = 0.80,
      targettype1 = 0.05,
      targetpce_h0 = 0.50,
      
      sustainn = 1L,
      
      integration = "quantile",
      quad_nodes = 16L,
      
      parallel = FALSE,
      parallel_backend = "none",
      ncores = 1L,
      
      returngrid = TRUE,
      returnmatrices = FALSE,
      
      progress = FALSE
    )
    
    expect_true(
      "pce_h0" %in% names(result$grid)
    )
    
    expect_true(
      "pce_h0_ok" %in% names(result$grid)
    )
    
    expect_true(
      "bayes_ok" %in% names(result$grid)
    )
    
    expect_equal(
      result$grid$pce_h0_ok,
      result$grid$pce_h0 >= 0.50
    )
    
    expect_true(
      all(
        result$grid$bayes_ok <=
          result$grid$pce_h0_ok
      )
    )
  }
)

test_that(
  "NULL targetpce_h0 does not constrain Bayesian feasibility",
  {
    prior_eq_test <- list(
      type = "baseline_difference",
      
      baseline = list(
        family = "beta",
        shape1 = 37,
        shape2 = 13
      ),
      
      difference = list(
        family = "normal",
        mean = 0,
        sd = 0.03,
        truncation = c(-0.10, 0.10)
      )
    )
    
    result <- design_twoarm_onestage_rope(
      nmin = 50L,
      nmax = 60L,
      nstep = 5L,
      
      delta = 0.10,
      gammaeq = 0.80,
      gammadiff = 0.975,
      
      analysis_prior_C = c(1, 1),
      analysis_prior_T = c(1, 1),
      
      design_prior_eq = prior_eq_test,
      design_prior_ne_C = c(37, 13),
      design_prior_ne_T = c(42, 8),
      
      calibration = "Bayesian",
      
      targetpower = 0.80,
      targettype1 = 0.05,
      targetpce_h0 = NULL,
      
      sustainn = 1L,
      
      integration = "quantile",
      quad_nodes = 16L,
      
      parallel = FALSE,
      parallel_backend = "none",
      ncores = 1L,
      
      returngrid = TRUE,
      returnmatrices = FALSE,
      
      progress = FALSE
    )
    
    expect_true(
      "pce_h0_ok" %in% names(result$grid)
    )
    
    expect_true(
      all(result$grid$pce_h0_ok)
    )
  }
)

test_that(
  "PSOCK and serial searches agree with targetpce_h0",
  {
    skip_on_cran()
    
    skip_if(
      parallel::detectCores(logical = FALSE) < 2L,
      "At least two logical cores are required for this PSOCK test."
    )
    
    prior_eq_test <- list(
      type = "baseline_difference",
      
      baseline = list(
        family = "beta",
        shape1 = 37,
        shape2 = 13
      ),
      
      difference = list(
        family = "normal",
        mean = 0,
        sd = 0.03,
        truncation = c(-0.10, 0.10)
      )
    )
    
    common_args <- list(
      nmin = 50L,
      nmax = 60L,
      nstep = 5L,
      
      delta = 0.10,
      gammaeq = 0.80,
      gammadiff = 0.975,
      
      analysis_prior_C = c(1, 1),
      analysis_prior_T = c(1, 1),
      
      design_prior_eq = prior_eq_test,
      design_prior_ne_C = c(37, 13),
      design_prior_ne_T = c(42, 8),
      
      calibration = "Bayesian",
      
      targetpower = 0.80,
      targettype1 = 0.05,
      targetpce_h0 = 0.50,
      
      sustainn = 1L,
      
      integration = "quantile",
      quad_nodes = 16L,
      
      returngrid = TRUE,
      returnmatrices = FALSE,
      
      progress = FALSE
    )
    
    serial <- do.call(
      design_twoarm_onestage_rope,
      c(
        common_args,
        list(
          parallel = FALSE,
          parallel_backend = "none",
          ncores = 1L
        )
      )
    )
    
    psock <- do.call(
      design_twoarm_onestage_rope,
      c(
        common_args,
        list(
          parallel = TRUE,
          parallel_backend = "psock",
          ncores = 2L
        )
      )
    )
    
    expect_equal(
      serial$grid,
      psock$grid,
      tolerance = 1e-12
    )
    
    expect_identical(
      serial$nstar,
      psock$nstar
    )
    
    expect_equal(
      serial$selected,
      psock$selected,
      tolerance = 1e-12
    )
  }
)
