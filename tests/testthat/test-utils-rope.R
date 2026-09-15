test_that("explicit equivalence direction reproduces the default", {
  expect_equal(
    posterior_rope_prob(
      y = 20,
      n = 60,
      p0 = 0.30,
      delta = 0.10,
      analysis_prior = c(1, 1)
    ),
    posterior_rope_prob(
      y = 20,
      n = 60,
      p0 = 0.30,
      delta = 0.10,
      analysis_prior = c(1, 1),
      direction = "equivalence"
    )
  )
})

test_that("noninferiority posterior probability equals beta upper-tail mass", {
  y <- 12
  n <- 40
  p0 <- 0.30
  delta <- 0.10
  analysis_prior <- c(1, 1)
  
  expected <- 1 - pbeta(
    p0 - delta,
    analysis_prior[1] + y,
    analysis_prior[2] + n - y
  )
  
  observed <- posterior_rope_prob(
    y = y,
    n = n,
    p0 = p0,
    delta = delta,
    analysis_prior = analysis_prior,
    direction = "noninferiority"
  )
  
  expect_equal(observed, expected)
})

test_that("superiority posterior probability equals beta upper-tail mass", {
  y <- 22
  n <- 40
  p0 <- 0.30
  delta <- 0.10
  analysis_prior <- c(1, 1)
  
  expected <- 1 - pbeta(
    p0 + delta,
    analysis_prior[1] + y,
    analysis_prior[2] + n - y
  )
  
  observed <- posterior_rope_prob(
    y = y,
    n = n,
    p0 = p0,
    delta = delta,
    analysis_prior = analysis_prior,
    direction = "superiority"
  )
  
  expect_equal(observed, expected)
})

test_that("one-sided posterior H1 support increases with response count", {
  ni_prob <- vapply(
    0:30,
    function(y) {
      posterior_rope_prob(
        y = y,
        n = 30,
        p0 = 0.30,
        delta = 0.10,
        analysis_prior = c(1, 1),
        direction = "noninferiority"
      )
    },
    numeric(1L)
  )
  
  sup_prob <- vapply(
    0:30,
    function(y) {
      posterior_rope_prob(
        y = y,
        n = 30,
        p0 = 0.30,
        delta = 0.10,
        analysis_prior = c(1, 1),
        direction = "superiority"
      )
    },
    numeric(1L)
  )
  
  expect_true(all(diff(ni_prob) >= 0))
  expect_true(all(diff(sup_prob) >= 0))
})

test_that("invalid direction-specific boundaries generate errors", {
  expect_error(
    design_singlearm_twostage_rope(
      p0 = 0.10,
      delta = 0.10,
      analysis_prior = c(1, 1),
      design_prior_h0 = c(5, 95),
      design_prior_h1 = c(30, 70),
      direction = "noninferiority"
    ),
    "p0 \\- delta"
  )
  
  expect_error(
    design_singlearm_twostage_rope(
      p0 = 0.90,
      delta = 0.10,
      analysis_prior = c(1, 1),
      design_prior_h0 = c(30, 70),
      design_prior_h1 = c(80, 20),
      direction = "superiority"
    ),
    "p0 \\+ delta"
  )
})

test_that("one-sided continuation region uses posterior support for H0 (noninferiority)", {
  n1 <- 30
  p0 <- 0.30
  delta <- 0.10
  gamma_1 <- 0.80
  analysis_prior <- c(1, 1)
  
  expected_ni <- (0:n1)[vapply(
    0:n1,
    function(y1) {
      posterior_rope_prob(
        y = y1,
        n = n1,
        p0 = p0,
        delta = delta,
        analysis_prior = analysis_prior,
        direction = "noninferiority"
      ) > 1 - gamma_1
    },
    logical(1L)
  )]
  
  observed_ni <- .continuation_region_twostage(
    n1 = n1,
    p0 = p0,
    delta = delta,
    gamma_1 = gamma_1,
    analysis_prior = analysis_prior,
    direction = "noninferiority"
  )
  
  expect_identical(observed_ni, expected_ni)
})

test_that("one-sided continuation region uses posterior support for H0 (superiority)", {
  n1 <- 30
  p0 <- 0.30
  delta <- 0.10
  gamma_1 <- 0.80
  analysis_prior <- c(1, 1)
  
  expected_ni <- (0:n1)[vapply(
    0:n1,
    function(y1) {
      posterior_rope_prob(
        y = y1,
        n = n1,
        p0 = p0,
        delta = delta,
        analysis_prior = analysis_prior,
        direction = "superiority"
      ) > 1 - gamma_1
    },
    logical(1L)
  )]
  
  observed_ni <- .continuation_region_twostage(
    n1 = n1,
    p0 = p0,
    delta = delta,
    gamma_1 = gamma_1,
    analysis_prior = analysis_prior,
    direction = "superiority"
  )
  
  expect_identical(observed_ni, expected_ni)
})
