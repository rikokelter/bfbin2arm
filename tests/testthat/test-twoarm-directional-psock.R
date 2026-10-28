test_that("PSOCK parallel search matches serial search for a SINGLE-inspired grid", {
  skip_unless_slow_tests()
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "Requires at least 2 cores.")
  
  design_prior_ni <- list(
    type = "baseline_difference",
    baseline = list(mean = 0.75, sd = 0.02, truncation = c(0, 1)),
    difference = list(mean = 0, sd = 1e-6, truncation = c(-0.10, 1))
  )
  
  design_prior_inf <- list(
    type = "baseline_difference",
    baseline = list(mean = 0.75, sd = 0.02, truncation = c(0, 1)),
    difference = list(mean = -0.15, sd = 1e-6, truncation = c(-1, -0.10))
  )
  
  common_args <- list(
    nmin = 40L, nmax = 400L, nstep = 20L,
    direction = "noninferiority",
    delta_ni = 0.10, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1),
    design_prior_ni = design_prior_ni,
    design_prior_inf = design_prior_inf,
    target_ni_power = 0.80,
    target_false_ni = 0.05,
    calibration = "Bayesian",
    quad_nodes = 64L,
    progress = FALSE,
    returngrid = TRUE
  )
  
  fit_serial <- do.call(
    design_twoarm_onestage_rope_directional,
    c(common_args, list(parallel = FALSE))
  )
  
  fit_parallel <- do.call(
    design_twoarm_onestage_rope_directional,
    c(common_args, list(parallel = TRUE, ncores = 2L, parallel_backend = "psock"))
  )
  
  expect_false(is.na(fit_serial$nstar))
  expect_equal(fit_serial$nstar, fit_parallel$nstar)
  
  merged <- merge(
    fit_serial$grid[, c("n", "ni_power", "false_ni")],
    fit_parallel$grid[, c("n", "ni_power", "false_ni")],
    by = "n", suffixes = c("_serial", "_parallel")
  )
  
  expect_equal(merged$ni_power_serial, merged$ni_power_parallel, tolerance = 1e-8)
  expect_equal(merged$false_ni_serial, merged$false_ni_parallel, tolerance = 1e-8)
})

test_that("PSOCK parallel search matches serial search for a nested design", {
  skip_unless_slow_tests()
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "Requires at least 2 cores.")
  
  design_prior_ni <- list(
    type = "baseline_difference",
    baseline = list(mean = 0.63, sd = 0.02, truncation = c(0, 1)),
    difference = list(mean = 0.028, sd = 1e-6, truncation = c(-0.0397, 1))
  )
  
  design_prior_sup <- list(
    type = "baseline_difference",
    baseline = list(mean = 0.63, sd = 0.02, truncation = c(0, 1)),
    difference = list(mean = 0.06, sd = 1e-6, truncation = c(0, 1))
  )
  
  common_args <- list(
    nmin = 50L, nmax = 800L, nstep = 25L,
    direction = "nested",
    delta_ni = 0.0397, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1),
    design_prior_ni = design_prior_ni,
    design_prior_sup = design_prior_sup,
    target_ni_power = 0.70,
    target_sup_power = 0.20,
    calibration = "Bayesian",
    quad_nodes = 64L,
    progress = FALSE,
    returngrid = TRUE
  )
  
  fit_serial <- do.call(
    design_twoarm_onestage_rope_directional,
    c(common_args, list(parallel = FALSE))
  )
  
  fit_parallel <- do.call(
    design_twoarm_onestage_rope_directional,
    c(common_args, list(parallel = TRUE, ncores = 2L, parallel_backend = "psock"))
  )
  
  expect_false(is.na(fit_serial$nstar))
  expect_equal(fit_serial$nstar, fit_parallel$nstar)
  
  merged <- merge(
    fit_serial$grid[, c("n", "ni_power", "sup_power")],
    fit_parallel$grid[, c("n", "ni_power", "sup_power")],
    by = "n", suffixes = c("_serial", "_parallel")
  )
  
  expect_equal(merged$ni_power_serial, merged$ni_power_parallel, tolerance = 1e-8)
  expect_equal(merged$sup_power_serial, merged$sup_power_parallel, tolerance = 1e-8)
})

