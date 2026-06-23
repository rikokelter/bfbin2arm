test_that("single-arm prior-predictive pmf sums to 1 and marginalizes correctly", {
  tol <- 1e-10
  n   <- 20
  
  # truncated Beta design prior
  pmf1 <- singlearm_priorpred_pmf(
    x  = 0:n,
    n  = n,
    da = 2, db = 3,
    dl = 0.1, du = 0.9,
    dp = NA_real_
  )
  expect_lt(abs(sum(pmf1) - 1), tol)
  
  # point design prior
  pmf2 <- singlearm_priorpred_pmf(
    x  = 0:n,
    n  = n,
    dp = 0.35
  )
  expect_lt(abs(sum(pmf2) - 1), tol)
  
  # joint truncated-beta predictive sums to 1 and marginalizes to pmf1 at n1
  n1 <- 10
  n2 <- 25
  
  joint_sum <- 0
  pmf_x     <- numeric(n1 + 1L)
  
  for (x in 0:n1) {
    inner <- numeric(n2 - n1 + 1L)
    for (z in 0:(n2 - n1)) {
      val <- dbinbin_truncbeta(
        x     = x,
        z     = z,
        n1    = n1,
        n2    = n2,
        a     = 2,
        b     = 3,
        lower = 0.1,
        upper = 0.9
      )
      joint_sum <- joint_sum + val
      inner[z + 1L] <- val
    }
    pmf_x[x + 1L] <- sum(inner)
  }
  
  expect_lt(abs(joint_sum - 1), tol)
  
  pmf_x_ref <- singlearm_priorpred_pmf(
    x  = 0:n1,
    n  = n1,
    da = 2, db = 3,
    dl = 0.1, du = 0.9,
    dp = NA_real_
  )
  
  expect_lt(max(abs(pmf_x - pmf_x_ref)), tol)
})

test_that("single-arm direction BF01 is monotone in x for superiority direction", {
  n   <- 20
  p0  <- 0.2
  bfs <- singlearm_bf01(
    x    = 0:n,
    n    = n,
    p0   = p0,
    a0   = 1,
    b0   = 1,
    a1   = 1,
    b1   = 1,
    type = "direction"
  )
  
  # For H0: p <= p0 vs H1: p > p0, BF01 should decrease with x (up to numerical noise)
  expect_true(all(diff(bfs) <= 1e-10))
})