test_that("powerbinbf01seq satisfies basic inequalities and bounds", {
  res <- powerbinbf01seq(
    n1 = 12,
    n2 = 28,
    k  = 1/10,
    kf = 3,
    p0 = 0.2,
    a0 = 1, b0 = 1,
    a1 = 1, b1 = 1,
    da0 = 1, db0 = 1,
    da1 = 1, db1 = 1,
    dp  = NA_real_,         # or set a specific design p if desired
    type = "direction",
    k_ce = 10
  )
  
  expect_lte(res$pfineff,  res$pnaive + 1e-12)
  expect_lte(res$pfineff0, res$pnaive0 + 1e-12)
  expect_gte(res$perased,  -1e-12)
  expect_gte(res$perased0, -1e-12)
  
  expect_gte(res$nexp,  res$n1 - 1e-12)
  expect_lte(res$nexp,  res$n2 + 1e-12)
  expect_gte(res$nexp0, res$n1 - 1e-12)
  expect_lte(res$nexp0, res$n2 + 1e-12)
  
  if (!is.na(res$pce0_naive) && !is.na(res$pce0_corr)) {
    expect_gte(res$pce0_corr, res$pce0_naive - 1e-12)
  }
})