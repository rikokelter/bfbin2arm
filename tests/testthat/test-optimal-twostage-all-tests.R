test_that("optimal_twostage_2arm_bf behaves sensibly for all supported Bayes factor tests", {
  skip_on_cran()
  
  test_opts <- c("BF01", "BF+0", "BF-0", "BF+-")
  
  common_args <- list(
    alpha = 0.10,
    beta = 0.20,
    k = 1/3,
    k_f = 3,
    n1_min = c(3, 3),
    n2_max = c(10, 10),
    alloc1 = 0.5,
    alloc2 = 0.5,
    power_cushion = 0,
    pceH0 = 0.50,
    interim_fraction = c(0.25, 0.75),
    grid_step = 1,
    progress = FALSE,
    max_iter = 100L,
    compute_freq_oc = FALSE,
    a_0_d = 1, b_0_d = 1,
    a_0_a = 1, b_0_a = 1,
    a_1_d = 1, b_1_d = 1,
    a_2_d = 1, b_2_d = 1,
    a_1_a = 1, b_1_a = 1,
    a_2_a = 1, b_2_a = 1,
    a_1_d_Hminus = 1, b_1_d_Hminus = 1,
    a_2_d_Hminus = 1, b_2_d_Hminus = 1,
    a_1_a_Hminus = 1, b_1_a_Hminus = 1,
    a_2_a_Hminus = 1, b_2_a_Hminus = 1
  )
  
  for (tt in test_opts) {
    args <- common_args
    args$test <- tt
    
    res <- do.call(optimal_twostage_2arm_bf, args)
    
    # Always expect a list with basic components
    expect_type(res, "list")
    expect_true(all(c("design", "occ", "naive_oc", "priors") %in% names(res)))
    
    # Check that the test label is carried through
    expect_true(is.list(res$priors))
    expect_identical(res$priors$test, tt)
    
    # If no feasible design was found, design or occ may contain NA/NULL.
    # In that case, just check that this is reflected by NA / NULL and
    # skip the stronger assertions below.
    if (anyNA(res$design) || is.null(res$occ) || anyNA(unlist(res$occ))) {
      next
    }
    
    # From here on, we assume a feasible design.
    expect_true(is.numeric(res$design))
    expect_length(res$design, 4)
    expect_false(anyNA(res$design))
    expect_true(all(res$design >= 1))
    expect_true(res$design[1] <= res$design[3])
    expect_true(res$design[2] <= res$design[4])
    
    # Occ joint vector
    expect_true(is.numeric(res$occ))
    expect_true(length(res$occ) >= 4)
    expect_false(anyNA(res$occ))
    
    occ_names <- names(res$occ)
    expect_true(!is.null(occ_names))
    
    get_occ <- function(x, aliases) {
      for (nm in aliases) {
        if (!is.na(match(nm, names(x)))) return(as.numeric(x[[nm]]))
      }
      NA_real_
    }
    
    power_val <- get_occ(res$occ, c("Power"))
    t1e_val   <- get_occ(res$occ, c("Type1_Error", "Type1Error"))
    ceh0_val  <- get_occ(res$occ, c("CE_H0", "CEH0"))
    enh0_val  <- get_occ(res$occ, c("E_H0_N", "EH0N"))
    fut_val   <- get_occ(res$occ, c("futility_prob", "futilityprob"))
    
    expect_true(is.finite(power_val))
    expect_true(is.finite(t1e_val))
    expect_true(is.finite(ceh0_val))
    expect_true(is.finite(enh0_val))
    
    expect_gte(power_val, 0); expect_lte(power_val, 1)
    expect_gte(t1e_val,   0); expect_lte(t1e_val,   1)
    expect_gte(ceh0_val,  0); expect_lte(ceh0_val,  1)
    
    expect_gte(enh0_val, sum(res$design[1:2]))
    expect_lte(enh0_val, sum(res$design[3:4]))
    
    if (is.finite(fut_val)) {
      expect_gte(fut_val, 0)
      expect_lte(fut_val, 1)
    }
    
    # Naive fixed-sample calibration
    expect_true(is.list(res$naive_oc))
    expect_true(all(c("n1", "n2", "power", "t1e", "pceH0") %in% names(res$naive_oc)))
    
    expect_true(is.finite(res$naive_oc$power))
    expect_true(is.finite(res$naive_oc$t1e))
    expect_true(is.finite(res$naive_oc$pceH0))
  }
})