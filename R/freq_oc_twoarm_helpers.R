#' Exact frequentist operating characteristics for a fixed two-stage two-arm BF design
#'
#' Internal helper. Computes the exact frequentist rejection probability of the
#' two-stage design with futility stopping at interim and efficacy decision at final.
#'
#' @keywords internal
freq_oc_twostage_twoarm_fixed <- function(
    n1_1, n1_2,
    n2_1, n2_2,
    k, k_f,
    test,
    p1, p2,
    a_0_a, b_0_a,
    a_1_a, b_1_a,
    a_2_a, b_2_a,
    a_1_a_Hminus, b_1_a_Hminus,
    a_2_a_Hminus, b_2_a_Hminus
) {
  test <- match.arg(test, c("BF01", "BF+0", "BF-0", "BF+-"))
  stopifnot(n1_1 <= n2_1, n1_2 <= n2_2)
  
  reject_prob   <- 0
  stop_fut_prob <- 0
  
  # Interim futility rule: stop for compelling evidence for null / against alternative
  bf_interim_futility <- function(y1, y2) {
    if (test == "BF01") {
      bf01 <- twoarmbinbf01(
        y1 = y1, y2 = y2, n1 = n1_1, n2 = n1_2,
        a_0_a = a_0_a, b_0_a = b_0_a,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a
      )
      return(bf01 >= k_f)
    }
    if (test == "BF+0") {
      bf01 <- twoarmbinbf01(
        y1 = y1, y2 = y2, n1 = n1_1, n2 = n1_2,
        a_0_a = a_0_a, b_0_a = b_0_a,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a
      )
      # BFplus1 has _d-named formals; pass analysis-prior values
      bf_plus_1 <- BFplus1(
        y1 = y1, y2 = y2, n1 = n1_1, n2 = n1_2,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a
      )
      bf_plus_0 <- BFplus0(BFplus1 = bf_plus_1, BF01 = bf01)
      return(bf_plus_0 <= 1 / k_f)
    }
    if (test == "BF-0") {
      bf01 <- twoarmbinbf01(
        y1 = y1, y2 = y2, n1 = n1_1, n2 = n1_2,
        a_0_a = a_0_a, b_0_a = b_0_a,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a
      )
      # BFminus1 has _d-named formals; pass analysis-prior Hminus values
      bf_minus_1 <- BFminus1(
        y1 = y1, y2 = y2, n1 = n1_1, n2 = n1_2,
        a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
        a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
      )
      bf_minus_0 <- BFminus0(BFminus1 = bf_minus_1, BF01 = bf01)
      return(bf_minus_0 >= k_f)
    }
    ## test == "BF+-"
    bf_plus_1 <- BFplus1(
      y1 = y1, y2 = y2, n1 = n1_1, n2 = n1_2,
      a_1_a = a_1_a,        b_1_a = b_1_a,
      a_2_a = a_2_a,        b_2_a = b_2_a
    )
    bf_minus_1 <- BFminus1(
      y1 = y1, y2 = y2, n1 = n1_1, n2 = n1_2,
      a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
      a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
    )
    bf_pm <- BFplusMinus(BFplus1 = bf_plus_1, BFminus1 = bf_minus_1)
    bf_pm >= k_f
  }
  
  # Final decision rule
  bf_final_reject <- function(y1, y2) {
    if (test == "BF01") {
      bf01 <- twoarmbinbf01(
        y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
        a_0_a = a_0_a, b_0_a = b_0_a,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a
      )
      return(bf01 < k)
    }
    if (test == "BF+0") {
      bf01 <- twoarmbinbf01(
        y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
        a_0_a = a_0_a, b_0_a = b_0_a,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a
      )
      bf_plus_1 <- BFplus1(
        y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a
      )
      bf_plus_0 <- BFplus0(BFplus1 = bf_plus_1, BF01 = bf01)
      return(bf_plus_0 > 1 / k)
    }
    if (test == "BF-0") {
      bf01 <- twoarmbinbf01(
        y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
        a_0_a = a_0_a, b_0_a = b_0_a,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a
      )
      bf_minus_1 <- BFminus1(
        y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
        a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
        a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
      )
      bf_minus_0 <- BFminus0(BFminus1 = bf_minus_1, BF01 = bf01)
      return(bf_minus_0 < k)
    }
    ## test == "BF+-"
    bf_plus_1 <- BFplus1(
      y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
      a_1_a = a_1_a,        b_1_a = b_1_a,
      a_2_a = a_2_a,        b_2_a = b_2_a
    )
    bf_minus_1 <- BFminus1(
      y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
      a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
      a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
    )
    bf_pm <- BFplusMinus(BFplus1 = bf_plus_1, BFminus1 = bf_minus_1)
    bf_pm < k
  }
  
  for (x1 in 0:n1_1) {
    for (x2 in 0:n1_2) {
      p_stage1 <- dbinom(x1, size = n1_1, prob = p1) *
        dbinom(x2, size = n1_2, prob = p2)
      
      stop_futility <- bf_interim_futility(x1, x2)
      if (stop_futility) {
        stop_fut_prob <- stop_fut_prob + p_stage1
        next
      }
      
      for (z1 in 0:(n2_1 - n1_1)) {
        for (z2 in 0:(n2_2 - n1_2)) {
          y1 <- x1 + z1
          y2 <- x2 + z2
          
          p_stage2 <- dbinom(z1, size = n2_1 - n1_1, prob = p1) *
            dbinom(z2, size = n2_2 - n1_2, prob = p2)
          
          if (bf_final_reject(y1, y2)) {
            reject_prob <- reject_prob + p_stage1 * p_stage2
          }
        }
      }
    }
  }
  
  c(
    Reject_Prob            = reject_prob,
    Stop_for_futility_prob = stop_fut_prob
  )
}


#' Two-stage frequentist type-I error supremum over a null grid
#'
#' Internal helper. For a given two-stage design (n1_1, n1_2, n2_1, n2_2),
#' computes the supremum of the two-stage frequentist type-I error over a grid
#' of null parameter configurations, using freq_oc_twostage_twoarm_fixed().
#'
#' @keywords internal
freq_t1e_twostage_twoarm_sup <- function(
    n1_1, n1_2, n2_1, n2_2,
    k, k_f, test,
    p_null_grid  = NULL,
    a_0_a, b_0_a,
    a_1_a, b_1_a,
    a_2_a, b_2_a,
    a_1_a_Hminus, b_1_a_Hminus,
    a_2_a_Hminus, b_2_a_Hminus,
    alpha_target = NULL,
    tol_excess   = 1e-4
) {
  test <- match.arg(test, c("BF01", "BF+0", "BF-0", "BF+-"))
  
  if (is.null(p_null_grid)) {
    p_grid_fine <- seq(0.05, 0.95, by = 0.05)
  } else {
    p_grid_fine <- sort(unique(p_null_grid))
  }
  
  p_grid_coarse <- sort(unique(
    c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)
  ))
  p_grid_coarse <- p_grid_coarse[
    p_grid_coarse >= min(p_grid_fine) & p_grid_coarse <= max(p_grid_fine)
  ]
  if (length(p_grid_coarse) == 0L) p_grid_coarse <- p_grid_fine
  
  # Helper: two-stage type-I error at a single (p1, p2) point
  t1e_two_at <- function(p1, p2) {
    freq_oc_twostage_twoarm_fixed(
      n1_1 = n1_1, n1_2 = n1_2,
      n2_1 = n2_1, n2_2 = n2_2,
      k = k, k_f = k_f, test = test,
      p1 = p1, p2 = p2,
      a_0_a = a_0_a, b_0_a = b_0_a,
      a_1_a = a_1_a, b_1_a = b_1_a,
      a_2_a = a_2_a, b_2_a = b_2_a,
      a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
      a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
    )[["Reject_Prob"]]
  }
  
  # ── BF01: supremum on the diagonal p1 = p2 ───────────────────────────────────
  if (test == "BF01") {
    coarse_vals <- vapply(p_grid_coarse, function(p0) t1e_two_at(p0, p0), numeric(1))
    
    max_coarse <- max(coarse_vals, na.rm = TRUE)
    if (!is.null(alpha_target) && max_coarse > alpha_target + tol_excess) {
      return(max_coarse)
    }
    
    max_idx        <- which(coarse_vals >= max_coarse - 1e-8)
    p0_star_coarse <- p_grid_coarse[max_idx]
    
    refine_points <- numeric(0L)
    for (p0_star in p0_star_coarse) {
      refine_points <- c(
        refine_points,
        p_grid_fine[p_grid_fine >= (p0_star - 0.05) &
                      p_grid_fine <= (p0_star + 0.05)]
      )
    }
    refine_points <- sort(unique(refine_points))
    if (length(refine_points) == 0L) refine_points <- p_grid_fine
    
    refine_vals <- vapply(refine_points, function(p0) t1e_two_at(p0, p0), numeric(1))
    return(max(refine_vals, na.rm = TRUE))
  }
  
  # ── Directional / two-sided: supremum over (p1, p2) grid ─────────────────────
  grid_pairs_coarse <- expand.grid(p1 = p_grid_coarse, p2 = p_grid_coarse)
  if (test == "BF+0") {
    grid_pairs_coarse <- grid_pairs_coarse[
      grid_pairs_coarse$p2 <= grid_pairs_coarse$p1, , drop = FALSE
    ]
  } else if (test == "BF-0") {
    grid_pairs_coarse <- grid_pairs_coarse[
      grid_pairs_coarse$p2 >= grid_pairs_coarse$p1, , drop = FALSE
    ]
  }
  
  coarse_vals <- vapply(
    seq_len(nrow(grid_pairs_coarse)),
    function(i) t1e_two_at(grid_pairs_coarse$p1[i], grid_pairs_coarse$p2[i]),
    numeric(1)
  )
  
  max_coarse <- max(coarse_vals, na.rm = TRUE)
  if (!is.null(alpha_target) && max_coarse > alpha_target + tol_excess) {
    return(max_coarse)
  }
  
  max_idx <- which(coarse_vals >= max_coarse - 1e-8)
  p1_star <- grid_pairs_coarse$p1[max_idx]
  p2_star <- grid_pairs_coarse$p2[max_idx]
  
  grid_pairs_fine <- expand.grid(p1 = p_grid_fine, p2 = p_grid_fine)
  if (test == "BF+0") {
    grid_pairs_fine <- grid_pairs_fine[
      grid_pairs_fine$p2 <= grid_pairs_fine$p1, , drop = FALSE
    ]
  } else if (test == "BF-0") {
    grid_pairs_fine <- grid_pairs_fine[
      grid_pairs_fine$p2 >= grid_pairs_fine$p1, , drop = FALSE
    ]
  }
  
  refine_rows <- rep(FALSE, nrow(grid_pairs_fine))
  for (j in seq_along(p1_star)) {
    refine_rows <- refine_rows |
      (grid_pairs_fine$p1 >= (p1_star[j] - 0.05) &
         grid_pairs_fine$p1 <= (p1_star[j] + 0.05) &
         grid_pairs_fine$p2 >= (p2_star[j] - 0.05) &
         grid_pairs_fine$p2 <= (p2_star[j] + 0.05))
  }
  if (!any(refine_rows)) refine_rows[] <- TRUE
  
  grid_pairs_refine <- grid_pairs_fine[refine_rows, , drop = FALSE]
  
  refine_vals <- vapply(
    seq_len(nrow(grid_pairs_refine)),
    function(i) t1e_two_at(grid_pairs_refine$p1[i], grid_pairs_refine$p2[i]),
    numeric(1)
  )
  
  max(refine_vals, na.rm = TRUE)
}


#' Exact frequentist OCs for one two-stage two-arm BF design
#'
#' @keywords internal
compute_freq_twostage_oc_2arm <- function(
    n1_1, n1_2, n2_1, n2_2,
    k, k_f, test,
    p1_power, p2_power,
    p_null_grid  = NULL,
    a_0_a, b_0_a,
    a_1_a, b_1_a,
    a_2_a, b_2_a,
    a_1_a_Hminus, b_1_a_Hminus,
    a_2_a_Hminus, b_2_a_Hminus
) {
  test <- match.arg(test, c("BF01", "BF+0", "BF-0", "BF+-"))
  stopifnot(n1_1 <= n2_1, n1_2 <= n2_2)
  
  final_n1 <- n2_1
  final_n2 <- n2_2
  
  # Fixed-sample frequentist power at (p1_power, p2_power)
  freq_fix_alt <- freq_oc_twoarm_fixed(
    n1 = final_n1, n2 = final_n2,
    k = k, k_f = k_f, test = test,
    p1_null = p1_power, p2_null = p2_power,
    p1_alt  = p1_power, p2_alt  = p2_power,
    a_0_a = a_0_a, b_0_a = b_0_a,
    a_1_a = a_1_a, b_1_a = b_1_a,
    a_2_a = a_2_a, b_2_a = b_2_a,
    a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
    a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
  )
  power_fix <- unname(freq_fix_alt[["Power_freq"]])
  
  # Fixed-sample supremum type-I error over p_null_grid
  t1e_fix <- freq_t1e_sup_fixed(
    n1 = final_n1, n2 = final_n2,
    k = k, k_f = k_f, test = test,
    p_null_grid = p_null_grid,
    a_0_a = a_0_a, b_0_a = b_0_a,
    a_1_a = a_1_a, b_1_a = b_1_a,
    a_2_a = a_2_a, b_2_a = b_2_a,
    a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
    a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
  )
  
  # Two-stage frequentist power at (p1_power, p2_power)
  freq_two <- freq_oc_twostage_twoarm_fixed(
    n1_1 = n1_1, n1_2 = n1_2,
    n2_1 = n2_1, n2_2 = n2_2,
    k = k, k_f = k_f, test = test,
    p1 = p1_power, p2 = p2_power,
    a_0_a = a_0_a, b_0_a = b_0_a,
    a_1_a = a_1_a, b_1_a = b_1_a,
    a_2_a = a_2_a, b_2_a = b_2_a,
    a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
    a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
  )
  power_two <- unname(freq_two[["Reject_Prob"]])
  
  # Two-stage supremum type-I error
  t1e_two <- freq_t1e_twostage_twoarm_sup(
    n1_1 = n1_1, n1_2 = n1_2,
    n2_1 = n2_1, n2_2 = n2_2,
    k = k, k_f = k_f, test = test,
    p_null_grid = p_null_grid,
    a_0_a = a_0_a, b_0_a = b_0_a,
    a_1_a = a_1_a, b_1_a = b_1_a,
    a_2_a = a_2_a, b_2_a = b_2_a,
    a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
    a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
  )
  
  c(
    Type1_Error_freq_fixed     = t1e_fix,
    Power_freq_fixed           = power_fix,
    Type1_Error_freq_two_stage = t1e_two,
    Power_freq_two_stage       = power_two
  )
}