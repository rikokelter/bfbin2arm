# Internal helper: frequentist OC for a fixed two-arm design
#
# Computes frequentist type-I error and power for a given fixed-sample
# two-arm Bayes-factor design (no interim), by exhaustive enumeration
# over all y1, y2.
#
# Bayes factors are always evaluated using analysis priors (_a / _a_Hminus).
# Design priors (_d) play no role here; they are only used in the Bayesian
# predictive-density computations (power_twostage_2arm.R).
#
# Arguments:
#   n1, n2         final sample sizes per arm
#   k, k_f         BF thresholds (efficacy and futility)
#   test           one of "BF01" / "BF+0" / "BF-0" / "BF+-"
#   p1_null        control response prob under H0 (frequentist)
#   p2_null        treatment response prob under H0 (frequentist)
#   p1_alt         control response prob under H1 (frequentist)
#   p2_alt         treatment response prob under H1 (frequentist)
#   a_0_a, b_0_a   analysis prior hyperparameters under H0 (common p)
#   a_1_a, b_1_a   analysis prior hyperparameters for arm 1 under H+ / H1
#   a_2_a, b_2_a   analysis prior hyperparameters for arm 2 under H+ / H1
#   a_1_a_Hminus, b_1_a_Hminus
#   a_2_a_Hminus, b_2_a_Hminus
#                  analysis prior hyperparameters under H- (arm 1 / arm 2)
#
# Returns a named numeric vector:
#   c(Type1_Error_freq = ..., Power_freq = ...)
#
freq_oc_twoarm_fixed <- function(
    n1, n2,
    k, k_f, test,
    p1_null, p2_null,
    p1_alt,  p2_alt,
    a_0_a, b_0_a,
    a_1_a, b_1_a,
    a_2_a, b_2_a,
    a_1_a_Hminus, b_1_a_Hminus,
    a_2_a_Hminus, b_2_a_Hminus
) {
  test <- match.arg(test, c("BF01", "BF+0", "BF-0", "BF+-"))
  
  reject_prob_null <- 0
  reject_prob_alt  <- 0
  
  for (y1 in 0:n1) {
    for (y2 in 0:n2) {
      
      if (test == "BF01") {
        
        bf01   <- twoarmbinbf01(
          y1 = y1, y2 = y2, n1 = n1, n2 = n2,
          a_0_a = a_0_a, b_0_a = b_0_a,
          a_1_a = a_1_a, b_1_a = b_1_a,
          a_2_a = a_2_a, b_2_a = b_2_a
        )
        reject <- (bf01 < k)
        
      } else if (test == "BF+0") {
        
        # twoarmbinbf_plus0_direct now takes _a parameters
        bf_plus_0 <- twoarmbinbf_plus0_direct(
          y1 = y1, y2 = y2, n1 = n1, n2 = n2,
          a_0_a = a_0_a, b_0_a = b_0_a,
          a_1_a = a_1_a, b_1_a = b_1_a,
          a_2_a = a_2_a, b_2_a = b_2_a
        )
        reject <- (bf_plus_0 > 1 / k)
        
      } else if (test == "BF-0") {
        
        bf01       <- twoarmbinbf01(
          y1 = y1, y2 = y2, n1 = n1, n2 = n2,
          a_0_a = a_0_a, b_0_a = b_0_a,
          a_1_a = a_1_a, b_1_a = b_1_a,
          a_2_a = a_2_a, b_2_a = b_2_a
        )
        # BFminus1 has _d-named formals; pass analysis-prior values into them
        bf_minus_1 <- BFminus1(
          y1 = y1, y2 = y2, n1 = n1, n2 = n2,
          a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
          a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
        )
        bf_minus_0 <- BFminus0(BFminus1 = bf_minus_1, BF01 = bf01)
        reject     <- (bf_minus_0 < k)
        
      } else if (test == "BF+-") {
        
        # BFplus1 / BFminus1 have _d-named formals; pass analysis-prior values
        bf_plus_1  <- BFplus1(
          y1 = y1, y2 = y2, n1 = n1, n2 = n2,
          a_1_a = a_1_a,        b_1_a = b_1_a,
          a_2_a = a_2_a,        b_2_a = b_2_a
        )
        bf_minus_1 <- BFminus1(
          y1 = y1, y2 = y2, n1 = n1, n2 = n2,
          a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
          a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
        )
        bf_pm  <- BFplusMinus(BFplus1 = bf_plus_1, BFminus1 = bf_minus_1)
        reject <- (bf_pm < k)
      }
      
      if (!reject) next
      
      p_y1_null <- dbinom(y1, size = n1, prob = p1_null)
      p_y2_null <- dbinom(y2, size = n2, prob = p2_null)
      p_y1_alt  <- dbinom(y1, size = n1, prob = p1_alt)
      p_y2_alt  <- dbinom(y2, size = n2, prob = p2_alt)
      
      reject_prob_null <- reject_prob_null + p_y1_null * p_y2_null
      reject_prob_alt  <- reject_prob_alt  + p_y1_alt  * p_y2_alt
    }
  }
  
  c(
    Type1_Error_freq = reject_prob_null,
    Power_freq       = reject_prob_alt
  )
}


#' Fixed-sample frequentist type-I error supremum over a null grid
#'
#' Internal helper. For a given fixed-sample design (n1, n2), computes
#' the supremum of the frequentist type-I error over a grid of null
#' parameter configurations, using \code{freq_oc_twoarm_fixed()}.
#'
#' @keywords internal
freq_t1e_sup_fixed <- function(
    n1, n2,
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
  
  # Single helper that calls freq_oc_twoarm_fixed with the unified signature
  t1e_at <- function(p1, p2) {
    freq_oc_twoarm_fixed(
      n1 = n1, n2 = n2,
      k = k, k_f = k_f, test = test,
      p1_null = p1, p2_null = p2,
      p1_alt  = p1, p2_alt  = p2,
      a_0_a = a_0_a, b_0_a = b_0_a,
      a_1_a = a_1_a, b_1_a = b_1_a,
      a_2_a = a_2_a, b_2_a = b_2_a,
      a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
      a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
    )[["Type1_Error_freq"]]
  }
  
  # ── BF01: supremum is on the diagonal p1 = p2 = p0 ──────────────────────────
  if (test == "BF01") {
    
    coarse_vals <- vapply(p_grid_coarse, function(p0) t1e_at(p0, p0), numeric(1))
    
    max_coarse <- max(coarse_vals, na.rm = TRUE)
    if (!is.null(alpha_target) && max_coarse > alpha_target + tol_excess) {
      return(max_coarse)
    }
    
    max_idx       <- which(coarse_vals >= max_coarse - 1e-8)
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
    
    refine_vals <- vapply(refine_points, function(p0) t1e_at(p0, p0), numeric(1))
    return(max(refine_vals, na.rm = TRUE))
  }
  
  # ── Directional / two-sided directional: supremum over (p1, p2) grid ────────
  
  # Build filtered coarse grid of (p1, p2) pairs
  grid_pairs_coarse <- expand.grid(p1 = p_grid_coarse, p2 = p_grid_coarse)
  if (test == "BF+0") {
    # H+ is rejected when p2 <= p1 (null region for H+)
    grid_pairs_coarse <- grid_pairs_coarse[
      grid_pairs_coarse$p2 <= grid_pairs_coarse$p1, , drop = FALSE
    ]
  } else if (test == "BF-0") {
    # H- is rejected when p2 >= p1 (null region for H-)
    grid_pairs_coarse <- grid_pairs_coarse[
      grid_pairs_coarse$p2 >= grid_pairs_coarse$p1, , drop = FALSE
    ]
  }
  # BF+-: use full grid (two-sided, any (p1, p2) can yield type-I error)
  
  coarse_vals <- vapply(
    seq_len(nrow(grid_pairs_coarse)),
    function(i) t1e_at(grid_pairs_coarse$p1[i], grid_pairs_coarse$p2[i]),
    numeric(1)
  )
  
  max_coarse <- max(coarse_vals, na.rm = TRUE)
  if (!is.null(alpha_target) && max_coarse > alpha_target + tol_excess) {
    return(max_coarse)
  }
  
  max_idx <- which(coarse_vals >= max_coarse - 1e-8)
  p1_star <- grid_pairs_coarse$p1[max_idx]
  p2_star <- grid_pairs_coarse$p2[max_idx]
  
  # Build fine grid and apply same directional filter
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
  
  # Restrict fine grid to neighbourhoods of the coarse maximisers
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
    function(i) t1e_at(grid_pairs_refine$p1[i], grid_pairs_refine$p2[i]),
    numeric(1)
  )
  
  max(refine_vals, na.rm = TRUE)
}