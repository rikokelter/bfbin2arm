#' Calibrate an optimal single-arm two-stage ROPE design
#'
#' Finds a single-arm two-stage Bayesian design based on the region of
#' practical equivalence (ROPE) for a binary endpoint, with a single interim
#' analysis allowing early stopping for futility. The design covers three
#' decision types via the \code{direction} argument:
#' \describe{
#' \item{\code{"equivalence"}}{Posterior mass inside the two-sided ROPE
#' \eqn{[p_0 - \delta,\, p_0 + \delta]} must exceed \eqn{\gamma_{\mathrm{eq}}}.}
#' \item{\code{"noninferiority"}}{Posterior probability
#' \eqn{\Pr(p \ge p_0 - \delta \mid Y)} must exceed
#' \eqn{\gamma_{\mathrm{eq}}}.}
#' \item{\code{"superiority"}}{Posterior probability
#' \eqn{\Pr(p > p_0 + \delta \mid Y)} must exceed
#' \eqn{\gamma_{\mathrm{eq}}}.}
#' }
#'
#' The search proceeds in two steps: (1) find the minimum fixed-sample size
#' \eqn{n^*} at which the one-stage constraints are satisfied; (2) enumerate
#' all two-stage splits \eqn{n_1 + n_2 = n^*} and retain those satisfying the
#' two-stage constraints. The optimal design minimises \eqn{\mathrm{EN}_0}
#' (or \eqn{n^*} under the minimax criterion) among all feasible splits.
#'
#' @param p0 Benchmark response probability.
#' @param delta ROPE half-width (\code{"equivalence"}), non-inferiority margin
#' (\code{"noninferiority"}), or superiority margin (\code{"superiority"}).
#' Must be a single positive number.
#' @param analysis_prior Numeric vector \code{c(a, b)} for the
#' \eqn{\mathrm{Beta}(a, b)} analysis prior on \eqn{p}. Defaults to
#' \code{c(1, 1)} (uniform).
#' @param design_prior_h0 Numeric vector \code{c(a, b)} for the null design
#' prior.
#' @param design_prior_h1 Numeric vector \code{c(a, b)} for the alternative
#' design prior.
#' @param gamma_1 Interim futility threshold in \eqn{(0, 1)} applied to the
#' posterior support for \eqn{H_0}. In the equivalence design, the trial stops
#' early for futility if
#' \eqn{\Pr(p \notin \mathcal{R}_p \mid Y_1) \ge \gamma_1},
#' equivalently if the interim posterior ROPE probability is at most
#' \eqn{1-\gamma_1}. Continuation to stage 2 occurs otherwise.
#' @param gamma_eq Final evidence threshold in \eqn{(0.5, 1)}: the appropriate
#' posterior ROPE probability must exceed \code{gamma_eq} to declare
#' equivalence, non-inferiority, or superiority.
#' @param gamma_diff Threshold for compelling evidence for \eqn{H_0}: the
#' complementary posterior ROPE probability must exceed \code{gamma_diff}.
#' In the two-stage design, compelling evidence for \eqn{H_0} may be obtained
#' either at the interim analysis or at the final analysis. Defaults to
#' \code{gamma_eq}.
#' @param alpha Target predictive type-I error level (upper bound).
#' @param power Target predictive power (lower bound).
#' @param pce Optional lower bound on predictive \code{PCE(H0)}.
#' @param alpha_freq Optional upper bound on the frequentist type-I error.
#' @param power_freq Optional lower bound on the frequentist power.
#' @param p_t1e Point at which the frequentist type-I error is evaluated.
#' @param p_power Point at which the frequentist power is evaluated.
#' @param nmax Upper bound on the fixed-sample size \eqn{n^*} searched in
#' step 1. An informative error is raised if no feasible size is found.
#' @param direction Character string specifying the decision type. One of
#' \code{"equivalence"} (default), \code{"noninferiority"}, or
#' \code{"superiority"}.
#' @param minimax Logical. If \code{TRUE}, minimise \eqn{n} (minimax
#' criterion); if \code{FALSE} (default), minimise \eqn{\mathrm{EN}_0}
#' (optimal criterion).
#' @param progress Logical. If \code{TRUE} (default), print progress messages.
#'
#' @return An object of class \code{"singlearm_rope_twostage_design"}.
#' @export
design_singlearm_twostage_rope <- function(
    p0,
    delta,
    analysis_prior = c(1, 1),
    design_prior_h0,
    design_prior_h1,
    gamma_1 = 0.50,
    gamma_eq = 0.90,
    gamma_diff = gamma_eq,
    alpha = 0.10,
    power = 0.80,
    pce = NULL,
    alpha_freq = NULL,
    power_freq = NULL,
    p_t1e = NULL,
    p_power = NULL,
    nmax = 300L,
    direction = c("equivalence", "noninferiority", "superiority"),
    minimax = FALSE,
    progress = TRUE
) {
  direction <- match.arg(direction)
  
  ## --- local fast helpers --------------------------------------------------
  .bbpmf_loc <- function(y, n, a, b) {
    exp(lchoose(n, y) + lbeta(a + y, b + n - y) - lbeta(a, b))
  }
  
  .post_prob <- switch(
    direction,
    equivalence = function(y, n) {
      lo <- max(0, p0 - delta)
      hi <- min(1, p0 + delta)
      aA <- analysis_prior[1]
      bA <- analysis_prior[2]
      pbeta(hi, aA + y, bA + n - y) - pbeta(lo, aA + y, bA + n - y)
    },
    noninferiority = function(y, n) {
      lo <- max(0, p0 - delta)
      aA <- analysis_prior[1]
      bA <- analysis_prior[2]
      1 - pbeta(lo, aA + y, bA + n - y)
    },
    superiority = function(y, n) {
      hi <- min(1, p0 + delta)
      aA <- analysis_prior[1]
      bA <- analysis_prior[2]
      1 - pbeta(hi, aA + y, bA + n - y)
    }
  )
  
  .freq_success_prob_onestage <- function(n, p_true) {
    y <- 0:n
    post <- .post_prob(y, n)
    sum(dbinom(y[post >= gamma_eq], size = n, prob = p_true))
  }
  
  .freq_success_prob_twostage <- function(n1, n2, p_true) {
    oc <- .evaluate_rope_twostage_oc(
      n1 = n1,
      n2 = n2,
      p0 = p0,
      delta = delta,
      gamma_1 = gamma_1,
      gamma_eq = gamma_eq,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior,
      design_prior_h0 = design_prior_h0,
      design_prior_h1 = design_prior_h1
    )
    
    cont <- oc$cont_region
    if (length(cont) == 0L) return(0)
    
    n <- n1 + n2
    out <- 0
    for (y1 in cont) {
      py1 <- dbinom(y1, size = n1, prob = p_true)
      for (y2 in 0:n2) {
        y <- y1 + y2
        py2 <- dbinom(y2, size = n2, prob = p_true)
        if (.post_prob(y, n) >= gamma_eq) {
          out <- out + py1 * py2
        }
      }
    }
    out
  }
  
  ## --- input validation ----------------------------------------------------
  .validate_probability(p0, "p0")
  
  if (!is.numeric(delta) || length(delta) != 1L ||
      !is.finite(delta) || delta <= 0) {
    stop("'delta' must be a single positive number.", call. = FALSE)
  }
  
  .validate_beta_prior(analysis_prior, "analysis_prior")
  .validate_beta_prior(design_prior_h0, "design_prior_h0")
  .validate_beta_prior(design_prior_h1, "design_prior_h1")
  
  if (!is.numeric(gamma_1) || length(gamma_1) != 1L ||
      gamma_1 <= 0 || gamma_1 >= 1) {
    stop("'gamma_1' must be a single number in (0, 1).", call. = FALSE)
  }
  
  if (!is.numeric(gamma_eq) || length(gamma_eq) != 1L ||
      gamma_eq <= 0.5 || gamma_eq > 1) {
    stop("'gamma_eq' must be a single number in (0.5, 1].", call. = FALSE)
  }
  
  if (!is.numeric(gamma_diff) || length(gamma_diff) != 1L ||
      gamma_diff <= 0.5 || gamma_diff > 1) {
    stop("'gamma_diff' must be a single number in (0.5, 1].", call. = FALSE)
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be in (0, 1).", call. = FALSE)
  }
  
  if (!is.numeric(power) || length(power) != 1L || power <= 0 || power >= 1) {
    stop("'power' must be in (0, 1).", call. = FALSE)
  }
  
  if (!is.null(pce)) {
    if (!is.numeric(pce) || length(pce) != 1L || pce <= 0 || pce >= 1) {
      stop("'pce' must be a single number in (0, 1).", call. = FALSE)
    }
  }
  
  if (!is.null(alpha_freq)) {
    if (!is.numeric(alpha_freq) || length(alpha_freq) != 1L ||
        alpha_freq <= 0 || alpha_freq >= 1) {
      stop("'alpha_freq' must be in (0, 1).", call. = FALSE)
    }
    if (is.null(p_t1e)) {
      stop("'p_t1e' must be supplied when 'alpha_freq' is non-NULL.", call. = FALSE)
    }
    .validate_probability(p_t1e, "p_t1e")
  } else if (!is.null(p_t1e)) {
    warning("'p_t1e' supplied but 'alpha_freq' is NULL; frequentist type-I constraint will not be imposed.")
  }
  
  if (!is.null(power_freq)) {
    if (!is.numeric(power_freq) || length(power_freq) != 1L ||
        power_freq <= 0 || power_freq >= 1) {
      stop("'power_freq' must be in (0, 1).", call. = FALSE)
    }
    if (is.null(p_power)) {
      stop("'p_power' must be supplied when 'power_freq' is non-NULL.", call. = FALSE)
    }
    .validate_probability(p_power, "p_power")
  } else if (!is.null(p_power)) {
    warning("'p_power' supplied but 'power_freq' is NULL; frequentist power constraint will not be imposed.")
  }
  
  nmax <- as.integer(nmax)
  if (nmax < 2L) {
    stop("'nmax' must be an integer >= 2.", call. = FALSE)
  }
  
  ## --- step 1: minimum feasible fixed-sample size n* ----------------------
  if (progress) {
    message("Step 1: searching for minimum feasible fixed-sample size n* ...")
  }
  
  n_star <- NA_integer_
  t1_star <- NA_real_
  pw_star <- NA_real_
  pce_star <- NA_real_
  ft1_star <- NA_real_
  fpw_star <- NA_real_
  
  for (n in 2L:nmax) {
    y <- 0:n
    pmf0 <- .bbpmf_loc(y, n, design_prior_h0[1], design_prior_h0[2])
    pmf1 <- .bbpmf_loc(y, n, design_prior_h1[1], design_prior_h1[2])
    pp <- .post_prob(y, n)
    
    t1 <- sum(pmf0[pp >= gamma_eq])
    pw <- sum(pmf1[pp >= gamma_eq])
    
    pce1 <- if (is.null(pce)) {
      NA_real_
    } else {
      .predictive_pce_onestage(
        n = n,
        p0 = p0,
        delta = delta,
        gamma_diff = gamma_diff,
        analysis_prior = analysis_prior,
        design_prior = design_prior_h0
      )
    }
    
    ft1 <- if (is.null(alpha_freq)) {
      NA_real_
    } else {
      .freq_success_prob_onestage(n = n, p_true = p_t1e)
    }
    
    fpw <- if (is.null(power_freq)) {
      NA_real_
    } else {
      .freq_success_prob_onestage(n = n, p_true = p_power)
    }
    
    cond_pred <- (t1 <= alpha && pw >= power)
    cond_pce  <- is.null(pce) || (pce1 >= pce)
    cond_ft1  <- is.null(alpha_freq) || (ft1 <= alpha_freq)
    cond_fpw  <- is.null(power_freq) || (fpw >= power_freq)
    
    if (cond_pred && cond_pce && cond_ft1 && cond_fpw) {
      n_star <- n
      t1_star <- t1
      pw_star <- pw
      pce_star <- pce1
      ft1_star <- ft1
      fpw_star <- fpw
      break
    }
  }
  
  if (is.na(n_star)) {
    msg <- paste0(
      "No feasible fixed-sample size found within nmax = ", nmax, ". ",
      "Consider increasing 'nmax', relaxing predictive constraints"
    )
    if (!is.null(pce)) {
      msg <- paste0(msg, ", relaxing 'pce'")
    }
    if (!is.null(alpha_freq)) {
      msg <- paste0(msg, ", relaxing 'alpha_freq'")
    }
    if (!is.null(power_freq)) {
      msg <- paste0(msg, ", relaxing 'power_freq'")
    }
    msg <- paste0(msg, ", or adjusting the design priors.")
    stop(msg, call. = FALSE)
  }
  
  if (progress) {
    msg <- sprintf(
      " => n* = %d (pred. type-I = %.4f, pred. power = %.4f",
      n_star, t1_star, pw_star
    )
    if (!is.null(pce)) {
      msg <- paste0(msg, sprintf(", pred. PCE(H0) = %.4f", pce_star))
    }
    if (!is.null(alpha_freq)) {
      msg <- paste0(msg, sprintf(", freq. type-I = %.4f", ft1_star))
    }
    if (!is.null(power_freq)) {
      msg <- paste0(msg, sprintf(", freq. power = %.4f", fpw_star))
    }
    msg <- paste0(msg, ")")
    message(msg)
  }
  
  ## --- step 2: evaluate all splits n1 + n2 = n* ---------------------------
  if (progress) {
    message(sprintf(
      "Step 2: evaluating all %d splits of n* = %d ...",
      n_star - 1L, n_star
    ))
  }
  
  candidates <- vector("list", n_star - 1L)
  k <- 0L
  
  for (n1 in 1L:(n_star - 1L)) {
    n2 <- n_star - n1
    
    oc <- .evaluate_rope_twostage_oc(
      n1 = n1,
      n2 = n2,
      p0 = p0,
      delta = delta,
      gamma_1 = gamma_1,
      gamma_eq = gamma_eq,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior,
      design_prior_h0 = design_prior_h0,
      design_prior_h1 = design_prior_h1
    )
    
    ft1_2st <- if (is.null(alpha_freq)) {
      NA_real_
    } else {
      .freq_success_prob_twostage(n1 = n1, n2 = n2, p_true = p_t1e)
    }
    
    fpw_2st <- if (is.null(power_freq)) {
      NA_real_
    } else {
      .freq_success_prob_twostage(n1 = n1, n2 = n2, p_true = p_power)
    }
    
    cond_pred <- (oc$type1_2st <= alpha && oc$power_2st >= power)
    cond_pce  <- is.null(pce) || (oc$pce_2st >= pce)
    cond_ft1  <- is.null(alpha_freq) || (ft1_2st <= alpha_freq)
    cond_fpw  <- is.null(power_freq) || (fpw_2st >= power_freq)
    
    if (cond_pred && cond_pce && cond_ft1 && cond_fpw) {
      k <- k + 1L
      candidates[[k]] <- c(
        unlist(oc[c(
          "n1", "n2", "n",
          "type1_1st", "power_1st", "pce_1st",
          "type1_2st", "power_2st", "pce_2st",
          "EN0", "EN1"
        )]),
        freq_type1_2st = ft1_2st,
        freq_power_2st = fpw_2st
      )
    }
  }
  
  if (progress) {
    message(sprintf(" => %d feasible two-stage design(s) found.", k))
  }
  
  if (k == 0L) {
    msg <- paste0(
      "No feasible two-stage ROPE design found among splits of n* = ", n_star,
      ". Consider relaxing 'gamma_1', 'alpha', or 'power'"
    )
    if (!is.null(pce)) {
      msg <- paste0(msg, ", or 'pce'")
    }
    if (!is.null(alpha_freq)) {
      msg <- paste0(msg, ", or 'alpha_freq'")
    }
    if (!is.null(power_freq)) {
      msg <- paste0(msg, ", or 'power_freq'")
    }
    msg <- paste0(msg, ".")
    stop(msg, call. = FALSE)
  }
  
  cand_df <- as.data.frame(do.call(rbind, candidates[seq_len(k)]))
  
  if (minimax) {
    cand_df <- cand_df[order(cand_df$n, cand_df$EN0, -cand_df$power_2st), , drop = FALSE]
  } else {
    cand_df <- cand_df[order(cand_df$EN0, cand_df$n, -cand_df$power_2st), , drop = FALSE]
  }
  
  rownames(cand_df) <- NULL
  best <- cand_df[1L, , drop = FALSE]
  
  if (progress) {
    message(sprintf(
      "Done. Optimal design: n1 = %d, n2 = %d, n = %d, EN0 = %.2f",
      as.integer(best$n1), as.integer(best$n2),
      as.integer(best$n), best$EN0
    ))
  }
  
  cont <- .continuation_region_twostage(
    n1 = as.integer(best$n1),
    p0 = p0,
    delta = delta,
    gamma_1 = gamma_1,
    analysis_prior = analysis_prior
  )
  
  structure(
    list(
      call = match.call(),
      p0 = p0,
      delta = delta,
      direction = direction,
      analysis_prior = analysis_prior,
      design_prior_h0 = design_prior_h0,
      design_prior_h1 = design_prior_h1,
      alpha = alpha,
      target_power = power,
      target_pce = pce,
      alpha_freq = alpha_freq,
      target_power_freq = power_freq,
      p_t1e = p_t1e,
      p_power = p_power,
      gamma_1 = gamma_1,
      gamma_eq = gamma_eq,
      gamma_diff = gamma_diff,
      optimality = if (minimax) "minimax" else "optimal",
      design = best,
      continuation_region = cont,
      candidates = cand_df
    ),
    class = "singlearm_rope_twostage_design"
  )
}