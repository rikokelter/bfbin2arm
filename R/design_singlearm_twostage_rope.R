#' Calibrate an optimal single-arm two-stage ROPE design
#'
#' Finds a single-arm two-stage Bayesian design based on the region of
#' practical equivalence (ROPE) for a binary endpoint, with a single interim
#' analysis allowing early stopping for futility. The design covers three
#' decision types via the \code{direction} argument:
#' \describe{
#'   \item{\code{"equivalence"}}{Posterior mass inside the two-sided ROPE
#'     \eqn{[p_0 - \delta,\, p_0 + \delta]} must exceed \eqn{\gamma_{\mathrm{eq}}}.}
#'   \item{\code{"noninferiority"}}{Posterior probability
#'     \eqn{\Pr(p \ge p_0 - \delta \mid Y)} must exceed
#'     \eqn{\gamma_{\mathrm{eq}}}.}
#'   \item{\code{"superiority"}}{Posterior probability
#'     \eqn{\Pr(p > p_0 + \delta \mid Y)} must exceed
#'     \eqn{\gamma_{\mathrm{eq}}}.}
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
#'   (\code{"noninferiority"}), or superiority margin (\code{"superiority"}).
#'   Must be a single positive number.
#' @param analysis_prior Numeric vector \code{c(a, b)} for the
#'   \eqn{\mathrm{Beta}(a, b)} analysis prior on \eqn{p}. Defaults to
#'   \code{c(1, 1)} (uniform).
#' @param design_prior_h0 Numeric vector \code{c(a, b)} for the null design
#'   prior.
#' @param design_prior_h1 Numeric vector \code{c(a, b)} for the alternative
#'   design prior.
#' @param gamma_1 Interim futility threshold in \eqn{(0, 1)}: continuation
#'   requires the interim posterior ROPE probability to exceed \code{gamma_1}.
#' @param gamma_eq Final evidence threshold in \eqn{(0.5, 1)}: the appropriate
#'   posterior ROPE probability must exceed \code{gamma_eq} to declare
#'   equivalence, non-inferiority, or superiority.
#' @param gamma_diff Threshold for compelling evidence for \eqn{H_0}: the
#'   complementary posterior ROPE probability must exceed \code{gamma_diff}.
#'   Defaults to \code{gamma_eq}.
#' @param alpha Target predictive type-I error level (upper bound).
#' @param power Target predictive power (lower bound).
#' @param nmax Upper bound on the fixed-sample size \eqn{n^*} searched in
#'   step 1. An informative error is raised if no feasible size is found.
#' @param direction Character string specifying the decision type. One of
#'   \code{"equivalence"} (default), \code{"noninferiority"}, or
#'   \code{"superiority"}.
#' @param minimax Logical. If \code{TRUE}, minimise \eqn{n} (minimax
#'   criterion); if \code{FALSE} (default), minimise \eqn{\mathrm{EN}_0}
#'   (optimal criterion).
#' @param progress Logical. If \code{TRUE} (default), print progress messages.
#'
#' @return An object of class \code{"singlearm_rope_twostage_design"} with
#'   components:
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{p0}, \code{delta}, \code{direction}}{Design parameters.}
#'     \item{\code{analysis_prior}, \code{design_prior_h0},
#'       \code{design_prior_h1}}{Prior specifications.}
#'     \item{\code{alpha}, \code{target_power}}{Constraint levels.}
#'     \item{\code{gamma_1}, \code{gamma_eq}, \code{gamma_diff}}{Evidence
#'       thresholds.}
#'     \item{\code{optimality}}{Either \code{"optimal"} or \code{"minimax"}.}
#'     \item{\code{design}}{A one-row data frame with the optimal design and
#'       its operating characteristics.}
#'     \item{\code{continuation_region}}{Integer vector of stage-1 response
#'       counts in \eqn{\mathcal{C}_1}.}
#'     \item{\code{candidates}}{Data frame of all feasible two-stage designs,
#'       sorted by the primary optimality criterion.}
#'   }
#'
#' @seealso \code{\link{print.singlearm_rope_twostage_design}},
#'   \code{\link{summary.singlearm_rope_twostage_design}},
#'   \code{\link{plot.singlearm_rope_twostage_design}}
#'
#' @export
design_singlearm_twostage_rope <- function(
    p0,
    delta,
    analysis_prior  = c(1, 1),
    design_prior_h0,
    design_prior_h1,
    gamma_1    = 0.50,
    gamma_eq   = 0.90,
    gamma_diff = gamma_eq,
    alpha      = 0.10,
    power      = 0.80,
    nmax       = 300L,
    direction  = c("equivalence", "noninferiority", "superiority"),
    minimax    = FALSE,
    progress   = TRUE
) {
  direction <- match.arg(direction)
  ## --- local fast helpers (no validation, vectorised) -----------------------
  .bbpmf_loc <- function(y, n, a, b)
    exp(lchoose(n, y) + lbeta(a + y, b + n - y) - lbeta(a, b))
  
  .post_prob <- switch(
    direction,
    equivalence = function(y, n) {
      lo <- max(0, p0 - delta); hi <- min(1, p0 + delta)
      aA <- analysis_prior[1];  bA <- analysis_prior[2]
      pbeta(hi, aA + y, bA + n - y) - pbeta(lo, aA + y, bA + n - y)
    },
    noninferiority = function(y, n) {
      lo <- max(0, p0 - delta)
      aA <- analysis_prior[1]; bA <- analysis_prior[2]
      1 - pbeta(lo, aA + y, bA + n - y)
    },
    superiority = function(y, n) {
      hi <- min(1, p0 + delta)
      aA <- analysis_prior[1]; bA <- analysis_prior[2]
      1 - pbeta(hi, aA + y, bA + n - y)
    }
  )
  
  ## --- input validation -----------------------------------------------------
  .validate_probability(p0, "p0")
  if (!is.numeric(delta) || length(delta) != 1L ||
      !is.finite(delta) || delta <= 0)
    stop("'delta' must be a single positive number.", call. = FALSE)
  .validate_beta_prior(analysis_prior,  "analysis_prior")
  .validate_beta_prior(design_prior_h0, "design_prior_h0")
  .validate_beta_prior(design_prior_h1, "design_prior_h1")
  if (!is.numeric(gamma_1) || length(gamma_1) != 1L ||
      gamma_1 <= 0 || gamma_1 >= 1)
    stop("'gamma_1' must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(gamma_eq) || length(gamma_eq) != 1L ||
      gamma_eq <= 0.5 || gamma_eq > 1)
    stop("'gamma_eq' must be a single number in (0.5, 1].", call. = FALSE)
  if (!is.numeric(gamma_diff) || length(gamma_diff) != 1L ||
      gamma_diff <= 0.5 || gamma_diff > 1)
    stop("'gamma_diff' must be a single number in (0.5, 1].", call. = FALSE)
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be in (0, 1).", call. = FALSE)
  if (!is.numeric(power) || power <= 0 || power >= 1)
    stop("'power' must be in (0, 1).", call. = FALSE)
  nmax <- as.integer(nmax)
  if (nmax < 2L)
    stop("'nmax' must be an integer >= 2.", call. = FALSE)
  
  ## --- direction-specific posterior probability function --------------------
  ## Returns Pr(H1 supported | y, n) for the given direction.
  .post_prob <- switch(
    direction,
    equivalence = function(y, n) {
      lo <- max(0, p0 - delta); hi <- min(1, p0 + delta)
      aA <- analysis_prior[1]; bA <- analysis_prior[2]
      pbeta(hi, aA + y, bA + n - y) - pbeta(lo, aA + y, bA + n - y)
    },
    noninferiority = function(y, n) {
      lo <- max(0, p0 - delta)
      aA <- analysis_prior[1]; bA <- analysis_prior[2]
      1 - pbeta(lo, aA + y, bA + n - y)
    },
    superiority = function(y, n) {
      hi <- min(1, p0 + delta)
      aA <- analysis_prior[1]; bA <- analysis_prior[2]
      1 - pbeta(hi, aA + y, bA + n - y)
    }
  )
  
  ## --- step 1: find minimum feasible fixed-sample size n* ------------------
  if (progress)
    message("Step 1: searching for minimum feasible fixed-sample size n* ...")
  
  n_star  <- NA_integer_
  t1_star <- NA_real_
  pw_star <- NA_real_
  
  for (n in 2L:nmax) {
    y   <- 0:n
    pmf0 <- .bbpmf_loc(y, n, design_prior_h0[1], design_prior_h0[2])
    pmf1 <- .bbpmf_loc(y, n, design_prior_h1[1], design_prior_h1[2])
    pp   <- .post_prob(y, n)
    t1   <- sum(pmf0[pp >= gamma_eq])
    pw   <- sum(pmf1[pp >= gamma_eq])
    if (t1 <= alpha && pw >= power) {
      n_star  <- n
      t1_star <- t1
      pw_star <- pw
      break
    }
  }
  
  if (is.na(n_star))
    stop(
      "No feasible fixed-sample size found within nmax = ", nmax, ". ",
      "Consider increasing 'nmax', relaxing 'alpha' or 'power', or ",
      "adjusting the design priors.",
      call. = FALSE
    )
  
  if (progress)
    message(sprintf(
      "  => n* = %d  (one-stage type-I = %.4f, power = %.4f)",
      n_star, t1_star, pw_star))
  
  ## --- step 2: evaluate all splits n1 + n2 = n* ----------------------------
  if (progress)
    message(sprintf(
      "Step 2: evaluating all %d splits of n* = %d ...",
      n_star - 1L, n_star))
  
  ## helper: continuation region for a given n1
  .cont <- function(n1) {
    y1 <- 0:n1
    y1[.post_prob(y1, n1) > gamma_1]
  }
  
  ## helper: two-stage predictive probability (accept or reject direction)
  .pred_2st <- function(n1, n2, threshold_fn, design_prior) {
    cont <- .cont(n1)
    if (length(cont) == 0L) return(0)
    n       <- n1 + n2
    aD      <- design_prior[1]; bD <- design_prior[2]
    y2_vals <- 0:n2
    out <- 0
    for (y1 in cont) {
      p_y1    <- .bbpmf_loc(y1, n1, aD, bD)
      y_tot   <- y1 + y2_vals
      p_y2    <- .bbpmf_loc(y2_vals, n2, aD + y1, bD + n1 - y1)
      pp      <- .post_prob(y_tot, n)
      out     <- out + p_y1 * sum(p_y2[threshold_fn(pp)])
    }
    out
  }
  
  ## helper: expected sample size
  .en <- function(n1, n2, design_prior) {
    cont   <- .cont(n1)
    p_cont <- if (length(cont) == 0L) 0 else
      sum(.bbpmf_loc(cont, n1, design_prior[1], design_prior[2]))
    n1 + n2 * p_cont
  }
  
  candidates <- vector("list", n_star - 1L)
  k <- 0L
  
  for (n1 in 1L:(n_star - 1L)) {
    n2 <- n_star - n1
    
    t1_2st <- .pred_2st(n1, n2, function(pp) pp >= gamma_eq,  design_prior_h0)
    pw_2st <- .pred_2st(n1, n2, function(pp) pp >= gamma_eq,  design_prior_h1)
    
    if (t1_2st <= alpha && pw_2st >= power) {
      k <- k + 1L
      
      ## one-stage OCs at same n (for reference)
      y    <- 0:(n_star)
      pp_n <- .post_prob(y, n_star)
      pm0  <- .bbpmf_loc(y, n_star, design_prior_h0[1], design_prior_h0[2])
      pm1  <- .bbpmf_loc(y, n_star, design_prior_h1[1], design_prior_h1[2])
      
      pce_2st <- .pred_2st(n1, n2,
                           function(pp) (1 - pp) >= gamma_diff,
                           design_prior_h0)
      
      candidates[[k]] <- c(
        n1        = n1,
        n2        = n2,
        n         = n_star,
        type1_1st = sum(pm0[pp_n >= gamma_eq]),
        power_1st = sum(pm1[pp_n >= gamma_eq]),
        pce_1st   = sum(pm0[(1 - pp_n) >= gamma_diff]),
        type1_2st = t1_2st,
        power_2st = pw_2st,
        pce_2st   = pce_2st,
        EN0       = .en(n1, n2, design_prior_h0),
        EN1       = .en(n1, n2, design_prior_h1)
      )
    }
  }
  
  if (progress)
    message(sprintf("  => %d feasible two-stage design(s) found.", k))
  
  if (k == 0L)
    stop(
      "No feasible two-stage ROPE design found among splits of n* = ",
      n_star, ". Consider relaxing 'gamma_1', 'alpha', or 'power'.",
      call. = FALSE
    )
  
  cand_df <- as.data.frame(do.call(rbind, candidates[seq_len(k)]))
  
  if (minimax) {
    cand_df <- cand_df[order(cand_df$n,   cand_df$EN0,  -cand_df$power_2st), ]
  } else {
    cand_df <- cand_df[order(cand_df$EN0, cand_df$n,    -cand_df$power_2st), ]
  }
  
  best <- cand_df[1L, , drop = FALSE]
  
  if (progress)
    message(sprintf(
      "Done. Optimal design: n1 = %d, n2 = %d, n = %d, EN0 = %.2f",
      as.integer(best$n1), as.integer(best$n2),
      as.integer(best$n),  best$EN0))
  
  cont <- .cont(best$n1)
  
  structure(
    list(
      call             = match.call(),
      p0               = p0,
      delta            = delta,
      direction        = direction,
      analysis_prior   = analysis_prior,
      design_prior_h0  = design_prior_h0,
      design_prior_h1  = design_prior_h1,
      alpha            = alpha,
      target_power     = power,
      gamma_1          = gamma_1,
      gamma_eq         = gamma_eq,
      gamma_diff       = gamma_diff,
      optimality       = if (minimax) "minimax" else "optimal",
      design           = best,
      continuation_region = cont,
      candidates       = cand_df
    ),
    class = "singlearm_rope_twostage_design"
  )
}