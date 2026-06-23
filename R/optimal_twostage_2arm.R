utils::globalVariables(c("n1_1", "n1_2"))
#' Optimal two-stage two-arm Bayes-factor design for binary endpoints
#'
#' Computes an optimal two-stage two-arm Bayes-factor design for binary
#' endpoints, minimizing the expected sample size under the null hypothesis
#' while correcting the operating characteristics for the possibility of
#' early stopping for futility.
#'
#' @inheritParams powertwoarmbinbf01
#' @param alpha Numeric scalar, Bayesian type-I-error target.
#' @param beta Numeric scalar, 1 minus the minimal Bayesian power target.
#' @param k Numeric scalar, efficacy threshold; evidence against the null
#'   hypothesis is declared when the corresponding Bayes factor is smaller
#'   than \code{k}.
#' @param k_f Numeric scalar, futility threshold; compelling evidence for the
#'   null hypothesis is declared when the corresponding Bayes factor is at
#'   least \code{k_f}.
#' @param n1_min Numeric vector of length 2, minimum interim sample sizes
#'   for arms 1 and 2.
#' @param n2_max Numeric vector of length 2, maximum final sample sizes
#'   for arms 1 and 2.
#' @param alloc1,alloc2 Positive numbers, allocation probabilities to arms 1
#'   and 2.
#' @param power_cushion Numeric scalar, optional extra power cushion used in
#'   the fixed-sample search of step 1.
#' @param pceH0 Optional numeric scalar in \code{[0,1]}. If specified,
#'   candidate two-stage designs must satisfy corrected \code{CE_H0 >= pceH0}.
#' @param interim_fraction Numeric vector of length 2 giving lower and upper
#'   bounds for the interim sample size in each arm as a fraction of the fixed
#'   sample size.
#' @param grid_step Positive integer giving the spacing of the interim design
#'   grid.
#' @param coarse_step Positive integer giving the spacing of the coarse
#'   fixed-sample search grid in step 1.
#' @param progress Logical; if \code{TRUE}, prints progress information.
#' @param max_iter Integer, maximum number of total fixed-sample sizes searched
#'   in step 1.
#' @param compute_freq_oc Logical or \code{NULL}. Controls whether frequentist
#'   operating characteristics are computed for candidate two-stage designs
#'   during the search.
#' @param ncores Integer; number of parallel worker processes to use in the
#'   calibration. Defaults to \code{getOption("bfbin2arm.ncores", 1L)}. In
#'   vignettes and examples, a conservative value (e.g. 1 or 2) is recommended
#'   for CRAN checks, whereas users can increase this to exploit all available
#'   cores on their own machines.
#' @param calibration_mode Character string specifying the calibration mode.
#'   Must be one of \code{"Bayesian"}, \code{"frequentist"}, or
#'   \code{"hybrid"}.
#' @param calibration_EN Character string or \code{NULL} specifying whether
#'   the design is ranked by Bayesian or frequentist expected sample size
#'   under the null hypothesis.
#' @param p1_EN_H0,p2_EN_H0 Numeric scalars specifying the null response
#'   probabilities in control and treatment arm used when
#'   \code{calibration_EN = "frequentist"}.
#' @param alpha_freq Numeric scalar, frequentist type-I error target.
#' @param beta_freq Numeric scalar, 1 minus the frequentist power target.
#' @param p1_power,p2_power Numeric scalars specifying the success probabilities
#'   in control and treatment arm used for the frequentist power calculation.
#' @param p_null_grid Optional numeric vector giving the grid of null response
#'   probabilities used for frequentist type-I-error maximization. If
#'   \code{NULL}, a default grid is used.
#' @param a_1_a_Hminus,b_1_a_Hminus Shape parameters of the analysis prior
#'   under the directional null hypothesis H0- for arm 1.
#' @param a_2_a_Hminus,b_2_a_Hminus Shape parameters of the analysis prior
#'   under the directional null hypothesis H0- for arm 2.
#'
#' @return A list with the following components:
#' \item{design}{Four-element integer vector containing the selected two-stage
#' design: interim sample sizes in arms 1 and 2 followed by final sample sizes
#' in arms 1 and 2.}
#' \item{naive_oc}{Named list of uncorrected fixed-sample operating
#' characteristics and fixed-sample sizes found in step 1.}
#' \item{occ}{Named numeric vector of corrected Bayesian operating
#' characteristics for the selected two-stage design.}
#' \item{priors}{List storing design hyperparameters and search settings.}
#' \item{freq_occ}{Named numeric vector with fixed-sample and two-stage
#' frequentist operating characteristics for the final design when frequentist
#' calibration or reporting is active; otherwise \code{NULL}.}
#' \item{conv}{Character string describing the search outcome. Typical values
#' include \code{"converged"}, \code{"no_feasible_fixed"},
#' \code{"no_interim_grid"}, and \code{"no_feasible_design"}. In frequentist
#' or hybrid calibration modes, additional informative status values may be
#' returned when the best available design is returned although all requested
#' constraints were not fully satisfied.}
#'
#' @examples
#' ## Fast Bayesian example with small search space
#' res <- optimal_twostage_2arm_bf(
#'   alpha = 0.10,
#'   beta = 0.20,
#'   k = 1 / 3,
#'   k_f = 3,
#'   n1_min = c(3, 3),
#'   n2_max = c(12, 12),
#'   alloc1 = 0.5,
#'   alloc2 = 0.5,
#'   power_cushion = 0,
#'   pceH0 = NULL,
#'   interim_fraction = c(0.25, 0.75),
#'   grid_step = 2L,
#'   coarse_step = 4L,
#'   progress = FALSE,
#'   max_iter = 24L,
#'   calibration_mode = "Bayesian",
#'   test = "BF01",
#'   a_0_d = 1, b_0_d = 1,
#'   a_0_a = 1, b_0_a = 1,
#'   a_1_d = 1, b_1_d = 1,
#'   a_2_d = 1, b_2_d = 1,
#'   a_1_a = 1, b_1_a = 1,
#'   a_2_a = 1, b_2_a = 1
#' )
#' res$design
#' res$occ
#'
#' \donttest{
#' res2 <- optimal_twostage_2arm_bf(
#'   alpha = 0.05,
#'   beta = 0.20,
#'   k = 1 / 3,
#'   k_f = 3,
#'   n1_min = c(5, 5),
#'   n2_max = c(20, 20),
#'   alloc1 = 0.5,
#'   alloc2 = 0.5,
#'   power_cushion = 0.02,
#'   pceH0 = 0.50,
#'   interim_fraction = c(0.25, 0.75),
#'   grid_step = 1L,
#'   coarse_step = 4L,
#'   progress = FALSE,
#'   max_iter = 40L,
#'   calibration_mode = "Bayesian",
#'   test = "BF+0",
#'   a_0_d = 1, b_0_d = 1,
#'   a_0_a = 1, b_0_a = 1,
#'   a_1_d = 1, b_1_d = 2,
#'   a_2_d = 2, b_2_d = 1,
#'   a_1_a = 1, b_1_a = 1,
#'   a_2_a = 1, b_2_a = 1,
#'   a_1_d_Hminus = 1, b_1_d_Hminus = 1,
#'   a_2_d_Hminus = 1, b_2_d_Hminus = 1,
#'   a_1_a_Hminus = 1, b_1_a_Hminus = 1,
#'   a_2_a_Hminus = 1, b_2_a_Hminus = 1
#' )
#' res2$design
#' res2$occ
#' }
#'
#' @import parallel
#' @importFrom utils flush.console
#' @export
optimal_twostage_2arm_bf <- function(
    alpha = 0.05,
    beta = 0.20,
    k = 1 / 3,
    k_f = 3,
    n1_min = c(5, 5),
    n2_max = c(50, 50),
    alloc1 = 0.5,
    alloc2 = 0.5,
    power_cushion = 0.0,
    pceH0 = NULL,
    interim_fraction = c(0, 1),
    grid_step = 1L,
    coarse_step = 10L,
    progress = TRUE,
    max_iter = 10000L,
    ncores = getOption("bfbin2arm.ncores", 1L),
    compute_freq_oc = NULL,
    calibration_mode = c("Bayesian", "frequentist", "hybrid"),
    calibration_EN = NULL,
    p1_EN_H0 = NULL,
    p2_EN_H0 = NULL,
    alpha_freq = alpha,
    beta_freq = beta,
    p1_power = NULL,
    p2_power = NULL,
    p_null_grid = NULL,
    
    test = "BF01", # one of c("BF01", "BF+0", "BF-0", "BF+-")
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
) {
  
  calibration_mode <- match.arg(calibration_mode)
  test <- match.arg(test, c("BF01", "BF+0", "BF-0", "BF+-"))
  
  # Resolve number of cores
  if (is.null(ncores) || is.na(ncores) || ncores < 1L) {
    ncores <- 1L
  } else {
    ncores <- as.integer(ncores)
  }
  
  ## Default EN criterion depends on calibration mode:
  ## - Bayesian -> Bayesian EN
  ## - frequentist -> frequentist EN
  ## - hybrid -> frequentist EN by default
  if (is.null(calibration_EN)) {
    calibration_EN <- switch(
      calibration_mode,
      Bayesian = "Bayesian",
      frequentist = "frequentist",
      hybrid = "frequentist"
    )
  } else {
    calibration_EN <- match.arg(calibration_EN, c("Bayesian", "frequentist"))
  }
  
  ## Default frequentist OC computation:
  ## - Bayesian: off by default, unless explicitly requested
  ## - frequentist / hybrid: on by default
  if (is.null(compute_freq_oc)) {
    compute_freq_oc <- calibration_mode %in% c("frequentist", "hybrid")
  }
  
  if (calibration_EN == "frequentist") {
    if (is.null(p1_EN_H0) || is.null(p2_EN_H0)) {
      stop("For calibration_EN='frequentist', please specify p1_EN_H0 and p2_EN_H0.")
    }
  }
  
  ## Internal: Bayesian corrected two-stage OCs
  compute_corrected_twostage_oc_2arm <- function(
    n1_1, n1_2,
    n2_1, n2_2,
    k, k_f, test,
    a_0_d, b_0_d,
    a_0_a, b_0_a,
    a_1_d, b_1_d, a_2_d, b_2_d,
    a_1_a, b_1_a, a_2_a, b_2_a,
    a_1_d_Hminus, b_1_d_Hminus, a_2_d_Hminus, b_2_d_Hminus,
    a_1_a_Hminus, b_1_a_Hminus, a_2_a_Hminus, b_2_a_Hminus
  ) {
    stopifnot(n1_1 <= n2_1, n1_2 <= n2_2)
    stopifnot(n2_1 - n1_1 >= 0, n2_2 - n1_2 >= 0)
    
    BFmat_interim <- matrix(0L, n1_1 + 1L, n1_2 + 1L)
    
    f1_0    <- matrix(0, n1_1 + 1L, n1_2 + 1L)
    f1_1    <- matrix(0, n1_1 + 1L, n1_2 + 1L)
    f1_plus <- matrix(0, n1_1 + 1L, n1_2 + 1L)
    f1_minus<- matrix(0, n1_1 + 1L, n1_2 + 1L)
    
    f2_0    <- matrix(0, n2_1 - n1_1 + 1L, n2_2 - n1_2 + 1L)
    f2_1    <- matrix(0, n2_1 - n1_1 + 1L, n2_2 - n1_2 + 1L)
    f2_plus <- matrix(0, n2_1 - n1_1 + 1L, n2_2 - n1_2 + 1L)
    f2_minus<- matrix(0, n2_1 - n1_1 + 1L, n2_2 - n1_2 + 1L)
    
    E2_mat     <- matrix(FALSE, n2_1 + 1L, n2_2 + 1L)
    CE2_H0_mat <- matrix(FALSE, n2_1 + 1L, n2_2 + 1L)
    
    ## Interim stage
    for (i1 in 1:(n1_1 + 1L)) {
      x1 <- i1 - 1L
      for (j1 in 1:(n1_2 + 1L)) {
        x2 <- j1 - 1L
        
        if (test == "BF01") {
          bf01_int <- twoarmbinbf01(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_0_a = a_0_a, b_0_a = b_0_a,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          BFmat_interim[i1, j1] <- as.integer(bf01_int >= k_f)
          
          f1_0[i1, j1] <- predictiveDensityH0(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_0_d = a_0_d, b_0_d = b_0_d
          )
          f1_1[i1, j1] <- predictiveDensityH1(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_d = a_1_d, b_1_d = b_1_d,
            a_2_d = a_2_d, b_2_d = b_2_d
          )
          
        } else if (test == "BF+0") {
          ## NOTE: this is the version from your current file; if you have
          ## switched to the direct BF+0 function, replace this block by
          ## the direct-BF implementation you already tested.
          bf01_int <- twoarmbinbf01(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_0_a = a_0_a, b_0_a = b_0_a,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          bf_plus_1_int <- BFplus1(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          bf_plus_0_int <- BFplus0(BFplus1 = bf_plus_1_int, BF01 = bf01_int)
          BFmat_interim[i1, j1] <- as.integer(bf_plus_0_int >= k_f)
          
          f1_0[i1, j1] <- predictiveDensityH0(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_0_d = a_0_d, b_0_d = b_0_d
          )
          f1_plus[i1, j1] <- predictiveDensityHplus_trunc(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_d = a_1_d, b_1_d = b_1_d,
            a_2_d = a_2_d, b_2_d = b_2_d
          )
          
        } else if (test == "BF-0") {
          bf01_int <- twoarmbinbf01(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_0_a = a_0_a, b_0_a = b_0_a,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          bf_minus_1_int <- BFminus1(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
            a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
          )
          bf_minus_0_int <- BFminus0(BFminus1 = bf_minus_1_int, BF01 = bf01_int)
          BFmat_interim[i1, j1] <- as.integer(bf_minus_0_int >= k_f)
          
          f1_0[i1, j1] <- predictiveDensityH0(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_0_d = a_0_d, b_0_d = b_0_d
          )
          f1_minus[i1, j1] <- predictiveDensityHminus_trunc(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_d = a_1_d_Hminus, b_1_d = b_1_d_Hminus,
            a_2_d = a_2_d_Hminus, b_2_d = b_2_d_Hminus
          )
          
        } else if (test == "BF+-") {
          bf_plus_1_int <- BFplus1(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          bf_minus_1_int <- BFminus1(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
            a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
          )
          bf_pm_int <- BFplusMinus(BFplus1 = bf_plus_1_int, BFminus1 = bf_minus_1_int)
          BFmat_interim[i1, j1] <- as.integer(bf_pm_int >= k_f)
          
          f1_0[i1, j1] <- predictiveDensityH0(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_0_d = a_0_d, b_0_d = b_0_d
          )
          f1_plus[i1, j1] <- predictiveDensityHplus_trunc(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_d = a_1_d, b_1_d = b_1_d,
            a_2_d = a_2_d, b_2_d = b_2_d
          )
          f1_minus[i1, j1] <- predictiveDensityHminus_trunc(
            y1 = x1, y2 = x2, n1 = n1_1, n2 = n1_2,
            a_1_d = a_1_d_Hminus, b_1_d = b_1_d_Hminus,
            a_2_d = a_2_d_Hminus, b_2_d = b_2_d_Hminus
          )
        }
      }
    }
    
    ## Final-stage decision regions
    for (i2 in 1:(n2_1 + 1L)) {
      y1 <- i2 - 1L
      for (j2 in 1:(n2_2 + 1L)) {
        y2 <- j2 - 1L
        
        if (test == "BF01") {
          bf01_fin <- twoarmbinbf01(
            y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
            a_0_a = a_0_a, b_0_a = b_0_a,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          E2_mat[i2, j2] <- (bf01_fin < k)
          CE2_H0_mat[i2, j2] <- (bf01_fin >= k_f)
          
        } else if (test == "BF+0") {
          bf01_fin <- twoarmbinbf01(
            y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
            a_0_a = a_0_a, b_0_a = b_0_a,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          bf_plus_1_fin <- BFplus1(
            y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          bf_plus_0_fin <- BFplus0(BFplus1 = bf_plus_1_fin, BF01 = bf01_fin)
          E2_mat[i2, j2] <- (bf_plus_0_fin < k)
          CE2_H0_mat[i2, j2] <- (bf_plus_0_fin >= k_f)
          
        } else if (test == "BF-0") {
          bf01_fin <- twoarmbinbf01(
            y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
            a_0_a = a_0_a, b_0_a = b_0_a,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          bf_minus_1_fin <- BFminus1(
            y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
            a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
            a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
          )
          bf_minus_0_fin <- BFminus0(BFminus1 = bf_minus_1_fin, BF01 = bf01_fin)
          E2_mat[i2, j2] <- (bf_minus_0_fin < k)
          CE2_H0_mat[i2, j2] <- (bf_minus_0_fin >= k_f)
          
        } else if (test == "BF+-") {
          bf_plus_1_fin <- BFplus1(
            y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a
          )
          bf_minus_1_fin <- BFminus1(
            y1 = y1, y2 = y2, n1 = n2_1, n2 = n2_2,
            a_1_a = a_1_a_Hminus, b_1_a = b_1_a_Hminus,
            a_2_a = a_2_a_Hminus, b_2_a = b_2_a_Hminus
          )
          bf_pm_fin <- BFplusMinus(BFplus1 = bf_plus_1_fin, BFminus1 = bf_minus_1_fin)
          E2_mat[i2, j2] <- (bf_pm_fin < k)
          CE2_H0_mat[i2, j2] <- (bf_pm_fin >= k_f)
        }
      }
    }
    
    ## Predictive densities for second stage
    for (i in 1:(n2_1 - n1_1 + 1L)) {
      z1 <- i - 1L
      for (j in 1:(n2_2 - n1_2 + 1L)) {
        z2 <- j - 1L
        
        if (test == "BF01") {
          f2_0[i, j] <- predictiveDensityH0(
            y1 = z1, y2 = z2,
            n1 = n2_1 - n1_1, n2 = n2_2 - n1_2,
            a_0_d = a_0_d, b_0_d = b_0_d
          )
          f2_1[i, j] <- predictiveDensityH1(
            y1 = z1, y2 = z2,
            n1 = n2_1 - n1_1, n2 = n2_2 - n1_2,
            a_1_d = a_1_d, b_1_d = b_1_d,
            a_2_d = a_2_d, b_2_d = b_2_d
          )
          
        } else if (test == "BF+0") {
          f2_0[i, j] <- predictiveDensityH0(
            y1 = z1, y2 = z2,
            n1 = n2_1 - n1_1, n2 = n2_2 - n1_2,
            a_0_d = a_0_d, b_0_d = b_0_d
          )
          f2_plus[i, j] <- predictiveDensityHplus_trunc(
            y1 = z1, y2 = z2,
            n1 = n2_1 - n1_1, n2 = n2_2 - n1_2,
            a_1_d = a_1_d, b_1_d = b_1_d,
            a_2_d = a_2_d, b_2_d = b_2_d
          )
          
        } else if (test == "BF-0") {
          f2_0[i, j] <- predictiveDensityH0(
            y1 = z1, y2 = z2,
            n1 = n2_1 - n1_1, n2 = n2_2 - n1_2,
            a_0_d = a_0_d, b_0_d = b_0_d
          )
          f2_minus[i, j] <- predictiveDensityHminus_trunc(
            y1 = z1, y2 = z2,
            n1 = n2_1 - n1_1, n2 = n2_2 - n1_2,
            a_1_d = a_1_d_Hminus, b_1_d = b_1_d_Hminus,
            a_2_d = a_2_d_Hminus, b_2_d = b_2_d_Hminus
          )
          
        } else if (test == "BF+-") {
          f2_plus[i, j] <- predictiveDensityHplus_trunc(
            y1 = z1, y2 = z2,
            n1 = n2_1 - n1_1, n2 = n2_2 - n1_2,
            a_1_d = a_1_d, b_1_d = b_1_d,
            a_2_d = a_2_d, b_2_d = b_2_d
          )
          f2_minus[i, j] <- predictiveDensityHminus_trunc(
            y1 = z1, y2 = z2,
            n1 = n2_1 - n1_1, n2 = n2_2 - n1_2,
            a_1_d = a_1_d_Hminus, b_1_d = b_1_d_Hminus,
            a_2_d = a_2_d_Hminus, b_2_d = b_2_d_Hminus
          )
        }
      }
    }
    
    ## Corrections
    Delta0 <- 0
    Delta1 <- 0
    DeltaCE0 <- 0
    
    for (i1 in 1:(n1_1 + 1L)) {
      x1 <- i1 - 1L
      for (j1 in 1:(n1_2 + 1L)) {
        x2 <- j1 - 1L
        
        if (BFmat_interim[i1, j1] == 1L) {
          p_final_ce_h0_given_x <- 0
          
          for (i2 in 1:(n2_1 - n1_1 + 1L)) {
            z1 <- i2 - 1L
            for (j2 in 1:(n2_2 - n1_2 + 1L)) {
              z2 <- j2 - 1L
              y1 <- x1 + z1
              y2 <- x2 + z2
              
              if (y1 > n2_1 || y2 > n2_2) next
              
              i_y <- y1 + 1L
              j_y <- y2 + 1L
              
              if (E2_mat[i_y, j_y]) {
                if (test == "BF01") {
                  Delta0 <- Delta0 + f1_0[i1, j1] * f2_0[i2, j2]
                  Delta1 <- Delta1 + f1_1[i1, j1] * f2_1[i2, j2]
                } else if (test == "BF+0") {
                  Delta0 <- Delta0 + f1_0[i1, j1] * f2_0[i2, j2]
                  Delta1 <- Delta1 + f1_plus[i1, j1] * f2_plus[i2, j2]
                } else if (test == "BF-0") {
                  Delta0 <- Delta0 + f1_0[i1, j1] * f2_0[i2, j2]
                  Delta1 <- Delta1 + f1_minus[i1, j1] * f2_minus[i2, j2]
                } else if (test == "BF+-") {
                  Delta0 <- Delta0 + f1_0[i1, j1] * f2_plus[i2, j2]
                  Delta1 <- Delta1 + f1_plus[i1, j1] * f2_plus[i2, j2]
                }
              }
              
              if (CE2_H0_mat[i_y, j_y]) {
                p_final_ce_h0_given_x <- p_final_ce_h0_given_x + f2_0[i2, j2]
              }
            }
          }
          
          DeltaCE0 <- DeltaCE0 + f1_0[i1, j1] * (1 - p_final_ce_h0_given_x)
        }
      }
    }
    
    ## Naive (fixed-sample) Bayesian OC at n2 and corrections
    oc <- powertwoarmbinbf01(
      n1 = n2_1, n2 = n2_2, k = k, k_f = k_f, test = test,
      a_0_d = a_0_d, b_0_d = b_0_d,
      a_0_a = a_0_a, b_0_a = b_0_a,
      a_1_d = a_1_d, b_1_d = b_1_d,
      a_2_d = a_2_d, b_2_d = b_2_d,
      a_1_a = a_1_a, b_1_a = b_1_a,
      a_2_a = a_2_a, b_2_a = b_2_a,
      a_1_d_Hminus = a_1_d_Hminus, b_1_d_Hminus = b_1_d_Hminus,
      a_2_d_Hminus = a_2_d_Hminus, b_2_d_Hminus = b_2_d_Hminus,
      output = "numeric"
    )
    
    naive_power <- unname(oc["Power"])
    naive_t1e   <- unname(oc["Type1_Error"])
    naive_pceH0 <- unname(oc["CE_H0"])
    
    # Corrections:
    # - Power, Type1_Error: subtract trajectories removed by interim futility
    # - CE_H0: add futility trajectories that would not be CE(H0) at final
    corr_power <- naive_power - Delta1
    corr_t1e   <- naive_t1e   - Delta0
    corr_pceH0 <- naive_pceH0 + DeltaCE0
    
    # Numerical safety: clamp to [0, 1]
    clamp01 <- function(x) pmin(pmax(x, 0), 1)
    
    corr_power <- clamp01(corr_power)
    corr_t1e   <- clamp01(corr_t1e)
    corr_pceH0 <- clamp01(corr_pceH0)
    
    futility_prob <- sum(f1_0 * (BFmat_interim == 1L))
    N1 <- n1_1 + n1_2
    N2 <- n2_1 + n2_2
    E_H0_N <- N1 * futility_prob + N2 * (1 - futility_prob)
    
    c(
      Power         = corr_power,
      Type1_Error   = corr_t1e,
      CE_H0         = corr_pceH0,
      futility_prob = futility_prob,
      E_H0_N        = E_H0_N
    )
  }
  
  ## Internal: unified feasibility check used both in Step 2 and post hoc
  check_design_feasibility <- function(
    bayes,
    freq_vec = NULL,
    alpha,
    beta,
    pceH0 = NULL,
    alpha_freq = NA_real_,
    beta_freq = NA_real_
  ) {
    bayes_ok <- TRUE
    
    if (!is.null(bayes)) {
      if (!is.na(alpha) && unname(bayes["Type1_Error"]) > alpha) {
        bayes_ok <- FALSE
      }
      if (!is.na(beta) && unname(bayes["Power"]) < (1 - beta)) {
        bayes_ok <- FALSE
      }
      if (!is.null(pceH0) && unname(bayes["CE_H0"]) < pceH0) {
        bayes_ok <- FALSE
      }
    }
    
    freq_ok <- TRUE
    
    if (!is.null(freq_vec)) {
      t1e_two   <- unname(freq_vec["Type1_Error_freq_two_stage"])
      power_two <- unname(freq_vec["Power_freq_two_stage"])
      
      if (!is.na(alpha_freq) && t1e_two > alpha_freq) {
        freq_ok <- FALSE
      }
      if (!is.na(beta_freq) && power_two < (1 - beta_freq)) {
        freq_ok <- FALSE
      }
    }
    
    list(
      bayes_ok = bayes_ok,
      freq_ok = freq_ok,
      hybrid_ok = bayes_ok && freq_ok
    )
  }
  
  ## Internal: frequentist EN under a single null point
  freq_EN_twoarm_null <- function(
    n1_1, n1_2,
    n2_1, n2_2,
    p1_null, p2_null,
    k, k_f, test,
    a_0_a, b_0_a,
    a_1_a, b_1_a,
    a_2_a, b_2_a,
    a_1_a_Hminus, b_1_a_Hminus,
    a_2_a_Hminus, b_2_a_Hminus
  ) {
    freq2 <- freq_oc_twostage_twoarm_fixed(
      n1_1 = n1_1, n1_2 = n1_2,
      n2_1 = n2_1, n2_2 = n2_2,
      k = k, k_f = k_f, test = test,
      p1 = p1_null, p2 = p2_null,
      a_0_a = a_0_a, b_0_a = b_0_a,
      a_1_a = a_1_a, b_1_a = b_1_a,
      a_2_a = a_2_a, b_2_a = b_2_a,
      a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
      a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
    )
    
    fut_prob <- unname(freq2["Stop_for_futility_prob"])
    if (is.na(fut_prob)) fut_prob <- 0
    
    N1 <- n1_1 + n1_2
    N2 <- n2_1 + n2_2
    N1 * fut_prob + N2 * (1 - fut_prob)
  }
  
  #####################
  #### Fixed-sample helper (unchanged)
  #####################
  
  find_fixed_sample_size <- function(
    calibration_mode,
    alpha, beta, k, k_f, test,
    n1_min, n2_max, alloc1, alloc2,
    power_cushion, pceH0,
    alpha_freq, beta_freq,
    p1_power, p2_power, p_null_grid,
    progress,
    max_iter,
    coarse_step,
    a_0_d, b_0_d,
    a_0_a, b_0_a,
    a_1_d, b_1_d, a_2_d, b_2_d,
    a_1_a, b_1_a, a_2_a, b_2_a,
    a_1_d_Hminus, b_1_d_Hminus,
    a_2_d_Hminus, b_2_d_Hminus
  ) {
    coarse_step <- coarse_step
    target_power <- (1 - beta) + power_cushion
    
    eval_fixed_candidate <- function(n_tot) {
      n_fix_1 <- as.integer(round(n_tot * alloc1))
      n_fix_2 <- as.integer(n_tot - n_fix_1)
      
      if (n_fix_1 < n1_min[1] || n_fix_2 < n1_min[2]) {
        return(NULL)
      }
      if (n_fix_1 > n2_max[1] || n_fix_2 > n2_max[2]) {
        return(NULL)
      }
      
      oc <- powertwoarmbinbf01(
        n1 = n_fix_1, n2 = n_fix_2, k = k, k_f = k_f, test = test,
        a_0_d = a_0_d, b_0_d = b_0_d,
        a_0_a = a_0_a, b_0_a = b_0_a,
        a_1_d = a_1_d, b_1_d = b_1_d,
        a_2_d = a_2_d, b_2_d = b_2_d,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a,
        a_1_d_Hminus = a_1_d_Hminus, b_1_d_Hminus = b_1_d_Hminus,
        a_2_d_Hminus = a_2_d_Hminus, b_2_d_Hminus = b_2_d_Hminus,
        output = "numeric"
      )
      
      if (is.na(oc["Power"]) || is.na(oc["Type1_Error"])) {
        return(NULL)
      }
      
      meets_bayes <- (oc["Power"] >= target_power) &&
        (oc["Type1_Error"] <= alpha)
      
      if (!is.null(pceH0)) {
        meets_bayes <- meets_bayes && (oc["CE_H0"] >= pceH0)
      }
      
      need_freq_step1 <- calibration_mode %in% c("frequentist", "hybrid")
      
      freq_power <- NA_real_
      freq_t1e <- NA_real_
      meets_freq <- NA
      
      if (need_freq_step1) {
        freq_power <- freq_oc_twoarm_fixed(
          n1 = n_fix_1, n2 = n_fix_2,
          k = k, k_f = k_f, test = test,
          p1_null = p1_power, p2_null = p2_power,
          p1_alt  = p1_power, p2_alt  = p2_power,
          a_0_a = a_0_a, b_0_a = b_0_a,
          a_1_a = a_1_a, b_1_a = b_1_a,
          a_2_a = a_2_a, b_2_a = b_2_a,
          a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
          a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
        )["Power_freq"]
        
        freq_t1e <- freq_t1e_sup_fixed(
          n1 = n_fix_1, n2 = n_fix_2,
          k = k, k_f = k_f, test = test,
          p_null_grid = p_null_grid,
          a_0_a = a_0_a, b_0_a = b_0_a,
          a_1_a = a_1_a, b_1_a = b_1_a,
          a_2_a = a_2_a, b_2_a = b_2_a,
          a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
          a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
        )
        
        meets_freq <- (!is.na(freq_power) && freq_power >= target_power) &&
          (!is.na(freq_t1e) && freq_t1e <= alpha_freq)
      }
      
      feasible <- switch(
        calibration_mode,
        Bayesian    = meets_bayes,
        frequentist = isTRUE(meets_freq),
        hybrid      = meets_bayes && isTRUE(meets_freq)
      )
      
      list(
        n_tot = n_tot,
        n1 = n_fix_1,
        n2 = n_fix_2,
        bayes_oc = oc,
        meets_bayes = meets_bayes,
        freq_power = freq_power,
        freq_t1e = freq_t1e,
        meets_freq = meets_freq,
        feasible = feasible
      )
    }
    
    ## Balanced total sample sizes only if alloc1 == alloc2
    min_tot <- sum(n1_min)
    max_tot <- min(sum(n2_max), max_iter)
    
    if (abs(alloc1 - alloc2) < .Machine$double.eps^0.5) {
      coarse_seq <- seq.int(min_tot, max_tot, by = coarse_step)
      coarse_seq <- unique(as.integer(2L * round(coarse_seq / 2L)))
    } else {
      coarse_seq <- seq.int(min_tot, max_tot, by = coarse_step)
    }
    
    if (tail(coarse_seq, 1) != max_tot) {
      coarse_seq <- unique(c(coarse_seq, max_tot))
    }
    
    if (progress) {
      cat("Step 1: coarse fixed-sample search...\n")
      utils::flush.console()
    }
    
    coarse_res <- vector("list", length(coarse_seq))
    first_feasible_idx <- NA_integer_
    
    for (i in seq_along(coarse_seq)) {
      res <- eval_fixed_candidate(coarse_seq[i])
      coarse_res[[i]] <- res
      
      if (!is.null(res) && progress) {
        msg <- sprintf(
          " Coarse grid[%3d]: n_tot=%3d | n1=%3d n2=%3d | Bayes Power=%.3f | Bayes T1E=%.3f | PCE(H0)=%.3f",
          i, res$n_tot, res$n1, res$n2,
          res$bayes_oc["Power"], res$bayes_oc["Type1_Error"], res$bayes_oc["CE_H0"]
        )
        if (calibration_mode %in% c("frequentist", "hybrid")) {
          msg <- paste0(
            msg,
            sprintf(" | Freq Power=%.3f | Freq supT1E=%.3f",
                    res$freq_power, res$freq_t1e)
          )
        }
        cat(msg, "\n")
        utils::flush.console()
      }
      
      if (!is.null(res) && isTRUE(res$feasible)) {
        first_feasible_idx <- i
        break
      }
    }
    
    if (is.na(first_feasible_idx)) {
      return(list(found = FALSE, best = NULL))
    }
    
    ## Refine interval
    left_tot <- if (first_feasible_idx == 1L) min_tot else coarse_seq[first_feasible_idx - 1L]
    right_tot <- coarse_seq[first_feasible_idx]
    
    if (progress) {
      cat(sprintf(
        "Refining fixed-sample search on [%d, %d]...\n",
        left_tot, right_tot
      ))
      utils::flush.console()
    }
    
    if (abs(alloc1 - alloc2) < .Machine$double.eps^0.5) {
      refine_seq <- seq.int(left_tot, right_tot, by = 2L)
      refine_seq <- unique(as.integer(2L * round(refine_seq / 2L)))
    } else {
      refine_seq <- seq.int(left_tot, right_tot, by = 1L)
    }
    
    best_res <- NULL
    
    for (n_tot in refine_seq) {
      res <- eval_fixed_candidate(n_tot)
      
      if (!is.null(res) && progress) {
        msg <- sprintf(
          " Refine n_tot=%3d | n1=%3d n2=%3d | Bayes Power=%.3f | Bayes T1E=%.3f | PCE(H0)=%.3f",
          res$n_tot, res$n1, res$n2,
          res$bayes_oc["Power"], res$bayes_oc["Type1_Error"], res$bayes_oc["CE_H0"]
        )
        if (calibration_mode %in% c("frequentist", "hybrid")) {
          msg <- paste0(
            msg,
            sprintf(" | Freq Power=%.3f | Freq supT1E=%.3f",
                    res$freq_power, res$freq_t1e)
          )
        }
        cat(msg, "\n")
        utils::flush.console()
      }
      
      if (!is.null(res) && isTRUE(res$feasible)) {
        best_res <- res
        break
      }
    }
    
    if (is.null(best_res)) {
      return(list(found = FALSE, best = NULL))
    }
    
    list(found = TRUE, best = best_res)
  }
  
  ## ---------------------------- Main body -------------------------------------
  
  test <- match.arg(test, c("BF01", "BF+0", "BF-0", "BF+-"))
  
  stopifnot(
    is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1,
    is.numeric(beta), length(beta) == 1L, beta > 0, beta < 1,
    is.numeric(power_cushion), length(power_cushion) == 1L, power_cushion >= 0,
    length(n1_min) == 2L,
    length(n2_max) == 2L,
    length(interim_fraction) == 2L,
    is.numeric(grid_step), length(grid_step) == 1L, grid_step >= 1
  )
  
  n1_min <- as.integer(n1_min)
  n2_max <- as.integer(n2_max)
  grid_step <- as.integer(grid_step)
  interim_fraction <- as.numeric(interim_fraction)
  max_iter <- as.integer(max_iter)
  
  stopifnot(
    all(n1_min >= 1L),
    all(n2_max >= n1_min),
    interim_fraction[1] >= 0,
    interim_fraction[2] <= 1,
    interim_fraction[1] <= interim_fraction[2]
  )
  
  if (!is.null(pceH0)) {
    stopifnot(is.numeric(pceH0), length(pceH0) == 1L, pceH0 >= 0, pceH0 <= 1)
  }
  
  need_freq_inputs <-
    calibration_mode %in% c("frequentist", "hybrid") || isTRUE(compute_freq_oc)
  
  if (need_freq_inputs) {
    if (is.null(p1_power) || is.null(p2_power)) {
      stop("When calibration_mode is 'frequentist' or 'hybrid', or when ",
           "'compute_freq_oc' is TRUE, both 'p1_power' and 'p2_power' must be specified.")
    }
  }
  
  if (!is.null(p_null_grid)) {
    stopifnot(is.numeric(p_null_grid), all(p_null_grid > 0), all(p_null_grid < 1))
  }
  
  alloc_sum <- alloc1 + alloc2
  if (alloc1 <= 0 || alloc2 <= 0) {
    stop("alloc1 and alloc2 must be positive.")
  }
  
  alloc1 <- alloc1 / alloc_sum
  alloc2 <- alloc2 / alloc_sum
  
  priors <- list(
    test = test,
    alpha = alpha,
    beta = beta,
    k = k,
    k_f = k_f,
    power_cushion = power_cushion,
    pceH0 = pceH0,
    interim_fraction = interim_fraction,
    grid_step = grid_step,
    n1_min = n1_min,
    n2_max = n2_max,
    alloc1 = alloc1,
    alloc2 = alloc2,
    compute_freq_oc = compute_freq_oc,
    alpha_freq = alpha_freq,
    beta_freq = beta_freq,
    p1_power = p1_power,
    p2_power = p2_power,
    p_null_grid = p_null_grid,
    calibration_mode = calibration_mode,
    calibration_EN   = calibration_EN,
    p1_EN_H0 = p1_EN_H0,
    p2_EN_H0 = p2_EN_H0,
    a_0_d = a_0_d, b_0_d = b_0_d,
    a_0_a = a_0_a, b_0_a = b_0_a,
    a_1_d = a_1_d, b_1_d = b_1_d,
    a_2_d = a_2_d, b_2_d = b_2_d,
    a_1_a = a_1_a, b_1_a = b_1_a,
    a_2_a = a_2_a, b_2_a = b_2_a,
    a_1_d_Hminus = a_1_d_Hminus, b_1_d_Hminus = b_1_d_Hminus,
    a_2_d_Hminus = a_2_d_Hminus, b_2_d_Hminus = b_2_d_Hminus,
    a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
    a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
  )
  
  ## Step 1: fixed-sample search
  
  if (progress) {
    cat(
      "Step 1: searching for fixed-sample sufficiency (alpha=",
      alpha, ", beta=", beta, ", cushion=", power_cushion, ")...\n", sep = ""
    )
    utils::flush.console()
  }
  
  fixed_search <- find_fixed_sample_size(
    calibration_mode = calibration_mode,
    alpha = alpha,
    beta = beta,
    k = k,
    k_f = k_f,
    test = test,
    n1_min = n1_min,
    n2_max = n2_max,
    alloc1 = alloc1,
    alloc2 = alloc2,
    power_cushion = power_cushion,
    pceH0 = pceH0,
    alpha_freq = alpha_freq,
    beta_freq = beta_freq,
    p1_power = p1_power,
    p2_power = p2_power,
    p_null_grid = p_null_grid,
    progress = progress,
    max_iter = max_iter,
    coarse_step = coarse_step,
    a_0_d = a_0_d, b_0_d = b_0_d,
    a_0_a = a_0_a, b_0_a = b_0_a,
    a_1_d = a_1_d, b_1_d = b_1_d,
    a_2_d = a_2_d, b_2_d = b_2_d,
    a_1_a = a_1_a, b_1_a = b_1_a,
    a_2_a = a_2_a, b_2_a = b_2_a,
    a_1_d_Hminus = a_1_d_Hminus, b_1_d_Hminus = b_1_d_Hminus,
    a_2_d_Hminus = a_2_d_Hminus, b_2_d_Hminus = b_2_d_Hminus
  )
  
  if (!fixed_search$found) {
    if (progress) {
      cat("Warning: no fixed-sample size reached target constraints.\n")
      utils::flush.console()
    }
    
    return(list(
      design = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
      naive_oc = list(
        n1 = NA_integer_,
        n2 = NA_integer_,
        power = NA_real_,
        t1e = NA_real_,
        pceH0 = NA_real_
      ),
      occ = NULL,
      priors = priors,
      conv = "no_feasible_fixed"
    ))
  }
  
  best_fixed <- fixed_search$best
  
  best_naive_oc <- list(
    n1 = best_fixed$n1,
    n2 = best_fixed$n2,
    power = unname(best_fixed$bayes_oc["Power"]),
    t1e = unname(best_fixed$bayes_oc["Type1_Error"]),
    pceH0 = unname(best_fixed$bayes_oc["CE_H0"])
  )
  
  final_n1 <- as.integer(best_fixed$n1)
  final_n2 <- as.integer(best_fixed$n2)
  
  if (progress) {
    cat(sprintf(
      " --> Fixed-sample size found: n_tot=%d (n1=%d, n2=%d, Power=%.3f, T1E=%.3f, PCE(H0)=%.3f)\n",
      best_fixed$n_tot, final_n1, final_n2,
      best_fixed$bayes_oc["Power"],
      best_fixed$bayes_oc["Type1_Error"],
      best_fixed$bayes_oc["CE_H0"]
    ))
    if (calibration_mode %in% c("frequentist", "hybrid")) {
      cat(sprintf(
        " Frequentist check passed: Power=%.3f, supT1E=%.3f\n",
        best_fixed$freq_power, best_fixed$freq_t1e
      ))
    }
    utils::flush.console()
  }
  
  ## Step 2: interim grid
  
  n1_1_min <- max(n1_min[1], ceiling(interim_fraction[1] * final_n1))
  n1_2_min <- max(n1_min[2], ceiling(interim_fraction[1] * final_n2))
  n1_1_max <- min(final_n1 - 1L, floor(interim_fraction[2] * final_n1))
  n1_2_max <- min(final_n2 - 1L, floor(interim_fraction[2] * final_n2))
  
  if (n1_1_min > n1_1_max || n1_2_min > n1_2_max) {
    if (progress) {
      cat("No interim design grid points; adjust n1_min, interim_fraction, or n2_max.\n")
      utils::flush.console()
    }
    
    return(list(
      design = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
      naive_oc = best_naive_oc,
      occ = NULL,
      priors = priors,
      conv = "no_interim_grid"
    ))
  }
  
  if (abs(alloc1 - alloc2) < .Machine$double.eps^0.5) {
    n1_bal_min <- max(n1_1_min, n1_2_min)
    n1_bal_max <- min(n1_1_max, n1_2_max)
    
    if (n1_bal_min > n1_bal_max) {
      if (progress) {
        cat("No balanced interim design grid points; adjust n1_min, interim_fraction, or n2_max.\n")
        utils::flush.console()
      }
      
      return(list(
        design = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
        naive_oc = best_naive_oc,
        occ = NULL,
        priors = priors,
        conv = "no_interim_grid"
      ))
    }
    
    n1_bal_grid <- seq.int(n1_bal_min, n1_bal_max, by = grid_step)
    n1_12_pairs <- data.frame(
      n1_1 = n1_bal_grid,
      n1_2 = n1_bal_grid
    )
    
  } else {
    n1_1_grid <- seq.int(n1_1_min, n1_1_max, by = grid_step)
    n1_2_grid <- seq.int(n1_2_min, n1_2_max, by = grid_step)
    n1_12_pairs <- expand.grid(n1_1 = n1_1_grid, n1_2 = n1_2_grid)
  }
  
  if (nrow(n1_12_pairs) == 0L) {
    if (progress) {
      cat("No interim design grid points; adjust grid_step or interim_fraction.\n")
      utils::flush.console()
    }
    
    return(list(
      design = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
      naive_oc = best_naive_oc,
      occ = NULL,
      priors = priors,
      conv = "no_interim_grid"
    ))
  }
  
  
  if (progress) {
    cat(" => Parallelizing over", nrow(n1_12_pairs),
        "interim designs using", ncores, "cores...\n")
    utils::flush.console()
  }
  
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, library(bfbin2arm))
  parallel::clusterExport(
    cl = cl,
    varlist = c(
      "compute_corrected_twostage_oc_2arm",
      "compute_freq_twostage_oc_2arm",
      "check_design_feasibility",
      "freq_EN_twoarm_null",
      "freq_oc_twostage_twoarm_fixed",
      "freq_t1e_twostage_twoarm_sup",
      "twoarmbinbf_plus0_direct",
      "final_n1", "final_n2",
      "k", "k_f", "test",
      "a_0_d", "b_0_d", "a_0_a", "b_0_a",
      "a_1_d", "b_1_d", "a_2_d", "b_2_d",
      "a_1_a", "b_1_a", "a_2_a", "b_2_a",
      "a_1_a_Hminus", "b_1_a_Hminus",
      "a_2_a_Hminus", "b_2_a_Hminus",
      "calibration_mode", "calibration_EN", "compute_freq_oc",
      "p1_EN_H0", "p2_EN_H0",
      "alpha", "beta", "pceH0",
      "alpha_freq", "beta_freq",
      "p1_power", "p2_power",
      "p_null_grid"
    ),
    envir = environment()
  )
  
  ## Build row list and initialize
  rows    <- asplit(n1_12_pairs, 1L)
  n_total <- length(rows)
  
  ## Choose a chunk size; 10 is a reasonable default
  chunk_size <- 10L
  n_chunks   <- ceiling(n_total / chunk_size)
  
  res_list <- vector("list", n_total)
  
  for (ch in seq_len(n_chunks)) {
    start_idx <- (ch - 1L) * chunk_size + 1L
    end_idx   <- min(ch * chunk_size, n_total)
    this_idx  <- start_idx:end_idx
    this_rows <- rows[this_idx]
    
    res_chunk <- parallel::parLapply(cl, this_rows, function(pair) {
      n1_1 <- as.integer(pair[["n1_1"]])
      n1_2 <- as.integer(pair[["n1_2"]])
      
      bayes <- tryCatch(
        compute_corrected_twostage_oc_2arm(
          n1_1 = n1_1,
          n1_2 = n1_2,
          n2_1 = final_n1,
          n2_2 = final_n2,
          k    = k,
          k_f  = k_f,
          test = test,
          a_0_d = a_0_d, b_0_d = b_0_d,
          a_0_a = a_0_a, b_0_a = b_0_a,
          a_1_d = a_1_d, b_1_d = b_1_d,
          a_2_d = a_2_d, b_2_d = b_2_d,
          a_1_a = a_1_a, b_1_a = b_1_a,
          a_2_a = a_2_a, b_2_a = b_2_a,
          a_1_d_Hminus = a_1_d_Hminus, b_1_d_Hminus = b_1_d_Hminus,
          a_2_d_Hminus = a_2_d_Hminus, b_2_d_Hminus = b_2_d_Hminus,
          a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
          a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
        ),
        error = function(e) NULL
      )
      if (is.null(bayes)) return(NULL)
      
      ## Expected sample size under H0 for ranking
      if (calibration_EN == "Bayesian" || is.null(p1_EN_H0) || is.null(p2_EN_H0)) {
        EN_val <- unname(bayes["E_H0_N"])
      } else {
        EN_val <- freq_EN_twoarm_null(
          n1_1 = n1_1, n1_2 = n1_2,
          n2_1 = final_n1, n2_2 = final_n2,
          p1_null = p1_EN_H0, p2_null = p2_EN_H0,
          k = k, k_f = k_f, test = test,
          a_0_a = a_0_a, b_0_a = b_0_a,
          a_1_a = a_1_a, b_1_a = b_1_a,
          a_2_a = a_2_a, b_2_a = b_2_a,
          a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
          a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
        )
      }
      
      ## Frequentist OCs for this candidate (only if needed)
      freq_vec <- NULL
      
      need_freq_candidate <-
        (calibration_mode %in% c("frequentist", "hybrid") || isTRUE(compute_freq_oc)) &&
        !is.null(p1_power) && !is.null(p2_power)
      
      if (need_freq_candidate) {
        freq_vec <- tryCatch(
          compute_freq_twostage_oc_2arm(
            n1_1 = n1_1, n1_2 = n1_2,
            n2_1 = final_n1, n2_2 = final_n2,
            k = k, k_f = k_f, test = test,
            p1_power = p1_power, p2_power = p2_power,
            p_null_grid = p_null_grid,
            a_0_a = a_0_a, b_0_a = b_0_a,
            a_1_a = a_1_a, b_1_a = b_1_a,
            a_2_a = a_2_a, b_2_a = b_2_a,
            a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
            a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
          ),
          error = function(e) NULL
        )
      }
      
      feas <- check_design_feasibility(
        bayes = bayes,
        freq_vec = freq_vec,
        alpha = alpha,
        beta = beta,
        pceH0 = pceH0,
        alpha_freq = alpha_freq,
        beta_freq = beta_freq
      )
      
      list(
        n1_1 = n1_1,
        n1_2 = n1_2,
        bayes = bayes,
        EN = EN_val,
        freq_oc = freq_vec,
        bayes_feasible = feas$bayes_ok,
        freq_feasible = feas$freq_ok,
        hybrid_feasible = feas$hybrid_ok
      )
    })
    
    ## Store chunk results back into the full list
    for (k_idx in seq_along(this_idx)) {
      res_list[[this_idx[k_idx]]] <- res_chunk[[k_idx]]
    }
    
    if (progress) {
      done <- end_idx
      cat(sprintf(
        "Step 2: evaluated %d / %d interim designs (%.1f%%)...\n",
        done, n_total, 100 * done / n_total
      ))
      utils::flush.console()
    }
  }
  
  ## Post-processing: mode-specific selection of the optimal two-stage design
  
  res_ok <- Filter(Negate(is.null), res_list)
  
  if (length(res_ok) == 0L) {
    if (progress) {
      cat("No two-stage design could be evaluated successfully.\n")
      utils::flush.console()
    }
    
    return(list(
      design = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
      naive_oc = best_naive_oc,
      occ = NULL,
      priors = priors,
      conv = "no_feasible_design"
    ))
  }
  
  bayes_ok  <- function(r) isTRUE(r$bayes_feasible)
  freq_ok   <- function(r) isTRUE(r$freq_feasible)
  hybrid_ok <- function(r) isTRUE(r$hybrid_feasible)
  
  if (calibration_mode == "Bayesian") {
    
    cand <- Filter(bayes_ok, res_ok)
    
    if (length(cand) > 0L) {
      ENs <- vapply(cand, function(r) r$EN, numeric(1))
      best_r <- cand[[which.min(ENs)]]
      conv_flag <- "converged"
    } else {
      ENs <- vapply(res_ok, function(r) r$EN, numeric(1))
      best_r <- res_ok[[which.min(ENs)]]
      conv_flag <- "no_fully_feasible_design_best_returned"
    }
    
  } else if (calibration_mode == "frequentist") {
    
    cand <- Filter(freq_ok, res_ok)
    
    if (length(cand) > 0L) {
      ENs <- vapply(cand, function(r) r$EN, numeric(1))
      best_r <- cand[[which.min(ENs)]]
      conv_flag <- "frequentist_EN_freq_constraints"
    } else {
      ENs <- vapply(res_ok, function(r) r$EN, numeric(1))
      best_r <- res_ok[[which.min(ENs)]]
      conv_flag <- "no_freq_feasible_design_best_returned"
    }
    
  } else if (calibration_mode == "hybrid") {
    
    cand <- Filter(hybrid_ok, res_ok)
    
    if (length(cand) > 0L) {
      ENs <- vapply(cand, function(r) r$EN, numeric(1))
      best_r <- cand[[which.min(ENs)]]
      conv_flag <- "hybrid_EN_bayesian_and_freq_constraints"
    } else {
      ENs <- vapply(res_ok, function(r) r$EN, numeric(1))
      best_r <- res_ok[[which.min(ENs)]]
      conv_flag <- "no_hybrid_feasible_design_best_returned"
    }
  }
  
  best_design <- c(best_r$n1_1, best_r$n1_2, final_n1, final_n2)
  best_occ <- best_r$bayes
  
  need_final_freq_oc <-
    calibration_mode %in% c("frequentist", "hybrid") ||
    isTRUE(compute_freq_oc)
  
  if (progress && need_final_freq_oc) {
    cat("Step 2 post-processing finished, calculating frequentist operating characteristics of the final design, this may take a while...\n")
    utils::flush.console()
  }
  
  out <- list(
    design = best_design,
    naive_oc = best_naive_oc,
    occ = best_occ,
    priors = priors,
    freq_occ = NULL,
    conv = conv_flag
  )
  
  freq_vec_final <- NULL
  
  if (need_final_freq_oc) {
    freq_vec <- tryCatch(
      compute_freq_twostage_oc_2arm(
        n1_1 = n1_1, n1_2 = n1_2,
        n2_1 = final_n1, n2_2 = final_n2,
        k = k, k_f = k_f, test = test,
        p1_power = p1_power, p2_power = p2_power,
        p_null_grid = p_null_grid,
        a_0_a = a_0_a, b_0_a = b_0_a,
        a_1_a = a_1_a, b_1_a = b_1_a,
        a_2_a = a_2_a, b_2_a = b_2_a,
        a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
        a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
      ),
      error = function(e) NULL
    )
    
    out$freq_occ <- freq_vec_final
    
    t1e_fix   <- unname(freq_vec_final[["Type1_Error_freq_fixed"]])
    power_fix <- unname(freq_vec_final[["Power_freq_fixed"]])
    t1e_two   <- unname(freq_vec_final[["Type1_Error_freq_two_stage"]])
    power_two <- unname(freq_vec_final[["Power_freq_two_stage"]])
    
    if (!is.na(t1e_fix) && !is.na(t1e_two) &&
        abs(t1e_two - t1e_fix) > 0.01) {
      warning(sprintf(
        paste0(
          "Exact two-stage frequentist type-I error (%.4f) differs from ",
          "fixed-sample (%.4f) by %.4f."
        ),
        t1e_two, t1e_fix, abs(t1e_two - t1e_fix)
      ))
    }
  }
  
  final_feas <- check_design_feasibility(
    bayes = out$occ,
    freq_vec = out$freq_occ,
    alpha = alpha,
    beta = beta,
    pceH0 = pceH0,
    alpha_freq = alpha_freq,
    beta_freq = beta_freq
  )
  
  if (!is.null(out$freq_occ)) {
    t1e_fix   <- unname(out$freq_occ[["Type1_Error_freq_fixed"]])
    power_fix <- unname(out$freq_occ[["Power_freq_fixed"]])
    t1e_two   <- unname(out$freq_occ[["Type1_Error_freq_two_stage"]])
    power_two <- unname(out$freq_occ[["Power_freq_two_stage"]])
    
    if (calibration_mode %in% c("frequentist", "hybrid") && !final_feas$freq_ok) {
      msg <- character(0)
      
      if (!is.na(alpha_freq) && t1e_two > alpha_freq) {
        msg <- c(msg, sprintf(
          "type-I error %.4f exceeds alpha_freq=%.4f",
          t1e_two, alpha_freq
        ))
      }
      
      if (!is.na(beta_freq) && power_two < (1 - beta_freq)) {
        msg <- c(msg, sprintf(
          "power %.4f is below target %.4f",
          power_two, 1 - beta_freq
        ))
      }
      
      if (length(msg) > 0L) {
        warning(paste(
          "Final exact two-stage frequentist constraints are not satisfied:",
          paste(msg, collapse = "; "),
          "."
        ))
      }
    }
    
    if (!is.na(t1e_fix) && !is.na(t1e_two) &&
        abs(t1e_two - t1e_fix) > 0.01) {
      warning(sprintf(
        paste0(
          "Exact two-stage frequentist type-I error (%.4f) differs from ",
          "fixed-sample (%.4f) by %.4f."
        ),
        t1e_two, t1e_fix, abs(t1e_two - t1e_fix)
      ))
    }
  }
  
  if (progress) {
    if (calibration_mode == "frequentist" && !final_feas$freq_ok) {
      cat(
        "Warning: no two-stage design satisfies the requested frequentist constraints under the final exact evaluation.\n",
        "Consider increasing n2_max, relaxing power_cushion, or loosening alpha_freq/beta_freq.\n",
        sep = ""
      )
      utils::flush.console()
    }
    
    if (calibration_mode == "hybrid" && !final_feas$hybrid_ok) {
      cat(
        "Warning: no two-stage design satisfies both Bayesian and frequentist constraints under the final exact evaluation.\n",
        "Consider increasing n2_max, relaxing power_cushion, or loosening alpha/beta and alpha_freq/beta_freq.\n",
        sep = ""
      )
      utils::flush.console()
    }
  }
  
  if (calibration_mode == "frequentist") {
    if (final_feas$freq_ok) {
      if (out$conv == "no_freq_feasible_design_best_returned") {
        out$conv <- "freq_constraints_satisfied_posthoc"
      }
    } else {
      if (out$conv == "frequentist_EN_freq_constraints") {
        out$conv <- "freq_constraints_violated"
      }
    }
  }
  
  if (calibration_mode == "hybrid") {
    if (final_feas$hybrid_ok) {
      if (out$conv == "no_hybrid_feasible_design_best_returned") {
        out$conv <- "hybrid_constraints_satisfied_posthoc"
      }
    } else {
      if (out$conv == "hybrid_EN_bayesian_and_freq_constraints") {
        out$conv <- "hybrid_constraints_violated"
      }
    }
  }
  
  return(out)
  
  out
}
