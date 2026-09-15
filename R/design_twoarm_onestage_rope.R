.format_twoarm_elapsed_time <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) {
    return("NA")
  }
  
  seconds <- round(seconds)
  
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  secs <- seconds %% 60L
  
  if (hours > 0L) {
    return(
      sprintf("%02d:%02d:%02d", hours, minutes, secs)
    )
  }
  
  sprintf("%02d:%02d", minutes, secs)
}

.twoarm_rope_evaluate_n_worker <- function(
    n,
    
    delta,
    gammaeq,
    gammadiff,
    
    analysis_prior_C,
    analysis_prior_T,
    
    design_prior_eq_C,
    design_prior_eq_T,
    design_prior_eq,
    
    design_prior_ne_C,
    design_prior_ne_T,
    
    z_eq,
    z_ne,
    
    integration,
    quad_nodes,
    rel.tol,
    
    compute_freq_power,
    freq_pC,
    freq_delta,
    
    compute_freq_type1,
    freq_grid_n,
    freq_pC_range,
    
    needs_bayes_power,
    needs_bayes_type1,
    needs_pce_h0,
    
    targetpower,
    targettype1,
    targetpce_h0,
    
    targetfreqpower,
    targetfreqtype1,
    
    calibration
) {
  fit <- .evaluate_twoarm_rope_design(
    nC = n,
    nT = n,
    
    delta = delta,
    gammaeq = gammaeq,
    gammadiff = gammadiff,
    
    analysis_prior_C = analysis_prior_C,
    analysis_prior_T = analysis_prior_T,
    
    design_prior_eq_C = design_prior_eq_C,
    design_prior_eq_T = design_prior_eq_T,
    design_prior_eq = design_prior_eq,
    
    design_prior_ne_C = design_prior_ne_C,
    design_prior_ne_T = design_prior_ne_T,
    
    z_eq = z_eq,
    z_ne = z_ne,
    
    integration = integration,
    quad_nodes = quad_nodes,
    rel.tol = rel.tol,
    
    compute_freq_power = compute_freq_power,
    freq_pC = freq_pC,
    freq_delta = freq_delta,
    
    compute_freq_type1 = compute_freq_type1,
    freq_grid_n = freq_grid_n,
    freq_pC_range = freq_pC_range,
    
    check_predictive = TRUE,
    return_matrices = FALSE
  )
  
  pce_h0_ok <- !needs_pce_h0 ||
    isTRUE(
      fit$p_diff_h0 >= targetpce_h0
    )
  
  bayes_ok <- (!needs_bayes_power ||
                 isTRUE(
                   fit$power >= targetpower
                 )) &&
    (!needs_bayes_type1 ||
       isTRUE(
         fit$type1 <= targettype1
       )) &&
    pce_h0_ok
  
  freq_ok <- (!compute_freq_power ||
                isTRUE(
                  fit$freqpower >= targetfreqpower
                )) &&
    (!compute_freq_type1 ||
       isTRUE(
         fit$freqtype1 <= targetfreqtype1
       ))
  
  hybrid_ok <- (!needs_bayes_power ||
                  isTRUE(
                    fit$power >= targetpower
                  )) &&
    (!compute_freq_type1 ||
       isTRUE(
         fit$freqtype1 <= targetfreqtype1
       )) &&
    pce_h0_ok
  
  feasible_pointwise <- switch(
    calibration,
    Bayesian = bayes_ok,
    frequentist = freq_ok,
    hybrid = hybrid_ok,
    full = bayes_ok && freq_ok
  )
  
  data.frame(
    nC = fit$nC,
    nT = fit$nT,
    n = n,
    N = fit$N,
    allocation = fit$allocation,
    
    power = fit$power,
    type1 = fit$type1,
    
    pce_h1 = fit$power,
    pce_h0 = fit$p_diff_h0,
    
    p_diff_h1 = fit$p_diff_h1,
    p_inc_h1 = fit$p_inc_h1,
    
    p_diff_h0 = fit$p_diff_h0,
    p_inc_h0 = fit$p_inc_h0,
    
    freqpower = fit$freqpower,
    freqdiff = fit$freqdiff,
    freqinc = fit$freqinc,
    
    freqtype1 = fit$freqtype1,
    freqtype1lower = fit$freqtype1lower,
    freqtype1upper = fit$freqtype1upper,
    
    freqtype1pClower = fit$freqtype1pClower,
    freqtype1pTlower = fit$freqtype1pTlower,
    
    freqtype1pCupper = fit$freqtype1pCupper,
    freqtype1pTupper = fit$freqtype1pTupper,
    
    z_eq = fit$z_eq,
    z_ne = fit$z_ne,
    
    predictive_sum_eq = fit$predictive_sum_eq,
    predictive_sum_ne = fit$predictive_sum_ne,
    
    pce_h0_ok = pce_h0_ok,
    
    bayes_ok = bayes_ok,
    freq_ok = freq_ok,
    hybrid_ok = hybrid_ok,
    
    feasible_pointwise = feasible_pointwise,
    
    stringsAsFactors = FALSE
  )
}

#' Calibrate a one-stage two-arm ROPE equivalence design for binary endpoints
#'
#' Calculates equal-allocation sample sizes for a parallel-group one-stage
#' equivalence study with binary endpoints. The treatment effect is the risk
#' difference `pT - pC`. Equivalence is declared when the posterior probability
#' of `abs(pT - pC) < delta` is at least `gammaeq`; meaningful non-equivalence
#' is declared when the posterior probability outside the ROPE is at least
#' `gammadiff`; otherwise the decision is inconclusive.
#'
#' @param nmin,nmax Minimum and maximum per-arm sample size.
#' @param nstep Positive per-arm sample-size increment.
#' @param delta ROPE half-width on the risk-difference scale.
#' @param gammaeq Posterior threshold for equivalence.
#' @param gammadiff Posterior threshold for meaningful non-equivalence.
#' @param analysis_prior_C,analysis_prior_T Length-two Beta analysis-prior shapes.
#' @param design_prior_eq_C,design_prior_eq_T Legacy equivalence-prior Beta kernels.
#' @param design_prior_eq Optional joint baseline-risk/risk-difference design prior.
#' @param design_prior_ne_C,design_prior_ne_T Non-equivalence-prior Beta kernels.
#' @param targetpower Target Bayesian predictive probability of declaring
#'   practical equivalence under the equivalence design distribution `G1`.
#' @param targettype1 Maximum Bayesian predictive probability of falsely
#'   declaring practical equivalence under the non-equivalence design
#'   distribution `G0`.
#' @param targetpce_h0 Optional minimum Bayesian predictive probability of
#'   compelling evidence for meaningful non-equivalence under `G0`. More
#'   precisely, this is the predictive probability of declaring meaningful
#'   difference, i.e. the posterior probability of `abs(pT - pC) >= delta`
#'   is at least `gammadiff`. If supplied, this criterion is enforced for
#'   `"Bayesian"`, `"hybrid"`, and `"full"` calibration.
#' @param sustainn Number of consecutive feasible evaluated sample sizes.
#' @param rel.tol Adaptive-integration relative tolerance.
#' @param integration Either `"adaptive"` or `"quantile"`.
#' @param quad_nodes Number of Gauss-Legendre nodes for quantile quadrature.
#' @param calibration One of `"Bayesian"`, `"frequentist"`, `"hybrid"`, or `"full"`.
#' @param freq_pC Control response probability for frequentist pointwise power.
#' @param freq_delta Treatment-control difference for frequentist pointwise power.
#' @param targetfreqpower Target pointwise frequentist equivalence power.
#' @param targetfreqtype1 Maximum boundary-grid frequentist type-I error.
#' @param freq_grid_n Boundary-grid size.
#' @param freq_pC_range Control-risk interval for boundary-grid evaluation.
#' @param parallel Logical; if `TRUE`, evaluate candidate per-arm sample sizes
#'   in parallel.
#' @param ncores Maximum number of worker processes used when
#'   `parallel = TRUE`.
#' @param parallel_backend Parallel backend. `"auto"` and `"psock"` use a
#'   PSOCK cluster. `"none"` evaluates sequentially. `"multicore"` is
#'   retained for compatibility but resolves to `"psock"`.
#' @param ncores Positive integer giving the maximum number of worker
#'   processes used when `parallel = TRUE`.
#' @param parallel_backend Parallel backend. `"auto"` and `"psock"` use a
#'   PSOCK cluster, which is supported on Windows, macOS, and Linux.
#'   `"none"` evaluates candidates sequentially. `"multicore"` is retained
#'   for compatibility but currently resolves to `"psock"` because forked
#'   processes can be unsafe with threaded numerical libraries.
#' @param ncores Number of workers if `parallel = TRUE`.
#' @param returngrid Return the full operating-characteristic grid.
#' @param returnmatrices Retain matrices for the selected design.
#' @param progress Logical; if `TRUE`, show a progress bar during sequential
#'   evaluation. Progress display is disabled for parallel evaluation.
#' @param progress_every Positive integer specifying the interval, in evaluated
#'   candidate sample sizes, at which sequential progress is updated.
#' @param returnmatrices Retain matrices for the selected design.
#'
#' @return An object of class `"bfbin2armrope2armdesign"`.
#' @export
design_twoarm_onestage_rope <- function(
    nmin,
    nmax,
    nstep = 1L,
    delta,
    gammaeq = 0.90,
    gammadiff = gammaeq,
    analysis_prior_C = c(1, 1),
    analysis_prior_T = c(1, 1),
    design_prior_eq_C = NULL,
    design_prior_eq_T = NULL,
    design_prior_eq = NULL,
    design_prior_ne_C,
    design_prior_ne_T,
    targetpower = NULL,
    targettype1 = NULL,
    targetpce_h0 = NULL,
    sustainn = 1L,
    rel.tol = 1e-10,
    integration = c("adaptive", "quantile"),
    quad_nodes = 128L,
    calibration = c("Bayesian", "frequentist", "hybrid", "full"),
    freq_pC = NULL,
    freq_delta = NULL,
    targetfreqpower = NULL,
    targetfreqtype1 = NULL,
    freq_grid_n = 101L,
    freq_pC_range = c(0, 1),
    parallel = FALSE,
    ncores = 1L,
    parallel_backend = c("auto", "psock", "multicore", "none"),
    returngrid = TRUE,
    returnmatrices = FALSE,
    progress = interactive(),
    progress_every = 1L
) {
  integration <- match.arg(integration)
  calibration <- match.arg(calibration)
  parallel_backend <- match.arg(parallel_backend)
  
  if (!isTRUE(parallel) || ncores <= 1L) {
    parallel_backend <- "none"
  }
  
  if (identical(parallel_backend, "auto")) {
    parallel_backend <- if (isTRUE(parallel)) "psock" else "none"
  }
  
  if (identical(parallel_backend, "multicore") &&
      .Platform$OS.type == "windows") {
    warning(
      "multicore parallelization is not available on Windows; using PSOCK.",
      call. = FALSE
    )
    
    parallel_backend <- "psock"
  }
  
  if (identical(parallel_backend, "multicore") &&
      !isTRUE(parallel)) {
    parallel_backend <- "none"
  }
  
  if (identical(parallel_backend, "psock") &&
      !isTRUE(parallel)) {
    parallel_backend <- "none"
  }
  
  if (identical(parallel_backend, "multicore")) {
    warning(
      "multicore backend is disabled for this design because forked BLAS ",
      "matrix multiplication may fail. Using PSOCK instead.",
      call. = FALSE
    )
    
    parallel_backend <- "psock"
  }
  
  .validate_count(nmin, "nmin")
  .validate_count(nmax, "nmax")
  .validate_count(nstep, "nstep")
  .validate_count(sustainn, "sustainn")
  .validate_count(freq_grid_n, "freq_grid_n")
  
  if (nmin < 1L) stop("nmin must be a positive integer.", call. = FALSE)
  if (nmax < nmin) stop("nmax must be greater than or equal to nmin.", call. = FALSE)
  if (nstep < 1L) stop("nstep must be a positive integer.", call. = FALSE)
  if (sustainn < 1L) stop("sustainn must be a positive integer.", call. = FALSE)
  if (freq_grid_n < 2L) stop("freq_grid_n must be at least 2.", call. = FALSE)
  
  if (!is.numeric(quad_nodes) || length(quad_nodes) != 1L ||
      !is.finite(quad_nodes) || quad_nodes < 2L ||
      quad_nodes != as.integer(quad_nodes)) {
    stop("quad_nodes must be an integer of at least 2.", call. = FALSE)
  }
  if (!is.numeric(rel.tol) || length(rel.tol) != 1L ||
      !is.finite(rel.tol) || rel.tol <= 0) {
    stop("rel.tol must be a positive finite number.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(ncores) || length(ncores) != 1L ||
      !is.finite(ncores) || ncores < 1L || ncores != as.integer(ncores)) {
    stop("ncores must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(freq_pC_range) || length(freq_pC_range) != 2L ||
      any(!is.finite(freq_pC_range)) || freq_pC_range[1L] < 0 ||
      freq_pC_range[2L] > 1 || freq_pC_range[1L] >= freq_pC_range[2L]) {
    stop("freq_pC_range must be c(lower, upper) with 0 <= lower < upper <= 1.", call. = FALSE)
  }
  
  if (!is.logical(progress) ||
      length(progress) != 1L ||
      is.na(progress)) {
    stop(
      "progress must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  .validate_count(progress_every, "progress_every")
  
  if (progress_every < 1L) {
    stop(
      "progress_every must be a positive integer.",
      call. = FALSE
    )
  }
  
  ncores <- as.integer(ncores)
  quad_nodes <- as.integer(quad_nodes)
  
  .validate_twoarm_rope_inputs(
    nC = nmin,
    nT = nmin,
    delta = delta,
    gammaeq = gammaeq,
    gammadiff = gammadiff,
    analysis_prior_C = analysis_prior_C,
    analysis_prior_T = analysis_prior_T
  )
  
  if (is.null(design_prior_eq)) {
    .validate_twoarm_design_prior(
      design_prior_eq_C,
      design_prior_eq_T,
      "design_prior_eq"
    )
  } else {
    .validate_baseline_difference_design_prior(
      design_prior_eq,
      delta,
      "design_prior_eq"
    )
  }
  
  .validate_twoarm_design_prior(
    design_prior_ne_C,
    design_prior_ne_T,
    "design_prior_ne"
  )
  
  needs_bayes_power <- calibration %in% c(
    "Bayesian",
    "hybrid",
    "full"
  )
  
  needs_bayes_type1 <- calibration %in% c(
    "Bayesian",
    "full"
  )
  
  needs_pce_h0 <- !is.null(targetpce_h0) &&
    calibration %in% c(
      "Bayesian",
      "hybrid",
      "full"
    )
  
  compute_freq_power <- calibration %in% c(
    "frequentist",
    "full"
  )
  
  compute_freq_type1 <- calibration %in% c(
    "frequentist",
    "hybrid",
    "full"
  )
  
  if (needs_bayes_power) .validate_probability(targetpower, "targetpower")
  if (needs_bayes_type1) .validate_probability(targettype1, "targettype1")
  
  if (!is.null(targetpce_h0)) {
    .validate_probability(
      targetpce_h0,
      "targetpce_h0"
    )
  }
  
  if (compute_freq_power) {
    .validate_probability(freq_pC, "freq_pC")
    .validate_probability(targetfreqpower, "targetfreqpower")
    
    if (!is.numeric(freq_delta) || length(freq_delta) != 1L ||
        !is.finite(freq_delta) || abs(freq_delta) >= delta ||
        freq_pC + freq_delta <= 0 || freq_pC + freq_delta >= 1) {
      stop(
        "freq_delta must define an interior equivalence scenario: abs(freq_delta) < delta and 0 < freq_pC + freq_delta < 1.",
        call. = FALSE
      )
    }
  }
  
  if (compute_freq_type1) {
    .validate_probability(targetfreqtype1, "targetfreqtype1")
  }
  
  if (isTRUE(parallel) && .Platform$OS.type == "windows") {
    warning(
      "parallel = TRUE is not supported through mclapply() on Windows; using sequential evaluation.",
      call. = FALSE
    )
    parallel <- FALSE
  }
  
  if (isTRUE(progress) && isTRUE(parallel)) {
    message(
      "Progress display is disabled when parallel = TRUE. ",
      "A completion summary will be printed when the search finishes."
    )
    
    progress <- FALSE
  }
  
  z_eq <- if (is.null(design_prior_eq)) {
    .trunc_beta_normconst_twoarm(
      delta, design_prior_eq_C, design_prior_eq_T,
      region = "equivalence", rel.tol = rel.tol
    )
  } else {
    NA_real_
  }
  
  z_ne <- .trunc_beta_normconst_twoarm(
    delta, design_prior_ne_C, design_prior_ne_T,
    region = "nonequivalence", rel.tol = rel.tol
  )
  
  nseq <- seq.int(nmin, nmax, by = nstep)
  if (tail(nseq, 1L) != nmax) nseq <- c(nseq, nmax)
  
  worker_args <- list(
    delta = delta,
    gammaeq = gammaeq,
    gammadiff = gammadiff,
    
    analysis_prior_C = analysis_prior_C,
    analysis_prior_T = analysis_prior_T,
    
    design_prior_eq_C = design_prior_eq_C,
    design_prior_eq_T = design_prior_eq_T,
    design_prior_eq = design_prior_eq,
    
    design_prior_ne_C = design_prior_ne_C,
    design_prior_ne_T = design_prior_ne_T,
    
    z_eq = z_eq,
    z_ne = z_ne,
    
    integration = integration,
    quad_nodes = quad_nodes,
    rel.tol = rel.tol,
    
    compute_freq_power = compute_freq_power,
    freq_pC = freq_pC,
    freq_delta = freq_delta,
    
    compute_freq_type1 = compute_freq_type1,
    freq_grid_n = freq_grid_n,
    freq_pC_range = freq_pC_range,
    
    needs_bayes_power = needs_bayes_power,
    needs_bayes_type1 = needs_bayes_type1,
    needs_pce_h0 = needs_pce_h0,
    
    targetpower = targetpower,
    targettype1 = targettype1,
    targetpce_h0 = targetpce_h0,
    
    targetfreqpower = targetfreqpower,
    targetfreqtype1 = targetfreqtype1,
    
    calibration = calibration
  )
  
  evaluate_n <- function(n) {
    do.call(
      .twoarm_rope_evaluate_n_worker,
      c(
        list(n = n),
        worker_args
      )
    )
  }
  
  search_start_time <- Sys.time()
  
  evaluate_one_n <- function(n) {
    do.call(
      .twoarm_rope_evaluate_n_worker,
      c(list(n = n), worker_args)
    )
  }
  
  if (identical(parallel_backend, "psock") &&
      length(nseq) > 1L) {
    
    n_workers <- min(ncores, length(nseq))
    
    message(
      "Evaluating ",
      length(nseq),
      " candidate per-arm sample sizes using ",
      n_workers,
      " PSOCK worker(s)..."
    )
    
    cl <- parallel::makeCluster(n_workers)
    
    on.exit(
      {
        if (!is.null(cl)) {
          parallel::stopCluster(cl)
        }
      },
      add = TRUE
    )
    
    # This code is needed only during package development with load_all().
    # In an installed package, the workers can use library(bfbin2arm).
    parallel::clusterExport(
      cl,
      varlist = c(
        ".twoarm_rope_evaluate_n_worker",
        
        ".evaluate_twoarm_rope_design",
        ".evaluate_twoarm_rope_design_base",
        ".normalize_twoarm_freq_outputs",
        
        ".validate_twoarm_rope_inputs",
        ".validate_twoarm_design_prior",
        ".validate_baseline_difference_design_prior",
        
        ".validate_count",
        ".validate_probability",
        ".validate_beta_prior",
        
        ".check_twoarm_counts",
        
        ".rope_twoarm_quad_rule",
        ".posterior_rope_prob_twoarm_quantile_row",
        ".posterior_rope_prob_twoarm",
        ".rope_prob_twoarm",
        ".rope_decision_matrix_twoarm",
        
        ".trunc_beta_normconst_twoarm",
        ".trunc_beta_predictive_twoarm_row_quantile",
        ".trunc_beta_predictive_twoarm",
        
        ".truncated_normal_quantile",
        ".baseline_difference_design_nodes",
        ".baseline_difference_predictive_matrix",
        
        ".pointwise_decision_probs_twoarm",
        ".freq_type1_boundary_twoarm",
        
        "%||%",
        
        "worker_args"
      ),
      envir = environment()
    )
    
    parallel::clusterEvalQ(
      cl,
      {
        if (!requireNamespace("statmod", quietly = TRUE)) {
          stop("Package 'statmod' is required on PSOCK workers.")
        }
        
        NULL
      }
    )
    
    rows <- parallel::parLapplyLB(
      cl,
      X = nseq,
      fun = function(n) {
        tryCatch(
          do.call(
            .twoarm_rope_evaluate_n_worker,
            c(list(n = n), worker_args)
          ),
          error = function(e) {
            structure(
              list(
                n = n,
                message = conditionMessage(e)
              ),
              class = "twoarm_rope_worker_error"
            )
          }
        )
      }
    )
    
  } else {
    
    if (isTRUE(progress)) {
      progress_bar <- utils::txtProgressBar(
        min = 0,
        max = length(nseq),
        style = 3
      )
      on.exit(close(progress_bar), add = TRUE)
    }
    
    rows <- vector("list", length(nseq))
    
    for (i in seq_along(nseq)) {
      rows[[i]] <- tryCatch(
        evaluate_one_n(nseq[i]),
        error = function(e) {
          structure(
            list(
              n = nseq[i],
              message = conditionMessage(e)
            ),
            class = "twoarm_rope_worker_error"
          )
        }
      )
      
      if (isTRUE(progress)) {
        utils::setTxtProgressBar(progress_bar, i)
      }
    }
  }
  
  search_elapsed_seconds <- as.numeric(
    difftime(Sys.time(), search_start_time, units = "secs")
  )
  
  if (isTRUE(progress) && identical(parallel_backend, "none")) {
    cat("\n")
  }
  
  message(
    "Sample-size search completed in ",
    .format_twoarm_elapsed_time(search_elapsed_seconds),
    "."
  )
  
  is_error <- vapply(
    rows,
    inherits,
    logical(1),
    what = "twoarm_rope_worker_error"
  )
  
  if (any(is_error)) {
    first_error <- rows[[which(is_error)[1L]]]
    
    stop(
      paste0(
        "Sample-size evaluation failed at n = ",
        first_error$n,
        ": ",
        first_error$message
      ),
      call. = FALSE
    )
  }
  
  if (!all(vapply(rows, is.data.frame, logical(1)))) {
    stop(
      "At least one sample-size evaluation did not return a data frame.",
      call. = FALSE
    )
  }
  
  grid <- do.call(rbind, rows)
  grid <- grid[order(grid$n), , drop = FALSE]
  rownames(grid) <- NULL
  
  grid$feasible <- FALSE
  for (i in seq_len(nrow(grid))) {
    j <- i + sustainn - 1L
    if (j <= nrow(grid)) {
      grid$feasible[i] <- all(grid$feasible_pointwise[i:j])
    }
  }
  
  if (any(grid$feasible)) {
    idx <- which(grid$feasible)[1L]
    nstar <- grid$n[idx]
    selected <- grid[idx, , drop = FALSE]
    
    selected_fit <- .evaluate_twoarm_rope_design(
      nC = nstar,
      nT = nstar,
      delta = delta,
      gammaeq = gammaeq,
      gammadiff = gammadiff,
      analysis_prior_C = analysis_prior_C,
      analysis_prior_T = analysis_prior_T,
      design_prior_eq_C = design_prior_eq_C,
      design_prior_eq_T = design_prior_eq_T,
      design_prior_eq = design_prior_eq,
      design_prior_ne_C = design_prior_ne_C,
      design_prior_ne_T = design_prior_ne_T,
      z_eq = z_eq,
      z_ne = z_ne,
      integration = integration,
      quad_nodes = quad_nodes,
      rel.tol = rel.tol,
      compute_freq_power = compute_freq_power,
      freq_pC = freq_pC,
      freq_delta = freq_delta,
      compute_freq_type1 = compute_freq_type1,
      freq_grid_n = freq_grid_n,
      freq_pC_range = freq_pC_range,
      check_predictive = TRUE,
      return_matrices = returnmatrices
    )
  } else {
    nstar <- NA_integer_
    selected <- NULL
    selected_fit <- NULL
  }
  
  out <- list(
    call = match.call(),
    inputs = list(
      nmin = nmin,
      nmax = nmax,
      nstep = nstep,
      delta = delta,
      gammaeq = gammaeq,
      gammadiff = gammadiff,
      analysis_prior_C = analysis_prior_C,
      analysis_prior_T = analysis_prior_T,
      design_prior_eq_C = design_prior_eq_C,
      design_prior_eq_T = design_prior_eq_T,
      design_prior_eq = design_prior_eq,
      design_prior_ne_C = design_prior_ne_C,
      design_prior_ne_T = design_prior_ne_T,
      calibration = calibration,
      targetpower = targetpower,
      targettype1 = targettype1,
      targetpce_h0 = targetpce_h0,
      freq_pC = freq_pC,
      freq_delta = freq_delta,
      targetfreqpower = targetfreqpower,
      targetfreqtype1 = targetfreqtype1,
      freq_grid_n = freq_grid_n,
      freq_pC_range = freq_pC_range,
      sustainn = sustainn,
      rel.tol = rel.tol,
      integration = integration,
      quad_nodes = quad_nodes,
      parallel = parallel,
      ncores = ncores,
      returngrid = returngrid,
      returnmatrices = returnmatrices,
      allocation = 1,
      progress = progress,
      progress_every = progress_every,
      parallel_backend = parallel_backend,
      parallel_workers = if (
        identical(parallel_backend, "psock")
      ) {
        min(ncores, length(nseq))
      } else {
        1L
      }
    ),
    nstar = nstar,
    Nstar = if (is.na(nstar)) NA_integer_ else 2L * nstar,
    selected = selected,
    selected_fit = selected_fit,
    grid = if (isTRUE(returngrid)) grid else NULL
  )
  
  class(out) <- "bfbin2armrope2armdesign"
  out
}
