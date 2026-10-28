.twoarm_directional_evaluate_n_worker <- function(
    n,
    
    delta_ni, delta_sup,
    gamma_ni, gamma_sup, gamma_inf,
    
    analysis_prior_C, analysis_prior_T,
    
    design_prior_inf, design_prior_ni,
    design_prior_sup, design_prior_not_sup,
    
    integration, quad_nodes, rel.tol,
    
    compute_freq_ni_power, freq_ni_pC, freq_ni_delta,
    compute_freq_sup_power, freq_sup_pC, freq_sup_delta,
    compute_freq_false_ni, compute_freq_false_sup,
    freq_grid_n, freq_pC_range,
    
    needs_ni_power, needs_sup_power,
    needs_false_ni, needs_false_sup,
    
    target_ni_power, target_sup_power,
    target_false_ni, target_false_sup,
    
    target_freq_ni_power, target_freq_sup_power,
    target_freq_false_ni, target_freq_false_sup,
    
    calibration
) {
  fit <- .evaluate_twoarm_directional_design(
    nC = n, nT = n,
    
    delta_ni = delta_ni, delta_sup = delta_sup,
    gamma_ni = gamma_ni, gamma_sup = gamma_sup, gamma_inf = gamma_inf,
    
    analysis_prior_C = analysis_prior_C, analysis_prior_T = analysis_prior_T,
    
    design_prior_inf = design_prior_inf, design_prior_ni = design_prior_ni,
    design_prior_sup = design_prior_sup, design_prior_not_sup = design_prior_not_sup,
    
    integration = integration, quad_nodes = quad_nodes, rel.tol = rel.tol,
    
    compute_freq_ni_power = compute_freq_ni_power,
    freq_ni_pC = freq_ni_pC, freq_ni_delta = freq_ni_delta,
    compute_freq_sup_power = compute_freq_sup_power,
    freq_sup_pC = freq_sup_pC, freq_sup_delta = freq_sup_delta,
    compute_freq_false_ni = compute_freq_false_ni,
    compute_freq_false_sup = compute_freq_false_sup,
    freq_grid_n = freq_grid_n, freq_pC_range = freq_pC_range,
    
    return_matrices = FALSE
  )
  
  bayes_ok <-
    (!needs_ni_power   || isTRUE(fit$ni_power   >= target_ni_power)) &&
    (!needs_sup_power  || isTRUE(fit$sup_power  >= target_sup_power)) &&
    (!needs_false_ni   || isTRUE(fit$false_ni   <= target_false_ni)) &&
    (!needs_false_sup  || isTRUE(fit$false_sup  <= target_false_sup))
  
  freq_ok <-
    (!compute_freq_ni_power  || isTRUE(fit$freq_ni_power  >= target_freq_ni_power)) &&
    (!compute_freq_sup_power || isTRUE(fit$freq_sup_power >= target_freq_sup_power)) &&
    (!compute_freq_false_ni  || isTRUE(fit$freq_false_ni  <= target_freq_false_ni)) &&
    (!compute_freq_false_sup || isTRUE(fit$freq_false_sup <= target_freq_false_sup))
  
  hybrid_ok <-
    (!needs_ni_power  || isTRUE(fit$ni_power  >= target_ni_power)) &&
    (!needs_sup_power || isTRUE(fit$sup_power >= target_sup_power)) &&
    (!compute_freq_false_ni  || isTRUE(fit$freq_false_ni  <= target_freq_false_ni)) &&
    (!compute_freq_false_sup || isTRUE(fit$freq_false_sup <= target_freq_false_sup))
  
  feasible_pointwise <- switch(
    calibration,
    Bayesian = bayes_ok,
    frequentist = freq_ok,
    hybrid = hybrid_ok,
    full = bayes_ok && freq_ok
  )
  
  data.frame(
    nC = fit$nC, nT = fit$nT, n = n, N = fit$N, allocation = fit$allocation,
    
    ni_power = fit$ni_power, sup_power = fit$sup_power,
    false_ni = fit$false_ni, false_sup = fit$false_sup,
    
    p_sup_g_inf = fit$p_sup_g_inf, p_ni_only_g_inf = fit$p_ni_only_g_inf,
    p_inf_g_inf = fit$p_inf_g_inf, p_inc_g_inf = fit$p_inc_g_inf,
    
    p_sup_g_ni = fit$p_sup_g_ni, p_ni_only_g_ni = fit$p_ni_only_g_ni,
    p_inf_g_ni = fit$p_inf_g_ni, p_inc_g_ni = fit$p_inc_g_ni,
    
    p_sup_g_sup = fit$p_sup_g_sup, p_ni_only_g_sup = fit$p_ni_only_g_sup,
    p_inf_g_sup = fit$p_inf_g_sup, p_inc_g_sup = fit$p_inc_g_sup,
    
    freq_ni_power = fit$freq_ni_power, freq_sup_power = fit$freq_sup_power,
    freq_false_ni = fit$freq_false_ni, freq_false_sup = fit$freq_false_sup,
    
    bayes_ok = bayes_ok, freq_ok = freq_ok, hybrid_ok = hybrid_ok,
    feasible_pointwise = feasible_pointwise,
    
    stringsAsFactors = FALSE
  )
}

#' Calibrate a one-stage two-arm directional ROPE design for binary endpoints
#'
#' Calculates equal-allocation sample sizes for a parallel-group one-stage
#' study with binary endpoints, where the treatment effect is the risk
#' difference `pT - pC`. Depending on `direction`, the design targets
#' non-inferiority, superiority, or a nested procedure in which superiority
#' is assessed only after non-inferiority has been established. Superiority
#' is declared when the posterior probability of `pT - pC > delta_sup` is at
#' least `gamma_sup`; non-inferiority is declared when the posterior
#' probability of `pT - pC > -delta_ni` is at least `gamma_ni`; harmful
#' inferiority is declared when the posterior probability of
#' `pT - pC <= -delta_ni` is at least `gamma_inf`; otherwise the decision is
#' inconclusive.
#'
#' @param nmin,nmax Minimum and maximum per-arm sample size.
#' @param nstep Positive per-arm sample-size increment.
#' @param direction One of `"noninferiority"`, `"superiority"`, or `"nested"`.
#' @param delta_ni Non-inferiority margin on the risk-difference scale.
#' @param delta_sup Superiority margin on the risk-difference scale, `>= 0`
#'   by convention, with `delta_sup > -delta_ni`.
#' @param gamma_ni Posterior threshold for non-inferiority.
#' @param gamma_sup Posterior threshold for superiority; must be at least
#'   `gamma_ni`.
#' @param gamma_inf Posterior threshold for harmful inferiority.
#' @param analysis_prior_C,analysis_prior_T Length-two Beta analysis-prior shapes.
#' @param design_prior_inf,design_prior_ni,design_prior_sup,design_prior_not_sup
#'   Baseline-risk/risk-difference design priors under the harmful-inferiority,
#'   non-inferiority, superiority, and not-superiority design distributions,
#'   respectively. Required depending on `direction` and the requested targets.
#' @param target_ni_power Target predictive probability of a non-inferiority
#'   or superiority conclusion, evaluated under `design_prior_ni`.
#' @param target_sup_power Target predictive probability of a superiority
#'   conclusion, evaluated under `design_prior_sup`.
#' @param target_false_ni Maximum predictive probability of a non-inferiority
#'   or superiority conclusion, evaluated under `design_prior_inf`.
#' @param target_false_sup Maximum predictive probability of a superiority
#'   conclusion, evaluated under `design_prior_not_sup`.
#' @param sustainn Number of consecutive feasible evaluated sample sizes.
#' @param rel.tol Adaptive-integration relative tolerance.
#' @param integration Either `"adaptive"` or `"quantile"`.
#' @param quad_nodes Number of Gauss-Legendre nodes for quantile quadrature.
#' @param calibration One of `"Bayesian"`, `"frequentist"`, `"hybrid"`, or `"full"`.
#' @param freq_ni_pC,freq_ni_delta Control response and true difference for
#'   pointwise frequentist non-inferiority power.
#' @param freq_sup_pC,freq_sup_delta Control response and true difference for
#'   pointwise frequentist superiority power.
#' @param target_freq_ni_power,target_freq_sup_power Target pointwise
#'   frequentist powers.
#' @param target_freq_false_ni,target_freq_false_sup Maximum boundary-grid
#'   frequentist error rates.
#' @param freq_grid_n Boundary-grid size.
#' @param freq_pC_range Control-risk interval for boundary-grid evaluation.
#' @param parallel Logical; if `TRUE`, evaluate candidate per-arm sample sizes
#'   in parallel.
#' @param ncores Positive integer giving the maximum number of worker
#'   processes used when `parallel = TRUE`.
#' @param parallel_backend Parallel backend. `"auto"` and `"psock"` use a
#'   PSOCK cluster. `"none"` evaluates sequentially.
#' @param returngrid Return the full operating-characteristic grid.
#' @param returnmatrices Retain matrices for the selected design.
#' @param progress Logical; if `TRUE`, show a progress bar during sequential
#'   evaluation.
#' @param progress_every Positive integer specifying the interval, in evaluated
#'   candidate sample sizes, at which sequential progress is updated.
#'
#' @return An object of class `"twoarm_rope_directional_onestage_design"`.
#' @export
design_twoarm_onestage_rope_directional <- function(
    nmin,
    nmax,
    nstep = 1L,
    
    direction = c("noninferiority", "superiority", "nested"),
    
    delta_ni,
    delta_sup = 0,
    
    gamma_ni = 0.975,
    gamma_sup = 0.975,
    gamma_inf = 0.975,
    
    analysis_prior_C = c(1, 1),
    analysis_prior_T = c(1, 1),
    
    design_prior_inf = NULL,
    design_prior_ni = NULL,
    design_prior_sup = NULL,
    design_prior_not_sup = NULL,
    
    target_ni_power = NULL,
    target_sup_power = NULL,
    target_false_ni = NULL,
    target_false_sup = NULL,
    
    sustainn = 1L,
    rel.tol = 1e-10,
    integration = c("adaptive", "quantile"),
    quad_nodes = 128L,
    
    calibration = c("Bayesian", "frequentist", "hybrid", "full"),
    
    freq_ni_pC = NULL,
    freq_ni_delta = NULL,
    freq_sup_pC = NULL,
    freq_sup_delta = NULL,
    
    target_freq_ni_power = NULL,
    target_freq_sup_power = NULL,
    target_freq_false_ni = NULL,
    target_freq_false_sup = NULL,
    
    freq_grid_n = 101L,
    freq_pC_range = c(0, 1),
    
    parallel = FALSE,
    ncores = 1L,
    parallel_backend = c("auto", "psock", "none"),
    
    returngrid = TRUE,
    returnmatrices = FALSE,
    
    progress = interactive(),
    progress_every = 1L
) {
  direction <- match.arg(direction)
  integration <- match.arg(integration)
  calibration <- match.arg(calibration)
  parallel_backend <- match.arg(parallel_backend)
  
  if (!isTRUE(parallel) || ncores <= 1L) parallel_backend <- "none"
  if (identical(parallel_backend, "auto")) {
    parallel_backend <- if (isTRUE(parallel)) "psock" else "none"
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
    stop("freq_pC_range must be c(lower, upper) with 0 <= lower < upper <= 1.",
         call. = FALSE)
  }
  
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("progress must be TRUE or FALSE.", call. = FALSE)
  }
  
  .validate_count(progress_every, "progress_every")
  if (progress_every < 1L) stop("progress_every must be a positive integer.", call. = FALSE)
  
  ncores <- as.integer(ncores)
  quad_nodes <- as.integer(quad_nodes)
  
  .validate_twoarm_directional_inputs(
    nC = nmin, nT = nmin,
    delta_ni = delta_ni, delta_sup = delta_sup,
    gamma_ni = gamma_ni, gamma_sup = gamma_sup, gamma_inf = gamma_inf,
    analysis_prior_C = analysis_prior_C, analysis_prior_T = analysis_prior_T
  )
  
  needs_ni_power  <- direction %in% c("noninferiority", "nested") && !is.null(target_ni_power)
  needs_sup_power <- direction %in% c("superiority", "nested") && !is.null(target_sup_power)
  needs_false_ni  <- direction %in% c("noninferiority", "nested") && !is.null(target_false_ni)
  needs_false_sup <- direction %in% c("superiority", "nested") && !is.null(target_false_sup)
  
  if (needs_ni_power) {
    .validate_twoarm_directional_design_prior(design_prior_ni, delta_ni, delta_sup, "design_prior_ni")
    .validate_probability(target_ni_power, "target_ni_power")
  }
  
  if (needs_sup_power) {
    .validate_twoarm_directional_design_prior(design_prior_sup, delta_ni, delta_sup, "design_prior_sup")
    .validate_probability(target_sup_power, "target_sup_power")
  }
  
  if (needs_false_ni) {
    .validate_twoarm_directional_design_prior(design_prior_inf, delta_ni, delta_sup, "design_prior_inf")
    .validate_probability(target_false_ni, "target_false_ni")
  }
  
  if (needs_false_sup) {
    .validate_twoarm_directional_design_prior(design_prior_not_sup, delta_ni, delta_sup, "design_prior_not_sup")
    .validate_probability(target_false_sup, "target_false_sup")
  }
  
  compute_freq_ni_power  <- !is.null(freq_ni_pC)
  compute_freq_sup_power <- !is.null(freq_sup_pC)
  compute_freq_false_ni  <- !is.null(target_freq_false_ni)
  compute_freq_false_sup <- !is.null(target_freq_false_sup)
  
  if (compute_freq_ni_power) {
    .validate_probability(freq_ni_pC, "freq_ni_pC")
    if (!is.null(target_freq_ni_power)) {
      .validate_probability(target_freq_ni_power, "target_freq_ni_power")
    }
    if (!is.numeric(freq_ni_delta) || length(freq_ni_delta) != 1L ||
        !is.finite(freq_ni_delta) || freq_ni_pC + freq_ni_delta <= 0 ||
        freq_ni_pC + freq_ni_delta >= 1) {
      stop("freq_ni_delta must define an admissible scenario for freq_ni_pC.",
           call. = FALSE)
    }
  }
  
  if (compute_freq_sup_power) {
    .validate_probability(freq_sup_pC, "freq_sup_pC")
    if (!is.null(target_freq_sup_power)) {
      .validate_probability(target_freq_sup_power, "target_freq_sup_power")
    }
    if (!is.numeric(freq_sup_delta) || length(freq_sup_delta) != 1L ||
        !is.finite(freq_sup_delta) || freq_sup_pC + freq_sup_delta <= 0 ||
        freq_sup_pC + freq_sup_delta >= 1) {
      stop("freq_sup_delta must define an admissible scenario for freq_sup_pC.",
           call. = FALSE)
    }
  }
  
  if (compute_freq_false_ni) .validate_probability(target_freq_false_ni, "target_freq_false_ni")
  if (compute_freq_false_sup) .validate_probability(target_freq_false_sup, "target_freq_false_sup")
  
  if (isTRUE(parallel) && .Platform$OS.type == "windows") {
    warning("parallel = TRUE is not supported on Windows; using sequential evaluation.",
            call. = FALSE)
    parallel <- FALSE
  }
  
  if (isTRUE(progress) && isTRUE(parallel)) {
    message("Progress display is disabled when parallel = TRUE.")
    progress <- FALSE
  }
  
  nseq <- seq.int(nmin, nmax, by = nstep)
  if (utils::tail(nseq, 1L) != nmax) nseq <- c(nseq, nmax)
  
  worker_args <- list(
    delta_ni = delta_ni, delta_sup = delta_sup,
    gamma_ni = gamma_ni, gamma_sup = gamma_sup, gamma_inf = gamma_inf,
    
    analysis_prior_C = analysis_prior_C, analysis_prior_T = analysis_prior_T,
    
    design_prior_inf = design_prior_inf, design_prior_ni = design_prior_ni,
    design_prior_sup = design_prior_sup, design_prior_not_sup = design_prior_not_sup,
    
    integration = integration, quad_nodes = quad_nodes, rel.tol = rel.tol,
    
    compute_freq_ni_power = compute_freq_ni_power,
    freq_ni_pC = freq_ni_pC, freq_ni_delta = freq_ni_delta,
    compute_freq_sup_power = compute_freq_sup_power,
    freq_sup_pC = freq_sup_pC, freq_sup_delta = freq_sup_delta,
    compute_freq_false_ni = compute_freq_false_ni,
    compute_freq_false_sup = compute_freq_false_sup,
    freq_grid_n = freq_grid_n, freq_pC_range = freq_pC_range,
    
    needs_ni_power = needs_ni_power, needs_sup_power = needs_sup_power,
    needs_false_ni = needs_false_ni, needs_false_sup = needs_false_sup,
    
    target_ni_power = target_ni_power, target_sup_power = target_sup_power,
    target_false_ni = target_false_ni, target_false_sup = target_false_sup,
    
    target_freq_ni_power = target_freq_ni_power,
    target_freq_sup_power = target_freq_sup_power,
    target_freq_false_ni = target_freq_false_ni,
    target_freq_false_sup = target_freq_false_sup,
    
    calibration = calibration
  )
  
  evaluate_one_n <- function(n) {
    do.call(.twoarm_directional_evaluate_n_worker, c(list(n = n), worker_args))
  }
  
  search_start_time <- Sys.time()
  
  if (identical(parallel_backend, "psock") && length(nseq) > 1L) {
    n_workers <- min(ncores, length(nseq))
    
    message("Evaluating ", length(nseq), " candidate per-arm sample sizes using ",
            n_workers, " PSOCK worker(s)...")
    
    cl <- parallel::makeCluster(n_workers)
    on.exit({ if (!is.null(cl)) parallel::stopCluster(cl) }, add = TRUE)
    
    parallel::clusterExport(
      cl,
      varlist = c(
        ".twoarm_directional_evaluate_n_worker",
        ".evaluate_twoarm_directional_design",
        
        ".validate_twoarm_directional_inputs",
        ".validate_twoarm_directional_design_prior",
        ".validate_count", ".validate_probability", ".validate_beta_prior",
        
        ".posterior_rd_cdf_twoarm",
        ".posterior_rd_cdf_twoarm_quantile_row",
        ".posterior_rd_cdf_twoarm_dispatch",
        ".posterior_directional_probs_twoarm",
        ".directional_decision_matrix_twoarm",
        
        ".rope_twoarm_quad_rule",
        ".truncated_normal_quantile",
        ".directional_design_nodes",
        ".directional_predictive_matrix_twoarm",
        
        ".pointwise_directional_probs_twoarm",
        ".freq_false_ni_boundary_twoarm",
        ".freq_false_sup_boundary_twoarm",
        
        "%||%",
        "worker_args"
      ),
      envir = environment()
    )
    
    parallel::clusterEvalQ(cl, {
      if (!requireNamespace("statmod", quietly = TRUE)) {
        stop("Package 'statmod' is required on PSOCK workers.")
      }
      NULL
    })
    
    rows <- parallel::parLapplyLB(
      cl, X = nseq,
      fun = function(n) {
        tryCatch(
          do.call(.twoarm_directional_evaluate_n_worker, c(list(n = n), worker_args)),
          error = function(e) {
            structure(list(n = n, message = conditionMessage(e)),
                      class = "twoarm_directional_worker_error")
          }
        )
      }
    )
    
    bad <- vapply(rows, function(r) inherits(r, "twoarm_directional_worker_error"), logical(1))
    if (any(bad)) {
      stop("Evaluation failed for n = ", rows[[which(bad)[1L]]]$n, ": ",
           rows[[which(bad)[1L]]]$message, call. = FALSE)
    }
    
  } else {
    
    if (isTRUE(progress)) {
      progress_bar <- utils::txtProgressBar(min = 0, max = length(nseq), style = 3)
    }
    
    rows <- vector("list", length(nseq))
    
    for (idx in seq_along(nseq)) {
      rows[[idx]] <- evaluate_one_n(nseq[idx])
      
      if (isTRUE(progress) && (idx %% progress_every == 0L || idx == length(nseq))) {
        utils::setTxtProgressBar(progress_bar, idx)
      }
    }
    
    if (isTRUE(progress)) close(progress_bar)
  }
  
  elapsed <- .format_twoarm_elapsed_time(as.numeric(difftime(Sys.time(), search_start_time, units = "secs")))
  
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
    
    selected_fit <- .evaluate_twoarm_directional_design(
      nC = nstar, nT = nstar,
      delta_ni = delta_ni, delta_sup = delta_sup,
      gamma_ni = gamma_ni, gamma_sup = gamma_sup, gamma_inf = gamma_inf,
      analysis_prior_C = analysis_prior_C, analysis_prior_T = analysis_prior_T,
      design_prior_inf = design_prior_inf, design_prior_ni = design_prior_ni,
      design_prior_sup = design_prior_sup, design_prior_not_sup = design_prior_not_sup,
      integration = integration, quad_nodes = quad_nodes, rel.tol = rel.tol,
      compute_freq_ni_power = compute_freq_ni_power,
      freq_ni_pC = freq_ni_pC, freq_ni_delta = freq_ni_delta,
      compute_freq_sup_power = compute_freq_sup_power,
      freq_sup_pC = freq_sup_pC, freq_sup_delta = freq_sup_delta,
      compute_freq_false_ni = compute_freq_false_ni,
      compute_freq_false_sup = compute_freq_false_sup,
      freq_grid_n = freq_grid_n, freq_pC_range = freq_pC_range,
      return_matrices = returnmatrices
    )
  } else {
    nstar <- NA_integer_
    selected <- NULL
    selected_fit <- NULL
    
    warning("No feasible per-arm sample size was found in [nmin, nmax].", call. = FALSE)
  }
  
  out <- list(
    inputs = list(
      nmin = nmin, nmax = nmax, nstep = nstep,
      direction = direction,
      delta_ni = delta_ni, delta_sup = delta_sup,
      gamma_ni = gamma_ni, gamma_sup = gamma_sup, gamma_inf = gamma_inf,
      analysis_prior_C = analysis_prior_C, analysis_prior_T = analysis_prior_T,
      design_prior_ni = design_prior_ni,
      design_prior_inf = design_prior_inf,
      design_prior_sup = design_prior_sup,
      design_prior_not_sup = design_prior_not_sup,
      target_ni_power = target_ni_power, target_sup_power = target_sup_power,
      target_false_ni = target_false_ni, target_false_sup = target_false_sup,
      target_freq_ni_power = target_freq_ni_power, target_freq_sup_power = target_freq_sup_power,
      target_freq_false_ni = target_freq_false_ni, target_freq_false_sup = target_freq_false_sup,
      calibration = calibration,
      elapsed = elapsed
    ),
    nstar = nstar,
    selected = selected,
    selected_fit = selected_fit,
    grid = if (isTRUE(returngrid)) grid else NULL
  )
  
  class(out) <- c("twoarm_rope_directional_onestage_design", "bfbin2arm_design")
  out
}