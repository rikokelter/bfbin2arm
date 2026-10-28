`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------------------
# Control settings
# -----------------------------------------------------------------------------

resolve_control_bf_ttest <- function(control = NULL) {
  defaults <- list(
    n_search = "two_stage",
    n_step_stage1 = 5L,
    rel_tol_m1 = 1e-8,
    subdivisions_m1 = 400L,
    grid_n = 401L,
    t_max = 10,
    t_max_cap = 40,
    root_tol = 1e-10,
    verify = TRUE,
    verify_rel_tol_m1 = 1e-10,
    verify_subdivisions_m1 = 800L,
    verify_grid_n = 801L,
    verify_t_max = 12,
    verify_root_tol = 1e-12,
    progress = interactive(),
    progress_style = 3L
  )
  
  if (is.null(control)) {
    return(defaults)
  }
  
  if (!is.list(control)) {
    stop("`control` must be a list or NULL.", call. = FALSE)
  }
  
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    warning(
      "Unknown `control` entries ignored: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  
  control <- control[names(control) %in% names(defaults)]
  out <- utils::modifyList(defaults, control)
  
  if (!out$n_search %in% c("two_stage", "full", "full_fast")) {
    stop("`control$n_search` must be 'two_stage', 'full', or 'full_fast'.", call. = FALSE)
  }
  
  out$n_step_stage1 <- max(1L, as.integer(out$n_step_stage1))
  out$grid_n <- max(3L, as.integer(out$grid_n))
  out$verify_grid_n <- max(3L, as.integer(out$verify_grid_n))
  out$subdivisions_m1 <- max(1L, as.integer(out$subdivisions_m1))
  out$verify_subdivisions_m1 <- max(1L, as.integer(out$verify_subdivisions_m1))
  out$t_max <- as.numeric(out$t_max)[1L]
  out$t_max_cap <- as.numeric(out$t_max_cap)[1L]
  out$root_tol <- as.numeric(out$root_tol)[1L]
  out$verify_t_max <- as.numeric(out$verify_t_max)[1L]
  out$verify_root_tol <- as.numeric(out$verify_root_tol)[1L]
  out$progress <- isTRUE(out$progress)
  out$progress_style <- as.integer(out$progress_style)[1L]
  
  if (!is.finite(out$t_max) || out$t_max <= 0) {
    stop("`control$t_max` must be one positive finite number.", call. = FALSE)
  }
  
  if (!is.finite(out$t_max_cap) || out$t_max_cap < out$t_max) {
    stop("`control$t_max_cap` must be finite and at least `control$t_max`.", call. = FALSE)
  }
  
  out
}

# -----------------------------------------------------------------------------
# Prior helpers
# -----------------------------------------------------------------------------

.make_design_prior <- function(spec) {
  if (!is.list(spec) || is.null(spec$type)) {
    stop(
      "Design prior specification must be a list with `type`.",
      call. = FALSE
    )
  }
  
  type <- tolower(as.character(spec$type)[1L])
  
  if (identical(type, "point")) {
    effect <- as.numeric(spec$effect)
    
    if (length(effect) != 1L || !is.finite(effect)) {
      stop(
        "Point design prior requires one finite `effect`.",
        call. = FALSE
      )
    }
    
    return(list(type = "point", effect = effect))
  }
  
  if (identical(type, "normal")) {
    mean <- as.numeric(spec$mean)
    sd <- as.numeric(spec$sd)
    
    if (length(mean) != 1L || !is.finite(mean)) {
      stop("`mean` must be one finite number.", call. = FALSE)
    }
    if (length(sd) != 1L || !is.finite(sd) || sd <= 0) {
      stop("`sd` must be one positive finite number.", call. = FALSE)
    }
    
    return(make_prior(
      family = "normal",
      location = mean,
      scale = sd,
      support = c(-Inf, Inf),
      name = sprintf(
        "N(%s, %s^2)",
        format(mean, trim = TRUE),
        format(sd, trim = TRUE)
      )
    ))
  }
  
  if (identical(type, "truncnormal")) {
    mean <- as.numeric(spec$mean)
    sd <- as.numeric(spec$sd)
    lower <- as.numeric(spec$lower)
    upper <- as.numeric(spec$upper)
    
    if (length(mean) != 1L || !is.finite(mean)) {
      stop("`mean` must be one finite number.", call. = FALSE)
    }
    if (length(sd) != 1L || !is.finite(sd) || sd <= 0) {
      stop("`sd` must be one positive finite number.", call. = FALSE)
    }
    if (length(lower) != 1L || length(upper) != 1L || lower >= upper) {
      stop("Need bounds with `lower < upper`.", call. = FALSE)
    }
    
    return(make_prior(
      family = "normal",
      location = mean,
      scale = sd,
      support = c(lower, upper),
      name = sprintf(
        "TN(%s, %s^2; [%s, %s])",
        format(mean, trim = TRUE),
        format(sd, trim = TRUE),
        format(lower, trim = TRUE),
        format(upper, trim = TRUE)
      )
    ))
  }
  
  if (identical(type, "cauchy")) {
    location <- as.numeric(spec$location %||% spec$mean)
    scale <- as.numeric(spec$scale)
    
    if (length(location) != 1L || !is.finite(location)) {
      stop(
        "Cauchy design prior requires one finite `location` (or `mean`).",
        call. = FALSE
      )
    }
    if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
      stop(
        "Cauchy design prior requires one positive finite `scale`.",
        call. = FALSE
      )
    }
    
    lower <- spec$lower %||% -Inf
    upper <- spec$upper %||% Inf
    lower <- as.numeric(lower)
    upper <- as.numeric(upper)
    
    if (length(lower) != 1L || length(upper) != 1L || lower >= upper) {
      stop("Need bounds with `lower < upper`.", call. = FALSE)
    }
    
    return(make_prior(
      family = "cauchy",
      location = location,
      scale = scale,
      support = c(lower, upper),
      name = sprintf(
        "Cauchy(%s, %s; [%s, %s])",
        format(location, trim = TRUE),
        format(scale, trim = TRUE),
        format(lower, trim = TRUE),
        format(upper, trim = TRUE)
      )
    ))
  }
  
  stop(
    "`type` must be one of: 'point', 'normal', 'truncnormal', or 'cauchy'.",
    call. = FALSE
  )
}

integrate_over_design_prior <- function(fun_delta, prior_spec, subdivisions = 400L) {
  prior <- .make_design_prior(prior_spec)
  
  if (identical(prior$type, "point")) {
    return(as.numeric(fun_delta(prior$effect)))
  }
  
  integrand <- function(x) fun_delta(x) * prior_pdf(x, prior)
  
  stats::integrate(
    integrand,
    lower = prior$support[1L],
    upper = prior$support[2L],
    subdivisions = subdivisions,
    rel.tol = 1e-8,
    stop.on.error = FALSE
  )$value
}

.make_bf_grid <- function(log_bf_fun, lower, upper, grid_n) {
  grid <- seq(lower, upper, length.out = as.integer(grid_n))
  values <- vapply(grid, log_bf_fun, numeric(1))
  keep <- is.finite(values)
  
  list(
    t = grid[keep],
    log_bf01 = values[keep]
  )
}

.find_roots_from_grid <- function(grid_t, grid_log_bf, threshold, log_bf_fun, root_tol = 1e-10) {
  target <- log(threshold)
  values <- grid_log_bf - target
  
  if (length(grid_t) < 2L) {
    return(numeric(0))
  }
  
  roots <- grid_t[abs(values) <= root_tol]
  idx <- which(values[-length(values)] * values[-1L] < 0)
  
  if (length(idx)) {
    g <- function(x) log_bf_fun(x) - target
    roots <- c(
      roots,
      vapply(
        idx,
        function(i) {
          stats::uniroot(g, interval = c(grid_t[i], grid_t[i + 1L]), tol = root_tol)$root
        },
        numeric(1)
      )
    )
  }
  
  roots <- sort(unique(roots))
  dedup_tol <- max(10 * root_tol, sqrt(.Machine$double.eps))
  if (length(roots) > 1L) {
    roots <- roots[c(TRUE, diff(roots) > dedup_tol)]
  }
  
  roots
}

find_bf_roots_general <- function(
    log_bf_fun,
    threshold,
    lower = -10,
    upper = 10,
    grid_n = 401L,
    root_tol = 1e-10
) {
  if (!is.function(log_bf_fun)) {
    stop("`log_bf_fun` must be a function.", call. = FALSE)
  }
  
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("`threshold` must be one positive finite number.", call. = FALSE)
  }
  
  if (length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop("`lower` and `upper` must be finite with `lower < upper`.", call. = FALSE)
  }
  
  bf_grid <- .make_bf_grid(log_bf_fun, lower = lower, upper = upper, grid_n = grid_n)
  
  .find_roots_from_grid(
    grid_t = bf_grid$t,
    grid_log_bf = bf_grid$log_bf01,
    threshold = threshold,
    log_bf_fun = log_bf_fun,
    root_tol = root_tol
  )
}

.classify_decision_intervals_bf_ttest <- function(
    log_bf_fun,
    threshold,
    roots,
    lower,
    upper,
    keep_when = c(">=", "<="),
    tail_probe = 1
) {
  keep_when <- match.arg(keep_when)
  
  if (!is.function(log_bf_fun)) {
    stop("`log_bf_fun` must be a function.", call. = FALSE)
  }
  
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("`threshold` must be one positive finite number.", call. = FALSE)
  }
  
  if (length(tail_probe) != 1L || !is.finite(tail_probe) || tail_probe <= 0) {
    stop("`tail_probe` must be one positive finite number.", call. = FALSE)
  }
  
  roots <- sort(unique(as.numeric(roots)))
  roots <- roots[is.finite(roots)]
  
  target <- log(threshold)
  cuts <- c(lower, roots, upper)
  kept <- list()
  
  for (i in seq_len(length(cuts) - 1L)) {
    a <- cuts[i]
    b <- cuts[i + 1L]
    
    probe <- if (is.infinite(a) && is.finite(b)) {
      b - tail_probe
    } else if (is.finite(a) && is.infinite(b)) {
      a + tail_probe
    } else if (is.finite(a) && is.finite(b)) {
      0.5 * (a + b)
    } else {
      next
    }
    
    val <- log_bf_fun(probe)
    if (!is.finite(val)) {
      next
    }
    
    keep <- if (identical(keep_when, ">=")) val >= target else val <= target
    if (keep) {
      kept[[length(kept) + 1L]] <- c(a, b)
    }
  }
  
  if (!length(kept)) {
    return(matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("lower", "upper"))))
  }
  
  out <- do.call(rbind, kept)
  colnames(out) <- c("lower", "upper")
  out
}

probability_over_intervals <- function(intervals, delta, df, n_eff) {
  if (is.null(intervals) || length(intervals) == 0L || nrow(intervals) == 0L) {
    return(rep.int(0, length(delta)))
  }
  
  if (!is.numeric(delta) || any(!is.finite(delta))) {
    stop("`delta` must contain only finite numeric values.", call. = FALSE)
  }
  
  cdf_nct <- function(x, df, ncp) {
    if (is.infinite(x) && x < 0) return(0)
    if (is.infinite(x) && x > 0) return(1)
    pnt_stable(x, df = df, ncp = ncp)
  }
  
  vapply(delta, function(d) {
    ncp <- sqrt(n_eff) * d
    sum(vapply(seq_len(nrow(intervals)), function(i) {
      lower <- as.numeric(intervals[i, 1L])
      upper <- as.numeric(intervals[i, 2L])
      
      cdf_nct(upper, df = df, ncp = ncp) -
        cdf_nct(lower, df = df, ncp = ncp)
    }, numeric(1)))
  }, numeric(1))
}

# -----------------------------------------------------------------------------
# Constraint handling
# -----------------------------------------------------------------------------

.normalize_freq_power_bf_ttest <- function(freq_power) {
  if (is.null(freq_power) || length(freq_power) == 0L) {
    return(list())
  }
  
  if (!is.list(freq_power)) {
    stop("`freq_power` must be a list or NULL.", call. = FALSE)
  }
  
  out <- lapply(seq_along(freq_power), function(i) {
    x <- freq_power[[i]]
    
    if (!is.list(x) || is.null(x$effect)) {
      stop(
        "Each `freq_power` specification must contain numeric `effect`; `target` is optional.",
        call. = FALSE
      )
    }
    
    if (length(x$effect) != 1L || !is.numeric(x$effect) || !is.finite(x$effect)) {
      stop("`freq_power[[i]]$effect` must be one finite number.", call. = FALSE)
    }
    
    target <- x$target
    if (!is.null(target) &&
        (length(target) != 1L || !is.numeric(target) || !is.finite(target) ||
         target < 0 || target > 1)) {
      stop("`freq_power[[i]]$target` must be NULL or a number in [0, 1].", call. = FALSE)
    }
    
    list(
      effect = as.numeric(x$effect),
      target = if (is.null(target)) NULL else as.numeric(target)
    )
  })
  
  nms <- names(freq_power)
  if (is.null(nms)) {
    nms <- rep("", length(out))
  }
  blank <- !nzchar(nms)
  nms[blank] <- paste0("power_", which(blank))
  names(out) <- nms
  
  out
}

.meets_constraints_bf_ttest <- function(
    oc,
    freq_t1e = NULL,
    freq_power = NULL,
    bayes_power = NULL,
    bayes_t1e = NULL,
    pce_h0 = NULL,
    tol = 1e-12
) {
  if (!isTRUE(oc$ok)) {
    return(FALSE)
  }
  
  freq_power <- .normalize_freq_power_bf_ttest(freq_power)
  
  if (!is.null(freq_t1e) && !(oc$freq_t1e <= freq_t1e + tol)) return(FALSE)
  if (!is.null(bayes_t1e) && !(oc$bayes_t1e <= bayes_t1e + tol)) return(FALSE)
  if (!is.null(bayes_power) && !(oc$bayes_power >= bayes_power - tol)) return(FALSE)
  if (!is.null(pce_h0) && !(oc$pce_h0 >= pce_h0 - tol)) return(FALSE)
  
  if (length(freq_power)) {
    if (is.null(oc$freq_power) || !length(oc$freq_power)) return(FALSE)
    for (nm in names(freq_power)) {
      tgt <- freq_power[[nm]]$target
      if (is.null(tgt)) next
      estimate <- oc$freq_power[[nm]]
      if (is.null(estimate) || is.na(estimate) || !(estimate >= tgt - tol)) {
        return(FALSE)
      }
    }
  }
  
  TRUE
}

.constraint_table_bf_ttest <- function(
    oc,
    freq_t1e = NULL,
    freq_power = NULL,
    bayes_power = NULL,
    bayes_t1e = NULL,
    pce_h0 = NULL,
    tol = 1e-12
) {
  rows <- list()
  
  add_row <- function(constraint, type, estimand, estimate, target, direction) {
    satisfied <- if (direction == "<=") estimate <= target + tol else estimate >= target - tol
    slack <- if (direction == "<=") target - estimate else estimate - target
    
    rows[[length(rows) + 1L]] <<- data.frame(
      constraint = constraint,
      type = type,
      estimand = estimand,
      estimate = unname(estimate),
      target = unname(target),
      direction = direction,
      slack = unname(slack),
      satisfied = unname(satisfied),
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.null(freq_t1e)) {
    add_row("freq_t1e", "frequentist", "P_delta=0(BF01 <= k1)", oc$freq_t1e, freq_t1e, "<=")
  }
  
  if (!is.null(bayes_t1e)) {
    add_row("bayes_t1e", "Bayesian", "P_design,H0(BF01 <= k1)", oc$bayes_t1e, bayes_t1e, "<=")
  }
  
  if (!is.null(bayes_power)) {
    add_row("bayes_power", "Bayesian", "P_design,H1(BF01 <= k1)", oc$bayes_power, bayes_power, ">=")
  }
  
  if (!is.null(pce_h0)) {
    add_row("pce_h0", "Bayesian", "P_design,H0(BF01 >= k0)", oc$pce_h0, pce_h0, ">=")
  }
  
  freq_power <- .normalize_freq_power_bf_ttest(freq_power)
  if (length(freq_power)) {
    for (nm in names(freq_power)) {
      if (is.null(freq_power[[nm]]$target)) next
      effect <- freq_power[[nm]]$effect
      add_row(
        paste0("freq_power_", nm),
        "frequentist",
        sprintf("P_delta=%s(BF01 <= k1)", format(effect, trim = TRUE)),
        oc$freq_power[[nm]],
        freq_power[[nm]]$target,
        ">="
      )
    }
  }
  
  if (!length(rows)) {
    return(data.frame())
  }
  
  do.call(rbind, rows)
}

.failure_score_bf_ttest <- function(
    oc,
    freq_t1e = NULL,
    freq_power = NULL,
    bayes_power = NULL,
    bayes_t1e = NULL,
    pce_h0 = NULL,
    tol = 1e-12
) {
  if (!isTRUE(oc$ok)) {
    return(Inf)
  }
  
  tab <- .constraint_table_bf_ttest(
    oc = oc,
    freq_t1e = freq_t1e,
    freq_power = freq_power,
    bayes_power = bayes_power,
    bayes_t1e = bayes_t1e,
    pce_h0 = pce_h0,
    tol = tol
  )
  
  if (!nrow(tab)) {
    return(0)
  }
  
  deficit <- ifelse(
    tab$direction == "<=",
    pmax(0, tab$estimate - tab$target),
    pmax(0, tab$target - tab$estimate)
  )
  
  max(deficit)
}

.best_candidate_bf_ttest <- function(
    results,
    grid,
    freq_t1e = NULL,
    freq_power = NULL,
    bayes_power = NULL,
    bayes_t1e = NULL,
    pce_h0 = NULL,
    tol = 1e-12
) {
  if (is.null(results) || !length(results)) {
    return(NULL)
  }
  
  nonnull_idx <- which(!vapply(results, is.null, logical(1)))
  if (!length(nonnull_idx)) {
    return(NULL)
  }
  
  ok_idx <- nonnull_idx[vapply(results[nonnull_idx], function(z) isTRUE(z$ok), logical(1))]
  if (length(ok_idx)) {
    scores <- vapply(
      ok_idx,
      function(j) .failure_score_bf_ttest(
        oc = results[[j]],
        freq_t1e = freq_t1e,
        freq_power = freq_power,
        bayes_power = bayes_power,
        bayes_t1e = bayes_t1e,
        pce_h0 = pce_h0,
        tol = tol
      ),
      numeric(1)
    )
    
    pick <- ok_idx[order(scores, grid[ok_idx])[1L]]
    oc <- results[[pick]]
    tab <- .constraint_table_bf_ttest(
      oc = oc,
      freq_t1e = freq_t1e,
      freq_power = freq_power,
      bayes_power = bayes_power,
      bayes_t1e = bayes_t1e,
      pce_h0 = pce_h0,
      tol = tol
    )
    
    return(list(
      n = grid[pick],
      oc = oc,
      constraint_table = tab,
      score = .failure_score_bf_ttest(
        oc = oc,
        freq_t1e = freq_t1e,
        freq_power = freq_power,
        bayes_power = bayes_power,
        bayes_t1e = bayes_t1e,
        pce_h0 = pce_h0,
        tol = tol
      )
    ))
  }
  
  pick <- nonnull_idx[1L]
  list(n = grid[pick], oc = results[[pick]], constraint_table = NULL, score = Inf)
}

# -----------------------------------------------------------------------------
# Bayes-factor model construction
# -----------------------------------------------------------------------------

make_log_bf01_for_hypothesis <- function(
    hypothesis,
    df,
    n_eff,
    prior_h1 = NULL,
    prior_h0 = NULL,
    rel.tol = 1e-9,
    subdivisions = 500L
) {
  hypothesis <- match.arg(
    hypothesis,
    c("two.sided.pointnull", "greater.pointnull", "greater.composite")
  )
  
  if (identical(hypothesis, "greater.composite")) {
    if (is.null(prior_h0) || is.null(prior_h1)) {
      stop("`greater.composite` requires both `prior_h0` and `prior_h1`.", call. = FALSE)
    }
    
    return(function(t) {
      log_m0 <- log_m1_t_generic(
        t = t,
        df = df,
        n_eff = n_eff,
        prior = prior_h0,
        rel.tol = rel.tol,
        subdivisions = subdivisions
      )
      
      log_m1 <- log_m1_t_generic(
        t = t,
        df = df,
        n_eff = n_eff,
        prior = prior_h1,
        rel.tol = rel.tol,
        subdivisions = subdivisions
      )
      
      log_m0 - log_m1
    })
  }
  
  if (is.null(prior_h1)) {
    stop("Point-null hypotheses require `prior_h1`.", call. = FALSE)
  }
  
  if (identical(hypothesis, "greater.pointnull")) {
    return(make_log_bf01_greater_pointnull_fun(
      df = df,
      n_eff = n_eff,
      prior_h1 = prior_h1,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    ))
  }
  
  function(t) {
    log_bf01_t_generic(
      t = t,
      df = df,
      n_eff = n_eff,
      prior = prior_h1,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    )
  }
}

.find_bracketing_roots_adaptive <- function(
    log_bf_fun,
    threshold,
    initial_t_max,
    t_max_cap,
    grid_n,
    root_tol,
    symmetric = FALSE
) {
  current_t_max <- initial_t_max
  
  repeat {
    if (symmetric) {
      bf_grid <- .make_bf_grid(
        log_bf_fun = log_bf_fun,
        lower = 0,
        upper = current_t_max,
        grid_n = grid_n
      )
      
      roots_pos <- .find_roots_from_grid(
        grid_t = bf_grid$t,
        grid_log_bf = bf_grid$log_bf01,
        threshold = threshold,
        log_bf_fun = log_bf_fun,
        root_tol = root_tol
      )
      
      if (length(roots_pos) >= 1L || current_t_max >= t_max_cap) {
        roots_pos <- roots_pos[roots_pos > 0]
        roots <- sort(c(-rev(roots_pos), roots_pos))
        return(list(
          roots = roots,
          t_max_used = current_t_max,
          bf_grid = bf_grid,
          symmetric = TRUE
        ))
      }
    } else {
      bf_grid <- .make_bf_grid(
        log_bf_fun = log_bf_fun,
        lower = -current_t_max,
        upper = current_t_max,
        grid_n = grid_n
      )
      
      roots <- .find_roots_from_grid(
        grid_t = bf_grid$t,
        grid_log_bf = bf_grid$log_bf01,
        threshold = threshold,
        log_bf_fun = log_bf_fun,
        root_tol = root_tol
      )
      
      if (length(roots) >= 1L || current_t_max >= t_max_cap) {
        return(list(
          roots = roots,
          t_max_used = current_t_max,
          bf_grid = bf_grid,
          symmetric = FALSE
        ))
      }
    }
    
    current_t_max <- min(t_max_cap, max(current_t_max + 2, ceiling(current_t_max * 1.5)))
  }
}

# -----------------------------------------------------------------------------
# Operating characteristics
# -----------------------------------------------------------------------------

compute_oc_one_sample <- function(
    n,
    k1,
    k0,
    analysis_prior = NULL,
    effect_null = 0,
    freq_power = list(),
    want_pce = FALSE,
    control = NULL,
    hypothesis = "two.sided.pointnull",
    design_prior_h0 = list(type = "point", effect = 0),
    design_prior_h1 = list(type = "point", effect = 0.4),
    analysis_prior_h0 = NULL,
    analysis_prior_h1 = NULL,
    analysis_prior_sd = 0.5
) {
  hypothesis <- match.arg(
    hypothesis,
    c("two.sided.pointnull", "greater.pointnull", "greater.composite")
  )
  
  control <- resolve_control_bf_ttest(control)
  
  if (length(n) != 1L || !is.finite(n) || n < 2L || n != as.integer(n)) {
    stop("`n` must be one integer greater than or equal to 2.", call. = FALSE)
  }
  
  if (length(k1) != 1L || !is.finite(k1) || k1 <= 0 || k1 >= 1) {
    stop("`k1` must be a scalar in (0, 1).", call. = FALSE)
  }
  
  if (length(k0) != 1L || !is.finite(k0) || k0 <= 1) {
    stop("`k0` must be a scalar greater than 1.", call. = FALSE)
  }
  
  if (length(effect_null) != 1L || !is.finite(effect_null)) {
    stop("`effect_null` must be one finite numeric value.", call. = FALSE)
  }
  
  if (
    hypothesis %in% c("two.sided.pointnull", "greater.pointnull") &&
    !isTRUE(all.equal(effect_null, 0))
  ) {
    stop(
      "For point-null hypotheses, the current Bayes-factor implementation assumes `effect_null = 0`.",
      call. = FALSE
    )
  }
  
  freq_power <- .normalize_freq_power_bf_ttest(freq_power)
  
  n <- as.integer(n)
  df <- n - 1L
  n_eff <- n
  
  if (identical(hypothesis, "two.sided.pointnull")) {
    resolved_analysis_prior_h0 <- make_point_prior(0)
    resolved_analysis_prior_h1 <- analysis_prior %||%
      get("make_analysis_prior_two_sided", envir = asNamespace("bfbin2arm"))(analysis_prior_sd)
  } else if (identical(hypothesis, "greater.pointnull")) {
    resolved_analysis_prior_h0 <- make_point_prior(0)
    resolved_analysis_prior_h1 <- analysis_prior %||%
      get("make_analysis_prior_greater", envir = asNamespace("bfbin2arm"))(analysis_prior_sd)
  } else {
    resolved_analysis_prior_h0 <- analysis_prior_h0 %||%
      get("make_analysis_prior_less_equal", envir = asNamespace("bfbin2arm"))(analysis_prior_sd)
    
    resolved_analysis_prior_h1 <- analysis_prior_h1 %||%
      get("make_analysis_prior_greater", envir = asNamespace("bfbin2arm"))(analysis_prior_sd)
  }
  
  .validate_prior_bf_ttest(resolved_analysis_prior_h1, allow_point = FALSE)
  
  if (identical(hypothesis, "greater.composite")) {
    .validate_prior_bf_ttest(resolved_analysis_prior_h0, allow_point = FALSE)
  }
  
  log_bf_fun <- make_log_bf01_for_hypothesis(
    hypothesis = hypothesis,
    df = df,
    n_eff = n_eff,
    prior_h0 = if (identical(hypothesis, "greater.composite")) resolved_analysis_prior_h0 else NULL,
    prior_h1 = resolved_analysis_prior_h1,
    rel.tol = control$rel_tol_m1,
    subdivisions = control$subdivisions_m1
  )
  
  symmetric_two_sided <- identical(hypothesis, "two.sided.pointnull")
  
  root_search_k1 <- .find_bracketing_roots_adaptive(
    log_bf_fun = log_bf_fun,
    threshold = k1,
    initial_t_max = control$t_max,
    t_max_cap = control$t_max_cap,
    grid_n = control$grid_n,
    root_tol = control$root_tol,
    symmetric = symmetric_two_sided
  )
  
  roots_k1 <- root_search_k1$roots
  t_max_used <- root_search_k1$t_max_used
  bf_grid <- root_search_k1$bf_grid
  
  if (identical(hypothesis, "two.sided.pointnull")) {
    if (length(roots_k1) != 2L) {
      return(list(
        ok = FALSE,
        reason = paste0(
          "Expected two symmetric k1 roots for `two.sided.pointnull`; found ",
          length(roots_k1), "."
        ),
        n = n,
        df = df,
        n_eff = n_eff,
        roots_k1 = roots_k1,
        roots_k0 = numeric(0),
        analysis_prior = resolved_analysis_prior_h1,
        analysis_prior_h0 = resolved_analysis_prior_h0,
        analysis_prior_h1 = resolved_analysis_prior_h1,
        t_max_used = t_max_used
      ))
    }
    
    h1_intervals <- matrix(
      c(-Inf, roots_k1[1L], roots_k1[2L], Inf),
      ncol = 2L,
      byrow = TRUE,
      dimnames = list(NULL, c("lower", "upper"))
    )
    c1 <- roots_k1[1L]
    c2 <- roots_k1[2L]
  } else if (identical(hypothesis, "greater.pointnull")) {
    positive_roots <- roots_k1[roots_k1 > 0]
    
    if (length(positive_roots) != 1L) {
      return(list(
        ok = FALSE,
        reason = paste0(
          "Expected exactly one positive k1 root for `greater.pointnull`; found ",
          length(positive_roots), "."
        ),
        n = n,
        df = df,
        n_eff = n_eff,
        roots_k1 = roots_k1,
        roots_k0 = numeric(0),
        analysis_prior = resolved_analysis_prior_h1,
        analysis_prior_h0 = resolved_analysis_prior_h0,
        analysis_prior_h1 = resolved_analysis_prior_h1,
        t_max_used = t_max_used
      ))
    }
    
    c1 <- NA_real_
    c2 <- positive_roots
    h1_intervals <- matrix(
      c(c2, Inf),
      nrow = 1L,
      dimnames = list(NULL, c("lower", "upper"))
    )
  } else {
    h1_intervals <- .classify_decision_intervals_bf_ttest(
      log_bf_fun = log_bf_fun,
      threshold = k1,
      roots = roots_k1,
      lower = -Inf,
      upper = Inf,
      keep_when = "<="
    )
    
    if (!nrow(h1_intervals)) {
      return(list(
        ok = FALSE,
        reason = "No evidence-for-H1 interval was found for BF01(t) <= k1.",
        n = n,
        df = df,
        n_eff = n_eff,
        roots_k1 = roots_k1,
        roots_k0 = numeric(0),
        h1_region = h1_intervals,
        h0_region = NULL,
        analysis_prior = resolved_analysis_prior_h1,
        analysis_prior_h0 = resolved_analysis_prior_h0,
        analysis_prior_h1 = resolved_analysis_prior_h1,
        t_max_used = t_max_used
      ))
    }
    
    c1 <- NA_real_
    c2 <- NA_real_
  }
  
  p_h1 <- function(delta) {
    probability_over_intervals(
      intervals = h1_intervals,
      delta = delta,
      df = df,
      n_eff = n_eff
    )
  }
  
  roots_k0 <- numeric(0)
  h0_intervals <- NULL
  a0 <- NA_real_
  b0 <- NA_real_
  
  p_h0 <- function(delta) {
    rep(NA_real_, length(delta))
  }
  
  if (isTRUE(want_pce)) {
    if (symmetric_two_sided) {
      roots_pos_k0 <- .find_roots_from_grid(
        grid_t = bf_grid$t,
        grid_log_bf = bf_grid$log_bf01,
        threshold = k0,
        log_bf_fun = log_bf_fun,
        root_tol = control$root_tol
      )
      roots_pos_k0 <- roots_pos_k0[roots_pos_k0 > 0]
      roots_k0 <- sort(c(-rev(roots_pos_k0), roots_pos_k0))
    } else {
      roots_k0 <- .find_roots_from_grid(
        grid_t = bf_grid$t,
        grid_log_bf = bf_grid$log_bf01,
        threshold = k0,
        log_bf_fun = log_bf_fun,
        root_tol = control$root_tol
      )
    }
    
    h0_intervals <- .classify_decision_intervals_bf_ttest(
      log_bf_fun = log_bf_fun,
      threshold = k0,
      roots = roots_k0,
      lower = -Inf,
      upper = Inf,
      keep_when = ">="
    )
    
    p_h0 <- function(delta) {
      probability_over_intervals(
        intervals = h0_intervals,
        delta = delta,
        df = df,
        n_eff = n_eff
      )
    }
    
    if (nrow(h0_intervals) == 1L) {
      a0 <- h0_intervals[1L, "lower"]
      b0 <- h0_intervals[1L, "upper"]
    }
  }
  
  power_effects <- vapply(freq_power, `[[`, numeric(1), "effect")
  
  freq_power_values <- if (length(power_effects)) {
    stats::setNames(
      vapply(power_effects, p_h1, numeric(1)),
      names(freq_power)
    )
  } else {
    numeric(0)
  }
  
  freq_t1e <- p_h1(effect_null)
  
  bayes_t1e <- integrate_over_design_prior(
    fun_delta = p_h1,
    prior_spec = design_prior_h0
  )
  
  bayes_power <- integrate_over_design_prior(
    fun_delta = p_h1,
    prior_spec = design_prior_h1
  )
  
  pce_h0 <- if (isTRUE(want_pce)) {
    integrate_over_design_prior(
      fun_delta = p_h0,
      prior_spec = design_prior_h0
    )
  } else {
    NA_real_
  }
  
  list(
    ok = TRUE,
    reason = NULL,
    n = n,
    df = df,
    n_eff = n_eff,
    k1 = k1,
    k0 = k0,
    c1 = c1,
    c2 = c2,
    a0 = a0,
    b0 = b0,
    roots_k1 = roots_k1,
    roots_k0 = roots_k0,
    h1_region = h1_intervals,
    h0_region = h0_intervals,
    freq_t1e = as.numeric(freq_t1e),
    bayes_t1e = as.numeric(bayes_t1e),
    freq_power = freq_power_values,
    bayes_power = as.numeric(bayes_power),
    pce_h0 = as.numeric(pce_h0),
    conditional_probability_h1 = p_h1,
    conditional_probability_h0 = p_h0,
    log_bf01_fun = log_bf_fun,
    bf_grid = bf_grid,
    t_max_used = t_max_used,
    analysis_prior = resolved_analysis_prior_h1,
    analysis_prior_h0 = resolved_analysis_prior_h0,
    analysis_prior_h1 = resolved_analysis_prior_h1
  )
}

verify_candidate_one_sample <- function(
    n,
    k1,
    k0,
    analysis_prior = NULL,
    effect_null = 0,
    freq_power = list(),
    want_pce = FALSE,
    control = NULL,
    hypothesis = "two.sided.pointnull",
    design_prior_h0 = list(type = "point", effect = 0),
    design_prior_h1 = list(type = "point", effect = 0.4),
    analysis_prior_h0 = NULL,
    analysis_prior_h1 = NULL,
    analysis_prior_sd = 0.5
) {
  control <- resolve_control_bf_ttest(control)
  
  ctl2 <- control
  ctl2$rel_tol_m1 <- control$verify_rel_tol_m1
  ctl2$subdivisions_m1 <- control$verify_subdivisions_m1
  ctl2$grid_n <- control$verify_grid_n
  ctl2$t_max <- control$verify_t_max
  ctl2$root_tol <- control$verify_root_tol
  
  compute_oc_one_sample(
    n = n,
    k1 = k1,
    k0 = k0,
    analysis_prior = analysis_prior,
    effect_null = effect_null,
    freq_power = freq_power,
    want_pce = want_pce,
    control = ctl2,
    hypothesis = hypothesis,
    design_prior_h0 = design_prior_h0,
    design_prior_h1 = design_prior_h1,
    analysis_prior_h0 = analysis_prior_h0,
    analysis_prior_h1 = analysis_prior_h1,
    analysis_prior_sd = analysis_prior_sd
  )
}

# -----------------------------------------------------------------------------
# Calibration
# -----------------------------------------------------------------------------

#' Calibrate a Bayesian t-test design
#'
#' Searches an integer sample-size range and returns the smallest design that
#' satisfies the requested frequentist and Bayesian operating-characteristic
#' constraints for a one-sample or paired Bayesian t-test.
#'
#' @param design Design type. One of `"one.sample"`, `"paired"`, or
#'   `"two.sample.pooled"`. The current implementation supports
#'   `"one.sample"` and `"paired"`.
#' @param hypothesis Hypothesis configuration. One of
#'   `"two.sided.pointnull"`, `"greater.pointnull"`, or
#'   `"greater.composite"`.
#' @param n_min Minimum candidate sample size.
#' @param n_max Maximum candidate sample size.
#' @param k1 Evidence threshold for evidence in favour of the alternative:
#'   evidence for H1 is declared when `BF01 <= k1`.
#' @param k0 Evidence threshold for evidence in favour of the null:
#'   evidence for H0 is declared when `BF01 >= k0`.
#' @param analysis_prior_sd Standard deviation of the default analysis prior.
#'   Used when a relevant explicit analysis prior is not supplied.
#' @param analysis_prior Analysis prior for point-null hypotheses. If `NULL`,
#'   a default prior determined by `hypothesis` and `analysis_prior_sd` is used.
#' @param freq_t1e Optional upper bound on the frequentist type-I error,
#'   calculated under the null effect.
#' @param freq_power Optional list of frequentist power specifications. Each
#'   element must contain `effect` and may contain a target `target`.
#' @param bayes_power Optional lower bound on Bayesian power, averaged under
#'   `design_prior_h1`.
#' @param bayes_t1e Optional upper bound on the Bayesian type-I error, averaged
#'   under `design_prior_h0`.
#' @param pce_h0 Optional lower bound on the probability of compelling evidence
#'   for H0, averaged under `design_prior_h0`.
#' @param design_prior_h1 Design-prior specification under H1.
#' @param design_prior_h0 Design-prior specification under H0.
#' @param analysis_prior_h0 Analysis prior under H0 for
#'   `"greater.composite"` hypotheses.
#' @param analysis_prior_h1 Analysis prior under H1 for
#'   `"greater.composite"` hypotheses.
#' @param control Optional named list controlling numerical integration,
#'   root-finding, search strategy, verification, and progress reporting.
#'
#' @return An object of class `"bf_ttest_design"`. On success, it contains the
#'   selected sample size (`n`), decision regions, operating characteristics,
#'   a constraint table, search results, and plot data. If no feasible design is
#'   found, it returns a `"bf_ttest_design"` object with `success = FALSE`.
#'
#' @examples
#' \dontrun{
#' design <- calibrate_bf_ttest(
#'   design = "one.sample",
#'   hypothesis = "greater.composite",
#'   n_min = 20,
#'   n_max = 120,
#'   k1 = 1 / 20,
#'   k0 = 3,
#'   freq_t1e = 0.05,
#'   freq_power = list(main = list(effect = 0.40, target = 0.80)),
#'   design_prior_h0 = list(
#'     type = "truncnormal",
#'     mean = -0.10,
#'     sd = 0.25,
#'     lower = -Inf,
#'     upper = 0
#'   ),
#'   design_prior_h1 = list(
#'     type = "truncnormal",
#'     mean = 0.40,
#'     sd = 0.20,
#'     lower = 0,
#'     upper = Inf
#'   ),
#'   control = list(progress = FALSE)
#' )
#' }
#'
#' @export
calibrate_bf_ttest <- function(
    design = c("one.sample", "paired", "two.sample.pooled"),
    hypothesis = c("two.sided.pointnull", "greater.pointnull", "greater.composite"),
    n_min,
    n_max,
    k1 = 1 / 10,
    k0 = 10,
    analysis_prior_sd = 0.5,
    analysis_prior = NULL,
    freq_t1e = NULL,
    freq_power = NULL,
    bayes_power = NULL,
    bayes_t1e = NULL,
    pce_h0 = NULL,
    design_prior_h1 = list(type = "point", effect = 0.4),
    design_prior_h0 = list(type = "point", effect = 0),
    analysis_prior_h0 = NULL,
    analysis_prior_h1 = NULL,
    control = NULL
) {
  design <- match.arg(design)
  hypothesis <- match.arg(hypothesis)
  control <- resolve_control_bf_ttest(control)
  freq_power <- .normalize_freq_power_bf_ttest(freq_power)
  
  if (!design %in% c("one.sample", "paired")) {
    stop(
      "This accelerated implementation currently supports `one.sample` and `paired` designs only.",
      call. = FALSE
    )
  }
  
  if (length(n_min) != 1L || length(n_max) != 1L ||
      !is.finite(n_min) || !is.finite(n_max) || n_min < 2 || n_max < n_min) {
    stop("`n_min` and `n_max` must be finite integers with 2 <= n_min <= n_max.", call. = FALSE)
  }
  
  n_min <- as.integer(n_min)
  n_max <- as.integer(n_max)
  want_pce <- !is.null(pce_h0)
  effect_null <- if (design_prior_h0$type %||% NA_character_ == "point") {
    as.numeric(design_prior_h0$effect)
  } else {
    0
  }
  
  coarse_grid <- if (identical(control$n_search, "two_stage")) {
    unique(c(seq.int(n_min, n_max, by = max(1L, control$n_step_stage1)), n_max))
  } else {
    seq.int(n_min, n_max, by = 1L)
  }
  
  coarse_res <- vector("list", length(coarse_grid))
  names(coarse_res) <- as.character(coarse_grid)
  
  make_pb <- function(total, label) {
    if (!isTRUE(control$progress) || total <= 0L) {
      return(NULL)
    }
    message(label)
    utils::txtProgressBar(min = 0, max = total, style = control$progress_style)
  }
  
  tick_pb <- function(pb, value) {
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, value)
    }
    invisible(NULL)
  }
  
  close_pb <- function(pb) {
    if (!is.null(pb)) {
      close(pb)
    }
    invisible(NULL)
  }
  
  feasible_hit <- NA_integer_
  pb1 <- make_pb(length(coarse_grid), "Coarse search:")
  on.exit(close_pb(pb1), add = TRUE)
  
  for (i in seq_along(coarse_grid)) {
    n <- coarse_grid[i]
    
    oc <- compute_oc_one_sample(
      n = n,
      k1 = k1,
      k0 = k0,
      analysis_prior = analysis_prior,
      effect_null = effect_null,
      freq_power = freq_power,
      want_pce = want_pce,
      control = control,
      hypothesis = hypothesis,
      design_prior_h0 = design_prior_h0,
      design_prior_h1 = design_prior_h1,
      analysis_prior_h0 = analysis_prior_h0,
      analysis_prior_h1 = analysis_prior_h1,
      analysis_prior_sd = analysis_prior_sd
    )
    
    coarse_res[[i]] <- oc
    tick_pb(pb1, i)
    
    if (isTRUE(oc$ok) && .meets_constraints_bf_ttest(
      oc,
      freq_t1e = freq_t1e,
      freq_power = freq_power,
      bayes_power = bayes_power,
      bayes_t1e = bayes_t1e,
      pce_h0 = pce_h0
    )) {
      feasible_hit <- n
      tick_pb(pb1, length(coarse_grid))
      break
    }
  }
  close_pb(pb1)
  
  if (is.na(feasible_hit)) {
    best <- .best_candidate_bf_ttest(
      results = coarse_res,
      grid = coarse_grid,
      freq_t1e = freq_t1e,
      freq_power = freq_power,
      bayes_power = bayes_power,
      bayes_t1e = bayes_t1e,
      pce_h0 = pce_h0
    )
    
    return(structure(list(
      success = FALSE,
      message = "No feasible design found in search range.",
      design_type = design,
      test_type = hypothesis,
      n = if (!is.null(best)) best$n else NA_integer_,
      df = if (!is.null(best)) best$oc$df else NA_real_,
      freq_t1e = if (!is.null(best)) best$oc$freq_t1e else NA_real_,
      bayes_t1e = if (!is.null(best)) best$oc$bayes_t1e else NA_real_,
      bayes_power = if (!is.null(best)) best$oc$bayes_power else NA_real_,
      pce_h0 = if (!is.null(best)) best$oc$pce_h0 else NA_real_,
      freq_power = if (!is.null(best)) best$oc$freq_power else NULL,
      constraint_table = if (!is.null(best)) best$constraint_table else NULL,
      best_candidate = best,
      inputs = list(
        n_min = n_min,
        n_max = n_max,
        k1 = k1,
        k0 = k0,
        freq_t1e_target = freq_t1e,
        freq_power_target = freq_power,
        bayes_power_target = bayes_power,
        bayes_t1e_target = bayes_t1e,
        pce_h0_target = pce_h0,
        analysis_prior_sd = analysis_prior_sd,
        analysis_prior = analysis_prior,
        design_prior_h1 = design_prior_h1,
        design_prior_h0 = design_prior_h0,
        analysis_prior_h0 = analysis_prior_h0,
        analysis_prior_h1 = analysis_prior_h1
      ),
      coarse_grid = coarse_grid,
      coarse_results = coarse_res,
      control = control
    ), class = "bf_ttest_design"))
  }
  
  refine_grid <- if (identical(control$n_search, "two_stage")) {
    seq.int(max(n_min, feasible_hit - control$n_step_stage1 + 1L), feasible_hit, by = 1L)
  } else {
    feasible_hit
  }
  
  refine_res <- vector("list", length(refine_grid))
  names(refine_res) <- as.character(refine_grid)
  
  n_star <- NA_integer_
  final_oc <- NULL
  
  pb2 <- make_pb(length(refine_grid), "Refinement:")
  on.exit(close_pb(pb2), add = TRUE)
  
  for (i in seq_along(refine_grid)) {
    n <- refine_grid[i]
    
    oc <- compute_oc_one_sample(
      n = n,
      k1 = k1,
      k0 = k0,
      analysis_prior = analysis_prior,
      effect_null = effect_null,
      freq_power = freq_power,
      want_pce = want_pce,
      control = control,
      hypothesis = hypothesis,
      design_prior_h0 = design_prior_h0,
      design_prior_h1 = design_prior_h1,
      analysis_prior_h0 = analysis_prior_h0,
      analysis_prior_h1 = analysis_prior_h1,
      analysis_prior_sd = analysis_prior_sd
    )
    
    refine_res[[i]] <- oc
    tick_pb(pb2, i)
    
    if (isTRUE(oc$ok) && .meets_constraints_bf_ttest(
      oc,
      freq_t1e = freq_t1e,
      freq_power = freq_power,
      bayes_power = bayes_power,
      bayes_t1e = bayes_t1e,
      pce_h0 = pce_h0
    )) {
      n_star <- n
      final_oc <- oc
      tick_pb(pb2, length(refine_grid))
      break
    }
  }
  close_pb(pb2)
  
  if (is.na(n_star)) {
    best <- .best_candidate_bf_ttest(
      results = refine_res,
      grid = refine_grid,
      freq_t1e = freq_t1e,
      freq_power = freq_power,
      bayes_power = bayes_power,
      bayes_t1e = bayes_t1e,
      pce_h0 = pce_h0
    )
    
    return(structure(list(
      success = FALSE,
      message = paste(
        "Coarse scan found a candidate region, but no feasible integer design survived refinement."
      ),
      design_type = design,
      test_type = hypothesis,
      n = if (!is.null(best)) best$n else NA_integer_,
      df = if (!is.null(best)) best$oc$df else NA_real_,
      freq_t1e = if (!is.null(best)) best$oc$freq_t1e else NA_real_,
      bayes_t1e = if (!is.null(best)) best$oc$bayes_t1e else NA_real_,
      bayes_power = if (!is.null(best)) best$oc$bayes_power else NA_real_,
      pce_h0 = if (!is.null(best)) best$oc$pce_h0 else NA_real_,
      freq_power = if (!is.null(best)) best$oc$freq_power else NULL,
      constraint_table = if (!is.null(best)) best$constraint_table else NULL,
      best_candidate = best,
      inputs = list(
        n_min = n_min,
        n_max = n_max,
        k1 = k1,
        k0 = k0,
        freq_t1e_target = freq_t1e,
        freq_power_target = freq_power,
        bayes_power_target = bayes_power,
        bayes_t1e_target = bayes_t1e,
        pce_h0_target = pce_h0,
        analysis_prior_sd = analysis_prior_sd,
        analysis_prior = analysis_prior,
        design_prior_h1 = design_prior_h1,
        design_prior_h0 = design_prior_h0,
        analysis_prior_h0 = analysis_prior_h0,
        analysis_prior_h1 = analysis_prior_h1
      ),
      coarse_grid = coarse_grid,
      coarse_results = coarse_res,
      refine_grid = refine_grid,
      refine_results = refine_res,
      control = control
    ), class = "bf_ttest_design"))
  }
  
  verified <- NULL
  if (isTRUE(control$verify)) {
    pb3 <- make_pb(1L, "Verification:")
    on.exit(close_pb(pb3), add = TRUE)
    
    verified <- verify_candidate_one_sample(
      n = n_star,
      k1 = k1,
      k0 = k0,
      analysis_prior = analysis_prior,
      effect_null = effect_null,
      freq_power = freq_power,
      want_pce = want_pce,
      control = control,
      hypothesis = hypothesis,
      design_prior_h0 = design_prior_h0,
      design_prior_h1 = design_prior_h1,
      analysis_prior_h0 = analysis_prior_h0,
      analysis_prior_h1 = analysis_prior_h1,
      analysis_prior_sd = analysis_prior_sd
    )
    
    tick_pb(pb3, 1L)
    close_pb(pb3)
    
    if (!isTRUE(verified$ok) || !.meets_constraints_bf_ttest(
      verified,
      freq_t1e = freq_t1e,
      freq_power = freq_power,
      bayes_power = bayes_power,
      bayes_t1e = bayes_t1e,
      pce_h0 = pce_h0,
      tol = 1e-10
    )) {
      warning(
        "Verified candidate did not satisfy all constraints at stricter tolerances. Consider increasing `n_max` or tightening exploratory controls.",
        call. = FALSE
      )
    } else {
      final_oc <- verified
    }
  }
  
  constraint_table <- .constraint_table_bf_ttest(
    oc = final_oc,
    freq_t1e = freq_t1e,
    freq_power = freq_power,
    bayes_power = bayes_power,
    bayes_t1e = bayes_t1e,
    pce_h0 = pce_h0
  )
  
  decision_regions <- list(
    evidence_h1 = list(
      intervals = final_oc$h1_region,
      roots = final_oc$roots_k1,
      threshold = k1,
      rule = "BF01 <= k1"
    ),
    evidence_h0 = list(
      intervals = final_oc$h0_region,
      roots = final_oc$roots_k0,
      threshold = k0,
      rule = "BF01 >= k0"
    )
  )
  
  pb4 <- make_pb(3L, "Building plot data:")
  on.exit(close_pb(pb4), add = TRUE)
  
  root_labels <- c(
    rep("k1", length(final_oc$roots_k1)),
    rep("k0", length(final_oc$roots_k0))
  )
  tick_pb(pb4, 1L)
  
  t_grid <- seq(
    -control$verify_t_max,
    control$verify_t_max,
    length.out = max(401L, control$verify_grid_n)
  )
  
  bf_fun_for_plot <- make_log_bf01_for_hypothesis(
    hypothesis = hypothesis,
    df = final_oc$df,
    n_eff = final_oc$n_eff,
    prior_h0 = if (identical(hypothesis, "greater.composite")) final_oc$analysis_prior_h0 else NULL,
    prior_h1 = final_oc$analysis_prior_h1,
    rel.tol = control$verify_rel_tol_m1,
    subdivisions = control$verify_subdivisions_m1
  )
  tick_pb(pb4, 2L)
  
  log_bf_vals <- vapply(t_grid, bf_fun_for_plot, numeric(1))
  
  plot_data <- list(
    bf_curve = data.frame(
      t = t_grid,
      log_bf01 = log_bf_vals,
      bf01 = exp(log_bf_vals),
      stringsAsFactors = FALSE
    ),
    t_decision_boundaries = data.frame(
      threshold = root_labels,
      t = c(final_oc$roots_k1, final_oc$roots_k0),
      stringsAsFactors = FALSE
    ),
    power_at_effect = data.frame(
      effect = c(effect_null, vapply(freq_power, `[[`, numeric(1), "effect")),
      probability_evidence_h1 = c(final_oc$freq_t1e, unname(final_oc$freq_power)),
      role = c("null", rep("alternative", length(freq_power))),
      stringsAsFactors = FALSE
    )
  )
  
  tick_pb(pb4, 3L)
  close_pb(pb4)
  
  out <- list(
    success = TRUE,
    design_type = design,
    test_type = hypothesis,
    n = n_star,
    n_total = n_star,
    n1 = NULL,
    n2 = NULL,
    df = final_oc$df,
    k1 = k1,
    k0 = k0,
    analysis_prior = final_oc$analysis_prior,
    analysis_prior_h0 = final_oc$analysis_prior_h0,
    analysis_prior_h1 = final_oc$analysis_prior_h1,
    design_priors = list(h0 = design_prior_h0, h1 = design_prior_h1),
    freq_t1e = final_oc$freq_t1e,
    bayes_t1e = final_oc$bayes_t1e,
    freq_power = final_oc$freq_power,
    power_at_effect = final_oc$freq_power,
    bayes_power = final_oc$bayes_power,
    pce_h0 = final_oc$pce_h0,
    roots_k1 = final_oc$roots_k1,
    roots_k0 = final_oc$roots_k0,
    c1 = final_oc$c1,
    c2 = final_oc$c2,
    decision_regions = decision_regions,
    operating_characteristics = final_oc,
    constraint_table = constraint_table,
    calibration = list(
      all_constraints_satisfied = all(constraint_table$satisfied),
      selected_n = n_star,
      objective = "minimum sample size",
      constraint_table = constraint_table
    ),
    plot_data = plot_data,
    inputs = list(
      n_min = n_min,
      n_max = n_max,
      analysis_prior_sd = analysis_prior_sd,
      analysis_prior = analysis_prior,
      freq_t1e_target = freq_t1e,
      freq_power_target = freq_power,
      bayes_power_target = bayes_power,
      bayes_t1e_target = bayes_t1e,
      pce_h0_target = pce_h0,
      design_prior_h1 = design_prior_h1,
      design_prior_h0 = design_prior_h0,
      analysis_prior_h0 = analysis_prior_h0,
      analysis_prior_h1 = analysis_prior_h1,
      k1 = k1,
      k0 = k0
    ),
    coarse_grid = coarse_grid,
    coarse_results = coarse_res,
    refine_grid = refine_grid,
    refine_results = refine_res,
    verified = verified,
    control = control
  )
  
  class(out) <- "bf_ttest_design"
  out
}