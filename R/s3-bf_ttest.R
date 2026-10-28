# =============================================================================
# s3-bf_ttest.R
# S3 methods for calibrated Bayesian t-test designs
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

.required_plot_components_bf_ttest <- c(
  "bf_curve",
  "t_decision_boundaries",
  "power_at_effect"
)

.has_complete_plot_data_bf_ttest <- function(x) {
  is.list(x$plot_data) &&
    all(.required_plot_components_bf_ttest %in% names(x$plot_data)) &&
    !is.null(x$plot_data$bf_curve) &&
    !is.null(x$plot_data$t_decision_boundaries) &&
    !is.null(x$plot_data$power_at_effect)
}

.get_plot_data_bf_ttest <- function(x, name) {
  plot_data <- x$plot_data %||% NULL
  if (!is.list(plot_data) || is.null(plot_data[[name]])) {
    stop("Required plot data component missing: `plot_data$", name, "`.", call. = FALSE)
  }
  plot_data[[name]]
}

.rebuild_plot_data_bf_ttest <- function(x) {
  oc <- x$operating_characteristics %||% NULL
  if (is.null(oc) || !is.list(oc)) return(NULL)
  
  control <- x$control %||% resolve_control_bf_ttest()
  hypothesis <- x$test_type %||% NULL
  k1 <- x$k1 %||% x$inputs$k1 %||% NULL
  k0 <- x$k0 %||% x$inputs$k0 %||% NULL
  
  if (is.null(hypothesis) || is.null(k1) || is.null(k0)) return(NULL)
  if (is.null(oc$df) || is.null(oc$n_eff)) return(NULL)
  
  log_bf_fun <- oc$log_bf01_fun %||% NULL
  
  if (is.null(log_bf_fun)) {
    prior_h1 <- x$analysis_prior_h1 %||% x$analysis_prior %||% oc$analysis_prior_h1 %||% oc$analysis_prior
    prior_h0 <- x$analysis_prior_h0 %||% oc$analysis_prior_h0 %||% NULL
    
    if (is.null(prior_h1)) return(NULL)
    
    log_bf_fun <- make_log_bf01_for_hypothesis(
      hypothesis = hypothesis,
      df = oc$df,
      n_eff = oc$n_eff,
      prior_h1 = prior_h1,
      prior_h0 = if (identical(hypothesis, "greater.composite")) prior_h0 else NULL,
      rel.tol = control$verify_rel_tol_m1 %||% control$rel_tol_m1,
      subdivisions = control$verify_subdivisions_m1 %||% control$subdivisions_m1
    )
  }
  
  t_max <- control$verify_t_max %||% control$t_max %||% 12
  grid_n <- max(401L, control$verify_grid_n %||% control$grid_n %||% 801L)
  t_grid <- seq(-t_max, t_max, length.out = grid_n)
  
  log_bf_vals <- vapply(t_grid, log_bf_fun, numeric(1))
  
  roots_k1 <- oc$roots_k1 %||% numeric(0)
  roots_k0 <- oc$roots_k0 %||% numeric(0)
  root_labels <- c(rep("k1", length(roots_k1)), rep("k0", length(roots_k0)))
  root_values <- c(roots_k1, roots_k0)
  
  freq_power <- oc$freq_power %||% numeric(0)
  freq_power_effects <- x$inputs$freq_power_target %||% list()
  
  alt_effects <- if (length(freq_power_effects)) {
    vapply(freq_power_effects, function(z) as.numeric(z$effect), numeric(1))
  } else if (length(freq_power)) {
    seq_along(freq_power) * NA_real_
  } else {
    numeric(0)
  }
  
  effect_null <- x$design_priors$h0$effect %||% 0
  
  power_df <- data.frame(
    effect = c(effect_null, alt_effects),
    probability_evidence_h1 = c(oc$freq_t1e %||% NA_real_, unname(freq_power)),
    role = c("null", rep("alternative", length(freq_power))),
    stringsAsFactors = FALSE
  )
  
  list(
    bf_curve = data.frame(
      t = t_grid,
      log_bf01 = log_bf_vals,
      bf01 = exp(log_bf_vals),
      stringsAsFactors = FALSE
    ),
    t_decision_boundaries = data.frame(
      threshold = root_labels,
      t = root_values,
      stringsAsFactors = FALSE
    ),
    power_at_effect = power_df
  )
}

.ensure_plot_data_bf_ttest <- function(x) {
  if (.has_complete_plot_data_bf_ttest(x)) return(x)
  
  rebuilt <- .rebuild_plot_data_bf_ttest(x)
  if (!is.null(rebuilt)) {
    x$plot_data <- rebuilt
  }
  
  x
}

#' @export
print.bf_ttest_design <- function(x, digits = 4L, ...) {
  k1 <- x$k1 %||% x$inputs$k1 %||% NA_real_
  k0 <- x$k0 %||% x$inputs$k0 %||% NA_real_
  
  cat("Calibrated Bayesian t-test design\n")
  cat("---------------------------------\n")
  
  if (isFALSE(x$success)) {
    cat("Status: failed\n")
    cat("Message:", x$message %||% "No message available.", "\n")
    cat("Design:", x$design_type %||% NA_character_, "\n")
    cat("Hypothesis:", x$test_type %||% NA_character_, "\n")
    
    if (!is.null(x$n) && is.finite(x$n)) {
      cat("Best evaluated n*:", x$n)
      if (!is.null(x$df) && is.finite(x$df)) {
        cat(" (df =", x$df, ")")
      }
      cat("\n")
    }
    
    if (!is.null(x$constraint_table) && nrow(x$constraint_table) > 0L) {
      cat("Nearest evaluated candidate constraints:\n")
      print(x$constraint_table, row.names = FALSE)
    }
    
    return(invisible(x))
  }
  
  cat("Design:", x$design_type %||% NA_character_, "\n")
  cat("Hypothesis:", x$test_type %||% NA_character_, "\n")
  cat("n*:", x$n %||% NA_integer_)
  if (!is.null(x$df) && is.finite(x$df)) {
    cat(" (df =", x$df, ")")
  }
  cat("\n")
  cat("Thresholds: k1 =", format(k1, digits = digits), ", k0 =", format(k0, digits = digits), "\n")
  
  cat("Analysis prior under H0:\n")
  print(x$analysis_prior_h0)
  
  cat("Analysis prior under H1:\n")
  print(x$analysis_prior_h1)
  
  cat("Evidence for H1 region(s) on t-scale:\n")
  if (!is.null(x$decision_regions$evidence_h1$intervals)) {
    print(x$decision_regions$evidence_h1$intervals)
  } else if (!is.null(x$operating_characteristics$h1_region)) {
    print(x$operating_characteristics$h1_region)
  } else {
    cat("none\n")
  }
  
  cat("Evidence for H0 region(s) on t-scale:\n")
  if (!is.null(x$decision_regions$evidence_h0$intervals)) {
    print(x$decision_regions$evidence_h0$intervals)
  } else if (!is.null(x$operating_characteristics$h0_region)) {
    print(x$operating_characteristics$h0_region)
  } else {
    cat("none\n")
  }
  
  if (!is.null(x$constraint_table) && nrow(x$constraint_table) > 0L) {
    cat("\nConstraints:\n")
    print(x$constraint_table, row.names = FALSE)
  }
  
  invisible(x)
}

#' @export
summary.bf_ttest_design <- function(object, ...) {
  print.bf_ttest_design(object, ...)
  invisible(object)
}

#' Plot method for calibrated Bayesian t-test designs
#'
#' @param x Object of class `bfttestdesign`.
#' @param which Character vector choosing panels among
#'   `"decision"`, `"analysisprior"`, `"designprior"`, `"effect"`,
#'   and `"ncurve"`.
#' @param annotate Logical; add global title/subtitle.
#' @param prior_points Number of grid points for prior density curves.
#' @param ... Ignored.
#'
#' @return A `ggplot` object or a patchwork object.
#' @export
plot.bf_ttest_design <- function(
    x,
    which = c("decision", "analysis_prior", "design_prior", "n_curve"),
    annotate = TRUE,
    prior_points = 801L,
    ...
) {
  if (!inherits(x, "bf_ttest_design")) {
    stop("`x` must inherit from class `bf_ttest_design`.", call. = FALSE)
  }
  
  if (!isTRUE(x$success)) {
    stop("Cannot plot a failed calibration result.", call. = FALSE)
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required.", call. = FALSE)
  }
  
  if (length(which) > 1L && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package `patchwork` is required for multi-panel plots.", call. = FALSE)
  }
  
  x <- .ensure_plot_data_bf_ttest(x)
  
  which <- match.arg(
    which,
    choices = c("decision", "analysis_prior", "design_prior", "effect", "n_curve"),
    several.ok = TRUE
  )
  
  theme_base <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  
  get_plot_data <- function(name) {
    .get_plot_data_bf_ttest(x, name)
  }
  
  has_intervals <- function(z) {
    is.matrix(z) && nrow(z) > 0L && ncol(z) == 2L
  }
  
  intervals_to_df <- function(intervals, region) {
    if (!has_intervals(intervals)) {
      return(data.frame(
        lower = numeric(0),
        upper = numeric(0),
        region = character(0),
        stringsAsFactors = FALSE
      ))
    }
    
    data.frame(
      lower = intervals[, 1L],
      upper = intervals[, 2L],
      region = region,
      stringsAsFactors = FALSE
    )
  }
  
  constraint_curve_df <- function() {
    extract_metric <- function(res_list, field) {
      vapply(res_list, function(z) {
        val <- z[[field]]
        if (is.null(val) || length(val) != 1L || !is.finite(val)) NA_real_ else val
      }, numeric(1))
    }
    
    if (is.null(x$coarse_results) || !length(x$coarse_results)) {
      return(data.frame())
    }
    
    n_vals <- as.integer(names(x$coarse_results))
    if (anyNA(n_vals)) n_vals <- seq_along(x$coarse_results)
    
    df <- data.frame(
      n = n_vals,
      freq_t1e = extract_metric(x$coarse_results, "freq_t1e"),
      bayes_t1e = extract_metric(x$coarse_results, "bayes_t1e"),
      bayes_power = extract_metric(x$coarse_results, "bayes_power"),
      pce_h0 = extract_metric(x$coarse_results, "pce_h0"),
      stringsAsFactors = FALSE
    )
    
    fp_targets <- x$inputs$freq_power_target
    if (!is.null(fp_targets) && length(fp_targets)) {
      for (nm in names(fp_targets)) {
        df[[paste0("freq_power:", nm)]] <- vapply(
          x$coarse_results,
          function(z) {
            val <- z$freq_power[[nm]]
            if (is.null(val) || length(val) != 1L || !is.finite(val)) NA_real_ else val
          },
          numeric(1)
        )
      }
    }
    
    keep <- setdiff(names(df), "n")
    keep <- keep[vapply(keep, function(nm) any(is.finite(df[[nm]])), logical(1))]
    
    if (!length(keep)) {
      return(data.frame())
    }
    
    out <- do.call(rbind, lapply(keep, function(nm) {
      data.frame(
        n = df$n,
        metric = nm,
        value = df[[nm]],
        stringsAsFactors = FALSE
      )
    }))
    rownames(out) <- NULL
    out
  }
  
  constraint_targets_df <- function() {
    rows <- list()
    
    add_row <- function(metric, target) {
      rows[[length(rows) + 1L]] <<- data.frame(
        metric = metric,
        target = target,
        stringsAsFactors = FALSE
      )
    }
    
    if (!is.null(x$inputs$freq_t1e_target)) add_row("freq_t1e", x$inputs$freq_t1e_target)
    if (!is.null(x$inputs$bayes_t1e_target)) add_row("bayes_t1e", x$inputs$bayes_t1e_target)
    if (!is.null(x$inputs$bayes_power_target)) add_row("bayes_power", x$inputs$bayes_power_target)
    if (!is.null(x$inputs$pce_h0_target)) add_row("pce_h0", x$inputs$pce_h0_target)
    
    fp_targets <- x$inputs$freq_power_target
    if (!is.null(fp_targets) && length(fp_targets)) {
      for (nm in names(fp_targets)) {
        tgt <- fp_targets[[nm]]$target
        if (!is.null(tgt)) add_row(paste0("freq_power:", nm), tgt)
      }
    }
    
    if (!length(rows)) {
      return(data.frame(metric = character(0), target = numeric(0)))
    }
    
    out <- do.call(rbind, rows)
    out[out$metric %in% unique(constraint_curve_df()$metric), , drop = FALSE]
  }
  
  prior_to_df <- function(prior, label, n = 801L) {
    if (is.null(prior)) {
      return(data.frame(
        delta = numeric(0),
        density = numeric(0),
        prior = character(0),
        prior_type = character(0),
        stringsAsFactors = FALSE
      ))
    }
    
    if (inherits(prior, "bf_ttest_point_prior")) {
      return(data.frame(
        delta = prior$effect,
        density = 1,
        prior = label,
        prior_type = "point",
        stringsAsFactors = FALSE
      ))
    }
    
    if (is.list(prior) && identical(prior$type, "point")) {
      return(data.frame(
        delta = prior$effect,
        density = 1,
        prior = label,
        prior_type = "point",
        stringsAsFactors = FALSE
      ))
    }
    
    if (!inherits(prior, "bf_ttest_prior")) {
      stop("Unsupported prior object encountered while plotting.", call. = FALSE)
    }
    
    lower <- prior$support[1L]
    upper <- prior$support[2L]
    
    if (!is.finite(lower)) lower <- prior$location - 4 * prior$scale
    if (!is.finite(upper)) upper <- prior$location + 4 * prior$scale
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      stop("Could not determine finite plotting support for prior.", call. = FALSE)
    }
    
    grid <- seq(lower, upper, length.out = max(201L, as.integer(n)))
    dens <- prior_pdf(grid, prior)
    
    data.frame(
      delta = grid,
      density = dens,
      prior = label,
      prior_type = "continuous",
      stringsAsFactors = FALSE
    )
  }
  
  plot_prior_panel <- function(df, title) {
    split_df <- split(df, df$prior)
    is_point <- vapply(split_df, function(z) unique(z$prior_type) == "point", logical(1))
    point_names <- names(is_point)[is_point]
    curve_names <- names(is_point)[!is_point]
    
    y_max_curve <- 0
    if (length(curve_names)) {
      y_max_curve <- max(df$density[df$prior %in% curve_names], na.rm = TRUE)
      if (!is.finite(y_max_curve)) y_max_curve <- 0
    }
    
    y_max <- max(1, 1.05 * y_max_curve)
    
    p <- ggplot2::ggplot()
    
    if (length(curve_names)) {
      p <- p +
        ggplot2::geom_line(
          data = df[df$prior %in% curve_names, , drop = FALSE],
          ggplot2::aes(x = .data$delta, y = .data$density, colour = .data$prior),
          linewidth = 0.9
        )
    }
    
    if (length(point_names)) {
      p <- p +
        ggplot2::geom_vline(
          data = df[df$prior %in% point_names, , drop = FALSE],
          ggplot2::aes(xintercept = .data$delta, colour = .data$prior),
          linewidth = 0.9,
          linetype = "dotted",
          show.legend = TRUE
        )
    }
    
    p +
      ggplot2::scale_y_continuous(
        limits = c(0, y_max),
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ) +
      ggplot2::labs(
        title = title,
        x = expression(delta),
        y = "Prior density",
        colour = NULL
      ) +
      theme_base +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2, byrow = TRUE))
  }
  
  make_top_annotation <- function(x) {
    title_grob <- grid::textGrob(
      "Calibrated Bayesian t-test design",
      x = 0.5, y = 0.82,
      just = "center",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
    
    line1 <- switch(
      x$test_type,
      "two.sided.pointnull" = bquote(H[0] * ":" ~ delta == 0 ~ "vs" ~ H[1] * ":" ~ delta != 0),
      "greater.pointnull" = bquote(H[0] * ":" ~ delta == 0 ~ "vs" ~ H[1] * ":" ~ delta > 0),
      "greater.composite" = bquote(H[0] * ":" ~ delta <= 0 ~ "vs" ~ H[1] * ":" ~ delta > 0),
      bquote(.(x$test_type))
    )
    
    line2 <- bquote(
      .(x$design_type) ~ "|" ~ n^"*" == .(x$n) ~ "|" ~
        k[1] == .(format(x$k1, trim = TRUE)) ~ "|" ~
        k[0] == .(format(x$k0, trim = TRUE))
    )
    
    line1_grob <- grid::textGrob(
      line1,
      x = 0.5, y = 0.48,
      just = "center",
      gp = grid::gpar(fontsize = 11)
    )
    
    line2_grob <- grid::textGrob(
      line2,
      x = 0.5, y = 0.18,
      just = "center",
      gp = grid::gpar(fontsize = 10)
    )
    
    patchwork::wrap_elements(
      full = grid::grobTree(title_grob, line1_grob, line2_grob)
    )
  }
  
  metric_label_map <- c(
    freq_t1e = "Frequentist type-I error",
    bayes_t1e = "Bayesian type-I error",
    bayes_power = "Bayesian power",
    pce_h0 = "PCE under H0",
    "freq_power:main" = "Frequentist power"
  )
  
  panels <- list()
  
  if ("decision" %in% which) {
    df_bf <- get_plot_data("bf_curve")
    h1_df <- intervals_to_df(x$decision_regions$evidence_h1$intervals, "Reject H0")
    h0_df <- intervals_to_df(x$decision_regions$evidence_h0$intervals, "Accept H0")
    
    x_rng <- range(df_bf$t, finite = TRUE)
    indecisive_df <- data.frame(
      lower = x_rng[1],
      upper = x_rng[2],
      region = "Indecisive",
      stringsAsFactors = FALSE
    )
    
    region_df <- indecisive_df
    if (nrow(h0_df)) region_df <- rbind(region_df, h0_df)
    if (nrow(h1_df)) region_df <- rbind(region_df, h1_df)
    
    h0_requested <- !is.null(x$inputs$pce_h0_target)
    
    p_decision <- ggplot2::ggplot(
      df_bf,
      ggplot2::aes(x = .data$t, y = .data$log_bf01)
    ) +
      ggplot2::geom_rect(
        data = region_df,
        ggplot2::aes(
          xmin = .data$lower,
          xmax = .data$upper,
          ymin = -Inf,
          ymax = Inf,
          fill = .data$region
        ),
        inherit.aes = FALSE,
        alpha = 0.14,
        colour = NA
      ) +
      ggplot2::geom_hline(
        yintercept = log(x$k1),
        colour = "#D55E00",
        linetype = "dashed",
        linewidth = 0.7
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        colour = "grey60",
        linewidth = 0.4
      ) +
      ggplot2::geom_line(linewidth = 0.9, colour = "black") +
      ggplot2::scale_fill_manual(
        values = c(
          "Reject H0" = "#D55E00",
          "Accept H0" = "#0072B2",
          "Indecisive" = "grey70"
        ),
        breaks = c("Reject H0", "Accept H0", "Indecisive"),
        name = NULL
      ) +
      ggplot2::labs(
        title = "Decision function",
        x = "Observed t statistic",
        y = "log(BF01)"
      ) +
      theme_base +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))
    
    if (h0_requested) {
      p_decision <- p_decision +
        ggplot2::geom_hline(
          yintercept = log(x$k0),
          colour = "#0072B2",
          linetype = "dashed",
          linewidth = 0.7
        )
    }
    
    if (!is.null(x$roots_k1) && length(x$roots_k1)) {
      p_decision <- p_decision +
        ggplot2::geom_vline(
          xintercept = x$roots_k1,
          colour = "#D55E00",
          linetype = "dotted",
          linewidth = 0.7,
          show.legend = FALSE
        )
    }
    
    if (!is.null(x$roots_k0) && length(x$roots_k0)) {
      p_decision <- p_decision +
        ggplot2::geom_vline(
          xintercept = x$roots_k0,
          colour = "#0072B2",
          linetype = "dotted",
          linewidth = 0.7,
          show.legend = FALSE
        )
    }
    
    panels$decision <- p_decision
  }
  
  if ("n_curve" %in% which) {
    df_long <- constraint_curve_df()
    if (!nrow(df_long)) {
      stop("Sample-size operating-characteristic panel requires `coarse_results`.", call. = FALSE)
    }
    
    target_df <- constraint_targets_df()
    
    df_long$metric_label <- ifelse(
      df_long$metric %in% names(metric_label_map),
      unname(metric_label_map[df_long$metric]),
      ifelse(substr(df_long$metric, 1L, 11L) == "freq_power:",
             "Frequentist power",
             df_long$metric)
    )
    
    df_long <- .deduplicate_constraint_curves_bf_ttest(df_long)
    
    if (nrow(target_df)) {
      target_df <- target_df[target_df$metric %in% unique(df_long$metric), , drop = FALSE]
      target_df$metric_label <- ifelse(
        target_df$metric %in% names(metric_label_map),
        unname(metric_label_map[target_df$metric]),
        ifelse(substr(target_df$metric, 1L, 11L) == "freq_power:",
               "Frequentist power",
               target_df$metric)
      )
    }
    
    present_breaks <- unique(df_long$metric_label)
    
    p_ncurve <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = .data$n, y = .data$value, colour = .data$metric_label)
    ) +
      ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
      ggplot2::geom_vline(
        xintercept = x$n,
        colour = "black",
        linetype = "dashed",
        linewidth = 0.7
      ) +
      ggplot2::scale_colour_discrete(breaks = present_breaks) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        expand = ggplot2::expansion(mult = c(0.02, 0.05))
      ) +
      ggplot2::labs(
        title = "Operating characteristics by sample size",
        x = "Sample size",
        y = "Probability",
        colour = NULL
      ) +
      theme_base +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2, byrow = TRUE))
    
    if (nrow(target_df)) {
      p_ncurve <- p_ncurve +
        ggplot2::geom_hline(
          data = target_df,
          ggplot2::aes(yintercept = .data$target, colour = .data$metric_label),
          linetype = "dotted",
          linewidth = 0.7,
          show.legend = FALSE
        )
    }
    
    panels$n_curve <- p_ncurve
  }
  
  if ("analysis_prior" %in% which) {
    df_analysis_prior <- rbind(
      prior_to_df(x$analysis_prior_h0, "H0 analysis prior", n = prior_points),
      prior_to_df(x$analysis_prior_h1 %||% x$analysis_prior, "H1 analysis prior", n = prior_points)
    )
    
    if (!nrow(df_analysis_prior)) {
      stop("Analysis priors are missing; cannot draw analysis prior panel.", call. = FALSE)
    }
    
    panels$analysis_prior <- plot_prior_panel(df_analysis_prior, "Analysis priors")
  }
  
  if ("design_prior" %in% which) {
    df_design_prior <- rbind(
      prior_to_df(.make_design_prior(x$design_priors$h0), "H0 design prior", n = prior_points),
      prior_to_df(.make_design_prior(x$design_priors$h1), "H1 design prior", n = prior_points)
    )
    
    if (!nrow(df_design_prior)) {
      stop("Design priors are missing; cannot draw design prior panel.", call. = FALSE)
    }
    
    panels$design_prior <- plot_prior_panel(df_design_prior, "Design priors")
  }
  
  if ("effect" %in% which) {
    df_oc <- get_plot_data("power_at_effect")
    
    p_effect <- ggplot2::ggplot(
      df_oc,
      ggplot2::aes(
        x = .data$effect,
        y = .data$probability_evidence_h1,
        colour = .data$role
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        expand = ggplot2::expansion(mult = c(0.02, 0.05))
      ) +
      ggplot2::labs(
        title = "Operating characteristics at selected effects",
        x = "True standardized effect",
        y = "P(BF01 <= k1)",
        colour = NULL
      ) +
      theme_base +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2, byrow = TRUE))
    
    panels$effect <- p_effect
  }
  
  panel_order <- c("decision", "n_curve", "analysis_prior", "design_prior")
  chosen_order <- panel_order[panel_order %in% names(panels)]
  
  if ("effect" %in% names(panels) && !"n_curve" %in% names(panels)) {
    panel_order <- c("decision", "effect", "analysis_prior", "design_prior")
    chosen_order <- panel_order[panel_order %in% names(panels)]
  }
  
  ordered_panels <- panels[chosen_order]
  
  if (length(ordered_panels) == 1L) {
    out <- ordered_panels[[1L]]
    if (annotate) {
      out <- make_top_annotation(x) / out + patchwork::plot_layout(heights = c(0.18, 1))
    }
    return(out)
  }
  
  core <- patchwork::wrap_plots(ordered_panels, ncol = 2, byrow = TRUE)
  
  if (annotate) {
    make_top_annotation(x) / core + patchwork::plot_layout(heights = c(0.14, 1))
  } else {
    core
  }
}

.constraint_curve_df_bf_ttest <- function(x) {
  extract_metric <- function(res_list, field) {
    vapply(res_list, function(z) {
      val <- z[[field]]
      if (is.null(val) || length(val) != 1L || !is.finite(val)) NA_real_ else val
    }, numeric(1))
  }
  
  if (is.null(x$coarse_results) || !length(x$coarse_results)) {
    return(data.frame())
  }
  
  n_vals <- as.integer(names(x$coarse_results))
  if (anyNA(n_vals)) {
    n_vals <- seq_along(x$coarse_results)
  }
  
  requested <- character(0)
  if (!is.null(x$inputs$freq_t1e_target)) requested <- c(requested, "freq_t1e")
  if (!is.null(x$inputs$bayes_t1e_target)) requested <- c(requested, "bayes_t1e")
  if (!is.null(x$inputs$bayes_power_target)) requested <- c(requested, "bayes_power")
  if (!is.null(x$inputs$pce_h0_target)) requested <- c(requested, "pce_h0")
  
  fp_targets <- x$inputs$freq_power_target
  requested_fp <- character(0)
  if (!is.null(fp_targets) && length(fp_targets)) {
    requested_fp <- names(fp_targets)[
      vapply(fp_targets, function(z) !is.null(z$target), logical(1))
    ]
    requested <- c(requested, paste0("freq_power:", requested_fp))
  }
  
  if (!length(requested)) {
    return(data.frame())
  }
  
  df <- data.frame(n = n_vals, stringsAsFactors = FALSE)
  
  if ("freq_t1e" %in% requested) {
    df$freq_t1e <- extract_metric(x$coarse_results, "freq_t1e")
  }
  if ("bayes_t1e" %in% requested) {
    df$bayes_t1e <- extract_metric(x$coarse_results, "bayes_t1e")
  }
  if ("bayes_power" %in% requested) {
    df$bayes_power <- extract_metric(x$coarse_results, "bayes_power")
  }
  if ("pce_h0" %in% requested) {
    df$pce_h0 <- extract_metric(x$coarse_results, "pce_h0")
  }
  
  if (length(requested_fp)) {
    for (nm in requested_fp) {
      df[[paste0("freq_power:", nm)]] <- vapply(
        x$coarse_results,
        function(z) {
          val <- z$freq_power[[nm]]
          if (is.null(val) || length(val) != 1L || !is.finite(val)) NA_real_ else val
        },
        numeric(1)
      )
    }
  }
  
  keep <- setdiff(names(df), "n")
  keep <- keep[vapply(keep, function(nm) any(is.finite(df[[nm]])), logical(1))]
  
  if (!length(keep)) {
    return(data.frame())
  }
  
  out <- do.call(rbind, lapply(keep, function(nm) {
    data.frame(
      n = df$n,
      metric = nm,
      value = df[[nm]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

.deduplicate_constraint_curves_bf_ttest <- function(df_long, tol = 1e-12) {
  if (!nrow(df_long)) return(df_long)
  
  split_df <- split(df_long, df_long$metric)
  
  curve_equal <- function(a, b) {
    xa <- split_df[[a]][order(split_df[[a]]$n), c("n", "value")]
    xb <- split_df[[b]][order(split_df[[b]]$n), c("n", "value")]
    
    if (nrow(xa) != nrow(xb)) return(FALSE)
    if (!identical(xa$n, xb$n)) return(FALSE)
    
    both_na <- is.na(xa$value) & is.na(xb$value)
    both_ok <- is.finite(xa$value) & is.finite(xb$value)
    
    if (!all(both_na | both_ok)) return(FALSE)
    if (any(abs(xa$value[both_ok] - xb$value[both_ok]) > tol)) return(FALSE)
    
    TRUE
  }
  
  metric_names <- names(split_df)
  drop_metrics <- character(0)
  
  if (all(c("freq_t1e", "bayes_t1e") %in% metric_names) &&
      curve_equal("freq_t1e", "bayes_t1e")) {
    drop_metrics <- c(drop_metrics, "bayes_t1e")
  }
  
  fp_metrics <- grep("^freq_power:", metric_names, value = TRUE)
  
  if ("bayes_power" %in% metric_names && length(fp_metrics)) {
    for (nm in fp_metrics) {
      if (curve_equal("bayes_power", nm)) {
        drop_metrics <- c(drop_metrics, "bayes_power")
        break
      }
    }
  }
  
  keep <- setdiff(metric_names, unique(drop_metrics))
  out <- do.call(rbind, split_df[keep])
  rownames(out) <- NULL
  out
}