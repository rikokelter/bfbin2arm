#' @export
print.twoarm_onestage_bf_design <- function(x, ...) {
  cat("\nOne-stage two-arm Bayes factor design\n")
  cat("------------------------------------\n")
  cat("Mode: ", x$mode, "\n", sep = "")
  cat("Status: ", x$status, "\n", sep = "")
  cat("Calibration: ", x$calibration, "\n", sep = "")
  
  report_freq_type1 <- isTRUE(x$inputs$report_freq_type1)
  cat(
    "Optional freq. Type-I reporting: ",
    if (report_freq_type1) "on" else "off",
    "\n",
    sep = ""
  )
  
  if (isTRUE(x$feasible) && !is.null(x$design)) {
    cat(
      "Design: n_total = ", x$design["n_total"],
      ", n1 = ", x$design["n1"],
      ", n2 = ", x$design["n2"], "\n",
      sep = ""
    )
  } else {
    cat("Design: no feasible design found\n")
  }
  
  oc <- x$operating_characteristics
  if (!is.null(oc)) {
    cat("\nOperating characteristics\n")
    cat(
      " Power = ", formatC(oc$power, digits = 4, format = "f"), "\n",
      " Type-I error = ", formatC(oc$type1, digits = 4, format = "f"), "\n",
      " CE(H0) = ", formatC(oc$ce_h0, digits = 4, format = "f"), "\n",
      sep = ""
    )
    
    if (!is.na(oc$freq_type1)) {
      cat(
        " Freq. Type-I = ", formatC(oc$freq_type1, digits = 4, format = "f"), "\n",
        sep = ""
      )
    }
    
    if (!is.na(oc$freq_power)) {
      cat(
        " Freq. Power = ", formatC(oc$freq_power, digits = 4, format = "f"), "\n",
        sep = ""
      )
    }
  }
  
  invisible(x)
}

#' @export
summary.twoarm_onestage_bf_design <- function(object, ...) {
  sr <- object$search_results
  
  out <- list(
    mode = object$mode,
    status = object$status,
    feasible = object$feasible,
    calibration = object$calibration,
    design = object$design,
    operating_characteristics = object$operating_characteristics,
    targets = list(
      target_power = object$inputs$target_power,
      target_type1 = object$inputs$target_type1,
      target_ce_h0 = object$inputs$target_ce_h0,
      target_freq_power = object$inputs$target_freq_power,
      target_freq_type1 = object$inputs$target_freq_type1,
      sustain_n = object$inputs$sustain_n,
      report_freq_type1 = object$inputs$report_freq_type1
    ),
    inputs = object$inputs
  )
  
  if (!is.null(sr)) {
    out$search_overview <- list(
      n_evaluated = nrow(sr),
      n_pointwise_feasible = sum(sr$feasible_pointwise, na.rm = TRUE),
      n_sustained_feasible = sum(sr$feasible, na.rm = TRUE),
      first_pointwise_n = if (any(sr$feasible_pointwise, na.rm = TRUE)) {
        sr$n_total[which(sr$feasible_pointwise)[1L]]
      } else {
        NA_integer_
      },
      first_sustained_n = if (any(sr$feasible, na.rm = TRUE)) {
        sr$n_total[which(sr$feasible)[1L]]
      } else {
        NA_integer_
      }
    )
  }
  
  class(out) <- "summary.twoarm_onestage_bf_design"
  out
}

#' @export
print.summary.twoarm_onestage_bf_design <- function(x, ...) {
  cat("\nSummary: One-stage two-arm Bayes factor design\n")
  cat("---------------------------------------------\n")
  cat("Mode:        ", x$mode, "\n", sep = "")
  cat("Status:      ", x$status, "\n", sep = "")
  cat("Calibration: ", x$calibration, "\n", sep = "")
  cat("Feasible:    ", if (isTRUE(x$feasible)) "yes" else "no", "\n", sep = "")
  
  if (!is.null(x$search_overview)) {
    cat("\nSearch overview\n")
    cat("  n evaluated          = ", x$search_overview$n_evaluated, "\n", sep = "")
    cat("  pointwise feasible   = ", x$search_overview$n_pointwise_feasible, "\n", sep = "")
    cat("  sustained feasible   = ", x$search_overview$n_sustained_feasible, "\n", sep = "")
    cat("  first pointwise n    = ", x$search_overview$first_pointwise_n, "\n", sep = "")
    cat("  first sustained n    = ", x$search_overview$first_sustained_n, "\n", sep = "")
  }
  
  if (!is.null(x$design) && !all(is.na(x$design))) {
    cat(
      "\nSelected design\n",
      "  n_total = ", x$design["n_total"],
      ", n1 = ", x$design["n1"],
      ", n2 = ", x$design["n2"], "\n",
      sep = ""
    )
  }
  
  invisible(x)
}

#' @importFrom rlang .data
#' @export
#' @importFrom rlang .data
#' @export
plot.twoarm_onestage_bf_design <- function(
    x,
    type = c("old", "oc", "feasibility"),
    ...
) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plotting.",
      call. = FALSE
    )
  }
  
  sr <- x$search_results
  
  if (is.null(sr) || nrow(sr) == 0L) {
    stop(
      "No search results available for plotting.",
      call. = FALSE
    )
  }
  
  # Standardize column aliases.
  if (!"n_total" %in% names(sr) && "n" %in% names(sr)) {
    sr$n_total <- sr$n
  }
  
  if (!"n" %in% names(sr) && "n_total" %in% names(sr)) {
    sr$n <- sr$n_total
  }
  
  if (!"type1" %in% names(sr) && "t1e" %in% names(sr)) {
    sr$type1 <- sr$t1e
  }
  
  if (!"t1e" %in% names(sr) && "type1" %in% names(sr)) {
    sr$t1e <- sr$type1
  }
  
  if (!"ce_h0" %in% names(sr) && "pceH0" %in% names(sr)) {
    sr$ce_h0 <- sr$pceH0
  }
  
  if (!"pceH0" %in% names(sr) && "ce_h0" %in% names(sr)) {
    sr$pceH0 <- sr$ce_h0
  }
  
  freq_power_available <-
    "freq_power" %in% names(sr) &&
    any(is.finite(sr$freq_power))
  
  freq_type1_available <-
    "freq_type1" %in% names(sr) &&
    any(is.finite(sr$freq_type1))
  
  freq_type1_colour <- "#D55E00"
  
  # ------------------------------------------------------------------
  # Compact operating-characteristic plot.
  # ------------------------------------------------------------------
  
  if (type == "oc") {
    df_long <- rbind(
      data.frame(
        n_total = sr$n_total,
        metric = "Power",
        value = sr$power
      ),
      data.frame(
        n_total = sr$n_total,
        metric = "Bayesian type-I error",
        value = sr$type1
      ),
      data.frame(
        n_total = sr$n_total,
        metric = "CE(H0)",
        value = sr$ce_h0
      )
    )
    
    if (freq_type1_available) {
      df_long <- rbind(
        df_long,
        data.frame(
          n_total = sr$n_total,
          metric = "Frequentist type-I error",
          value = sr$freq_type1
        )
      )
    }
    
    p <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(
        x = .data$n_total,
        y = .data$value,
        colour = .data$metric
      )
    ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_colour_manual(
        values = c(
          "Power" = "black",
          "Bayesian type-I error" = "blue",
          "CE(H0)" = "red",
          "Frequentist type-I error" = freq_type1_colour
        )
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_minimal(base_size = 12)
    
    return(p)
  }
  
  # ------------------------------------------------------------------
  # Feasibility plot.
  # ------------------------------------------------------------------
  
  if (type == "feasibility") {
    df_long <- rbind(
      data.frame(
        n_total = sr$n_total,
        status = "Pointwise feasible",
        value = as.numeric(sr$feasible_pointwise)
      ),
      data.frame(
        n_total = sr$n_total,
        status = "Sustained feasible",
        value = as.numeric(sr$feasible)
      )
    )
    
    p <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(
        x = .data$n_total,
        y = .data$value,
        colour = .data$status
      )
    ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_y_continuous(
        breaks = c(0, 1),
        limits = c(-0.05, 1.05)
      ) +
      ggplot2::theme_minimal(base_size = 12)
    
    return(p)
  }
  
  # ------------------------------------------------------------------
  # Original multi-panel layout.
  # ------------------------------------------------------------------
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required for type = 'old'.",
      call. = FALSE
    )
  }
  
  df <- sr
  ns <- df$n_total
  inp <- x$inputs
  
  powervec_k <- df$power
  t1evec_k <- df$t1e
  pceH0_vec_k_f <- df$pceH0
  
  power_target <- inp$target_power
  alpha_target <- inp$target_type1
  pce_target <- inp$target_ce_h0
  freq_power_target <- inp$target_freq_power
  freq_type1_target <- inp$target_freq_type1
  sustain_n <- inp$sustain_n
  
  valid_target <- function(z) {
    is.numeric(z) &&
      length(z) == 1L &&
      is.finite(z)
  }
  
  if (
    is.null(sustain_n) ||
    length(sustain_n) != 1L ||
    !is.finite(sustain_n)
  ) {
    sustain_n <- 0L
  }
  
  sustain_n <- max(0L, as.integer(sustain_n))
  
  power_target_available <- valid_target(power_target)
  alpha_target_available <- valid_target(alpha_target)
  pce_target_available <- valid_target(pce_target)
  freq_type1_target_available <- valid_target(freq_type1_target)
  
  # Preserve the original forward-window convention.
  first_sustained_index <- function(ok_vec, sustain_n) {
    n <- length(ok_vec)
    
    if (n == 0L) {
      return(NA_integer_)
    }
    
    for (i in seq_len(n)) {
      j <- min(n, i + sustain_n)
      window_ok <- ok_vec[i:j]
      
      if (
        length(window_ok) > 0L &&
        all(!is.na(window_ok)) &&
        all(window_ok)
      ) {
        return(i)
      }
    }
    
    NA_integer_
  }
  
  attainment_n <- function(values, target, direction) {
    if (!valid_target(target) || is.null(values)) {
      return(NA_real_)
    }
    
    ok <- if (direction == ">=") {
      is.finite(values) & values >= target
    } else {
      is.finite(values) & values <= target
    }
    
    idx <- first_sustained_index(ok, sustain_n)
    
    if (is.na(idx)) {
      return(NA_real_)
    }
    
    ns[idx]
  }
  
  n_power <- attainment_n(
    powervec_k,
    power_target,
    ">="
  )
  
  n_t1e <- attainment_n(
    t1evec_k,
    alpha_target,
    "<="
  )
  
  n_pceH0 <- if (
    pce_target_available &&
    pce_target > 0
  ) {
    attainment_n(
      pceH0_vec_k_f,
      pce_target,
      ">="
    )
  } else {
    NA_real_
  }
  
  n_freq_power <- if (freq_power_available) {
    attainment_n(
      df$freq_power,
      freq_power_target,
      ">="
    )
  } else {
    NA_real_
  }
  
  n_freq_type1 <- if (freq_type1_available) {
    attainment_n(
      df$freq_type1,
      freq_type1_target,
      "<="
    )
  } else {
    NA_real_
  }
  
  # Middle-panel labels contain only n = ...
  format_n <- function(n) {
    if (length(n) == 1L && is.finite(n)) {
      sprintf("n = %d", as.integer(n))
    } else {
      "n = NA"
    }
  }
  
  n_power_text <- format_n(n_power)
  n_t1e_text <- format_n(n_t1e)
  n_freq_power_text <- format_n(n_freq_power)
  n_freq_type1_text <- format_n(n_freq_type1)
  
  # Preserve the previous bottom-panel annotation.
  n_pce_text <- if (is.finite(n_pceH0)) {
    sprintf(
      "Null-evidence probability only: n=%d",
      as.integer(n_pceH0)
    )
  } else {
    "Null-evidence probability only: not reached"
  }
  
  x_annot <- min(ns) + (max(ns) - min(ns)) * 0.02
  
  # Increased spacing for middle-panel annotations.
  y_base <- 0.90
  y_step <- -0.14
  
  # ------------------------------------------------------------------
  # Middle panel: power and type-I error.
  # ------------------------------------------------------------------
  
  p1_plot <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$n_total)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = .data$power,
        color = "Bayes Power"
      ),
      linewidth = 0.9
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = .data$t1e,
        color = "Bayesian type-I error"
      ),
      linewidth = 0.9
    )
  
  if (freq_power_available) {
    p1_plot <- p1_plot +
      ggplot2::geom_line(
        ggplot2::aes(
          y = .data$freq_power,
          color = "Frequentist power"
        ),
        linewidth = 0.9
      )
  }
  
  if (freq_type1_available) {
    p1_plot <- p1_plot +
      ggplot2::geom_line(
        ggplot2::aes(
          y = .data$freq_type1,
          color = "Frequentist type-I error"
        ),
        linewidth = 0.9
      )
  }
  
  p1_plot <- p1_plot +
    ggplot2::scale_color_manual(
      values = c(
        "Bayes Power" = "black",
        "Bayesian type-I error" = "blue",
        "Frequentist power" = "green",
        "Frequentist type-I error" = freq_type1_colour
      )
    )
  
  if (power_target_available) {
    p1_plot <- p1_plot +
      ggplot2::geom_hline(
        yintercept = power_target,
        linetype = "dashed",
        color = "black",
        alpha = 0.7
      )
  }
  
  if (alpha_target_available) {
    p1_plot <- p1_plot +
      ggplot2::geom_hline(
        yintercept = alpha_target,
        linetype = "dashed",
        color = "blue",
        alpha = 0.7
      )
  }
  
  if (freq_type1_available && freq_type1_target_available) {
    p1_plot <- p1_plot +
      ggplot2::geom_hline(
        yintercept = freq_type1_target,
        linetype = "dotdash",
        color = freq_type1_colour,
        alpha = 0.8
      )
  }
  
  if (is.finite(n_power)) {
    p1_plot <- p1_plot +
      ggplot2::geom_vline(
        xintercept = n_power,
        color = "black",
        linewidth = 1.1
      )
  }
  
  if (is.finite(n_t1e)) {
    p1_plot <- p1_plot +
      ggplot2::geom_vline(
        xintercept = n_t1e,
        color = "blue",
        linewidth = 1.1
      )
  }
  
  if (is.finite(n_freq_power)) {
    p1_plot <- p1_plot +
      ggplot2::geom_vline(
        xintercept = n_freq_power,
        color = "green",
        linewidth = 1.1
      )
  }
  
  if (is.finite(n_freq_type1)) {
    p1_plot <- p1_plot +
      ggplot2::geom_vline(
        xintercept = n_freq_type1,
        color = freq_type1_colour,
        linewidth = 1.1
      )
  }
  
  labels_to_draw <- list(
    list(
      text = n_power_text,
      color = "black"
    ),
    list(
      text = n_t1e_text,
      color = "blue"
    )
  )
  
  if (freq_power_available) {
    labels_to_draw <- append(
      labels_to_draw,
      list(
        list(
          text = n_freq_power_text,
          color = "green"
        )
      )
    )
  }
  
  if (freq_type1_available) {
    labels_to_draw <- append(
      labels_to_draw,
      list(
        list(
          text = n_freq_type1_text,
          color = freq_type1_colour
        )
      )
    )
  }
  
  for (i in seq_along(labels_to_draw)) {
    p1_plot <- p1_plot +
      ggplot2::annotate(
        "text",
        x = x_annot,
        y = y_base + (i - 1L) * y_step,
        label = labels_to_draw[[i]]$text,
        color = labels_to_draw[[i]]$color,
        size = 4,
        fontface = "bold",
        hjust = 0
      )
  }
  
  p1_plot <- p1_plot +
    ggplot2::labs(
      title = "Power and Type-I Error Rate",
      x = "Total sample size",
      y = "Probability"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )
  
  # ------------------------------------------------------------------
  # Bottom panel: probability of compelling evidence.
  # ------------------------------------------------------------------
  
  bottom_title_expr <- if (inp$test == "BF+-") {
    expression(
      "Probability of compelling evidence for " * H["-"]
    )
  } else {
    expression("Probability of compelling evidence for H0")
  }
  
  p2_plot <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$n_total,
      y = .data$pceH0
    )
  ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::geom_line(
      color = "red",
      linewidth = 0.9
    )
  
  if (pce_target_available && pce_target > 0) {
    p2_plot <- p2_plot +
      ggplot2::geom_hline(
        yintercept = pce_target,
        linetype = "dashed",
        color = "red",
        alpha = 0.7
      )
    
    if (is.finite(n_pceH0)) {
      p2_plot <- p2_plot +
        ggplot2::geom_vline(
          xintercept = n_pceH0,
          color = "red",
          linewidth = 1.1
        ) +
        ggplot2::annotate(
          "text",
          x = x_annot,
          y = 0.9,
          label = n_pce_text,
          color = "red",
          size = 4,
          fontface = "bold",
          hjust = 0
        )
    }
  }
  
  p2_plot <- p2_plot +
    ggplot2::labs(
      title = bottom_title_expr,
      x = "Total sample size",
      y = "Probability"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 14
      )
    )
  
  # ------------------------------------------------------------------
  # Prior panels: original beta basis-density calculations.
  # ------------------------------------------------------------------
  
  xseq <- seq(0, 1, 0.01)
  
  make_prior_df <- function(
    a1d, b1d, a1a, b1a,
    a2d, b2d, a2a, b2a,
    hypothesis) {
    data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, a1d, b1d),
        stats::dbeta(xseq, a1a, b1a),
        stats::dbeta(xseq, a2d, b2d),
        stats::dbeta(xseq, a2a, b2a)
      ),
      prior_type = rep(
        c("Design", "Analysis", "Design", "Analysis"),
        each = length(xseq)
      ),
      param = rep(
        c("p[1]", "p[1]", "p[2]", "p[2]"),
        each = length(xseq)
      ),
      hypothesis = hypothesis
    )
  }
  
  if (inp$test == "BF01") {
    prior_df <- make_prior_df(
      a1d = inp$a_1_d,
      b1d = inp$b_1_d,
      a1a = inp$a_1_a,
      b1a = inp$b_1_a,
      a2d = inp$a_2_d,
      b2d = inp$b_2_d,
      a2a = inp$a_2_a,
      b2a = inp$b_2_a,
      hypothesis = "H[1]"
    )
  } else if (inp$test == "BF+0") {
    prior_df <- make_prior_df(
      a1d = inp$a_1_d,
      b1d = inp$b_1_d,
      a1a = inp$a_1_a,
      b1a = inp$b_1_a,
      a2d = inp$a_2_d,
      b2d = inp$b_2_d,
      a2a = inp$a_2_a,
      b2a = inp$b_2_a,
      hypothesis = "H[+]"
    )
  } else if (inp$test == "BF-0") {
    prior_df <- make_prior_df(
      a1d = inp$a_1_d_Hminus,
      b1d = inp$b_1_d_Hminus,
      a1a = inp$a_1_a_Hminus,
      b1a = inp$b_1_a_Hminus,
      a2d = inp$a_2_d_Hminus,
      b2d = inp$b_2_d_Hminus,
      a2a = inp$a_2_a_Hminus,
      b2a = inp$b_2_a_Hminus,
      hypothesis = "H[-]"
    )
  } else {
    prior_df_Hplus <- make_prior_df(
      a1d = inp$a_1_d,
      b1d = inp$b_1_d,
      a1a = inp$a_1_a,
      b1a = inp$b_1_a,
      a2d = inp$a_2_d,
      b2d = inp$b_2_d,
      a2a = inp$a_2_a,
      b2a = inp$b_2_a,
      hypothesis = "H[+]"
    )
    
    prior_df_Hminus <- make_prior_df(
      a1d = inp$a_1_d_Hminus,
      b1d = inp$b_1_d_Hminus,
      a1a = inp$a_1_a_Hminus,
      b1a = inp$b_1_a_Hminus,
      a2d = inp$a_2_d_Hminus,
      b2d = inp$b_2_d_Hminus,
      a2a = inp$a_2_a_Hminus,
      b2a = inp$b_2_a_Hminus,
      hypothesis = "H[-]"
    )
    
    prior_df <- rbind(
      prior_df_Hplus,
      prior_df_Hminus
    )
  }
  
  prior_title_expr <- switch(
    inp$test,
    "BF01" = expression(
      "Design and analysis priors under " *
        H[1] * ": " * p[1] != p[2]
    ),
    "BF+0" = expression(
      "Design and analysis priors under " *
        H["+"] * ": " * p[2] > p[1]
    ),
    "BF-0" = expression(
      "Design and analysis priors under " *
        H["-"] * ": " * p[1] > p[2]
    ),
    "BF+-" = expression(
      "Design and analysis priors under " *
        H["+"] * " and " * H["-"]
    )
  )
  
  p_priors <- ggplot2::ggplot(
    prior_df,
    ggplot2::aes(
      x = .data$x,
      y = .data$density,
      linetype = .data$prior_type
    )
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::facet_grid(
      hypothesis ~ param,
      scales = "free_y",
      labeller = ggplot2::labeller(
        param = ggplot2::label_parsed
      )
    ) +
    ggplot2::scale_linetype_manual(
      breaks = c("Design", "Analysis"),
      values = c(
        "Design" = "solid",
        "Analysis" = "dotted"
      ),
      labels = c(
        "Design prior (solid)",
        "Analysis prior (dotted)"
      )
    ) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(
        override.aes = list(
          linetype = c("solid", "dotted"),
          color = c("black", "black")
        )
      )
    ) +
    ggplot2::labs(
      title = prior_title_expr,
      linetype = "Prior type",
      y = "Density",
      x = "Probability"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 14
      ),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(
        size = 12,
        face = "bold"
      ),
      legend.key.width = grid::unit(1.5, "cm")
    )
  
  prior_height <- ifelse(
    inp$test == "BF+-",
    2.0,
    1.4
  )
  
  p_priors / p1_plot / p2_plot +
    patchwork::plot_layout(
      heights = c(prior_height, 2.5, 1.3)
    )
}