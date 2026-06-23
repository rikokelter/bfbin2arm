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
plot.twoarm_onestage_bf_design <- function(
    x,
    type = c("old", "oc", "feasibility"),
    ...
) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
  
  sr <- x$search_results
  if (is.null(sr) || nrow(sr) == 0L) {
    stop("No search results available for plotting.", call. = FALSE)
  }
  
  if (type == "oc") {
    df_long <- rbind(
      data.frame(n_total = sr$n_total, metric = "Power", value = sr$power),
      data.frame(n_total = sr$n_total, metric = "Type-I error", value = sr$type1),
      data.frame(n_total = sr$n_total, metric = "CE(H0)", value = sr$ce_h0)
    )
    
    p <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = .data$n_total, y = .data$value, colour = .data$metric)
    ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_colour_manual(
        values = c("Power" = "black", "Type-I error" = "blue", "CE(H0)" = "red")
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_minimal(base_size = 12)
    
    return(p)
  }
  
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
      ggplot2::aes(x = .data$n_total, y = .data$value, colour = .data$status)
    ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_y_continuous(breaks = c(0, 1), limits = c(-0.05, 1.05)) +
      ggplot2::theme_minimal(base_size = 12)
    
    return(p)
  }
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for type = 'old'.", call. = FALSE)
  }
  
  df <- sr
  
  if (!"n_total" %in% names(df) && "n" %in% names(df)) {
    df$n_total <- df$n
  }
  if (!"n" %in% names(df) && "n_total" %in% names(df)) {
    df$n <- df$n_total
  }
  if (!"t1e" %in% names(df) && "type1" %in% names(df)) {
    df$t1e <- df$type1
  }
  if (!"pceH0" %in% names(df) && "ce_h0" %in% names(df)) {
    df$pceH0 <- df$ce_h0
  }
  
  ns <- df$n_total
  powervec_k <- df$power
  t1evec_k <- df$t1e
  pceH0_vec_k_f <- df$pceH0
  
  power_target <- x$inputs$target_power
  alpha_target <- x$inputs$target_type1
  pce_target <- x$inputs$target_ce_h0
  freq_power_target <- x$inputs$target_freq_power
  sustain_n <- x$inputs$sustain_n
  
  first_sustained_index <- function(ok_vec, sustain_n) {
    n <- length(ok_vec)
    if (n == 0L) return(NA_integer_)
    
    for (i in seq_len(n)) {
      j <- min(n, i + sustain_n)
      window_ok <- ok_vec[i:j]
      if (length(window_ok) > 0L &&
          all(!is.na(window_ok)) &&
          all(window_ok)) {
        return(i)
      }
    }
    
    NA_integer_
  }
  
  power_ok <- !is.na(powervec_k) & (powervec_k >= power_target)
  t1e_ok <- !is.na(t1evec_k) & (t1evec_k <= alpha_target)
  pce_ok <- !is.na(pceH0_vec_k_f) & (pceH0_vec_k_f >= pce_target)
  
  i_power <- first_sustained_index(power_ok, sustain_n)
  i_t1e <- first_sustained_index(t1e_ok, sustain_n)
  i_pce <- if (pce_target > 0) first_sustained_index(pce_ok, sustain_n) else NA_integer_
  
  n_power <- if (!is.na(i_power)) ns[i_power] else NA_integer_
  n_t1e <- if (!is.na(i_t1e)) ns[i_t1e] else NA_integer_
  n_pceH0 <- if (!is.na(i_pce)) ns[i_pce] else NA_integer_
  
  freq_power_available <- "freq_power" %in% names(df) && any(!is.na(df$freq_power))
  if (freq_power_available) {
    freq_power_ok <- !is.na(df$freq_power) & (df$freq_power >= freq_power_target)
    i_freq_power <- first_sustained_index(freq_power_ok, sustain_n)
    n_freq_power <- if (!is.na(i_freq_power)) ns[i_freq_power] else NA_integer_
  } else {
    i_freq_power <- NA_integer_
    n_freq_power <- NA_integer_
  }
  
  n_power_text <- if (!is.na(n_power)) sprintf("n=%d", n_power) else "not reached"
  n_t1e_text <- if (!is.na(n_t1e)) sprintf("n=%d", n_t1e) else "not reached"
  n_freq_power_text <- if (!is.na(n_freq_power)) sprintf("n=%d", n_freq_power) else "not reached"
  n_pce_text <- if (!is.na(n_pceH0)) sprintf("n=%d", n_pceH0) else "not reached"
  
  x_annot <- min(ns) + (max(ns) - min(ns)) * 0.02
  y_base <- 0.9
  y_step <- -0.08
  
  p1_plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n_total)) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$power, color = "Bayes Power"),
      linewidth = 0.9
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$t1e, color = "Type I Error"),
      linewidth = 0.9
    )
  
  if (freq_power_available) {
    p1_plot <- p1_plot +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$freq_power, color = "Frequentist power"),
        linewidth = 0.9
      )
  }
  
  p1_plot <- p1_plot +
    ggplot2::scale_color_manual(
      values = c(
        "Bayes Power" = "black",
        "Type I Error" = "blue",
        "Frequentist power" = "green"
      )
    ) +
    ggplot2::geom_hline(
      yintercept = power_target,
      linetype = "dashed",
      color = "black",
      alpha = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = alpha_target,
      linetype = "dashed",
      color = "blue",
      alpha = 0.7
    )
  
  if (!is.na(n_power)) {
    p1_plot <- p1_plot +
      ggplot2::geom_vline(
        xintercept = n_power,
        color = "black",
        linewidth = 1.1
      )
  }
  
  if (!is.na(n_t1e)) {
    p1_plot <- p1_plot +
      ggplot2::geom_vline(
        xintercept = n_t1e,
        color = "blue",
        linewidth = 1.1
      )
  }
  
  if (!is.na(n_freq_power)) {
    p1_plot <- p1_plot +
      ggplot2::geom_vline(
        xintercept = n_freq_power,
        color = "green",
        linewidth = 1.1
      )
  }
  
  labels_to_draw <- list(
    list(text = n_power_text, color = "black"),
    list(text = n_t1e_text, color = "blue")
  )
  
  if (freq_power_available) {
    labels_to_draw <- append(
      labels_to_draw,
      list(list(text = n_freq_power_text, color = "green"))
    )
  }
  
  for (i in seq_along(labels_to_draw)) {
    p1_plot <- p1_plot +
      ggplot2::annotate(
        "text",
        x = x_annot,
        y = y_base + (i - 1) * y_step,
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
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )
  
  bottom_title_expr <- if (x$inputs$test == "BF+-") {
    expression("Probability of compelling evidence for " * H["-"])
  } else {
    expression("Probability of compelling evidence for H0")
  }
  
  p2_plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n_total, y = .data$pceH0)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::geom_line(color = "red", linewidth = 0.9)
  
  if (pce_target > 0) {
    p2_plot <- p2_plot +
      ggplot2::geom_hline(
        yintercept = pce_target,
        linetype = "dashed",
        color = "red",
        alpha = 0.7
      )
    
    if (!is.na(n_pceH0)) {
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
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
    )
  
  xseq <- seq(0, 1, 0.01)
  inp <- x$inputs
  
  if (inp$test == "BF01") {
    prior_df <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, inp$a_1_d, inp$b_1_d),
        stats::dbeta(xseq, inp$a_1_a, inp$b_1_a),
        stats::dbeta(xseq, inp$a_2_d, inp$b_2_d),
        stats::dbeta(xseq, inp$a_2_a, inp$b_2_a)
      ),
      prior_type = rep(c("Design", "Analysis", "Design", "Analysis"), each = length(xseq)),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"), each = length(xseq)),
      hypothesis = "H[1]"
    )
  } else if (inp$test == "BF+0") {
    prior_df <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, inp$a_1_d, inp$b_1_d),
        stats::dbeta(xseq, inp$a_1_a, inp$b_1_a),
        stats::dbeta(xseq, inp$a_2_d, inp$b_2_d),
        stats::dbeta(xseq, inp$a_2_a, inp$b_2_a)
      ),
      prior_type = rep(c("Design", "Analysis", "Design", "Analysis"), each = length(xseq)),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"), each = length(xseq)),
      hypothesis = "H[+]"
    )
  } else if (inp$test == "BF-0") {
    prior_df <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, inp$a_1_d_Hminus, inp$b_1_d_Hminus),
        stats::dbeta(xseq, inp$a_1_a_Hminus, inp$b_1_a_Hminus),
        stats::dbeta(xseq, inp$a_2_d_Hminus, inp$b_2_d_Hminus),
        stats::dbeta(xseq, inp$a_2_a_Hminus, inp$b_2_a_Hminus)
      ),
      prior_type = rep(c("Design", "Analysis", "Design", "Analysis"), each = length(xseq)),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"), each = length(xseq)),
      hypothesis = "H[-]"
    )
  } else {
    prior_df_Hplus <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, inp$a_1_d, inp$b_1_d),
        stats::dbeta(xseq, inp$a_1_a, inp$b_1_a),
        stats::dbeta(xseq, inp$a_2_d, inp$b_2_d),
        stats::dbeta(xseq, inp$a_2_a, inp$b_2_a)
      ),
      prior_type = rep(c("Design", "Analysis", "Design", "Analysis"), each = length(xseq)),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"), each = length(xseq)),
      hypothesis = "H[+]"
    )
    
    prior_df_Hminus <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, inp$a_1_d_Hminus, inp$b_1_d_Hminus),
        stats::dbeta(xseq, inp$a_1_a_Hminus, inp$b_1_a_Hminus),
        stats::dbeta(xseq, inp$a_2_d_Hminus, inp$b_2_d_Hminus),
        stats::dbeta(xseq, inp$a_2_a_Hminus, inp$b_2_a_Hminus)
      ),
      prior_type = rep(c("Design", "Analysis", "Design", "Analysis"), each = length(xseq)),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"), each = length(xseq)),
      hypothesis = "H[-]"
    )
    
    prior_df <- rbind(prior_df_Hplus, prior_df_Hminus)
  }
  
  prior_title_expr <- switch(
    inp$test,
    "BF01" = expression("Design and analysis priors under " * H[1] * ": " * p[1] != p[2]),
    "BF+0" = expression("Design and analysis priors under " * H["+"] * ": " * p[2] > p[1]),
    "BF-0" = expression("Design and analysis priors under " * H["-"] * ": " * p[1] > p[2]),
    "BF+-" = expression("Design and analysis priors under " * H["+"] * " and " * H["-"])
  )
  
  p_priors <- ggplot2::ggplot(
    prior_df,
    ggplot2::aes(x = x, y = density, linetype = prior_type)
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::facet_grid(
      hypothesis ~ param,
      scales = "free_y",
      labeller = ggplot2::labeller(param = ggplot2::label_parsed)
    ) +
    ggplot2::scale_linetype_manual(
      breaks = c("Design", "Analysis"),
      values = c("Design" = "solid", "Analysis" = "dotted"),
      labels = c("Design prior (solid)", "Analysis prior (dotted)")
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
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      legend.key.width = grid::unit(1.5, "cm")
    )
  
  prior_height <- ifelse(inp$test == "BF+-", 2.0, 1.4)
  
  p_priors / p1_plot / p2_plot +
    patchwork::plot_layout(heights = c(prior_height, 2.5, 1.3))
}