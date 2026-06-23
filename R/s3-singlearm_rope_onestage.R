#' @keywords internal
#' @noRd
.rope_compact_runs <- function(y) {
  if (length(y) == 0L) return("none")
  y <- sort(unique(as.integer(y)))
  
  starts <- c(y[1], y[which(diff(y) > 1) + 1])
  ends   <- c(y[which(diff(y) > 1)], y[length(y)])
  
  parts <- mapply(
    function(a, b) {
      if (a == b) {
        as.character(a)
      } else {
        paste0(a, "-", b)
      }
    },
    starts, ends,
    USE.NAMES = FALSE
  )
  
  paste(parts, collapse = ", ")
}

#' @keywords internal
#' @noRd
.rope_selected_regions <- function(x) {
  if (is.null(x$selected) || is.na(x$n_star)) {
    return(list(y_acc = integer(0), y_h0 = integer(0)))
  }
  
  n <- x$n_star
  p0 <- x$inputs$p0
  delta <- x$inputs$delta
  gamma_eq <- x$inputs$gamma_eq
  gamma_diff <- if (!is.null(x$inputs$gamma_diff)) x$inputs$gamma_diff else x$inputs$gamma_eq
  direction <- if (!is.null(x$inputs$direction)) x$inputs$direction else "equivalence"
  analysis_prior <- c(x$inputs$a, x$inputs$b)
  
  if (direction == "equivalence") {
    y_acc_all <- equivalence_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_eq = gamma_eq,
      analysis_prior = analysis_prior
    )
    y_h0_all <- nonequivalence_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior
    )
    y_acc <- y_acc_all
    y_h0  <- setdiff(y_h0_all, y_acc_all)
    
  } else if (direction == "noninferiority") {
    y_acc <- noninferiority_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_eq = gamma_eq,
      analysis_prior = analysis_prior
    )
    y_h0 <- inferiority_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior
    )
    
  } else { # superiority
    y_acc <- superiority_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_eq = gamma_eq,
      analysis_prior = analysis_prior
    )
    y_h0 <- nonsuperiority_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior
    )
  }
  
  list(y_acc = y_acc, y_h0 = y_h0)
}

#' @export
print.bfbin2arm_rope_design <- function(x, ...) {
  cat("One-stage single-arm ROPE design\n")
  if (!is.null(x$inputs$direction)) {
    cat("Direction:", x$inputs$direction, "\n")
  }
  cat("Calibration:", x$inputs$calibration, "\n")
  cat("Search range n:", x$inputs$n_min, "to", x$inputs$n_max, "\n")
  cat("Null probability p0:", x$inputs$p0, "\n")
  cat("Margin delta:", x$inputs$delta, "\n")
  cat("Probability threshold gamma_eq:", x$inputs$gamma_eq, "\n")
  if (!is.null(x$inputs$gamma_diff)) {
    cat("Probability threshold gamma_diff:", x$inputs$gamma_diff, "\n")
  }
  cat("Analysis prior: Beta(", x$inputs$a, ", ", x$inputs$b, ")\n", sep = "")
  cat("Design prior (H0): Beta(", x$inputs$da0, ", ", x$inputs$db0, ")\n", sep = "")
  cat("Design prior (H1): Beta(", x$inputs$da1, ", ", x$inputs$db1, ")\n", sep = "")
  
  if (!is.null(x$inputs$target_power)) {
    cat("Target Bayesian power:", x$inputs$target_power, "\n")
  }
  if (!is.null(x$inputs$target_type1)) {
    cat("Target Bayesian type-I error:", x$inputs$target_type1, "\n")
  }
  if (!is.null(x$inputs$target_pce_h0)) {
    cat("Target PCE(H0):", x$inputs$target_pce_h0, "\n")
  }
  if (!is.null(x$inputs$dp)) {
    cat("Frequentist power point dp:", x$inputs$dp, "\n")
  }
  if (!is.null(x$inputs$target_freq_power)) {
    cat("Target frequentist power:", x$inputs$target_freq_power, "\n")
  }
  if (!is.null(x$inputs$target_freq_type1)) {
    cat("Target frequentist type-I error:", x$inputs$target_freq_type1, "\n")
  }
  
  cat("Sustain n:", x$inputs$sustain_n, "\n")
  
  if (is.null(x$selected)) {
    cat("No feasible design found in the search range.\n")
    return(invisible(x))
  }
  
  cat("Selected sample size n*:", x$n_star, "\n")
  cat("Bayesian power(n*):", formatC(x$selected$power, digits = 4, format = "f"), "\n")
  cat("Bayesian type-I(n*):", formatC(x$selected$type1, digits = 4, format = "f"), "\n")
  
  if (!is.null(x$selected$pce_h0)) {
    cat("PCE(H0)(n*):", formatC(x$selected$pce_h0, digits = 4, format = "f"), "\n")
  }
  
  if (isTRUE(x$inputs$compute_freq_power)) {
    cat("Frequentist power(n*):", formatC(x$selected$freq_power, digits = 4, format = "f"), "\n")
  }
  
  if (isTRUE(x$inputs$compute_freq_type1)) {
    cat("Frequentist type-I(n*):", formatC(x$selected$freq_type1, digits = 4, format = "f"), "\n")
    if (!is.null(x$selected$freq_type1_lower) && !is.na(x$selected$freq_type1_lower)) {
      cat(" at p0 - delta:", formatC(x$selected$freq_type1_lower, digits = 4, format = "f"), "\n")
    }
    if (!is.null(x$selected$freq_type1_upper) && !is.na(x$selected$freq_type1_upper)) {
      cat(" at p0 + delta:", formatC(x$selected$freq_type1_upper, digits = 4, format = "f"), "\n")
    }
  }
  
  regs <- .rope_selected_regions(x)
  
  region_name <- if (!is.null(x$inputs$direction) && x$inputs$direction == "equivalence") {
    "Equivalence region"
  } else if (!is.null(x$inputs$direction) && x$inputs$direction == "noninferiority") {
    "Non-inferiority region"
  } else {
    "Superiority region"
  }
  
  h0_name <- if (!is.null(x$inputs$direction) && x$inputs$direction == "equivalence") {
    "Compelling evidence for non-equivalence region"
  } else if (!is.null(x$inputs$direction) && x$inputs$direction == "noninferiority") {
    "Compelling evidence for inferiority region"
  } else {
    "Compelling evidence for non-superiority region"
  }
  
  cat(region_name, ": {", .rope_compact_runs(regs$y_acc), "}\n", sep = "")
  cat(h0_name, ": {", .rope_compact_runs(regs$y_h0), "}\n", sep = "")
  
  invisible(x)
}

#' @export
summary.bfbin2arm_rope_design <- function(object, ...) {
  out <- list(
    inputs = object$inputs,
    selected = object$selected,
    head = if (!is.null(object$grid)) utils::head(object$grid, 10) else NULL,
    tail = if (!is.null(object$grid)) utils::tail(object$grid, 10) else NULL
  )
  class(out) <- "summary.bfbin2arm_rope_design"
  out
}

#' @export
print.summary.bfbin2arm_rope_design <- function(x, ...) {
  cat("Summary of one-stage single-arm ROPE design\n")
  cat("Direction:", x$inputs$direction, "\n")
  cat("Calibration:", x$inputs$calibration, "\n")
  
  if (!is.null(x$selected)) {
    cat("Selected design point:\n")
    print(x$selected)
  } else {
    cat("No feasible design found.\n")
  }
  
  if (!is.null(x$head)) {
    cat("\nFirst rows of evaluation grid:\n")
    print(x$head)
    cat("\nLast rows of evaluation grid:\n")
    print(x$tail)
  }
  
  invisible(x)
}

.rope_shaded_region <- function(p0, delta, direction) {
  if (direction == "equivalence") {
    c(max(0, p0 - delta), min(1, p0 + delta))
  } else if (direction == "noninferiority") {
    c(max(0, p0 - delta), 1)
  } else {
    c(min(1, p0 + delta), 1)
  }
}

.rope_beta_panel <- function(shape1, shape2, main, p0, delta, direction, col = 4,
                             legend_pos = "topright") {
  xx <- seq(0, 1, length.out = 1000)
  yy <- dbeta(xx, shape1 = shape1, shape2 = shape2)
  rg <- .rope_shaded_region(p0, delta, direction)
  shade_col <- grDevices::adjustcolor("gray80", alpha.f = 0.5)
  yref_col <- grDevices::adjustcolor("gray60", alpha.f = 0.8)
  ymax <- max(yy)
  
  plot(xx, yy, type = "l", lwd = 2, col = col,
       xlab = "p", ylab = "Density", main = main,
       ylim = c(0, ymax * 1.05))
  rect(rg[1], 0, rg[2], ymax * 1.05, border = NA, col = shade_col)
  lines(xx, yy, lwd = 2, col = col)
  abline(v = p0, lty = 3, lwd = 1.5)
  
  legend(legend_pos,
         legend = c("Prior density", .rope_region_label(direction), expression(p[0])),
         col = c(col, yref_col, 1),
         lwd = c(2, 6, 1.5),
         lty = c(1, 1, 3),
         bty = "n")
}

.rope_beta_panel_two <- function(a1, b1, a2, b2, p0, delta, direction,
                                 main = "Design priors",
                                 col1 = 2, col2 = 4,
                                 legend_pos = "topright") {
  xx <- seq(0, 1, length.out = 1000)
  yy1 <- dbeta(xx, shape1 = a1, shape2 = b1)
  yy2 <- dbeta(xx, shape1 = a2, shape2 = b2)
  rg <- .rope_shaded_region(p0, delta, direction)
  shade_col <- grDevices::adjustcolor("gray80", alpha.f = 0.5)
  yref_col <- grDevices::adjustcolor("gray60", alpha.f = 0.8)
  ymax <- max(c(yy1, yy2))
  
  plot(xx, yy1, type = "l", lwd = 2, col = col1,
       xlab = "p", ylab = "Density", ylim = c(0, ymax * 1.05),
       main = main)
  rect(rg[1], 0, rg[2], ymax * 1.05, border = NA, col = shade_col)
  lines(xx, yy1, lwd = 2, col = col1)
  lines(xx, yy2, lwd = 2, col = col2)
  abline(v = p0, lty = 3, lwd = 1.5)
  
  legend(legend_pos,
         legend = c(
           paste0("Design prior (H0): Beta(", a1, ", ", b1, ")"),
           paste0("Design prior (H1): Beta(", a2, ", ", b2, ")"),
           .rope_region_label(direction),
           expression(p[0])
         ),
         col = c(col1, col2, yref_col, 1),
         lwd = c(2, 2, 6, 1.5),
         lty = c(1, 1, 1, 3),
         bty = "n")
}

#' @importFrom graphics segments
.rope_summary_panel <- function(x, cex_label = 0.88, cex_value = 0.88) {
  plot.new()
  title("Design summary")
  
  fmt4 <- function(z) {
    if (is.null(z) || is.na(z)) "NA" else sprintf("%.4f", z)
  }
  
  region_label <- .rope_region_label(x$inputs$direction)
  
  labels <- list(
    "Direction",
    "Calibration",
    expression("Null probability " * p[0]),
    expression("Margin " * delta),
    expression("Threshold " * gamma[eq]),
    expression("Threshold " * gamma[diff]),
    "Analysis prior",
    expression("Design prior (" * H[0] * ")"),
    expression("Design prior (" * H[1] * ")"),
    "Target Bayesian power",
    "Target Bayesian type-I error",
    "Target PCE(H0)",
    expression("Point alternative " * d[p]),
    "Target frequentist power",
    "Target frequentist type-I error",
    "Sustain n",
    expression("Selected sample size " * n^"*"),
    expression("Bayesian power(" * n^"*" * ")"),
    expression("Bayesian type-I(" * n^"*" * ")"),
    "PCE(H0)",
    expression("Frequentist power(" * n^"*" * ")"),
    expression("Frequentist type-I(" * n^"*" * ")"),
    paste0(tools::toTitleCase(region_label), " region")
  )
  
  values <- list(
    x$inputs$direction,
    x$inputs$calibration,
    x$inputs$p0,
    x$inputs$delta,
    x$inputs$gamma_eq,
    x$inputs$gamma_diff,
    paste0("Beta(", x$inputs$a, ", ", x$inputs$b, ")"),
    paste0("Beta(", x$inputs$da0, ", ", x$inputs$db0, ")"),
    paste0("Beta(", x$inputs$da1, ", ", x$inputs$db1, ")"),
    if (!is.null(x$inputs$target_power)) fmt4(x$inputs$target_power) else "-",
    if (!is.null(x$inputs$target_type1)) fmt4(x$inputs$target_type1) else "-",
    if (!is.null(x$inputs$target_pce_h0)) fmt4(x$inputs$target_pce_h0) else "-",
    if (!is.null(x$inputs$dp)) x$inputs$dp else "-",
    if (!is.null(x$inputs$target_freq_power)) fmt4(x$inputs$target_freq_power) else "-",
    if (!is.null(x$inputs$target_freq_type1)) fmt4(x$inputs$target_freq_type1) else "-",
    x$inputs$sustain_n,
    if (!is.na(x$n_star)) x$n_star else "none",
    if (!is.null(x$selected)) fmt4(x$selected$power) else "NA",
    if (!is.null(x$selected)) fmt4(x$selected$type1) else "NA",
    if (!is.null(x$selected)) fmt4(x$selected$pce_h0) else "NA",
    if (!is.null(x$selected) && isTRUE(x$inputs$compute_freq_power)) fmt4(x$selected$freq_power) else "-",
    if (!is.null(x$selected) && isTRUE(x$inputs$compute_freq_type1)) fmt4(x$selected$freq_type1) else "-",
    if (!is.null(x$selected)) paste0("[", x$selected$y_acc_min, ", ", x$selected$y_acc_max, "]") else "NA"
  )
  
  yloc <- seq(0.96, 0.05, length.out = length(labels))
  x_label <- 0.03
  x_value <- 0.72
  
  for (i in seq_along(labels)) {
    text(x = x_label, y = yloc[i], labels = labels[[i]], adj = c(0, 0.5), cex = cex_label)
    text(x = x_value, y = yloc[i], labels = as.character(values[[i]]), adj = c(0, 0.5), cex = cex_value)
  }
  
  segments(x0 = 0.68, y0 = 0.04, x1 = 0.68, y1 = 0.97, col = "gray70")
}

#' @export
plot.bfbin2arm_rope_design <- function(
    x,
    what = c("overview", "operating_characteristics", "decision_region"),
    ...
) {
  what <- match.arg(what)
  if (is.null(x$grid)) stop("Plotting requires return_grid = TRUE.")
  
  region_label <- .rope_region_label(x$inputs$direction)
  
  col_bayes_power <- "#0072B2"
  col_bayes_type1 <- "#D55E00"
  col_pce_h0 <- "#6A3D9A"
  col_freq_power <- "#009E73"
  col_freq_type1 <- "#CC79A7"
  col_ref <- "gray55"
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (what == "operating_characteristics") {
    ylim_all <- c(x$grid$power, x$grid$type1, x$grid$pce_h0)
    if (isTRUE(x$inputs$compute_freq_power)) {
      ylim_all <- c(ylim_all, x$grid$freq_power)
    }
    if (isTRUE(x$inputs$compute_freq_type1)) {
      ylim_all <- c(ylim_all, x$grid$freq_type1)
    }
    
    plot(
      x$grid$n, x$grid$power,
      type = "l", lwd = 2, lty = 1,
      ylim = range(ylim_all, na.rm = TRUE),
      xlab = "Sample size n", ylab = "Probability",
      col = col_bayes_power
    )
    lines(x$grid$n, x$grid$type1, lwd = 2, lty = 1, col = col_bayes_type1)
    lines(x$grid$n, x$grid$pce_h0, lwd = 2, lty = 1, col = col_pce_h0)
    
    if (isTRUE(x$inputs$compute_freq_power)) {
      lines(x$grid$n, x$grid$freq_power, lwd = 2, lty = 1, col = col_freq_power)
    }
    if (isTRUE(x$inputs$compute_freq_type1)) {
      lines(x$grid$n, x$grid$freq_type1, lwd = 2, lty = 1, col = col_freq_type1)
    }
    
    if (!is.null(x$inputs$target_power)) abline(h = x$inputs$target_power, lty = 1, lwd = 1, col = col_ref)
    if (!is.null(x$inputs$target_type1)) abline(h = x$inputs$target_type1, lty = 1, lwd = 1, col = col_ref)
    if (!is.null(x$inputs$target_pce_h0)) abline(h = x$inputs$target_pce_h0, lty = 1, lwd = 1, col = col_ref)
    if (!is.null(x$inputs$target_freq_power)) abline(h = x$inputs$target_freq_power, lty = 1, lwd = 1, col = col_ref)
    if (!is.null(x$inputs$target_freq_type1)) abline(h = x$inputs$target_freq_type1, lty = 1, lwd = 1, col = col_ref)
    if (!is.na(x$n_star)) abline(v = x$n_star, lty = 1, lwd = 1, col = "gray70")
    
    leg <- c("Bayesian power under H1", "Bayesian type-I under H0", "PCE(H0) under H0")
    cols <- c(col_bayes_power, col_bayes_type1, col_pce_h0)
    lwds <- c(2, 2, 2)
    
    if (isTRUE(x$inputs$compute_freq_power)) {
      leg <- c(leg, "Frequentist power at dp")
      cols <- c(cols, col_freq_power)
      lwds <- c(lwds, 2)
    }
    if (isTRUE(x$inputs$compute_freq_type1)) {
      leg <- c(leg, "Frequentist type-I (boundary worst-case)")
      cols <- c(cols, col_freq_type1)
      lwds <- c(lwds, 2)
    }
    
    legend("right", legend = leg, col = cols, lty = 1, lwd = lwds, bty = "n")
    
  } else if (what == "decision_region") {
    ok <- !is.na(x$grid$y_acc_min)
    
    plot(
      x$grid$n[ok], x$grid$y_acc_min[ok],
      type = "l", lwd = 2, lty = 1,
      xlab = "Sample size n", ylab = "Responder count y",
      col = col_bayes_power,
      ylim = range(c(x$grid$y_acc_min[ok], x$grid$y_acc_max[ok]), na.rm = TRUE)
    )
    lines(x$grid$n[ok], x$grid$y_acc_max[ok], lwd = 2, lty = 1, col = col_bayes_type1)
    
    legend(
      "topleft",
      legend = c(
        paste0("Min y for ", region_label),
        paste0("Max y for ", region_label)
      ),
      col = c(col_bayes_power, col_bayes_type1),
      lwd = 2, lty = 1, bty = "n"
    )
    
  } else {
    layout(
      matrix(c(1, 1, 2,
               3, 3, 4), nrow = 2, byrow = TRUE),
      widths = c(1, 1, 1),
      heights = c(1, 1)
    )
    
    par(mar = c(4, 4, 3, 1))
    ylim_all <- c(x$grid$power, x$grid$type1, x$grid$pce_h0)
    if (isTRUE(x$inputs$compute_freq_power)) {
      ylim_all <- c(ylim_all, x$grid$freq_power)
    }
    if (isTRUE(x$inputs$compute_freq_type1)) {
      ylim_all <- c(ylim_all, x$grid$freq_type1)
    }
    
    plot(
      x$grid$n, x$grid$power,
      type = "l", lwd = 2, lty = 1,
      ylim = range(ylim_all, na.rm = TRUE),
      xlab = "Sample size n", ylab = "Probability",
      col = col_bayes_power,
      main = "Operating characteristics"
    )
    lines(x$grid$n, x$grid$type1, lwd = 2, lty = 1, col = col_bayes_type1)
    lines(x$grid$n, x$grid$pce_h0, lwd = 2, lty = 1, col = col_pce_h0)
    
    if (isTRUE(x$inputs$compute_freq_power)) {
      lines(x$grid$n, x$grid$freq_power, lwd = 2, lty = 1, col = col_freq_power)
    }
    if (isTRUE(x$inputs$compute_freq_type1)) {
      lines(x$grid$n, x$grid$freq_type1, lwd = 2, lty = 1, col = col_freq_type1)
    }
    
    if (!is.null(x$inputs$target_power)) abline(h = x$inputs$target_power, lty = 1, lwd = 1, col = col_ref)
    if (!is.null(x$inputs$target_type1)) abline(h = x$inputs$target_type1, lty = 1, lwd = 1, col = col_ref)
    if (!is.null(x$inputs$target_pce_h0)) abline(h = x$inputs$target_pce_h0, lty = 1, lwd = 1, col = col_ref)
    if (!is.null(x$inputs$target_freq_power)) abline(h = x$inputs$target_freq_power, lty = 1, lwd = 1, col = col_ref)
    if (!is.null(x$inputs$target_freq_type1)) abline(h = x$inputs$target_freq_type1, lty = 1, lwd = 1, col = col_ref)
    if (!is.na(x$n_star)) abline(v = x$n_star, lty = 1, lwd = 1, col = "gray70")
    
    leg <- c("Bayesian power under H1", "Bayesian type-I under H0", "PCE(H0) under H0")
    cols <- c(col_bayes_power, col_bayes_type1, col_pce_h0)
    lwds <- c(2, 2, 2)
    
    if (isTRUE(x$inputs$compute_freq_power)) {
      leg <- c(leg, "Frequentist power at dp")
      cols <- c(cols, col_freq_power)
      lwds <- c(lwds, 2)
    }
    if (isTRUE(x$inputs$compute_freq_type1)) {
      leg <- c(leg, "Frequentist type-I (boundary worst-case)")
      cols <- c(cols, col_freq_type1)
      lwds <- c(lwds, 2)
    }
    
    legend("right", legend = leg, col = cols, lty = 1, lwd = lwds, bty = "n")
    
    par(mar = c(1, 1, 3, 1))
    .rope_summary_panel(x)
    
    par(mar = c(4, 4, 3, 1))
    .rope_beta_panel_two(
      x$inputs$da0, x$inputs$db0,
      x$inputs$da1, x$inputs$db1,
      p0 = x$inputs$p0,
      delta = x$inputs$delta,
      direction = x$inputs$direction,
      main = "Design priors",
      col1 = col_bayes_type1,
      col2 = col_bayes_power,
      legend_pos = "topright"
    )
    
    par(mar = c(4, 4, 3, 1))
    .rope_beta_panel(
      x$inputs$a, x$inputs$b,
      p0 = x$inputs$p0,
      delta = x$inputs$delta,
      direction = x$inputs$direction,
      main = "Analysis prior",
      col = "black",
      legend_pos = "bottomright"
    )
  }
  
  invisible(x)
}