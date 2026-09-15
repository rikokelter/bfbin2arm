# S3 methods for two-arm one-stage ROPE equivalence designs.
# Place this file in R/s3-twoarm_rope_onestage.R.

.format_rope2arm_number <- function(x, digits = 4L) {
  if (is.null(x) ||
      length(x) != 1L ||
      !is.finite(x)) {
    return("NA")
  }
  
  formatC(
    x,
    digits = digits,
    format = "f"
  )
}


.rope2arm_prior_label <- function(prior) {
  paste0("Beta(", prior[1L], ", ", prior[2L], ")")
}


.rope2arm_selected_fit <- function(x) {
  if (!is.null(x$selected_fit)) {
    return(x$selected_fit)
  }
  
  NULL
}


#' @export
print.bfbin2armrope2armdesign <- function(x, ...) {
  cat("One-stage two-arm ROPE equivalence design\n")
  cat("Equal allocation: nC = nT\n")
  cat(
    "Search range per arm:",
    x$inputs$nmin,
    "to",
    x$inputs$nmax,
    "(step",
    x$inputs$nstep,
    ")\n"
  )
  cat("ROPE half-width delta:", x$inputs$delta, "\n")
  cat(
    "Thresholds: gammaeq =",
    x$inputs$gammaeq,
    ", gammadiff =",
    x$inputs$gammadiff,
    "\n"
  )
  
  cat("Integration:", x$inputs$integration)
  if (identical(x$inputs$integration, "quantile")) {
    cat("(quantile quadrature; nodes =", x$inputs$quad_nodes, ")")
  }
  cat("\n")
  
  cat("Calibration:", x$inputs$calibration, "\n")
  
  if (!is.null(x$inputs$targetpower) &&
      is.finite(x$inputs$targetpower)) {
    cat(
      "Target PCE(H1):",
      .format_rope2arm_number(x$inputs$targetpower),
      "\n"
    )
  }
  
  if (!is.null(x$inputs$targettype1) &&
      is.finite(x$inputs$targettype1)) {
    cat(
      "Maximum predictive false-equivalence probability:",
      .format_rope2arm_number(x$inputs$targettype1),
      "\n"
    )
  }
  
  if (!is.null(x$inputs$targetpce_h0) &&
      is.finite(x$inputs$targetpce_h0)) {
    cat(
      "Target PCE(H0):",
      .format_rope2arm_number(x$inputs$targetpce_h0),
      "\n"
    )
  }
  
  if (!is.null(x$inputs$targetfreqpower) &&
      is.finite(x$inputs$targetfreqpower)) {
    cat(
      "Target frequentist power:",
      .format_rope2arm_number(x$inputs$targetfreqpower),
      "\n"
    )
  }
  
  if (!is.null(x$inputs$targetfreqtype1) &&
      is.finite(x$inputs$targetfreqtype1)) {
    cat(
      "Maximum frequentist type-I error:",
      .format_rope2arm_number(x$inputs$targetfreqtype1),
      "\n"
    )
  }
  
  cat(
    "Sustain n:",
    x$inputs$sustainn,
    "\n"
  )
  
  if (is.null(x$selected)) {
    cat("\nNo sustained feasible design found in the search range.\n")
    return(invisible(x))
  }
  
  fit <- .rope2arm_selected_fit(x)
  
  cat("\nSelected design\n")
  cat(
    "nC =",
    fit$nC,
    ", nT =",
    fit$nT,
    ", total N =",
    fit$N,
    "\n"
  )
  
  cat(
    "PCE(H1), practical equivalence:",
    .format_rope2arm_number(fit$power),
    "\n"
  )
  
  cat(
    "Predictive false-equivalence probability:",
    .format_rope2arm_number(fit$type1),
    "\n"
  )
  
  cat(
    "P(meaningful difference | G1):",
    .format_rope2arm_number(fit$p_diff_h1),
    "\n"
  )
  
  cat(
    "P(inconclusive | G1):",
    .format_rope2arm_number(fit$p_inc_h1),
    "\n"
  )
  
  cat(
    "PCE(H0), meaningful non-equivalence:",
    .format_rope2arm_number(fit$p_diff_h0),
    "\n"
  )
  
  cat(
    "P(inconclusive | G0):",
    .format_rope2arm_number(fit$p_inc_h0),
    "\n"
  )
  
  if (!is.null(x$inputs$targetpce_h0) &&
      is.finite(x$inputs$targetpce_h0)) {
    cat(
      "PCE(H0) target satisfied:",
      isTRUE(
        fit$p_diff_h0 >= x$inputs$targetpce_h0
      ),
      "\n"
    )
  }
  
  if (!is.null(fit$freqpower) &&
      is.finite(fit$freqpower)) {
    cat(
      "Pointwise frequentist equivalence probability:",
      .format_rope2arm_number(fit$freqpower),
      "\n"
    )
  }
  
  if (!is.null(fit$freqtype1) &&
      is.finite(fit$freqtype1)) {
    cat(
      "Boundary-grid frequentist type-I error:",
      .format_rope2arm_number(fit$freqtype1),
      "\n"
    )
  }
  
  invisible(x)
}


#' @export
summary.bfbin2armrope2armdesign <- function(object, ...) {
  out <- list(
    inputs = object$inputs,
    nstar = object$nstar,
    Nstar = object$Nstar,
    selected = object$selected,
    selected_fit = object$selected_fit,
    
    first_pointwise_feasible = if (!is.null(object$grid) &&
                                   any(object$grid$feasible_pointwise)) {
      object$grid$n[which(object$grid$feasible_pointwise)[1L]]
    } else {
      NA_integer_
    },
    
    first_sustained_feasible = if (!is.null(object$grid) &&
                                   any(object$grid$feasible)) {
      object$grid$n[which(object$grid$feasible)[1L]]
    } else {
      NA_integer_
    }
  )
  
  class(out) <- "summary.bfbin2armrope2armdesign"
  out
}


#' @export
print.summary.bfbin2armrope2armdesign <- function(x, ...) {
  cat("Summary: one-stage two-arm ROPE equivalence design\n")
  
  cat(
    "First pointwise feasible per-arm n:",
    x$first_pointwise_feasible,
    "\n"
  )
  
  cat(
    "First sustained feasible per-arm n:",
    x$first_sustained_feasible,
    "\n"
  )
  
  if (is.null(x$selected)) {
    cat("No sustained feasible design found.\n")
  } else {
    print(x$selected)
  }
  
  invisible(x)
}


.rope2arm_plot_operating_characteristics <- function(x) {
  grid <- x$grid
  
  if (is.null(grid) || nrow(grid) == 0L) {
    graphics::plot.new()
    graphics::title("Operating characteristics")
    
    graphics::text(
      0.5,
      0.5,
      "No grid is available.\nSet returngrid = TRUE.",
      cex = 0.9
    )
    
    return(invisible(NULL))
  }
  
  cols <- c(
    bayes_power = "#0072B2",
    bayes_type1 = "#D55E00",
    bayes_difference_h1 = "#009E73",
    bayes_inconclusive_h1 = "#56B4E9",
    bayes_difference_h0 = "#CC79A7",
    bayes_inconclusive_h0 = "#999999",
    freq_power = "#000000",
    freq_type1 = "#E69F00"
  )
  
  has_freq_power <- "freqpower" %in% names(grid) &&
    any(is.finite(grid$freqpower))
  
  has_freq_type1 <- "freqtype1" %in% names(grid) &&
    any(is.finite(grid$freqtype1))
  
  y_values <- c(
    grid$power,
    grid$type1,
    grid$p_diff_h1,
    grid$p_inc_h1,
    grid$p_diff_h0,
    grid$p_inc_h0
  )
  
  if (has_freq_power) {
    y_values <- c(
      y_values,
      grid$freqpower
    )
  }
  
  if (has_freq_type1) {
    y_values <- c(
      y_values,
      grid$freqtype1
    )
  }
  
  y_values <- y_values[is.finite(y_values)]
  
  ylim <- if (length(y_values) == 0L) {
    c(0, 1)
  } else {
    c(0, max(1, max(y_values)))
  }
  
  graphics::plot(
    grid$n,
    grid$power,
    type = "l",
    lwd = 2,
    lty = 1,
    col = cols["bayes_power"],
    ylim = ylim,
    xlab = "Per-arm sample size",
    ylab = "Probability",
    main = "Operating characteristics"
  )
  
  graphics::lines(
    grid$n,
    grid$type1,
    lwd = 2,
    lty = 1,
    col = cols["bayes_type1"]
  )
  
  graphics::lines(
    grid$n,
    grid$p_diff_h1,
    lwd = 1.8,
    lty = 1,
    col = cols["bayes_difference_h1"]
  )
  
  graphics::lines(
    grid$n,
    grid$p_inc_h1,
    lwd = 1.8,
    lty = 1,
    col = cols["bayes_inconclusive_h1"]
  )
  
  graphics::lines(
    grid$n,
    grid$p_diff_h0,
    lwd = 1.8,
    lty = 1,
    col = cols["bayes_difference_h0"]
  )
  
  graphics::lines(
    grid$n,
    grid$p_inc_h0,
    lwd = 1.8,
    lty = 1,
    col = cols["bayes_inconclusive_h0"]
  )
  
  if (has_freq_power) {
    graphics::lines(
      grid$n,
      grid$freqpower,
      lwd = 3,
      lty = 1,
      col = cols["freq_power"]
    )
  }
  
  if (has_freq_type1) {
    graphics::lines(
      grid$n,
      grid$freqtype1,
      lwd = 3,
      lty = 1,
      col = cols["freq_type1"]
    )
  }
  
  if (!is.null(x$inputs$targetpower) &&
      length(x$inputs$targetpower) == 1L &&
      is.finite(x$inputs$targetpower)) {
    graphics::abline(
      h = x$inputs$targetpower,
      lty = 3,
      lwd = 1,
      col = cols["bayes_power"]
    )
  }
  
  if (!is.null(x$inputs$targettype1) &&
      length(x$inputs$targettype1) == 1L &&
      is.finite(x$inputs$targettype1)) {
    graphics::abline(
      h = x$inputs$targettype1,
      lty = 3,
      lwd = 1,
      col = cols["bayes_type1"]
    )
  }
  
  if (!is.null(x$inputs$targetpce_h0) &&
      length(x$inputs$targetpce_h0) == 1L &&
      is.finite(x$inputs$targetpce_h0)) {
    graphics::abline(
      h = x$inputs$targetpce_h0,
      lty = 3,
      lwd = 1.2,
      col = cols["bayes_difference_h0"]
    )
  }
  
  if (!is.null(x$inputs$targetfreqpower) &&
      length(x$inputs$targetfreqpower) == 1L &&
      is.finite(x$inputs$targetfreqpower)) {
    graphics::abline(
      h = x$inputs$targetfreqpower,
      lty = 3,
      lwd = 1.2,
      col = cols["freq_power"]
    )
  }
  
  if (!is.null(x$inputs$targetfreqtype1) &&
      length(x$inputs$targetfreqtype1) == 1L &&
      is.finite(x$inputs$targetfreqtype1)) {
    graphics::abline(
      h = x$inputs$targetfreqtype1,
      lty = 3,
      lwd = 1.2,
      col = cols["freq_type1"]
    )
  }
  
  if (!is.na(x$nstar)) {
    graphics::abline(
      v = x$nstar,
      lty = 3,
      lwd = 1,
      col = "grey35"
    )
  }
  
  legend_labels <- c(
    "Bayesian predictive power",
    "Bayesian predictive type-I error",
    "P(difference | equivalence)",
    "P(inconclusive | equivalence)",
    "PCE(H0): compelling non-equivalence",
    "P(inconclusive | non-equivalence)"
  )
  
  legend_cols <- c(
    cols["bayes_power"],
    cols["bayes_type1"],
    cols["bayes_difference_h1"],
    cols["bayes_inconclusive_h1"],
    cols["bayes_difference_h0"],
    cols["bayes_inconclusive_h0"]
  )
  
  legend_lwd <- c(
    2,
    2,
    1.8,
    1.8,
    1.8,
    1.8
  )
  
  if (has_freq_power) {
    legend_labels <- c(
      legend_labels,
      "Frequentist power"
    )
    
    legend_cols <- c(
      legend_cols,
      cols["freq_power"]
    )
    
    legend_lwd <- c(
      legend_lwd,
      3
    )
  }
  
  if (has_freq_type1) {
    legend_labels <- c(
      legend_labels,
      "Frequentist type-I error"
    )
    
    legend_cols <- c(
      legend_cols,
      cols["freq_type1"]
    )
    
    legend_lwd <- c(
      legend_lwd,
      3
    )
  }
  
  graphics::legend(
    "right",
    legend = legend_labels,
    col = legend_cols,
    lty = 1,
    lwd = legend_lwd,
    cex = 0.62,
    bty = "n"
  )
  
  invisible(x)
}


.rope2arm_plot_joint_design_prior <- function(
    prior_C,
    prior_T,
    delta,
    region = c("H0", "H1"),
    title_text,
    n_grid = 151L
) {
  region <- match.arg(region)
  
  pC <- seq(0, 1, length.out = n_grid)
  pT <- seq(0, 1, length.out = n_grid)
  
  in_rope <- outer(
    pC,
    pT,
    FUN = function(pc, pt) abs(pt - pc) < delta
  )
  
  base_density <- outer(
    pC,
    pT,
    FUN = function(pc, pt) {
      dbeta(pc, prior_C[1L], prior_C[2L]) *
        dbeta(pt, prior_T[1L], prior_T[2L])
    }
  )
  
  support <- if (region == "H1") in_rope else !in_rope
  
  density <- base_density
  density[!support] <- NA_real_
  
  col_h1 <- grDevices::adjustcolor("#BBDDF0", alpha.f = 0.65)
  col_h0 <- grDevices::adjustcolor("#F3C0C0", alpha.f = 0.65)
  
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    asp = 1,
    xlab = expression(p[C]),
    ylab = expression(p[T]),
    main = title_text
  )
  
  background <- ifelse(in_rope, 1L, 2L)
  
  image(
    x = pC,
    y = pT,
    z = background,
    add = TRUE,
    axes = FALSE,
    col = c(col_h1, col_h0),
    breaks = c(0.5, 1.5, 2.5),
    useRaster = TRUE
  )
  
  max_density <- max(density, na.rm = TRUE)
  
  if (is.finite(max_density) && max_density > 0) {
    image(
      x = pC,
      y = pT,
      z = density,
      add = TRUE,
      axes = FALSE,
      col = hcl.colors(64L, palette = "YlOrRd", rev = TRUE),
      breaks = seq(0, max_density, length.out = 65L),
      useRaster = TRUE
    )
  }
  
  abline(a = delta, b = 1, lty = 2, lwd = 1.4, col = "grey20")
  abline(a = -delta, b = 1, lty = 2, lwd = 1.4, col = "grey20")
  abline(a = 0, b = 1, lty = 3, lwd = 1.0, col = "grey35")
  
  legend(
    "topright",
    legend = c(
      expression(H[1]~":"~abs(p[T] - p[C]) < delta),
      expression(H[0]~":"~abs(p[T] - p[C]) >= delta),
      "ROPE boundaries",
      expression(p[T] == p[C]),
      "Product-beta kernel on design-prior support"
    ),
    fill = c(col_h1, col_h0, NA, NA, NA),
    border = NA,
    col = c(NA, NA, "grey20", "grey35", "#8C510A"),
    lty = c(NA, NA, 2, 3, 1),
    lwd = c(NA, NA, 1.4, 1.0, 5),
    cex = 0.64,
    bty = "n"
  )
  
  mtext(
    paste0(
      "Control kernel: ",
      .rope2arm_prior_label(prior_C),
      "   Treatment kernel: ",
      .rope2arm_prior_label(prior_T)
    ),
    side = 3,
    line = 0.25,
    cex = 0.66
  )
  
  invisible(NULL)
}


.rope2arm_plot_baseline_difference_prior <- function(
    design_prior_eq,
    delta,
    main = "Design prior under H1: equivalence"
) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  graphics::layout(matrix(c(1, 2), nrow = 1L))
  
  baseline <- design_prior_eq$baseline
  difference <- design_prior_eq$difference
  
  p <- seq(0, 1, length.out = 1000L)
  
  d_baseline <- dbeta(
    p,
    shape1 = baseline$shape1,
    shape2 = baseline$shape2
  )
  
  graphics::par(mar = c(4.2, 4.2, 3.2, 1.0))
  
  plot(
    p,
    d_baseline,
    type = "l",
    lwd = 2,
    col = "#0072B2",
    xlab = expression(p[C]),
    ylab = "Density",
    main = "Baseline response prior"
  )
  
  abline(
    v = baseline$shape1 /
      (baseline$shape1 + baseline$shape2),
    lty = 3,
    col = "grey35"
  )
  
  legend(
    "topright",
    legend = c(
      paste0(
        "Beta(",
        baseline$shape1,
        ", ",
        baseline$shape2,
        ")"
      ),
      paste0(
        "Mean = ",
        formatC(
          baseline$shape1 /
            (baseline$shape1 + baseline$shape2),
          digits = 3,
          format = "f"
        )
      )
    ),
    col = c("#0072B2", "grey35"),
    lty = c(1, 3),
    lwd = c(2, 1),
    bty = "n",
    cex = 0.72
  )
  
  truncation <- difference$truncation
  
  d_min <- max(-0.20, truncation[1L] - 0.02)
  d_max <- min(0.20, truncation[2L] + 0.02)
  
  d <- seq(d_min, d_max, length.out = 1500L)
  
  lower_cdf <- pnorm(
    truncation[1L],
    mean = difference$mean,
    sd = difference$sd
  )
  
  upper_cdf <- pnorm(
    truncation[2L],
    mean = difference$mean,
    sd = difference$sd
  )
  
  d_difference <- dnorm(
    d,
    mean = difference$mean,
    sd = difference$sd
  ) / (upper_cdf - lower_cdf)
  
  d_difference[d < truncation[1L] |
                 d > truncation[2L]] <- 0
  
  graphics::par(mar = c(4.2, 4.2, 3.2, 1.0))
  
  plot(
    d,
    d_difference,
    type = "l",
    lwd = 2,
    col = "#009E73",
    xlab = expression(Delta == p[T] - p[C]),
    ylab = "Density",
    main = "Treatment-difference prior"
  )
  
  polygon(
    x = c(-delta, delta, delta, -delta),
    y = c(0, 0, max(d_difference) * 1.05,
          max(d_difference) * 1.05),
    col = grDevices::adjustcolor("#BBDDF0", alpha.f = 0.35),
    border = NA
  )
  
  lines(
    d,
    d_difference,
    lwd = 2,
    col = "#009E73"
  )
  
  abline(
    v = c(-delta, delta),
    lty = 2,
    lwd = 1.2,
    col = "grey25"
  )
  
  abline(
    v = 0,
    lty = 3,
    lwd = 1,
    col = "grey35"
  )
  
  legend(
    "topright",
    legend = c(
      paste0(
        "Truncated N(",
        difference$mean,
        ", ",
        difference$sd,
        "^2)"
      ),
      "ROPE",
      expression(Delta == 0)
    ),
    col = c("#009E73", "#BBDDF0", "grey35"),
    lty = c(1, NA, 3),
    lwd = c(2, NA, 1),
    pch = c(NA, 15, NA),
    pt.cex = c(NA, 2, NA),
    bty = "n",
    cex = 0.72
  )
  
  mtext(
    main,
    side = 3,
    outer = TRUE,
    line = -1.2,
    cex = 0.90
  )
  
  invisible(NULL)
}

.rope2arm_plot_baseline_difference_prior_overview <- function(
    design_prior_eq,
    delta,
    main = "Design prior under H1: equivalence"
) {
  baseline <- design_prior_eq$baseline
  difference <- design_prior_eq$difference
  
  p <- seq(0, 1, length.out = 1000L)
  
  d_baseline <- dbeta(
    p,
    shape1 = baseline$shape1,
    shape2 = baseline$shape2
  )
  
  plot(
    p,
    d_baseline,
    type = "l",
    lwd = 2,
    col = "#0072B2",
    xlab = expression(p[C]),
    ylab = "Density",
    main = main
  )
  
  baseline_mean <- baseline$shape1 /
    (baseline$shape1 + baseline$shape2)
  
  abline(
    v = baseline_mean,
    lty = 3,
    lwd = 1,
    col = "grey35"
  )
  
  truncation <- difference$truncation
  
  delta_grid <- seq(
    truncation[1L],
    truncation[2L],
    length.out = 1000L
  )
  
  lower_cdf <- pnorm(
    truncation[1L],
    mean = difference$mean,
    sd = difference$sd
  )
  
  upper_cdf <- pnorm(
    truncation[2L],
    mean = difference$mean,
    sd = difference$sd
  )
  
  d_delta <- dnorm(
    delta_grid,
    mean = difference$mean,
    sd = difference$sd
  ) / (upper_cdf - lower_cdf)
  
  delta_scale <- max(d_delta) / max(d_baseline)
  
  scaled_delta <- d_delta / delta_scale
  
  lines(
    baseline_mean + delta_grid,
    scaled_delta,
    lwd = 2,
    col = "#009E73"
  )
  
  abline(
    v = baseline_mean - delta,
    lty = 2,
    lwd = 1,
    col = "grey65"
  )
  
  abline(
    v = baseline_mean + delta,
    lty = 2,
    lwd = 1,
    col = "grey65"
  )
  
  legend(
    "topright",
    legend = c(
      paste0(
        expression(p[C] %~% Beta),
        "(",
        baseline$shape1,
        ", ",
        baseline$shape2,
        ")"
      ),
      expression(Delta %~% N(0, 0.03^2)),
      "Vertical line: E[pC]",
      "Dashed lines: E[pC] +/- delta"
    ),
    col = c("#0072B2", "#009E73", "grey35", "grey65"),
    lty = c(1, 1, 3, 2),
    lwd = c(2, 2, 1, 1),
    bty = "n",
    cex = 0.58
  )
  
  mtext(
    paste0(
      "Difference prior: N(",
      difference$mean,
      ", ",
      difference$sd,
      "^2), truncated to [",
      truncation[1L],
      ", ",
      truncation[2L],
      "]"
    ),
    side = 3,
    line = 0.20,
    cex = 0.58
  )
  
  invisible(NULL)
}


.rope2arm_plot_summary <- function(x) {
  fit <- .rope2arm_selected_fit(x)
  
  graphics::plot.new()
  graphics::title("Selected design summary")
  
  if (is.null(fit)) {
    graphics::text(
      0.5,
      0.5,
      "No sustained feasible design found",
      cex = 1
    )
    
    return(invisible(NULL))
  }
  
  finite_input <- function(name) {
    value <- x$inputs[[name]]
    
    !is.null(value) &&
      length(value) == 1L &&
      is.finite(value)
  }
  
  bayesian_target_text <- {
    pieces <- character(0)
    
    if (finite_input("targetpower")) {
      pieces <- c(
        pieces,
        paste0(
          "PCE(H1) >= ",
          .format_rope2arm_number(x$inputs$targetpower)
        )
      )
    }
    
    if (finite_input("targettype1")) {
      pieces <- c(
        pieces,
        paste0(
          "False equivalence <= ",
          .format_rope2arm_number(x$inputs$targettype1)
        )
      )
    }
    
    if (finite_input("targetpce_h0")) {
      pieces <- c(
        pieces,
        paste0(
          "PCE(H0) >= ",
          .format_rope2arm_number(x$inputs$targetpce_h0)
        )
      )
    }
    
    if (length(pieces) == 0L) {
      "Not used for feasibility"
    } else {
      paste(pieces, collapse = "; ")
    }
  }
  
  frequentist_target_text <- {
    pieces <- character(0)
    
    if (finite_input("targetfreqpower")) {
      pieces <- c(
        pieces,
        paste0(
          "Power >= ",
          .format_rope2arm_number(x$inputs$targetfreqpower)
        )
      )
    }
    
    if (finite_input("targetfreqtype1")) {
      pieces <- c(
        pieces,
        paste0(
          "Type-I <= ",
          .format_rope2arm_number(x$inputs$targetfreqtype1)
        )
      )
    }
    
    if (length(pieces) == 0L) {
      "Not used for feasibility"
    } else {
      paste(pieces, collapse = "; ")
    }
  }
  
  has_freq_power <- !is.null(fit$freqpower) &&
    length(fit$freqpower) == 1L &&
    is.finite(fit$freqpower)
  
  has_freq_type1 <- !is.null(fit$freqtype1) &&
    length(fit$freqtype1) == 1L &&
    is.finite(fit$freqtype1)
  
  labels <- c(
    "Design",
    "Calibration",
    "ROPE half-width",
    "Thresholds",
    "Bayesian targets",
    "Frequentist targets",
    "Selected nC, nT",
    "Total N",
    "PCE(H1): practical equivalence",
    "Predictive false equivalence",
    "PCE(H0): non-equivalence",
    "P(inconclusive | G1)",
    "P(inconclusive | G0)",
    "Frequentist power",
    "Frequentist type-I error"
  )
  
  values <- c(
    "One-stage two-arm equivalence",
    
    x$inputs$calibration,
    
    .format_rope2arm_number(x$inputs$delta),
    
    paste0(
      "eq = ",
      .format_rope2arm_number(x$inputs$gammaeq),
      "; diff = ",
      .format_rope2arm_number(x$inputs$gammadiff)
    ),
    
    bayesian_target_text,
    frequentist_target_text,
    
    paste(fit$nC, fit$nT, sep = ", "),
    fit$N,
    
    .format_rope2arm_number(fit$power),
    .format_rope2arm_number(fit$type1),
    .format_rope2arm_number(fit$p_diff_h0),
    .format_rope2arm_number(fit$p_inc_h1),
    .format_rope2arm_number(fit$p_inc_h0),
    
    if (has_freq_power) {
      .format_rope2arm_number(fit$freqpower)
    } else {
      "Not evaluated"
    },
    
    if (has_freq_type1) {
      .format_rope2arm_number(fit$freqtype1)
    } else {
      "Not evaluated"
    }
  )
  
  yy <- seq(
    from = 0.95,
    to = 0.05,
    length.out = length(labels)
  )
  
  for (i in seq_along(labels)) {
    graphics::text(
      0.03,
      yy[i],
      labels[i],
      adj = c(0, 0.5),
      cex = 0.66
    )
    
    graphics::text(
      0.53,
      yy[i],
      values[i],
      adj = c(0, 0.5),
      cex = 0.66
    )
  }
  
  graphics::segments(
    0.49,
    0.02,
    0.49,
    0.98,
    col = "grey75"
  )
  
  invisible(x)
}


.rope2arm_plot_prior_panel <- function(
    prior_C,
    prior_T,
    delta,
    title_text,
    region = c("analysis", "H0", "H1"),
    colour_C = "#0072B2",
    colour_T = "#D55E00"
) {
  region <- match.arg(region)
  
  p <- seq(0, 1, length.out = 1000L)
  
  dC <- dbeta(p, prior_C[1L], prior_C[2L])
  dT <- dbeta(p, prior_T[1L], prior_T[2L])
  
  ymax <- max(dC, dT) * 1.08
  
  plot(
    p,
    dC,
    type = "l",
    lwd = 2,
    col = colour_C,
    ylim = c(0, ymax),
    xlab = "Response probability",
    ylab = "Density",
    main = title_text
  )
  
  lines(p, dT, lwd = 2, col = colour_T)
  
  if (region == "H1") {
    polygon(
      x = c(0, 1, 1, 0),
      y = c(0, 0, ymax, ymax),
      col = grDevices::adjustcolor("#BBDDF0", alpha.f = 0.18),
      border = NA
    )
    
    lines(p, dC, lwd = 2, col = colour_C)
    lines(p, dT, lwd = 2, col = colour_T)
    
    text(
      0.5,
      ymax * 0.94,
      "Marginal beta kernels; joint prior truncated to H1",
      cex = 0.64
    )
  }
  
  if (region == "H0") {
    polygon(
      x = c(0, 1, 1, 0),
      y = c(0, 0, ymax, ymax),
      col = grDevices::adjustcolor("#F3C0C0", alpha.f = 0.18),
      border = NA
    )
    
    lines(p, dC, lwd = 2, col = colour_C)
    lines(p, dT, lwd = 2, col = colour_T)
    
    text(
      0.5,
      ymax * 0.94,
      "Marginal beta kernels; joint prior truncated to H0",
      cex = 0.64
    )
  }
  
  legend(
    "topright",
    legend = c(
      paste0("Control: ", .rope2arm_prior_label(prior_C)),
      paste0("Treatment: ", .rope2arm_prior_label(prior_T))
    ),
    col = c(colour_C, colour_T),
    lwd = 2,
    bty = "n",
    cex = 0.68
  )
  
  invisible(NULL)
}


.rope2arm_plot_decision_heatmap <- function(x) {
  fit <- .rope2arm_selected_fit(x)
  
  if (is.null(fit) || is.null(fit$posterior_rope)) {
    plot.new()
    title("Decision heatmap")
    
    text(
      0.5,
      0.5,
      "Set returnmatrices = TRUE\nfor the selected design",
      cex = 0.90
    )
    
    return(invisible(NULL))
  }
  
  z <- fit$posterior_rope
  
  image(
    x = 0:fit$nC,
    y = 0:fit$nT,
    z = z,
    col = hcl.colors(64L, palette = "YlOrRd", rev = TRUE),
    xlab = expression(y[C]),
    ylab = expression(y[T]),
    main = "Posterior ROPE probability",
    zlim = c(0, 1)
  )
  
  contour(
    x = 0:fit$nC,
    y = 0:fit$nT,
    z = z,
    levels = c(
      1 - x$inputs$gammadiff,
      x$inputs$gammaeq
    ),
    add = TRUE,
    drawlabels = TRUE,
    col = c("#D55E00", "#0072B2"),
    lwd = 2
  )
  
  legend(
    "topleft",
    legend = c(
      paste0(
        "Difference threshold: ",
        .format_rope2arm_number(1 - x$inputs$gammadiff)
      ),
      paste0(
        "Equivalence threshold: ",
        .format_rope2arm_number(x$inputs$gammaeq)
      )
    ),
    col = c("#D55E00", "#0072B2"),
    lwd = 2,
    bty = "n",
    cex = 0.72
  )
  
  invisible(NULL)
}


.rope2arm_plot_priors <- function(x) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  graphics::layout(matrix(c(1, 2), nrow = 1L))
  
  graphics::par(mar = c(4.2, 4.2, 3.8, 1.2))
  
  .rope2arm_plot_joint_design_prior(
    prior_C = x$inputs$design_prior_ne_C,
    prior_T = x$inputs$design_prior_ne_T,
    delta = x$inputs$delta,
    region = "H0",
    title_text = "Design prior under H0: non-equivalence"
  )
  
  graphics::par(mar = c(4.2, 4.2, 3.8, 1.2))
  
  if (!is.null(x$inputs$design_prior_eq) &&
      identical(
        x$inputs$design_prior_eq$type,
        "baseline_difference"
      )) {
    .rope2arm_plot_baseline_difference_prior(
      design_prior_eq = x$inputs$design_prior_eq,
      delta = x$inputs$delta,
      main = "Design prior under H1: equivalence"
    )
  } else {
    .rope2arm_plot_joint_design_prior(
      prior_C = x$inputs$design_prior_eq_C,
      prior_T = x$inputs$design_prior_eq_T,
      delta = x$inputs$delta,
      region = "H1",
      title_text = "Design prior under H1: equivalence"
    )
  }
  
  invisible(x)
}


#' Plot a two-arm one-stage ROPE equivalence design
#'
#' @param x A `bfbin2armrope2armdesign` object.
#' @param type Plot type.
#' @param ... Currently unused.
#'
#' @return Invisibly returns `x`.
#' @export
plot.bfbin2armrope2armdesign <- function(
    x,
    type = c("overview", "operatingcharacteristics", "heatmap", "prior"),
    ...
) {
  type <- match.arg(type)
  
  if (is.null(x$grid)) {
    stop("Plotting requires returngrid = TRUE.", call. = FALSE)
  }
  
  if (type == "prior") {
    return(.rope2arm_plot_priors(x))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  if (type == "operatingcharacteristics") {
    .rope2arm_plot_operating_characteristics(x)
    return(invisible(x))
  }
  
  if (type == "heatmap") {
    .rope2arm_plot_decision_heatmap(x)
    return(invisible(x))
  }
  
  graphics::layout(matrix(1:6, nrow = 2L, ncol = 3L, byrow = TRUE))
  
  graphics::par(mar = c(4.0, 4.0, 3.0, 1.0))
  .rope2arm_plot_operating_characteristics(x)
  
  graphics::par(mar = c(4.0, 4.0, 3.0, 1.0))
  .rope2arm_plot_decision_heatmap(x)
  
  graphics::par(mar = c(1.0, 1.0, 3.0, 1.0))
  .rope2arm_plot_summary(x)
  
  graphics::par(mar = c(4.0, 4.0, 3.0, 1.0))
  .rope2arm_plot_prior_panel(
    prior_C = x$inputs$analysis_prior_C,
    prior_T = x$inputs$analysis_prior_T,
    delta = x$inputs$delta,
    title_text = "Analysis priors",
    region = "analysis"
  )
  
  graphics::par(mar = c(4.0, 4.0, 3.0, 1.0))
  .rope2arm_plot_prior_panel(
    prior_C = x$inputs$design_prior_ne_C,
    prior_T = x$inputs$design_prior_ne_T,
    delta = x$inputs$delta,
    title_text = "Design-prior kernels under H0",
    region = "H0"
  )
  
  graphics::par(mar = c(4.0, 4.0, 3.0, 1.0))
  
  if (!is.null(x$inputs$design_prior_eq) &&
      identical(
        x$inputs$design_prior_eq$type,
        "baseline_difference"
      )) {
    
    .rope2arm_plot_baseline_difference_prior_overview(
      design_prior_eq = x$inputs$design_prior_eq,
      delta = x$inputs$delta,
      main = "Design prior under H1: equivalence"
    )
    
  } else {
    
    .rope2arm_plot_prior_panel(
      prior_C = x$inputs$design_prior_eq_C,
      prior_T = x$inputs$design_prior_eq_T,
      delta = x$inputs$delta,
      title_text = "Design-prior kernels under H1",
      region = "H1"
    )
  }
  
  invisible(x)
}