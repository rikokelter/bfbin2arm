#' @export
print.twoarm_twostage_bf_design <- function(x, ...) {
  cat("\nOptimal two-stage two-arm Bayes factor design\n")
  cat("------------------------------------\n")
  cat("Mode: ", x$mode, "\n", sep = "")
  cat("Status: ", x$status, "\n", sep = "")
  cat("Calibration: ", x$calibration, "\n", sep = "")
  cat("Convergence flag: ", x$optimizer$conv, "\n", sep = "")
  
  if (isTRUE(x$feasible) && !is.null(x$design) && !all(is.na(x$design))) {
    cat(
      "Design: n1 = (", x$design["n1_1"], ", ", x$design["n1_2"], ")",
      ", n2 = (", x$design["n2_1"], ", ", x$design["n2_2"], ")\n",
      sep = ""
    )
  } else {
    cat("Design: no feasible design found\n")
  }
  
  if (!is.null(x$fixed_design) && !all(is.na(x$fixed_design))) {
    cat(
      "Fixed-sample anchor: n = (", x$fixed_design["n_fixed_1"], ", ",
      x$fixed_design["n_fixed_2"], ")\n",
      sep = ""
    )
  }
  
  oc <- x$operating_characteristics
  if (!is.null(oc)) {
    cat("\nCorrected operating characteristics\n")
    cat(
      " Power = ", formatC(oc$power, digits = 4, format = "f"), "\n",
      " Type-I error = ", formatC(oc$type1, digits = 4, format = "f"), "\n",
      " CE(H0) = ", formatC(oc$ce_h0, digits = 4, format = "f"), "\n",
      " EN (Bayesian) = ", formatC(oc$en_bayes, digits = 2, format = "f"), "\n",
      sep = ""
    )
    
    if (!is.na(oc$freq_type1)) {
      cat(" Freq. Type-I = ", formatC(oc$freq_type1, digits = 4, format = "f"), "\n", sep = "")
    }
    if (!is.na(oc$freq_power)) {
      cat(" Freq. Power = ", formatC(oc$freq_power, digits = 4, format = "f"), "\n", sep = "")
    }
    if (!is.na(oc$en_freq)) {
      cat(" EN (Frequentist) = ", formatC(oc$en_freq, digits = 2, format = "f"), "\n", sep = "")
    }
  }
  
  invisible(x)
}

#' @export
summary.twoarm_twostage_bf_design <- function(object, ...) {
  sr <- object$search_results
  
  out <- list(
    mode = object$mode,
    status = object$status,
    feasible = object$feasible,
    calibration = object$calibration,
    conv = object$optimizer$conv,
    design = object$design,
    fixed_design = object$fixed_design,
    operating_characteristics = object$operating_characteristics,
    fixed_operating_characteristics = object$fixed_operating_characteristics,
    targets = list(
      target_power = object$inputs$target_power,
      target_type1 = object$inputs$target_type1,
      target_ce_h0 = object$inputs$target_ce_h0,
      target_freq_power = object$inputs$target_freq_power,
      target_freq_type1 = object$inputs$target_freq_type1
    ),
    inputs = object$inputs
  )
  
  if (!is.null(sr) && nrow(sr) > 0L) {
    active_feasible <- switch(
      object$calibration,
      Bayesian = sr$bayes_feasible,
      frequentist = sr$freq_feasible,
      hybrid = sr$hybrid_feasible,
      full = sr$bayes_feasible & sr$freq_feasible
    )
    
    out$search_overview <- list(
      n_evaluated = nrow(sr),
      n_bayes_feasible = sum(sr$bayes_feasible, na.rm = TRUE),
      n_freq_feasible = sum(sr$freq_feasible, na.rm = TRUE),
      n_hybrid_feasible = sum(sr$hybrid_feasible, na.rm = TRUE),
      n_active_feasible = sum(active_feasible, na.rm = TRUE),
      first_active_feasible = if (any(active_feasible, na.rm = TRUE)) {
        sr[which(active_feasible)[1L], c("n1_1", "n1_2", "n2_1", "n2_2")]
      } else {
        NULL
      }
    )
  }
  
  class(out) <- "summary.twoarm_twostage_bf_design"
  out
}

#' @export
print.summary.twoarm_twostage_bf_design <- function(x, ...) {
  cat("\nSummary: Optimal two-stage two-arm Bayes factor design\n")
  cat("---------------------------------------------\n")
  cat("Mode: ", x$mode, "\n", sep = "")
  cat("Status: ", x$status, "\n", sep = "")
  cat("Calibration: ", x$calibration, "\n", sep = "")
  cat("Convergence flag: ", x$conv, "\n", sep = "")
  cat("Feasible: ", if (isTRUE(x$feasible)) "yes" else "no", "\n", sep = "")
  
  if (!is.null(x$search_overview)) {
    cat("\nSearch overview\n")
    cat(" candidates evaluated = ", x$search_overview$n_evaluated, "\n", sep = "")
    cat(" Bayesian-feasible = ", x$search_overview$n_bayes_feasible, "\n", sep = "")
    cat(" Frequentist-feasible = ", x$search_overview$n_freq_feasible, "\n", sep = "")
    cat(" Hybrid-feasible = ", x$search_overview$n_hybrid_feasible, "\n", sep = "")
    cat(" Active-feasible = ", x$search_overview$n_active_feasible, "\n", sep = "")
  }
  
  if (!is.null(x$design) && !all(is.na(x$design))) {
    cat(
      "\nSelected design\n",
      " n1 = (", x$design["n1_1"], ", ", x$design["n1_2"], ")",
      ", n2 = (", x$design["n2_1"], ", ", x$design["n2_2"], ")\n",
      sep = ""
    )
  }
  
  invisible(x)
}

#' @export
plot.twoarm_twostage_bf_design <- function(
    x,
    type = c("old"),
    main = "Optimal two-stage two-arm Bayes factor design",
    ...
) {
  type <- match.arg(type)
  
  if (type != "old") {
    stop("Currently only type = 'old' is implemented for two-stage designs.", call. = FALSE)
  }
  
  ## Reconstruct the old engine-style list -------------------------------
  eng <- x$engine_output
  if (is.null(eng)) {
    stop("No stored engine_output available; cannot reproduce old plot.", call. = FALSE)
  }
  
  if (is.null(eng$design) || length(eng$design) != 4L || anyNA(eng$design)) {
    stop("engine_output$design must be a valid four-element design vector.", call. = FALSE)
  }
  if (is.null(eng$priors)) stop("engine_output$priors is missing.")
  if (is.null(eng$occ))    stop("engine_output$occ is missing.")
  
  res <- list(
    design = eng$design,
    priors = eng$priors,
    occ    = eng$occ
  )
  
  ## Inline copy of the original plot_twostage_2arm_bf -------------------
  if (is.null(res$design) || length(res$design) != 4L || anyNA(res$design)) {
    stop("res$design must be a valid four-element design vector.")
  }
  if (is.null(res$priors)) {
    stop("res$priors is missing.")
  }
  if (is.null(res$occ)) {
    stop("res$occ is missing.")
  }
  
  pri <- res$priors
  des <- as.integer(res$design)
  occ <- res$occ
  
  n11 <- des[1]
  n12 <- des[2]
  n21 <- des[3]
  n22 <- des[4]
  
  test <- pri$test
  k    <- pri$k
  k_f  <- pri$k_f
  
  a_0_d <- pri$a_0_d; b_0_d <- pri$b_0_d
  a_0_a <- pri$a_0_a; b_0_a <- pri$b_0_a
  
  a_1_d <- pri$a_1_d; b_1_d <- pri$b_1_d
  a_2_d <- pri$a_2_d; b_2_d <- pri$b_2_d
  a_1_a <- pri$a_1_a; b_1_a <- pri$b_1_a
  a_2_a <- pri$a_2_a; b_2_a <- pri$b_2_a
  
  a_1_d_Hminus <- pri$a_1_d_Hminus; b_1_d_Hminus <- pri$b_1_d_Hminus
  a_2_d_Hminus <- pri$a_2_d_Hminus; b_2_d_Hminus <- pri$b_2_d_Hminus
  a_1_a_Hminus <- pri$a_1_a_Hminus; b_1_a_Hminus <- pri$b_1_a_Hminus
  a_2_a_Hminus <- pri$a_2_a_Hminus; b_2_a_Hminus <- pri$b_2_a_Hminus
  
  get_occ <- function(occ, aliases) {
    for (nm in aliases) {
      if (!is.null(occ[[nm]]) && length(occ[[nm]]) == 1L && is.finite(occ[[nm]])) {
        return(unname(as.numeric(occ[[nm]])))
      }
      if (!is.na(match(nm, names(occ)))) {
        val <- occ[[nm]]
        if (length(val) == 1L && is.finite(val)) {
          return(unname(as.numeric(val)))
        }
      }
    }
    NA_real_
  }
  
  power_val    <- get_occ(occ, c("Power"))
  t1e_val      <- get_occ(occ, c("Type1_Error", "Type1Error"))
  ceh0_val     <- get_occ(occ, c("CE_H0", "CEH0"))
  futility_val <- get_occ(occ, c("futility_prob", "futilityprob"))
  enh0_val     <- get_occ(occ, c("E_H0_N", "EH0N"))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  layout(
    matrix(c(1, 1, 2, 2,
             3, 4, 5, 6), nrow = 2, byrow = TRUE),
    widths = c(1, 1, 1, 1),
    heights = c(1, 1)
  )
  par(oma = c(0, 0, 3, 0))
  
  p_grid <- seq(0, 1, length.out = 1000)
  
  line_cols <- c("royalblue4", "firebrick3")
  line_ltys <- c(1, 2)
  line_lwds <- c(2, 2)
  legend_expr <- c(expression(p[1]), expression(p[2]))
  
  ## 1. Design schematic
  par(mar = c(3.5, 3.5, 3, 1))
  plot(
    NA, NA,
    xlim = c(0, 10), ylim = c(0, 10),
    axes = FALSE, xlab = "", ylab = "",
    main = "Design schematic"
  )
  
  rect(0.7, 6.2, 4.0, 8.7, col = "grey95", border = "grey35", lwd = 1.5)
  text(2.35, 8.05, "Interim", cex = 1.0, font = 2)
  text(2.35, 7.20, bquote(n[1]^1 == .(n11) ~ "," ~ n[1]^2 == .(n12)), cex = 0.95)
  
  rect(6.0, 6.2, 9.3, 8.7, col = "grey95", border = "grey35", lwd = 1.5)
  text(7.65, 8.05, "Final", cex = 1.0, font = 2)
  text(7.65, 7.20, bquote(n[2]^1 == .(n21) ~ "," ~ n[2]^2 == .(n22)), cex = 0.95)
  
  arrows(4.1, 7.45, 5.9, 7.45, length = 0.08, lwd = 1.5)
  text(5.0, 9.05, "Continue if no futility stop", cex = 0.85)
  
  text(
    5.0, 5.2,
    bquote(
      k == .(formatC(k, format = "f", digits = 3)) ~ "," ~
        k[f] == .(formatC(k_f, format = "f", digits = 3))
    ),
    cex = 0.95
  )
  if (!is.na(enh0_val)) {
    text(
      5.0, 4.1,
      bquote(E[H[0]](N) == .(formatC(enh0_val, format = "f", digits = 2))),
      cex = 0.95
    )
  }
  text(5.0, 2.9, paste("Test:", test), cex = 0.9)
  
  ## 2. Operating characteristics
  par(mar = c(4.5, 4, 3, 1))
  oc_vals <- c(power_val, t1e_val, ceh0_val, futility_val)
  oc_labs <- c("Power", "T1E", "CE(H0)", "Pr(Futility)")
  oc_cols <- c("#8FB9E3", "#C97B84", "#5E8C61", "#B39B6B")
  
  ylim_top <- max(c(oc_vals, 0), na.rm = TRUE)
  if (!is.finite(ylim_top) || ylim_top <= 0) ylim_top <- 1
  ylim_top <- max(1, 1.15 * ylim_top)
  
  mids <- barplot(
    height = oc_vals,
    names.arg = oc_labs,
    col = oc_cols,
    border = "grey35",
    ylim = c(0, ylim_top),
    ylab = "Probability",
    main = "Operating characteristics",
    las = 1,
    cex.names = 0.9
  )
  
  valid_idx <- which(is.finite(oc_vals))
  if (length(valid_idx) > 0L) {
    text(
      x = mids[valid_idx],
      y = oc_vals[valid_idx],
      labels = formatC(oc_vals[valid_idx], format = "f", digits = 2),
      pos = 3,
      cex = 0.9
    )
  }
  
  ## 3. Design prior under H0
  par(mar = c(4, 4, 3, 1))
  d0 <- dbeta(p_grid, a_0_d, b_0_d)
  plot(
    p_grid, d0, type = "l", lwd = 2, col = "black",
    xlab = "Success probability", ylab = "Density",
    main = expression("Design prior under " * H[0])
  )
  legend("topright", legend = "common", col = "black", lty = 1, lwd = 2, bty = "n")
  
  ## 4. Design prior under alternative
  if (identical(test, "BF01")) {
    d1 <- dbeta(p_grid, a_1_d, b_1_d)
    d2 <- dbeta(p_grid, a_2_d, b_2_d)
    ylim_d <- range(c(d1, d2), finite = TRUE)
    
    plot(
      p_grid, d1, type = "l",
      lwd = line_lwds[1], lty = line_ltys[1], col = line_cols[1],
      xlab = "Success probability", ylab = "Density",
      ylim = ylim_d,
      main = expression("Design prior under " * H[1])
    )
    lines(p_grid, d2, lwd = line_lwds[2], lty = line_ltys[2], col = line_cols[2])
    legend("topright", legend = legend_expr, col = line_cols,
           lty = line_ltys, lwd = line_lwds, bty = "n")
    
  } else if (identical(test, "BF+0")) {
    d1 <- dbeta(p_grid, a_1_d, b_1_d)
    d2 <- dbeta(p_grid, a_2_d, b_2_d)
    ylim_d <- range(c(d1, d2), finite = TRUE)
    
    plot(
      p_grid, d1, type = "l",
      lwd = line_lwds[1], lty = line_ltys[1], col = line_cols[1],
      xlab = "Success probability", ylab = "Density",
      ylim = ylim_d,
      main = expression("Design prior under " * H["+"])
    )
    lines(p_grid, d2, lwd = line_lwds[2], lty = line_ltys[2], col = line_cols[2])
    legend("topright", legend = legend_expr, col = line_cols,
           lty = line_ltys, lwd = line_lwds, bty = "n")
    
  } else if (identical(test, "BF-0")) {
    d1 <- dbeta(p_grid, a_1_d_Hminus, b_1_d_Hminus)
    d2 <- dbeta(p_grid, a_2_d_Hminus, b_2_d_Hminus)
    ylim_d <- range(c(d1, d2), finite = TRUE)
    
    plot(
      p_grid, d1, type = "l",
      lwd = line_lwds[1], lty = line_ltys[1], col = line_cols[1],
      xlab = "Success probability", ylab = "Density",
      ylim = ylim_d,
      main = expression("Design prior under " * H[0]^"-")
    )
    lines(p_grid, d2, lwd = line_lwds[2], lty = line_ltys[2], col = line_cols[2])
    legend("topright", legend = legend_expr, col = line_cols,
           lty = line_ltys, lwd = line_lwds, bty = "n")
    
  } else if (identical(test, "BF+-")) {
    d1  <- dbeta(p_grid, a_1_d,        b_1_d)
    d2  <- dbeta(p_grid, a_2_d,        b_2_d)
    d1m <- dbeta(p_grid, a_1_d_Hminus, b_1_d_Hminus)
    d2m <- dbeta(p_grid, a_2_d_Hminus, b_2_d_Hminus)
    ylim_d <- range(c(d1, d2, d1m, d2m), finite = TRUE)
    
    plot(
      p_grid, d1, type = "l",
      lwd = 2, lty = 1, col = line_cols[1],
      xlab = "Success probability", ylab = "Density",
      ylim = ylim_d,
      main = expression("Design priors under " * H["+"] * " and " * H[0]^"-")
    )
    lines(p_grid, d2,  lwd = 2,   lty = 2, col = line_cols[2])
    lines(p_grid, d1m, lwd = 1.5, lty = 3, col = line_cols[1])
    lines(p_grid, d2m, lwd = 1.5, lty = 4, col = line_cols[2])
    legend(
      "topright",
      legend = c(
        expression(p[1] * " (" * H["+"] * ")"),
        expression(p[2] * " (" * H["+"] * ")"),
        expression(p[1] * " (" * H[0]^"-" * ")"),
        expression(p[2] * " (" * H[0]^"-" * ")")
      ),
      col = c(line_cols[1], line_cols[2], line_cols[1], line_cols[2]),
      lty = c(1, 2, 3, 4),
      lwd = c(2, 2, 1.5, 1.5),
      bty = "n"
    )
  } else {
    plot.new()
    title(main = "Design prior")
    text(0.5, 0.5, "Unsupported test setting.")
  }
  
  ## 5. Analysis prior under H0
  par(mar = c(4, 4, 3, 1))
  a0 <- dbeta(p_grid, a_0_a, b_0_a)
  plot(
    p_grid, a0, type = "l", lwd = 2, col = "black",
    xlab = "Success probability", ylab = "Density",
    main = expression("Analysis prior under " * H[0])
  )
  legend("topright", legend = "common", col = "black", lty = 1, lwd = 2, bty = "n")
  
  ## 6. Analysis prior under alternative
  if (identical(test, "BF01")) {
    a1 <- dbeta(p_grid, a_1_a, b_1_a)
    a2 <- dbeta(p_grid, a_2_a, b_2_a)
    ylim_a <- range(c(a1, a2), finite = TRUE)
    
    plot(
      p_grid, a1, type = "l",
      lwd = line_lwds[1], lty = line_ltys[1], col = line_cols[1],
      xlab = "Success probability", ylab = "Density",
      ylim = ylim_a,
      main = expression("Analysis prior under " * H[1])
    )
    lines(p_grid, a2, lwd = line_lwds[2], lty = line_ltys[2], col = line_cols[2])
    legend("topright", legend = legend_expr, col = line_cols,
           lty = line_ltys, lwd = line_lwds, bty = "n")
    
  } else if (identical(test, "BF+0")) {
    a1 <- dbeta(p_grid, a_1_a, b_1_a)
    a2 <- dbeta(p_grid, a_2_a, b_2_a)
    ylim_a <- range(c(a1, a2), finite = TRUE)
    
    plot(
      p_grid, a1, type = "l",
      lwd = line_lwds[1], lty = line_ltys[1], col = line_cols[1],
      xlab = "Success probability", ylab = "Density",
      ylim = ylim_a,
      main = expression("Analysis prior under " * H["+"])
    )
    lines(p_grid, a2, lwd = line_lwds[2], lty = line_ltys[2], col = line_cols[2])
    legend("topright", legend = legend_expr, col = line_cols,
           lty = line_ltys, lwd = line_lwds, bty = "n")
    
  } else if (identical(test, "BF-0")) {
    a1 <- dbeta(p_grid, a_1_a_Hminus, b_1_a_Hminus)
    a2 <- dbeta(p_grid, a_2_a_Hminus, b_2_a_Hminus)
    ylim_a <- range(c(a1, a2), finite = TRUE)
    
    plot(
      p_grid, a1, type = "l",
      lwd = line_lwds[1], lty = line_ltys[1], col = line_cols[1],
      xlab = "Success probability", ylab = "Density",
      ylim = ylim_a,
      main = expression("Analysis prior under " * H[0]^"-")
    )
    lines(p_grid, a2, lwd = line_lwds[2], lty = line_ltys[2], col = line_cols[2])
    legend("topright", legend = legend_expr, col = line_cols,
           lty = line_ltys, lwd = line_lwds, bty = "n")
    
  } else if (identical(test, "BF+-")) {
    a1p <- dbeta(p_grid, a_1_a, b_1_a)
    a2p <- dbeta(p_grid, a_2_a, b_2_a)
    a1m <- dbeta(p_grid, a_1_a_Hminus, b_1_a_Hminus)
    a2m <- dbeta(p_grid, a_2_a_Hminus, b_2_a_Hminus)
    ylim_a <- range(c(a1p, a2p, a1m, a2m), finite = TRUE)
    
    plot(
      p_grid, a1p, type = "l",
      lwd = 2, lty = 1, col = line_cols[1],
      xlab = "Success probability", ylab = "Density",
      ylim = ylim_a,
      main = expression("Analysis priors under " * H["+"] * " and " * H[0]^"-")
    )
    lines(p_grid, a2p, lwd = 2, lty = 2, col = line_cols[2])
    lines(p_grid, a1m, lwd = 1.5, lty = 3, col = line_cols[1])
    lines(p_grid, a2m, lwd = 1.5, lty = 4, col = line_cols[2])
    legend(
      "topright",
      legend = c(
        expression(p[1] * " (" * H["+"] * ")"),
        expression(p[2] * " (" * H["+"] * ")"),
        expression(p[1] * " (" * H[0]^"-" * ")"),
        expression(p[2] * " (" * H[0]^"-" * ")")
      ),
      col = c(line_cols[1], line_cols[2], line_cols[1], line_cols[2]),
      lty = c(1, 2, 3, 4),
      lwd = c(2, 2, 1.5, 1.5),
      bty = "n"
    )
  } else {
    plot.new()
    title(main = "Analysis prior")
    text(0.5, 0.5, "Unsupported test setting.")
  }
  
  mtext(main, outer = TRUE, cex = 1.2)
  invisible(NULL)
}