#' Summary for single-arm BF designs
#'
#' @param object An object of class "singlearm_bf_design".
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class "summary.singlearm_bf_design".
#' @export
summary.singlearm_bf_design <- function(object, ...) {
  oc <- object$operating_characteristics
  
  out <- list(
    design = object$design,
    feasible = isTRUE(object$feasible),
    operating_characteristics = oc,
    message = if (!is.null(object$message)) object$message else object$status,
    status = object$status,
    calibration = object$calibration,
    inputs = object$inputs,
    fixed_anchor = object$fixed_anchor
  )
  
  class(out) <- "summary.singlearm_bf_design"
  out
}


#' Print method for summary.singlearm_bf_design
#'
#' @param x An object of class "summary.singlearm_bf_design".
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return The input object \code{x}, invisibly.
#' @export
print.summary.singlearm_bf_design <- function(x, ...) {
  cat("Summary: Single-arm two-stage Bayes factor design\n")
  cat("---------------------------------------------------------\n")
  
  cat("Feasible: ", if (isTRUE(x$feasible)) "TRUE" else "FALSE", "\n", sep = "")
  
  if (!is.null(x$calibration)) {
    cat("Calibration: ", x$calibration, "\n", sep = "")
  }
    
    ## H0: directional -> truncated; point -> point mass at p0
    if (!is.null(x$inputs) && !is.null(x$inputs$prior)) {
      pr <- x$inputs$prior
      
      if (!is.null(pr$da0) && !is.null(pr$db0)) {
        if (!is.null(x$inputs$type) && identical(x$inputs$type, "direction")) {
          cat(sprintf(
            "Design prior under H0: Beta(%s, %s) truncated to [0, p0]\n",
            format(pr$da0), format(pr$db0)
          ))
        } else {
          cat("Design prior under H0: point mass at p0\n")
        }
      }
      
      if (!is.null(pr$da1) && !is.null(pr$db1)) {
        if (!is.null(x$inputs$type) && identical(x$inputs$type, "direction")) {
          cat(sprintf(
            "Design prior under H1: Beta(%s, %s) truncated to (p0, 1]\n",
            format(pr$da1), format(pr$db1)
          ))
        } else {
          cat(sprintf(
            "Design prior under H1: Beta(%s, %s)\n",
            format(pr$da1), format(pr$db1)
          ))
        }
      }
    }
  
  cat("\n")
  
  if (isTRUE(x$feasible)) {
    cat(
      "Selected design: n1 = ",
      x$design["n1"],
      ", n2 = ",
      x$design["n2"],
      "\n",
      sep = ""
    )
  } else {
    if (!is.null(x$status)) {
      cat(x$status, "\n")
    }    
    if (!is.null(x$fixed_anchor) && !is.na(x$fixed_anchor)) {
      cat(
        "The reported fixed-sample anchor from step 1 is n2 = ",
        x$fixed_anchor,
        ".\n",
        sep = ""
      )
    }
  }
  
  oc <- x$operating_characteristics
  if (!is.null(oc) && isTRUE(x$feasible)) {
    cat("\nBayesian operating characteristics\n")
    cat(sprintf("  Power: %.4f\n", oc$power))
    cat(sprintf("  Type-I: %.4f\n", oc$type1))
    cat(sprintf("  CE H0: %s\n", ifelse(is.na(oc$ce_h0), "NA", sprintf("%.4f", oc$ce_h0))))
    cat(sprintf("  EN H0: %.2f\n", oc$en_h0))
    cat(sprintf("  EN H1: %.2f\n", oc$en_h1))
    
    cat("\nFrequentist operating characteristics\n")
    cat(sprintf("  Power: %.4f\n", oc$freq_power))
    cat(sprintf("  Type-I: %.4f\n", oc$freq_type1))
    cat(sprintf("  EN H0: %.2f\n", oc$freq_en_h0))
    cat(sprintf("  EN H1: %.2f\n", oc$freq_en_h1))
  }
  
  invisible(x)
}


#' Plot a single-arm Bayes factor design
#'
#' Produces a diagnostic plot for a fitted single-arm two-stage Bayes factor
#' design. Depending on the available information in the object, the plot shows
#' the interim-search results, selected operating characteristics, and the
#' design and analysis priors under \eqn{H_0} and \eqn{H_1}.
#'
#' @param x An object of class `"singlearm_bf_design"`.
#' @param ... Currently unused.
#'
#' @return
#' The input object `x`, invisibly.
#'
#' @seealso [summary.singlearm_bf_design()], [design_singlearm_bf()],
#'   [optimal_twostage_singlearm_bf()]
#'
#' @export
plot.singlearm_bf_design <- function(x, ...) {
  sr <- x$search_results
  if (is.null(sr) || nrow(sr) == 0L) {
    stop("No search results available to plot.", call. = FALSE)
  }
  calibration <- if (!is.null(x$inputs$calibration)) x$inputs$calibration else "Bayesian"
  
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  inp <- x$inputs
  if (is.null(inp)) inp <- x$args
  if (is.null(inp)) inp <- list()
  
  des <- x$design
  oc  <- x$operating_characteristics
  if (is.null(oc)) oc <- x$operatingcharacteristics
  if (is.null(oc) && !is.null(des)) oc <- des
  
  get_design_value <- function(des, name) {
    if (is.null(des)) return(NA_real_)
    if (is.data.frame(des)) {
      if (name %in% names(des) && nrow(des) >= 1L) return(des[[name]][1])
    }
    if (is.list(des)) {
      if (!is.null(des[[name]])) return(des[[name]][1])
    }
    if (is.atomic(des) && !is.null(names(des)) && name %in% names(des)) {
      return(unname(des[name]))
    }
    NA_real_
  }
  
  get_oc_value <- function(oc, name, fallback = NA_real_) {
    if (is.null(oc)) return(fallback)
    if (is.data.frame(oc)) {
      if (name %in% names(oc) && nrow(oc) >= 1L) return(oc[[name]][1])
    }
    if (is.list(oc)) {
      if (!is.null(oc[[name]])) return(oc[[name]][1])
    }
    if (is.atomic(oc) && !is.null(names(oc)) && name %in% names(oc)) {
      return(unname(oc[name]))
    }
    fallback
  }
  
  get_first_finite <- function(...) {
    vals <- list(...)
    for (v in vals) {
      if (length(v) == 1L && is.finite(v)) return(v)
    }
    NA_real_
  }
  
  n1_opt <- get_design_value(des, "n1")
  n2_opt <- get_design_value(des, "n2")
  k_opt  <- get_design_value(des, "k")
  kf_opt <- get_design_value(des, "k_f")
  if (!is.finite(kf_opt)) kf_opt <- get_design_value(des, "kf")
  
  bayes_power <- get_first_finite(
    get_oc_value(oc, "power", NA_real_),
    get_design_value(des, "power")
  )
  bayes_t1e <- get_first_finite(
    get_oc_value(oc, "type1", NA_real_),
    get_design_value(des, "type1")
  )
  bayes_pce <- get_first_finite(
    get_oc_value(oc, "ce_h0", NA_real_),
    get_oc_value(oc, "ceh0", NA_real_),
    get_design_value(des, "ce_h0"),
    get_design_value(des, "ceh0")
  )
  bayes_en_h0 <- get_first_finite(
    get_oc_value(oc, "ess_h0", NA_real_),
    get_oc_value(oc, "enh0", NA_real_),
    get_oc_value(oc, "en_h0", NA_real_),
    get_design_value(des, "ess_h0"),
    get_design_value(des, "enh0"),
    get_design_value(des, "en_h0")
  )
  bayes_en_h1 <- get_first_finite(
    get_oc_value(oc, "ess_h1", NA_real_),
    get_oc_value(oc, "enh1", NA_real_),
    get_oc_value(oc, "en_h1", NA_real_),
    get_design_value(des, "ess_h1"),
    get_design_value(des, "enh1"),
    get_design_value(des, "en_h1")
  )
  
  freq_power <- get_first_finite(
    get_oc_value(oc, "freq_power", NA_real_),
    get_oc_value(oc, "power_freq", NA_real_),
    get_design_value(des, "freq_power"),
    get_design_value(des, "power_freq")
  )
  freq_t1e <- get_first_finite(
    get_oc_value(oc, "freq_type1", NA_real_),
    get_oc_value(oc, "type1_freq", NA_real_),
    get_oc_value(oc, "freq_t1e", NA_real_),
    get_design_value(des, "freq_type1"),
    get_design_value(des, "type1_freq"),
    get_design_value(des, "freq_t1e")
  )
  freq_en_h0 <- get_first_finite(
    get_oc_value(oc, "freq_ess_h0", NA_real_),
    get_oc_value(oc, "freq_en_h0", NA_real_),
    get_oc_value(oc, "freq_en0", NA_real_),
    get_oc_value(oc, "en_h0_freq", NA_real_),
    get_oc_value(oc, "enh0_freq", NA_real_),
    get_design_value(des, "freq_ess_h0"),
    get_design_value(des, "freq_en_h0"),
    get_design_value(des, "freq_en0"),
    get_design_value(des, "en_h0_freq"),
    get_design_value(des, "enh0_freq")
  )
  freq_en_h1 <- get_first_finite(
    get_oc_value(oc, "freq_ess_h1", NA_real_),
    get_oc_value(oc, "freq_en_h1", NA_real_),
    get_oc_value(oc, "freq_en1", NA_real_),
    get_oc_value(oc, "en_h1_freq", NA_real_),
    get_oc_value(oc, "enh1_freq", NA_real_),
    get_design_value(des, "freq_ess_h1"),
    get_design_value(des, "freq_en_h1"),
    get_design_value(des, "freq_en1"),
    get_design_value(des, "en_h1_freq"),
    get_design_value(des, "enh1_freq")
  )
  
  p_fut_h0_bayes <- get_first_finite(
    get_oc_value(oc, "p_fut_h0", NA_real_),
    get_oc_value(oc, "pfut_h0", NA_real_),
    get_oc_value(oc, "prob_futility_h0", NA_real_),
    get_oc_value(oc, "prob_stop_futility_h0", NA_real_),
    get_design_value(des, "p_fut_h0"),
    get_design_value(des, "pfut_h0")
  )
  if (!is.finite(p_fut_h0_bayes) &&
      is.finite(bayes_en_h0) && is.finite(n1_opt) && is.finite(n2_opt) &&
      (n2_opt - n1_opt) != 0) {
    p_fut_h0_bayes <- (n2_opt - bayes_en_h0) / (n2_opt - n1_opt)
  }
  
  p_fut_h1_bayes <- get_first_finite(
    get_oc_value(oc, "p_fut_h1", NA_real_),
    get_oc_value(oc, "pfut_h1", NA_real_),
    get_oc_value(oc, "prob_futility_h1", NA_real_),
    get_oc_value(oc, "prob_stop_futility_h1", NA_real_),
    get_design_value(des, "p_fut_h1"),
    get_design_value(des, "pfut_h1")
  )
  if (!is.finite(p_fut_h1_bayes) &&
      is.finite(bayes_en_h1) && is.finite(n1_opt) && is.finite(n2_opt) &&
      (n2_opt - n1_opt) != 0) {
    p_fut_h1_bayes <- (n2_opt - bayes_en_h1) / (n2_opt - n1_opt)
  }
  
  p_fut_h0_freq <- get_first_finite(
    get_oc_value(oc, "freq_p_fut_h0", NA_real_),
    get_oc_value(oc, "freq_pfut_h0", NA_real_),
    get_oc_value(oc, "p_fut_h0_freq", NA_real_),
    get_oc_value(oc, "pfut_h0_freq", NA_real_),
    get_oc_value(oc, "prob_stop_futility_h0_freq", NA_real_),
    get_design_value(des, "freq_p_fut_h0"),
    get_design_value(des, "freq_pfut_h0")
  )
  if (!is.finite(p_fut_h0_freq) &&
      is.finite(freq_en_h0) && is.finite(n1_opt) && is.finite(n2_opt) &&
      (n2_opt - n1_opt) != 0) {
    p_fut_h0_freq <- (n2_opt - freq_en_h0) / (n2_opt - n1_opt)
  }
  
  p_fut_h1_freq <- get_first_finite(
    get_oc_value(oc, "freq_p_fut_h1", NA_real_),
    get_oc_value(oc, "freq_pfut_h1", NA_real_),
    get_oc_value(oc, "p_fut_h1_freq", NA_real_),
    get_oc_value(oc, "pfut_h1_freq", NA_real_),
    get_oc_value(oc, "prob_stop_futility_h1_freq", NA_real_),
    get_design_value(des, "freq_p_fut_h1"),
    get_design_value(des, "freq_pfut_h1")
  )
  if (!is.finite(p_fut_h1_freq) &&
      is.finite(freq_en_h1) && is.finite(n1_opt) && is.finite(n2_opt) &&
      (n2_opt - n1_opt) != 0) {
    p_fut_h1_freq <- (n2_opt - freq_en_h1) / (n2_opt - n1_opt)
  }
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(oldpar)
  }, add = TRUE)
  
  graphics::layout(matrix(1:4, nrow = 2, byrow = TRUE), heights = c(1.9, 1.1))
  graphics::par(mar = c(4.2, 4.2, 3.2, 1.2), oma = c(0, 0, 1.6, 0))
  
  ## --------------------------------------------------
  ## upper left: interim sample size plot
  ## --------------------------------------------------
  cols <- c(
    bayes_power = "#D55E00",
    bayes_t1e   = "#0072B2",
    freq_power  = "#009E73",
    pce_h0      = "#CC79A7"
  )
  
  yvals <- c(
    sr$power,
    sr$type1,
    if ("freq_power" %in% names(sr)) sr$freq_power else NA_real_,
    if ("ce_h0" %in% names(sr)) sr$ce_h0 else NA_real_,
    inp$target_power %||% NA_real_,
    inp$target_type1 %||% NA_real_,
    inp$target_ce_h0 %||% NA_real_
  )
  ymax <- suppressWarnings(max(yvals, na.rm = TRUE))
  if (!is.finite(ymax)) ymax <- 1
  ymax <- max(1, ymax)
  
  graphics::plot(
    sr$n1, sr$power,
    type = "b", pch = 19, lty = 1, lwd = 1.3, col = cols["bayes_power"],
    ylim = c(0, ymax),
    xlab = "Interim sample size n1",
    ylab = "Operating characteristic",
    main = "Interim sample size search"
  )
  
  if ("type1" %in% names(sr) && any(!is.na(sr$type1))) {
    graphics::lines(sr$n1, sr$type1, type = "b", pch = 17, lty = 1, lwd = 1.3,
                    col = cols["bayes_t1e"])
  }
  if ("freq_power" %in% names(sr) && any(!is.na(sr$freq_power))) {
    graphics::lines(sr$n1, sr$freq_power, type = "b", pch = 15, lty = 1, lwd = 1.3,
                    col = cols["freq_power"])
  }
  if ("ce_h0" %in% names(sr) && any(!is.na(sr$ce_h0))) {
    graphics::lines(sr$n1, sr$ce_h0, type = "b", pch = 18, lty = 1, lwd = 1.3,
                    col = cols["pce_h0"])
  }
  
  if (is.finite(inp$target_power %||% NA_real_)) {
    graphics::abline(h = inp$target_power, col = cols["bayes_power"], lty = 3, lwd = 1.1)
  }
  if (is.finite(inp$target_type1 %||% NA_real_)) {
    graphics::abline(h = inp$target_type1, col = cols["bayes_t1e"], lty = 2, lwd = 1.1)
  }
  if (is.finite(inp$target_ce_h0 %||% NA_real_) &&
      (inp$target_ce_h0 %||% 0) > 0 &&
      "ce_h0" %in% names(sr) && any(!is.na(sr$ce_h0))) {
    graphics::abline(h = inp$target_ce_h0, col = cols["pce_h0"], lty = 2, lwd = 1.1)
  }
  if (is.finite(n1_opt)) {
    graphics::abline(v = n1_opt, col = "grey40", lwd = 1.3)
  }
  
  leg <- c("Bayesian power", "Bayesian type-I-error")
  leg_col <- c(cols["bayes_power"], cols["bayes_t1e"])
  leg_pch <- c(19, 17)
  if ("freq_power" %in% names(sr) && any(!is.na(sr$freq_power))) {
    leg <- c(leg, "Frequentist power")
    leg_col <- c(leg_col, cols["freq_power"])
    leg_pch <- c(leg_pch, 15)
  }
  if ("ce_h0" %in% names(sr) && any(!is.na(sr$ce_h0))) {
    leg <- c(leg, "PCE(H0)")
    leg_col <- c(leg_col, cols["pce_h0"])
    leg_pch <- c(leg_pch, 18)
  }
  graphics::legend("right", legend = leg, col = leg_col, pch = leg_pch, lty = 1,
                   bty = "n", cex = 0.85)
  
  ## --------------------------------------------------
  ## upper right: text operating characteristics panel
  ## --------------------------------------------------
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::title(main = "Operating characteristics")
  
  y <- 0.99
  line_gap <- 0.050
  section_gap <- 0.03
  
  add_line <- function(txt, cex = 1.05, font = 1, x = 0.02) {
    graphics::text(x, y, labels = txt, adj = c(0, 1), cex = cex, font = font)
    y <<- y - line_gap
  }
  
  add_blank <- function(gap = section_gap) {
    y <<- y - gap
  }
  
  add_line("Optimal design:", font = 2)
  add_line(sprintf("Interim sample size: n1 = %s",
                   if (is.finite(n1_opt)) format(n1_opt, trim = TRUE) else "NA"),
           font = 1)
  add_line(sprintf("Final sample size: n2 = %s",
                   if (is.finite(n2_opt)) format(n2_opt, trim = TRUE) else "NA"),
           font = 1)
  
  if (is.finite(k_opt) || is.finite(kf_opt)) {
    add_line(sprintf("Thresholds: k = %s, kf = %s",
                     if (is.finite(k_opt)) format(round(k_opt, 4), nsmall = 4) else "NA",
                     if (is.finite(kf_opt)) format(round(kf_opt, 4), nsmall = 4) else "NA"))
  }
  
  add_blank()
  
  add_line("Bayesian operating characteristics:", font = 2)
  add_line(sprintf("Bayesian power: %s",
                   if (is.finite(bayes_power)) format(round(bayes_power, 4), nsmall = 4) else "NA"))
  add_line(sprintf("Bayesian type-I-error rate: %s",
                   if (is.finite(bayes_t1e)) format(round(bayes_t1e, 4), nsmall = 4) else "NA"))
  add_line(sprintf("PCE(H0): %s",
                   if (is.finite(bayes_pce)) format(round(bayes_pce, 4), nsmall = 4) else "NA"))
  add_line(sprintf("Expected sample size under H0: %s",
                   if (is.finite(bayes_en_h0)) format(round(bayes_en_h0, 2), nsmall = 2) else "NA"))
  add_line(sprintf("Expected sample size under H1: %s",
                   if (is.finite(bayes_en_h1)) format(round(bayes_en_h1, 2), nsmall = 2) else "NA"))
  add_line(sprintf("Pr(stop early for futility | H0): %s",
                   if (is.finite(p_fut_h0_bayes)) format(round(p_fut_h0_bayes, 4), nsmall = 4) else "NA"))
  add_line(sprintf("Pr(stop early for futility | H1): %s",
                   if (is.finite(p_fut_h1_bayes)) format(round(p_fut_h1_bayes, 4), nsmall = 4) else "NA"))
  
  add_blank()
  
  add_line("Frequentist operating characteristics:", font = 2)
  add_line(sprintf("Frequentist power: %s",
                   if (is.finite(freq_power)) format(round(freq_power, 4), nsmall = 4) else "NA"))
  add_line(sprintf("Frequentist type-I-error rate: %s",
                   if (is.finite(freq_t1e)) format(round(freq_t1e, 4), nsmall = 4) else "NA"))
  add_line(sprintf("Expected sample size under H0: %s",
                   if (is.finite(freq_en_h0)) format(round(freq_en_h0, 2), nsmall = 2) else "NA"))
  add_line(sprintf("Expected sample size under H1: %s",
                   if (is.finite(freq_en_h1)) format(round(freq_en_h1, 2), nsmall = 2) else "NA"))
  add_line(sprintf("Pr(stop early for futility | H0): %s",
                   if (is.finite(p_fut_h0_freq)) format(round(p_fut_h0_freq, 4), nsmall = 4) else "NA"))
  add_line(sprintf("Pr(stop early for futility | H1): %s",
                   if (is.finite(p_fut_h1_freq)) format(round(p_fut_h1_freq, 4), nsmall = 4) else "NA"))
  
  ## --------------------------------------------------
  ## lower left: design priors
  ## --------------------------------------------------
  grid <- seq(0, 1, length.out = 1000)
  p0 <- inp$p0 %||% NA_real_
  type <- inp$type %||% "point"
  
  d_h0 <- rep(NA_real_, length(grid))
  d_h1 <- rep(NA_real_, length(grid))
  
  if (identical(type, "direction")) {
    da0 <- inp$da0 %||% 1
    db0 <- inp$db0 %||% 1
    da1 <- inp$da1 %||% inp$da %||% 1
    db1 <- inp$db1 %||% inp$db %||% 1
    
    if (is.finite(p0) && p0 > 0) {
      idx0 <- grid <= p0
      z0 <- grid[idx0] / p0
      d_h0[idx0] <- stats::dbeta(z0, shape1 = da0, shape2 = db0) / p0
    }
    if (is.finite(p0) && p0 < 1) {
      idx1 <- grid >= p0
      z1 <- (grid[idx1] - p0) / (1 - p0)
      d_h1[idx1] <- stats::dbeta(z1, shape1 = da1, shape2 = db1) / (1 - p0)
    }
    
    ymax_d <- suppressWarnings(max(c(d_h0, d_h1), na.rm = TRUE))
    if (!is.finite(ymax_d)) ymax_d <- 1
    
    graphics::plot(grid, d_h0, type = "l", lwd = 2, lty = 1, col = "black",
                   xlab = "p", ylab = "Density", ylim = c(0, ymax_d),
                   main = "Design priors under H0 / H1")
    graphics::lines(grid, d_h1, lwd = 2, lty = 3, col = "black")
  } else {
    da1 <- inp$da1 %||% inp$da %||% 1
    db1 <- inp$db1 %||% inp$db %||% 1
    d_h1 <- stats::dbeta(grid, shape1 = da1, shape2 = db1)
    
    ymax_d <- suppressWarnings(max(d_h1, na.rm = TRUE))
    if (!is.finite(ymax_d) || ymax_d <= 0) ymax_d <- 1
    ymax_d <- max(1, ymax_d)
    
    graphics::plot(grid, d_h1, type = "l", lwd = 2, lty = 3, col = "black",
                   xlab = "p", ylab = "Density", ylim = c(0, ymax_d),
                   main = "Design priors under H0 / H1")
    
    if (is.finite(p0)) {
      graphics::segments(x0 = p0, y0 = 0, x1 = p0, y1 = 0.92 * ymax_d,
                         lwd = 2, lty = 1, col = "black")
      graphics::points(p0, 0.92 * ymax_d, pch = 16, col = "black")
      graphics::text(p0, 0.97 * ymax_d, labels = expression(delta[p[0]]),
                     pos = 3, cex = 0.9)
    }
  }
  
  if (is.finite(p0)) graphics::abline(v = p0, col = "grey40", lty = 2, lwd = 1.1)
  graphics::legend("topright",
                   legend = c("Design prior under H0", "Design prior under H1"),
                   lty = c(1, 3), lwd = 2, col = "black",
                   bty = "n", cex = 0.9)
  
  ## --------------------------------------------------
  ## lower right: analysis priors
  ## --------------------------------------------------
  a0 <- inp$a0 %||% inp$a_null %||% 1
  b0 <- inp$b0 %||% inp$b_null %||% 1
  a1 <- inp$a1 %||% inp$a %||% 1
  b1 <- inp$b1 %||% inp$b %||% 1
  
  ad_h0 <- rep(NA_real_, length(grid))
  ad_h1 <- rep(NA_real_, length(grid))
  
  if (identical(type, "direction")) {
    if (is.finite(p0) && p0 > 0) {
      idx0 <- grid <= p0
      z0 <- grid[idx0] / p0
      ad_h0[idx0] <- stats::dbeta(z0, shape1 = a0, shape2 = b0) / p0
    }
    if (is.finite(p0) && p0 < 1) {
      idx1 <- grid >= p0
      z1 <- (grid[idx1] - p0) / (1 - p0)
      ad_h1[idx1] <- stats::dbeta(z1, shape1 = a1, shape2 = b1) / (1 - p0)
    }
    
    ymax_a <- suppressWarnings(max(c(ad_h0, ad_h1), na.rm = TRUE))
    if (!is.finite(ymax_a)) ymax_a <- 1
    
    graphics::plot(grid, ad_h0, type = "l", lwd = 2, lty = 1, col = "black",
                   xlab = "p", ylab = "Density", ylim = c(0, ymax_a),
                   main = "Analysis priors under H0 / H1")
    graphics::lines(grid, ad_h1, lwd = 2, lty = 3, col = "black")
  } else {
    ad_h1 <- stats::dbeta(grid, shape1 = a1, shape2 = b1)
    
    ymax_a <- suppressWarnings(max(ad_h1, na.rm = TRUE))
    if (!is.finite(ymax_a) || ymax_a <= 0) ymax_a <- 1
    ymax_a <- max(1, ymax_a)
    
    graphics::plot(grid, ad_h1, type = "l", lwd = 2, lty = 3, col = "black",
                   xlab = "p", ylab = "Density", ylim = c(0, ymax_a),
                   main = "Analysis priors under H0 / H1")
    
    if (is.finite(p0)) {
      graphics::segments(x0 = p0, y0 = 0, x1 = p0, y1 = 0.92 * ymax_a,
                         lwd = 2, lty = 1, col = "black")
      graphics::points(p0, 0.92 * ymax_a, pch = 16, col = "black")
      graphics::text(p0, 0.97 * ymax_a, labels = expression(delta[p[0]]),
                     pos = 3, cex = 0.9)
    }
  }
  
  if (is.finite(p0)) graphics::abline(v = p0, col = "grey40", lty = 2, lwd = 1.1)
  graphics::legend("topright",
                   legend = c("Analysis prior under H0", "Analysis prior under H1"),
                   lty = c(1, 3), lwd = 2, col = "black",
                   bty = "n", cex = 0.9)
  
  graphics::mtext("Optimal single-arm two-stage Bayes factor design",
                  outer = TRUE, cex = 1.05, font = 2)
  invisible(x)
}