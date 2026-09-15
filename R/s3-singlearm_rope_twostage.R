# =============================================================================
# S3 methods for "singlearm_rope_twostage_design" objects
# =============================================================================


# -----------------------------------------------------------------------------
# Direction-specific labels for user-facing output
# -----------------------------------------------------------------------------
.rope_twostage_direction_labels <- function(
    direction = c("equivalence", "noninferiority", "superiority")
) {
  direction <- match.arg(direction)
  
  switch(
    direction,
    
    equivalence = list(
      title = "equivalence",
      h1 = "equivalence",
      h0 = "non-equivalence",
      interim_h0 = "non-equivalence",
      success = "Declare equivalence",
      h0_evidence = "Compelling evidence for non-equivalence",
      reject = "No equivalence",
      delta = "ROPE half-width delta"
    ),
    
    noninferiority = list(
      title = "non-inferiority",
      h1 = "non-inferiority",
      h0 = "clinically relevant inferiority",
      interim_h0 = "clinically relevant inferiority",
      success = "Declare non-inferiority",
      h0_evidence = "Compelling evidence for inferiority",
      reject = "No non-inferiority",
      delta = "NI margin delta"
    ),
    
    superiority = list(
      title = "superiority",
      h1 = "superiority",
      h0 = "non-superiority",
      interim_h0 = "non-superiority",
      success = "Declare superiority",
      h0_evidence = "Compelling evidence for non-superiority",
      reject = "No superiority",
      delta = "Superiority margin delta"
    )
  )
}


# -----------------------------------------------------------------------------
#' Print a single-arm two-stage ROPE design
#'
#' @param x An object of class \code{"singlearm_rope_twostage_design"}.
#' @param ... Further arguments (currently ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.singlearm_rope_twostage_design <- function(x, ...) {
  d <- x$design
  direction <- x$direction
  labels <- .rope_twostage_direction_labels(direction)
  
  title_label <- labels$title
  
  rope_label <- switch(
    direction,
    equivalence = sprintf(
      "ROPE interval : [%.4f, %.4f]",
      max(0, x$p0 - x$delta), min(1, x$p0 + x$delta)
    ),
    noninferiority = sprintf(
      "NI boundary (p >=) : %.4f [%.4f, 1)",
      max(0, x$p0 - x$delta), max(0, x$p0 - x$delta)
    ),
    superiority = sprintf(
      "Sup. boundary (p >) : %.4f (%.4f, 1]",
      min(1, x$p0 + x$delta), min(1, x$p0 + x$delta)
    )
  )
  
  delta_label <- labels$delta
  
  cont_str <- if (length(x$continuation_region) <= 20L) {
    paste(x$continuation_region, collapse = ", ")
  } else {
    paste0(
      paste(head(x$continuation_region, 10L), collapse = ", "),
      " ... ",
      paste(tail(x$continuation_region, 5L), collapse = ", ")
    )
  }
  
  crit <- switch(
    direction,
    equivalence = sprintf(
      "Pr(p in [%.4f, %.4f] | Y) >= %.2f",
      max(0, x$p0 - x$delta),
      min(1, x$p0 + x$delta),
      x$gamma_eq
    ),
    noninferiority = sprintf(
      "Pr(p >= %.4f | Y) >= %.2f",
      max(0, x$p0 - x$delta), x$gamma_eq
    ),
    superiority = sprintf(
      "Pr(p > %.4f | Y) >= %.2f",
      min(1, x$p0 + x$delta), x$gamma_eq
    )
  )
  
  cat(sprintf("\nSingle-arm two-stage ROPE %s design\n", title_label))
  cat(strrep("=", nchar(sprintf(
    "Single-arm two-stage ROPE %s design", title_label
  ))), "\n")
  cat(sprintf(" Benchmark p0 : %.4f\n", x$p0))
  cat(sprintf(" %s : %.4f\n", delta_label, x$delta))
  cat(sprintf(" %s\n", rope_label))
  cat(sprintf(" Analysis prior : Beta(%g, %g)\n",
              x$analysis_prior[1], x$analysis_prior[2]))
  cat(sprintf(" Null design prior : Beta(%g, %g)\n",
              x$design_prior_h0[1], x$design_prior_h0[2]))
  cat(sprintf(" Alt. design prior : Beta(%g, %g)\n",
              x$design_prior_h1[1], x$design_prior_h1[2]))
  cat(sprintf(" Interim threshold gamma_1 : %.4f\n", x$gamma_1))
  cat(sprintf(" Final threshold gamma_eq : %.4f\n", x$gamma_eq))
  cat(sprintf(" Diff. threshold gamma_diff : %.4f\n", x$gamma_diff))
  cat(sprintf(" Target predictive alpha : %.4f\n", x$alpha))
  cat(sprintf(" Target predictive power : %.4f\n", x$target_power))
  if (!is.null(x$target_pce)) {
    cat(sprintf(" Target predictive PCE(H0) : %.4f\n", x$target_pce))
  }
  if (!is.null(x$alpha_freq)) {
    cat(sprintf(" Target frequentist type-I error : %.4f at p = %.4f\n",
                x$alpha_freq, x$p_t1e))
  }
  if (!is.null(x$target_power_freq)) {
    cat(sprintf(" Target frequentist power : %.4f at p = %.4f\n",
                x$target_power_freq, x$p_power))
  }
  cat(sprintf(" Optimality criterion : %s\n", x$optimality))
  cat(sprintf(" Direction : %s\n\n", title_label))
  
  cat("Optimal design\n")
  cat(strrep("-", 30L), "\n")
  cat(sprintf(" n1 (stage 1) : %d\n", as.integer(d$n1)))
  cat(sprintf(" n2 (stage 2) : %d\n", as.integer(d$n2)))
  cat(sprintf(" n (maximum) : %d\n", as.integer(d$n)))
  cat(sprintf(" Continuation region C1 : {%s}\n\n", cont_str))
  
  cat("Operating characteristics\n")
  cat(strrep("-", 47L), "\n")
  cat(sprintf("%-40s %7s %10s\n", "", "1-stage", "2-stage"))
  cat(sprintf(
    " %-38s %7.4f %10.4f\n",
    sprintf("Predictive type-I error (%s)", labels$h0),
    d$type1_1st,
    d$type1_2st
  ))
  
  cat(sprintf(
    " %-38s %7.4f %10.4f\n",
    sprintf("Predictive power (%s)", labels$h1),
    d$power_1st,
    d$power_2st
  ))
  
  cat(sprintf(
    " %-38s %7.4f %10.4f\n",
    sprintf("PCE(%s), any stage", labels$h0),
    d$pce_1st,
    d$pce_2st
  ))
  cat("\n")
  cat(sprintf(" EN under H0 prior : %.4f\n", d$EN0))
  cat(sprintf(" EN under H1 prior : %.4f\n", d$EN1))
  
  if ("freq_type1_2st" %in% names(d) && !is.na(d$freq_type1_2st)) {
    cat(sprintf(" Frequentist type-I error : %.4f", d$freq_type1_2st))
    if (!is.null(x$alpha_freq)) {
      cat(sprintf("  (target <= %.4f at p = %.4f)", x$alpha_freq, x$p_t1e))
    }
    cat("\n")
  }
  
  if ("freq_power_2st" %in% names(d) && !is.na(d$freq_power_2st)) {
    cat(sprintf(" Frequentist power : %.4f", d$freq_power_2st))
    if (!is.null(x$target_power_freq)) {
      cat(sprintf("  (target >= %.4f at p = %.4f)",
                  x$target_power_freq, x$p_power))
    }
    cat("\n")
  }
  
  cat("\n")
  cat(sprintf(" %s : %s\n\n", labels$success, crit))
  
  invisible(x)
}

# -----------------------------------------------------------------------------
#' Summary for a single-arm two-stage ROPE design
#'
#' @param object An object of class \code{"singlearm_rope_twostage_design"}.
#' @param ... Further arguments (currently ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.singlearm_rope_twostage_design <- function(object, ...) {
  direction <- object$direction
  labels <- .rope_twostage_direction_labels(direction)
  
  title_label <- labels$title
  h0_label <- labels$h0
  h1_label <- labels$h1
  
  cat(sprintf("\nSummary: Single-arm two-stage ROPE %s design\n", title_label))
  cat(strrep("=", 55L), "\n\n")
  
  cat("Hypotheses\n")
  cat(strrep("-", 30L), "\n")
  switch(
    direction,
    equivalence = {
      lo <- max(0, object$p0 - object$delta)
      hi <- min(1, object$p0 + object$delta)
      cat(sprintf(" H0 : p outside [%.4f, %.4f] (%s)\n", lo, hi, h0_label))
      cat(sprintf(" H1 : p inside  [%.4f, %.4f] (%s)\n\n", lo, hi, h1_label))
    },
    noninferiority = {
      bnd <- max(0, object$p0 - object$delta)
      cat(sprintf(" H0 : p < %.4f (%s)\n", bnd, h0_label))
      cat(sprintf(" H1 : p >= %.4f (%s)\n\n", bnd, h1_label))
    },
    superiority = {
      bnd <- min(1, object$p0 + object$delta)
      cat(sprintf(" H0 : p <= %.4f (%s)\n", bnd, h0_label))
      cat(sprintf(" H1 : p > %.4f (%s)\n\n", bnd, h1_label))
    }
  )
  
  cat("Decision rules\n")
  cat(strrep("-", 30L), "\n")
  
  crit_interim <- switch(
    direction,
    
    equivalence = sprintf(
      "Continue if Pr(p in ROPE | y1, n1) > %.2f",
      1 - object$gamma_1
    ),
    
    noninferiority = sprintf(
      "Continue if Pr(p >= %.4f | y1, n1) > %.2f",
      object$p0 - object$delta,
      1 - object$gamma_1
    ),
    
    superiority = sprintf(
      "Continue if Pr(p > %.4f | y1, n1) > %.2f",
      object$p0 + object$delta,
      1 - object$gamma_1
    )
  )
  
  crit_final <- switch(
    direction,
    equivalence = sprintf(
      "Declare equivalence if Pr(p in ROPE | y, n) >= %.2f",
      object$gamma_eq
    ),
    noninferiority = sprintf(
      "Declare non-inferiority if Pr(p >= %.4f | y, n) >= %.2f",
      max(0, object$p0 - object$delta), object$gamma_eq
    ),
    superiority = sprintf(
      "Declare superiority if Pr(p > %.4f | y, n) >= %.2f",
      min(1, object$p0 + object$delta), object$gamma_eq
    )
  )
  
  crit_h0 <- switch(
    direction,
    equivalence = sprintf(
      "Compelling evidence for H0 if 1 - Pr(p in ROPE | y_j, n_j) >= %.2f",
      object$gamma_diff
    ),
    noninferiority = sprintf(
      "Compelling evidence for H0 if 1 - Pr(p >= %.4f | y_j, n_j) >= %.2f",
      max(0, object$p0 - object$delta), object$gamma_diff
    ),
    superiority = sprintf(
      "Compelling evidence for H0 if 1 - Pr(p > %.4f | y_j, n_j) >= %.2f",
      min(1, object$p0 + object$delta), object$gamma_diff
    )
  )
  
  cat(sprintf(" Interim : %s\n", crit_interim))
  cat(sprintf(" Final   : %s\n", crit_final))
  cat(sprintf(" %s : %s\n\n", labels$h0_evidence, crit_h0))
  
  cat("Calibration targets\n")
  cat(strrep("-", 30L), "\n")
  cat(sprintf(" Predictive type-I error <= %.4f\n", object$alpha))
  cat(sprintf(" Predictive power >= %.4f\n", object$target_power))
  if (!is.null(object$target_pce)) {
    cat(sprintf(" Predictive PCE(H0) >= %.4f\n", object$target_pce))
  }
  if (!is.null(object$alpha_freq)) {
    cat(sprintf(" Frequentist type-I error <= %.4f at p = %.4f\n",
                object$alpha_freq, object$p_t1e))
  }
  if (!is.null(object$target_power_freq)) {
    cat(sprintf(" Frequentist power >= %.4f at p = %.4f\n",
                object$target_power_freq, object$p_power))
  }
  cat("\n")
  
  print(object, ...)
  invisible(object)
}

# -----------------------------------------------------------------------------
#' Print a compact summary of a single-arm two-stage ROPE design
#'
#' @param x An object of class \code{"summary.singlearm_rope_twostage_design"}.
#' @param digits Number of digits used to display operating characteristics.
#' @param ... Further arguments (currently ignored).
#' @export
print.summary.singlearm_rope_twostage_design <- function(x, digits = 4L, ...) {
  d <- x$design
  
  cat("\n")
  cat("Summary: single-arm two-stage ROPE design\n")
  cat("-----------------------------------------\n")
  cat(sprintf(" Benchmark p0 : %.4f\n", x$p0))
  cat(sprintf(" Delta : %.4f\n", x$delta))
  cat(sprintf(" Target predictive alpha : %.4f\n", x$alpha))
  cat(sprintf(" Target predictive power : %.4f\n", x$target_power))
  if (!is.null(x$target_pce)) {
    cat(sprintf(" Target predictive PCE(H0) : %.4f\n", x$target_pce))
  }
  if (!is.null(x$alpha_freq)) {
    cat(sprintf(" Target frequentist type-I error : %.4f at p = %.4f\n",
                x$alpha_freq, x$p_t1e))
  }
  if (!is.null(x$target_power_freq)) {
    cat(sprintf(" Target frequentist power : %.4f at p = %.4f\n",
                x$target_power_freq, x$p_power))
  }
  cat(sprintf(" Feasible designs : %d\n", nrow(x$candidates)))
  cat("\n")
  
  out <- data.frame(
    n1 = as.integer(d$n1),
    n2 = as.integer(d$n2),
    n = as.integer(d$n),
    type1_2st = round(d$type1_2st, digits),
    power_2st = round(d$power_2st, digits),
    pce_2st = round(d$pce_2st, digits),
    EN0 = round(d$EN0, digits),
    EN1 = round(d$EN1, digits)
  )
  
  if ("freq_type1_2st" %in% names(d)) {
    out$freq_type1_2st <- round(d$freq_type1_2st, digits)
  }
  if ("freq_power_2st" %in% names(d)) {
    out$freq_power_2st <- round(d$freq_power_2st, digits)
  }
  
  print(out, row.names = FALSE)
  invisible(x)
}

# -----------------------------------------------------------------------------
#' Extract design coefficients
#'
#' @param object An object of class \code{"singlearm_rope_twostage_design"}.
#' @param ... Further arguments (currently ignored).
#' @return A named numeric vector of key design quantities.
#' @export
coef.singlearm_rope_twostage_design <- function(object, ...) {
  d <- object$design
  c(n1        = d$n1,
    n2        = d$n2,
    n         = d$n,
    type1_2st = d$type1_2st,
    power_2st = d$power_2st,
    pce_2st   = d$pce_2st,
    EN0       = d$EN0,
    EN1       = d$EN1)
}

# -----------------------------------------------------------------------------
#' Convert candidates to a data frame
#'
#' @param x An object of class \code{"singlearm_rope_twostage_design"}.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Further arguments (currently ignored).
#' @return The candidates data frame.
#' @export
as.data.frame.singlearm_rope_twostage_design <- function(
    x, row.names = NULL, optional = FALSE, ...) {
  x$candidates
}

#' Plot a single-arm two-stage ROPE design
#'
#' Produces a six-panel figure in a 2-row by 3-column layout:
#' \describe{
#' \item{Top-left}{Predictive type-I error and power as functions of the
#' stage-1 sample size \eqn{n_1} for the fixed maximum \eqn{n^*}, with
#' the optimal \eqn{n_1} marked. Optional dashed frequentist curves are added
#' when corresponding calibration targets are present.}
#' \item{Top-centre}{Predictive PCE under \eqn{H_0} as a function of
#' \eqn{n_1} for the fixed maximum \eqn{n^*}, with the optimal \eqn{n_1}
#' marked.}
#' \item{Top-right}{Compact textual summary of the optimal design and its
#' operating characteristics.}
#' \item{Bottom-left}{Null design prior density with the
#' direction-appropriate ROPE region shaded.}
#' \item{Bottom-centre}{Alternative design prior density with the
#' direction-appropriate ROPE region shaded.}
#' \item{Bottom-right}{Analysis prior density with the
#' direction-appropriate ROPE region shaded.}
#' }
#'
#' Additional plot types are available via the \code{type} argument:
#' \describe{
#' \item{\code{"default"}}{The 2x3 summary layout described above (default).}
#' \item{\code{"interim"}}{Interim posterior probability vs. stage-1
#' responses, showing the continuation region \eqn{C_1}.}
#' \item{\code{"final"}}{Final posterior probability vs. total
#' responses, showing the acceptance region.}
#' }
#'
#' @param x An object of class \code{"singlearm_rope_twostage_design"}.
#' @param type Character string selecting the plot type. One of
#' \code{"default"} (default), \code{"interim"}, or \code{"final"}.
#' @param ... Further graphical parameters passed to \code{par()} for the
#' \code{"interim"} and \code{"final"} types.
#' @return Invisibly returns \code{x}.
#' @export
plot.singlearm_rope_twostage_design <- function(
    x,
    type = c("default", "interim", "final"),
    ...
) {
  type <- match.arg(type)
  
  labels <- .rope_twostage_direction_labels(x$direction)
  
  .bbpmf_loc <- function(y, n, a, b)
    exp(lchoose(n, y) + lbeta(a + y, b + n - y) - lbeta(a, b))
  
  .post_prob <- function(y, n) {
    posterior_rope_prob(
      y = y,
      n = n,
      p0 = x$p0,
      delta = x$delta,
      analysis_prior = x$analysis_prior,
      direction = x$direction
    )
  }
  
  .cont_region <- function(n1) {
    y1 <- 0:n1
    y1[.post_prob(y1, n1) > (1 - x$gamma_1)]
  }
  
  d <- x$design
  n_opt <- as.integer(d$n)
  n1_opt <- as.integer(d$n1)
  n2_opt <- as.integer(d$n2)
  cont <- x$continuation_region
  aA <- x$analysis_prior[1]; bA <- x$analysis_prior[2]
  aH0 <- x$design_prior_h0[1]; bH0 <- x$design_prior_h0[2]
  aH1 <- x$design_prior_h1[1]; bH1 <- x$design_prior_h1[2]
  lo_rope <- max(0, x$p0 - x$delta)
  hi_rope <- min(1, x$p0 + x$delta)
  p_seq <- seq(0, 1, length.out = 500)
  
  direction_label <- tools::toTitleCase(labels$title)
  h1_label <- labels$h1
  h0_label <- labels$h0
  
  crit_str <- switch(
    x$direction,
    equivalence = sprintf("Pr(p in [%.3f, %.3f] | Y) >= %.2f",
                          lo_rope, hi_rope, x$gamma_eq),
    noninferiority = sprintf("Pr(p >= %.4f | Y) >= %.2f",
                             lo_rope, x$gamma_eq),
    superiority = sprintf("Pr(p > %.4f | Y) >= %.2f",
                          hi_rope, x$gamma_eq)
  )
  
  ylab_post <- switch(
    x$direction,
    equivalence = expression(Pr(p %in% R[p] ~ "|" ~ y, n)),
    noninferiority = expression(Pr(p >= p[0] - delta ~ "|" ~ y, n)),
    superiority = expression(Pr(p > p[0] + delta ~ "|" ~ y, n))
  )
  
  ylab_post1 <- switch(
    x$direction,
    equivalence = expression(Pr(p %in% R[p] ~ "|" ~ y[1], n[1])),
    noninferiority = expression(Pr(p >= p[0] - delta ~ "|" ~ y[1], n[1])),
    superiority = expression(Pr(p > p[0] + delta ~ "|" ~ y[1], n[1]))
  )
  
  accept_label <- labels$success
  reject_label <- labels$reject
  
  ## palette
  col_cont <- "#2E86AB"
  col_stop <- "#E84855"
  col_bnd  <- "#F18F01"
  col_grid <- "#DDDDDD"
  col_h0   <- "#E84855"
  col_h1   <- "#2E86AB"
  col_ap   <- "#6A0572"
  col_rope <- grDevices::adjustcolor("#F18F01", alpha.f = 0.15)
  col_opt  <- "black"
  
  ## Bayesian vs frequentist colours (colour-blind friendlier)
  col_bayes_power <- "#1F77B4"   # blue
  col_bayes_t1e   <- "#D55E00"   # vermillion
  col_freq_power  <- "#6A3D9A"   # purple
  col_freq_t1e    <- "#E69F00"   # orange
  col_pce         <- "#009E73"   # bluish green
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  .prior_panel <- function(a, b, col_line, title) {
    dens <- dbeta(p_seq, a, b)
    plot(p_seq, dens, type = "l", col = col_line, lwd = 2,
         xlab = "p", ylab = "Density", main = title,
         las = 1L, bty = "l",
         ylim = c(0, max(dens) * 1.08))
    
    idx <- switch(
      x$direction,
      equivalence = p_seq >= lo_rope & p_seq <= hi_rope,
      noninferiority = p_seq >= lo_rope,
      superiority = p_seq >= hi_rope
    )
    
    polygon(c(p_seq[idx], rev(p_seq[idx])),
            c(dens[idx], rep(0, sum(idx))),
            col = col_rope, border = NA)
    
    if (x$direction == "equivalence") {
      abline(v = c(lo_rope, hi_rope), col = col_bnd, lty = 2L, lwd = 1.2)
    } else if (x$direction == "noninferiority") {
      abline(v = lo_rope, col = col_bnd, lty = 2L, lwd = 1.2)
    } else {
      abline(v = hi_rope, col = col_bnd, lty = 2L, lwd = 1.2)
    }
    
    abline(v = x$p0, col = "grey40", lty = 3L, lwd = 1.0)
  }
  
  if (type == "default") {
    
    n1_vals <- seq_len(n_opt - 1L)
    message(sprintf(
      "Computing OC curves over %d interim positions (n* = %d) ...",
      length(n1_vals), n_opt))
    
    t1e_vec <- pwr_vec <- pce_vec <- numeric(length(n1_vals))
    freq_t1e_vec <- freq_pwr_vec <- rep(NA_real_, length(n1_vals))
    
    has_freq_t1e <- !is.null(x$p_t1e)
    has_freq_pwr <- !is.null(x$p_power)
    
    for (i in seq_along(n1_vals)) {
      n1i <- n1_vals[i]
      n2i <- n_opt - n1i
      
      oc <- .evaluate_rope_twostage_oc(
        n1 = n1i,
        n2 = n2i,
        p0 = x$p0,
        delta = x$delta,
        gamma_1 = x$gamma_1,
        gamma_eq = x$gamma_eq,
        gamma_diff = x$gamma_diff,
        analysis_prior = x$analysis_prior,
        design_prior_h0 = x$design_prior_h0,
        design_prior_h1 = x$design_prior_h1,
        p_t1e = if (has_freq_t1e) x$p_t1e else NULL,
        p_power = if (has_freq_pwr) x$p_power else NULL,
        direction = x$direction
      )
      
      t1e_vec[i] <- oc$type1_2st
      pwr_vec[i] <- oc$power_2st
      pce_vec[i] <- oc$pce_2st
      
      if (has_freq_t1e && !is.null(oc$freq_type1_2st)) {
        freq_t1e_vec[i] <- oc$freq_type1_2st
      }
      if (has_freq_pwr && !is.null(oc$freq_power_2st)) {
        freq_pwr_vec[i] <- oc$freq_power_2st
      }
    }
    
    layout(matrix(c(1, 2, 3,
                    4, 5, 6),
                  nrow = 2L, byrow = TRUE),
           heights = c(0.54, 0.46))
    
    ## Panel 1
    par(mar = c(4.2, 4.8, 3, 1.5))
    y_max <- max(c(t1e_vec, pwr_vec, freq_t1e_vec, freq_pwr_vec), na.rm = TRUE)
    plot(n1_vals, pwr_vec, type = "l", col = col_bayes_power, lwd = 2.4,
         xlim = c(1, n_opt - 1L),
         ylim = c(0, min(1, y_max * 1.08)),
         xlab = expression(n[1] ~ "(stage-1 sample size)"),
         ylab = "Operating characteristic",
         main = bquote("Power & type-I error vs." ~ n[1] ~ (n == .(n_opt))),
         las = 1L, bty = "l")
    abline(h = seq(0, 1, by = 0.1), col = col_grid, lty = 3L)
    
    lines(n1_vals, t1e_vec, col = col_bayes_t1e, lwd = 2.4)
    
    if (has_freq_pwr && any(is.finite(freq_pwr_vec))) {
      lines(n1_vals, freq_pwr_vec, col = col_freq_power, lwd = 2.2)
    }
    if (has_freq_t1e && any(is.finite(freq_t1e_vec))) {
      lines(n1_vals, freq_t1e_vec, col = col_freq_t1e, lwd = 2.2)
    }
    
    abline(h = x$alpha, col = col_bayes_t1e, lty = 3L, lwd = 1.2)
    abline(h = x$target_power, col = col_bayes_power, lty = 3L, lwd = 1.2)
    if (!is.null(x$alpha_freq)) {
      abline(h = x$alpha_freq, col = col_freq_t1e, lty = 3L, lwd = 1.2)
    }
    if (!is.null(x$target_power_freq)) {
      abline(h = x$target_power_freq, col = col_freq_power, lty = 3L, lwd = 1.2)
    }
    
    abline(v = n1_opt, col = col_opt, lty = 2L, lwd = 1.5)
    points(n1_opt, pwr_vec[n1_opt], pch = 19L, col = col_bayes_power, cex = 1.3)
    points(n1_opt, t1e_vec[n1_opt], pch = 19L, col = col_bayes_t1e, cex = 1.3)
    
    legend_items <- c(
      sprintf("Bayesian power (H1: %s)", h1_label),
      sprintf("Bayesian type-I error (H0: %s)", h0_label)
    )
    legend_cols <- c(col_bayes_power, col_bayes_t1e)
    legend_lty  <- c(1L, 1L)
    legend_lwd  <- c(2.4, 2.4)
    
    if (has_freq_pwr && any(is.finite(freq_pwr_vec))) {
      legend_items <- c(legend_items, "Frequentist power")
      legend_cols <- c(legend_cols, col_freq_power)
      legend_lty  <- c(legend_lty, 1L)
      legend_lwd  <- c(legend_lwd, 2.2)
    }
    if (has_freq_t1e && any(is.finite(freq_t1e_vec))) {
      legend_items <- c(legend_items, "Frequentist type-I error")
      legend_cols <- c(legend_cols, col_freq_t1e)
      legend_lty  <- c(legend_lty, 1L)
      legend_lwd  <- c(legend_lwd, 2.2)
    }
    
    legend_items <- c(legend_items, bquote("optimal" ~ n[1] == .(n1_opt)))
    legend_cols <- c(legend_cols, col_opt)
    legend_lty  <- c(legend_lty, 2L)
    legend_lwd  <- c(legend_lwd, 1.5)
    
    legend("right", bty = "n", cex = 0.78,
           legend = legend_items,
           col = legend_cols, lty = legend_lty, lwd = legend_lwd)
    
    ## Panel 2
    par(mar = c(4.2, 4.5, 3, 1.5))
    plot(n1_vals, pce_vec, type = "l", col = col_pce, lwd = 2.2,
         xlim = c(1, n_opt - 1L),
         ylim = c(0, max(pce_vec, na.rm = TRUE) * 1.1),
         xlab = expression(n[1] ~ "(stage-1 sample size)"),
         ylab = sprintf("Bayesian PCE(%s) under H0", h0_label),
         main = bquote("PCE(H0) vs." ~ n[1] ~ (n == .(n_opt))),
         las = 1L, bty = "l")
    abline(h = seq(0, max(pce_vec) * 1.1, length.out = 6),
           col = col_grid, lty = 3L)
    abline(v = n1_opt, col = col_opt, lty = 2L, lwd = 1.5)
    points(n1_opt, pce_vec[n1_opt], pch = 19L, col = col_pce, cex = 1.3)
    legend("topright", bty = "n", cex = 0.82,
           legend = c("Bayesian PCE(H0)",
                      bquote("optimal" ~ n[1] == .(n1_opt))),
           col = c(col_pce, col_opt),
           lty = c(1L, 2L),
           lwd = c(2.2, 1.5))
    
    ## Panel 3
    par(mar = c(0.3, 0.5, 0.6, 0.5))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    rect(0, 0, 1, 1, border = col_grid, col = NA)
    
    ap <- x$analysis_prior
    bnd_str <- switch(
      x$direction,
      equivalence = sprintf("ROPE: [%.3f, %.3f]", lo_rope, hi_rope),
      noninferiority = sprintf("NI boundary: p >= %.4f", lo_rope),
      superiority = sprintf("Sup. boundary: p > %.4f", hi_rope)
    )
    
    txt <- c(
      sprintf("Optimal design (%s)", direction_label),
      "",
      sprintf("p0 = %.3f, delta = %.3f", x$p0, x$delta),
      bnd_str,
      sprintf("Analysis prior: Beta(%g, %g)", ap[1], ap[2]),
      sprintf("Null design prior: Beta(%g, %g)", aH0, bH0),
      sprintf("Alt. design prior: Beta(%g, %g)", aH1, bH1),
      "",
      sprintf("gamma_1 = %.2f", x$gamma_1),
      sprintf("gamma_eq = %.2f", x$gamma_eq),
      sprintf("gamma_diff = %.2f", x$gamma_diff),
      "",
      sprintf("n1 = %d, n2 = %d, n = %d", n1_opt, n2_opt, n_opt),
      sprintf("Bayesian type-I error : %.4f (<= %.2f)", d$type1_2st, x$alpha),
      sprintf("Bayesian power : %.4f (>= %.2f)", d$power_2st, x$target_power),
      sprintf("Bayesian PCE(H0) : %.4f", d$pce_2st),
      sprintf("EN(H0 prior) : %.2f", d$EN0),
      sprintf("EN(H1 prior) : %.2f", d$EN1),
      if (!is.null(d$freq_type1_2st))
        sprintf("Freq. type-I err.: %.4f", d$freq_type1_2st) else NULL,
      if (!is.null(d$freq_power_2st))
        sprintf("Freq. power      : %.4f", d$freq_power_2st) else NULL,
      "",
      sprintf("Optimality : %s", x$optimality),
      sprintf("Feasible designs: %d", nrow(x$candidates)),
      "",
      sprintf("Success: %s", crit_str),
      "",
      paste0("C1: {",
             if (length(cont) <= 12L)
               paste(cont, collapse = ", ")
             else
               paste0(paste(head(cont, 6L), collapse = ", "),
                      " ... ",
                      paste(tail(cont, 3L), collapse = ", ")),
             "}")
    )
    txt <- txt[!vapply(txt, is.null, logical(1L))]
    
    y_pos <- seq(0.975, 0.03, length.out = length(txt))
    text(0.05, y_pos, txt, adj = c(0, 1), cex = 0.70, family = "mono",
         col = ifelse(seq_along(txt) == 1L, "black", "grey20"),
         font = ifelse(seq_along(txt) == 1L, 2L, 1L))
    
    ## Panel 4
    par(mar = c(4, 4.2, 3, 1.5))
    .prior_panel(aH0, bH0, col_h0,
                 bquote("Null design prior" ~
                          "Beta(" * .(aH0) * "," ~ .(bH0) * ")"))
    leg_pos <- if (aH0 <= 2) "topright" else "topleft"
    legend(leg_pos, bty = "n", cex = 0.78,
           legend = c(sprintf("%s region", h1_label), bquote(p[0])),
           fill = c(col_rope, NA), border = c("grey60", NA),
           lty = c(NA, 3L), col = c(NA, "grey40"), merge = FALSE)
    
    ## Panel 5
    par(mar = c(4, 4.2, 3, 1.5))
    .prior_panel(aH1, bH1, col_h1,
                 bquote("Alt. design prior" ~
                          "Beta(" * .(aH1) * "," ~ .(bH1) * ")"))
    
    ## Panel 6
    par(mar = c(4, 4.2, 3, 1.5))
    .prior_panel(aA, bA, col_ap,
                 bquote("Analysis prior" ~
                          "Beta(" * .(aA) * "," ~ .(bA) * ")"))
    
  } else if (type == "interim") {
    
    par(mar = c(4.5, 4.5, 3, 3.5), mfrow = c(1L, 1L), ...)
    y1_all <- 0:n1_opt
    post1 <- .post_prob(y1_all, n1_opt)
    in_cont <- y1_all %in% cont
    col_pts <- ifelse(in_cont, col_cont, col_stop)
    pch_pts <- ifelse(in_cont, 19L, 4L)
    
    plot(y1_all, post1, type = "n",
         xlim = c(-0.5, n1_opt + 0.5), ylim = c(0, 1),
         xlab = expression(y[1] ~ "(stage-1 responses)"),
         ylab = ylab_post1,
         main = bquote("Stage 1 interim decision" ~
                         (n[1] == .(n1_opt) * "," ~
                            gamma[1] == .(x$gamma_1) * "," ~
                            "continue if posterior H1 support > " * .(1 - x$gamma_1))),
         las = 1L, bty = "l")
    abline(h = seq(0, 1, by = 0.1), col = col_grid, lty = 3L)
    abline(h = 1 - x$gamma_1, col = col_bnd, lty = 2L, lwd = 1.8)
    mtext(bquote(1 - gamma[1] == .(1 - x$gamma_1)),
          side = 4, at = 1 - x$gamma_1, las = 1, cex = 0.78,
          col = col_bnd, line = 0.5)
    if (length(cont) > 0L)
      rect(min(cont) - 0.45, -0.01, max(cont) + 0.45, 1.01,
           col = grDevices::adjustcolor(col_cont, alpha.f = 0.10), border = NA)
    lines(y1_all, post1, col = "grey50", lwd = 1)
    points(y1_all, post1, col = col_pts, pch = pch_pts, cex = 1.3)
    legend("topright", bty = "n", cex = 0.85,
           legend = c("Continue to stage 2", "Stop for futility",
                      bquote(gamma[1] ~ "boundary")),
           col = c(col_cont, col_stop, col_bnd),
           pch = c(19L, 4L, NA),
           lty = c(NA, NA, 2L),
           lwd = c(NA, NA, 1.8))
    
  } else {
    
    par(mar = c(4.5, 4.5, 3, 3.5), mfrow = c(1L, 1L), ...)
    y_all <- 0:n_opt
    post_f <- .post_prob(y_all, n_opt)
    accept <- post_f >= x$gamma_eq
    
    reachable <- vapply(y_all, function(y) {
      any(cont <= y & (y - cont) <= n2_opt)
    }, logical(1L))
    
    col_pts <- ifelse(accept, col_cont, col_stop)
    pch_pts <- ifelse(reachable, 19L, 1L)
    
    plot(y_all, post_f, type = "n",
         xlim = c(-0.5, n_opt + 0.5), ylim = c(0, 1),
         xlab = expression(y ~ "(total responses)"),
         ylab = ylab_post,
         main = bquote("Final decision" ~
                         (n == .(n_opt) * "," ~ gamma[eq] == .(x$gamma_eq))),
         las = 1L, bty = "l")
    abline(h = seq(0, 1, by = 0.1), col = col_grid, lty = 3L)
    abline(h = x$gamma_eq, col = col_bnd, lty = 2L, lwd = 1.8)
    mtext(bquote(gamma[eq] == .(x$gamma_eq)),
          side = 4, at = x$gamma_eq, las = 1, cex = 0.78,
          col = col_bnd, line = 0.5)
    acc_y <- y_all[accept]
    if (length(acc_y) > 0L)
      rect(min(acc_y) - 0.45, -0.01, max(acc_y) + 0.45, 1.01,
           col = grDevices::adjustcolor(col_cont, alpha.f = 0.10), border = NA)
    lines(y_all, post_f, col = "grey50", lwd = 1)
    points(y_all, post_f, col = col_pts, pch = pch_pts, cex = 1.3)
    legend("bottomright", bty = "n", cex = 0.85,
           legend = c(accept_label, reject_label,
                      "Reachable via C1", "Unreachable",
                      bquote(gamma[eq] ~ "boundary")),
           col = c(col_cont, col_stop, "black", "black", col_bnd),
           pch = c(19L, 4L, 19L, 1L, NA),
           lty = c(NA, NA, NA, NA, 2L),
           lwd = c(NA, NA, NA, NA, 1.8))
  }
  
  invisible(x)
}