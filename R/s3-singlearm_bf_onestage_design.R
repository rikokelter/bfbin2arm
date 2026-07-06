#' Print method for one-stage single-arm BF designs
#'
#' @param x An object of class \code{"singlearm_onestage_bf_design"}.
#' @param ... Currently unused.
#'
#' @return The input object \code{x}, invisibly.
#' @export
print.singlearm_onestage_bf_design <- function(x, ...) {
  cat("\nOne-stage single-arm Bayes factor design\n")
  cat("---------------------------------------\n")
  cat("Calibration mode:", x$calibration, "\n")
  cat("Sustain: ", x$inputs$sustain_n, " future n\n")
  cat("Status :", x$status, "\n")
  
  if (isTRUE(x$feasible)) {
    cat("Selected sample size :", x$design["n"], "\n")
    cat("Evidence threshold k (efficacy) :", x$inputs$k, "\n")
    if (!is.null(x$inputs$k_ce)) {
      cat("Evidence threshold k_f (futility) :", x$inputs$k_ce, "\n")
    }
  }
  
  invisible(x)
}


#' Summarize a one-stage single-arm BF design
#'
#' @param object An object of class \code{"singlearm_onestage_bf_design"}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.singlearm_onestage_bf_design"}.
#' @export
summary.singlearm_onestage_bf_design <- function(object, ...) {
  out <- list(
    call                      = object$call,
    feasible                  = object$feasible,
    status                    = object$status,
    calibration               = object$calibration,
    design                    = object$design,
    inputs                    = object$inputs,
    operating_characteristics = object$operating_characteristics
  )
  class(out) <- "summary.singlearm_onestage_bf_design"
  out
}


#' Print method for summaries of one-stage single-arm BF designs
#'
#' @param x An object of class \code{"summary.singlearm_onestage_bf_design"}.
#' @param digits Number of digits to print.
#' @param ... Currently unused.
#'
#' @return The input object \code{x}, invisibly.
#' @export
print.summary.singlearm_onestage_bf_design <- function(x, digits = 3, ...) {
  cat("Summary: One-stage single-arm Bayes factor design\n")
  cat("------------------------------------------------\n")
  cat("Calibration:", x$calibration, "\n")
  cat("Sustain: ", x$inputs$sustain_n, " future n\n")
  cat("Feasible :", x$feasible, "\n")
  cat("Status :", x$status, "\n")
  
  if (isTRUE(x$feasible)) {
    oc <- x$operating_characteristics
    cat("\nSelected design\n")
    cat(" n    :", x$design["n"], "\n")
    cat(" k    :", formatC(x$inputs$k, digits = digits, format = "fg"), "\n")
    if (!is.null(x$inputs$k_ce)) {
      cat(" k_ce :", formatC(x$inputs$k_ce, digits = digits, format = "fg"), "\n")
    }
    
    cat("\nOperating characteristics\n")
    cat(" Bayes power  :", formatC(oc$pfineff,       digits = digits, format = "f"), "\n")
    cat(" Bayes type-I :", formatC(oc$pfineff0,      digits = digits, format = "f"), "\n")
    cat(" CE(H0)       :", formatC(oc$pce0_corr,     digits = digits, format = "f"), "\n")
    cat(" Freq power   :", formatC(oc$pfineff_freq,  digits = digits, format = "f"), "\n")
    cat(" Freq type-I  :", formatC(oc$pfineff_freq0, digits = digits, format = "f"), "\n")
  }
  
  invisible(x)
}


## internal helper: truncated Beta(a, b) density on [lower, upper]
.trunc_beta_density <- function(a, b, lower, upper, n_grid = 512) {
  norm <- pbeta(upper, a, b) - pbeta(lower, a, b)
  if (!is.finite(norm) || norm <= 0) return(NULL)
  
  p_grid <- seq(lower, upper, length.out = n_grid)
  p_grid <- pmin(pmax(p_grid, 1e-9), 1 - 1e-9)
  dens   <- dbeta(p_grid, a, b) / norm
  
  list(x = p_grid, y = dens)
}


#' Plot a one-stage single-arm BF design
#'
#' Produces a 2x2 figure:
#' top-left: operating characteristic curves;
#' top-right: table-like summary;
#' bottom-left: design priors under H0 and H1;
#' bottom-right: analysis priors under H0 and H1.
#'
#' @param x An object of class \code{"singlearm_onestage_bf_design"}.
#' @param what Character string; one of \code{"all"} or \code{"oc"}.
#' @param legend_pos Position passed to \code{legend()}.
#' @param legend_inset Numeric inset for \code{legend()}.
#' @param col_h0 Colour for H0 priors in bottom panels.
#' @param col_h1 Colour for H1 priors in bottom panels.
#' @param prior_lwd Line width for prior density curves.
#' @param ... Currently unused.
#'
#' @return Invisibly returns \code{x}.
#' @importFrom graphics abline plot.window polygon
#' @importFrom grDevices adjustcolor
#' @importFrom utils head
#' @export
plot.singlearm_onestage_bf_design <- function(
    x,
    what         = c("all", "oc"),
    legend_pos   = "right",
    legend_inset = 0,
    col_h0       = "#0072B2",
    col_h1       = "#D55E00",
    prior_lwd    = 2,
    ...
) {
  what <- match.arg(what)
  sr   <- x$search_results
  
  if (is.null(sr) || nrow(sr) == 0L) {
    stop("No search results available for plotting.")
  }
  
  show_ce         <- isTRUE(x$inputs$target_ce_h0 > 0)
  show_freq_power <- !all(is.na(sr$freq_power))
  show_freq_type1 <- !all(is.na(sr$freq_type1))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  layout(matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE),
         widths  = c(1.0, 1.1),
         heights = c(1.0, 1.0))
  
  ## -----------------------------------------------------------------------
  ## Panel 1 (top-left): OC curves
  ## -----------------------------------------------------------------------
  par(mar = c(4, 4, 3, 0.5))
  
  y_candidates <- c(
    sr$power,
    sr$type1,
    if (show_ce)         sr$ce_h0      else NULL,
    if (show_freq_power) sr$freq_power else NULL,
    if (show_freq_type1) sr$freq_type1 else NULL,
    x$inputs$target_power + x$inputs$power_cushion,
    x$inputs$target_type1,
    if (show_ce)              x$inputs$target_ce_h0                           else NULL,
    if (!is.na(x$inputs$dp)) x$inputs$target_freq_power + x$inputs$power_cushion else NULL,
    if (!is.na(x$inputs$dp)) x$inputs$target_freq_type1                       else NULL
  )
  y_candidates <- y_candidates[is.finite(y_candidates)]
  
  ylim_use <- if (length(y_candidates)) {
    c(max(0, min(y_candidates, na.rm = TRUE) - 0.05),
      min(1, max(y_candidates, na.rm = TRUE) + 0.05))
  } else {
    c(0, 1)
  }
  
  plot(sr$n, sr$power,
       type = "l", lwd = 2, lty = 1, col = "#D55E00",
       ylim = ylim_use,
       xlab = "Total sample size n",
       ylab = "Operating characteristic",
       main = "Operating characteristics")
  
  lines(sr$n, sr$type1, lwd = 2, lty = 1, col = "#0072B2")
  
  if (show_ce)         lines(sr$n, sr$ce_h0,      lwd = 2, lty = 1, col = "#CC79A7")
  if (show_freq_power) lines(sr$n, sr$freq_power, lwd = 2, lty = 1, col = "#E69F00")
  if (show_freq_type1) lines(sr$n, sr$freq_type1, lwd = 2, lty = 1, col = "#56B4E9")
  
  abline(h = x$inputs$target_power + x$inputs$power_cushion, lty = 3, col = "#D55E00")
  abline(h = x$inputs$target_type1,                          lty = 3, col = "#0072B2")
  if (show_ce)              abline(h = x$inputs$target_ce_h0,                                 lty = 3, col = "#CC79A7")
  if (!is.na(x$inputs$dp)) abline(h = x$inputs$target_freq_power + x$inputs$power_cushion,   lty = 3, col = "#E69F00")
  if (!is.na(x$inputs$dp)) abline(h = x$inputs$target_freq_type1,                            lty = 3, col = "#56B4E9")
  if (!is.na(x$design["n"])) abline(v = x$design["n"], lty = 3, col = "grey40")
  
  leg_labels <- c("Bayesian power", "Bayesian type-I")
  leg_cols   <- c("#D55E00",        "#0072B2")
  leg_lty    <- c(1, 1)
  if (show_ce)         { leg_labels <- c(leg_labels, "Compelling evidence for H0"); leg_cols <- c(leg_cols, "#CC79A7"); leg_lty <- c(leg_lty, 1) }
  if (show_freq_power) { leg_labels <- c(leg_labels, "Frequentist power");         leg_cols <- c(leg_cols, "#E69F00"); leg_lty <- c(leg_lty, 1) }
  if (show_freq_type1) { leg_labels <- c(leg_labels, "Frequentist type-I");        leg_cols <- c(leg_cols, "#56B4E9"); leg_lty <- c(leg_lty, 1) }
  
  leg_args <- list(legend = leg_labels, col = leg_cols, lty = leg_lty,
                   lwd = 2, bty = "n", cex = 0.9, inset = legend_inset)
  if (is.character(legend_pos) && length(legend_pos) == 1L) {
    do.call(graphics::legend, c(list(x = legend_pos), leg_args))
  } else if (is.numeric(legend_pos) && length(legend_pos) == 2L) {
    do.call(graphics::legend, c(list(x = legend_pos[1], y = legend_pos[2]), leg_args))
  } else {
    stop("'legend_pos' must be a single character keyword or a numeric vector of length 2.")
  }
  
  ## -----------------------------------------------------------------------
  ## Panel 2 (top-right): compact table-like summary + hypotheses
  ## -----------------------------------------------------------------------
  par(mar = c(4, 1, 3, 1))
  plot.new()
  title("Selected design")
  
  oc         <- x$operating_characteristics
  test_label <- if (identical(x$inputs$type, "point")) "two-sided" else "directional"
  
  ## Hypotheses as plotmath
  h0_expr <- if (identical(x$inputs$type, "point")) {
    substitute(H[0] * ": p == " * p0v, list(p0v = x$inputs$p0))
  } else {
    substitute(H[0] * ": p <= " * p0v, list(p0v = x$inputs$p0))
  }
  h1_expr <- if (identical(x$inputs$type, "point")) {
    substitute(H[1] * ": p != " * p0v, list(p0v = x$inputs$p0))
  } else {
    substitute(H[1] * ": p > " * p0v, list(p0v = x$inputs$p0))
  }
  
  freq_lines <- if (is.na(x$inputs$dp)) {
    c("Frequentist power calculations carried out",
      "under: not requested")
  } else {
    c("Frequentist power calculations carried out",
      paste0("under: dp = ", formatC(x$inputs$dp, digits = 3, format = "fg")))
  }
  
  ## Build rows: hypotheses as dedicated rows, rendered via plotmath in left column
  rows <- list(
    c("Calibration mode:", x$calibration),
    c("Sustain:",     paste0(x$inputs$sustain_n, " future n")),
    c("Status:",      x$status),
    c("",            ""),
    c("Selected n:",  as.character(x$design["n"])),
    c("Efficacy threshold k:", formatC(x$inputs$k, digits = 3, format = "fg")),
    if (!is.null(x$inputs$k_ce) && show_ce)
      c("Futility threshold k_f:", formatC(x$inputs$k_ce, digits = 3, format = "fg")) else NULL,
    c("",            ""),
    c("Bayesian power:",   ifelse(is.null(oc), "NA",
                                 formatC(oc$pfineff,  digits = 3, format = "f"))),
    c("Bayesian type-I-error:",  ifelse(is.null(oc), "NA",
                                 formatC(oc$pfineff0, digits = 3, format = "f"))),
    if (show_ce)
      c("Compelling evidence for H0:", ifelse(is.null(oc) || is.na(oc$pce0_corr), "NA",
                                             formatC(oc$pce0_corr, digits = 3, format = "f"))) else NULL,
    c("Frequentist power:",  ifelse(is.null(oc) || is.na(oc$pfineff_freq),  "NA",
                                   formatC(oc$pfineff_freq,  digits = 3, format = "f"))),
    c("Frequentist type-I-error:", ifelse(is.null(oc) || is.na(oc$pfineff_freq0), "NA",
                                   formatC(oc$pfineff_freq0, digits = 3, format = "f"))),
    c("",              ""),
    ## placeholder rows for H0 and H1 (label handled via plotmath)
    c("H0_row",        ""),
    c("H1_row",        ""),
    c("",              ""),
    c(freq_lines[1],   ""),
    c(freq_lines[2],   "")
  )
  rows <- rows[!vapply(rows, is.null, logical(1))]
  
  n_rows <- length(rows)
  y_vals <- seq(0.95, 0.05, length.out = n_rows)
  
  ## Tight spacing between label and value columns
  x_lab <- 0.06
  x_val <- 0.45
  
  ## Draw all regular label–value pairs; skip the hypothesis placeholder rows
  for (i in seq_len(n_rows)) {
    lab <- rows[[i]][1]
    val <- rows[[i]][2]
    y   <- y_vals[i]
    
    if (!lab %in% c("H0_row", "H1_row")) {
      if (nzchar(lab)) {
        text(x = x_lab, y = y, labels = lab, adj = c(0, 1), cex = 0.90)
      }
      if (nzchar(val)) {
        text(x = x_val, y = y, labels = val, adj = c(0, 1), cex = 0.90)
      }
    }
  }
  
  ## Draw hypotheses with plotmath in the left column at their dedicated rows
  idx_h0 <- which(vapply(rows, function(r) r[1] == "H0_row", logical(1)))
  idx_h1 <- which(vapply(rows, function(r) r[1] == "H1_row", logical(1)))
  
  if (length(idx_h0) == 1L) {
    y_h0 <- y_vals[idx_h0]
    text(x = x_lab, y = y_h0, labels = h0_expr, adj = c(0, 1), cex = 0.90)
  }
  if (length(idx_h1) == 1L) {
    y_h1 <- y_vals[idx_h1]
    text(x = x_lab, y = y_h1, labels = h1_expr, adj = c(0, 1), cex = 0.90)
  }
  
  ## -----------------------------------------------------------------------
  ## Helper: prior panels
  ## -----------------------------------------------------------------------
  .draw_prior_panel <- function(a0, b0, a1, b1, p0, type, main_title) {
    if (identical(type, "direction")) {
      d0   <- .trunc_beta_density(a0, b0, lower = 0,  upper = p0)
      d1   <- .trunc_beta_density(a1, b1, lower = p0, upper = 1)
      lab0 <- paste0("H0: TrBeta(", a0, ", ", b0, ") on [0, ",   p0, "]")
      lab1 <- paste0("H1: TrBeta(", a1, ", ", b1, ") on (", p0, ", 1]")
    } else {
      d0   <- .trunc_beta_density(a0, b0, lower = 0, upper = 1)
      d1   <- .trunc_beta_density(a1, b1, lower = 0, upper = 1)
      lab0 <- paste0("H0: Beta(", a0, ", ", b0, ")")
      lab1 <- paste0("H1: Beta(", a1, ", ", b1, ")")
    }
    
    all_y <- c(if (!is.null(d0)) d0$y else NULL,
               if (!is.null(d1)) d1$y else NULL)
    ylim  <- if (length(all_y)) c(0, max(all_y) * 1.12) else c(0, 1)
    
    par(mar = c(4, 4, 3, 0.5))
    plot(NULL,
         xlim = c(0, 1), ylim = ylim,
         xlab = "Response probability p",
         ylab = "Prior density",
         main = main_title)
    
    if (!is.null(d0)) {
      polygon(c(d0$x[1], d0$x, d0$x[length(d0$x)]),
              c(0,        d0$y, 0),
              col = adjustcolor(col_h0, alpha.f = 0.20), border = NA)
      lines(d0$x, d0$y, col = col_h0, lwd = prior_lwd)
    }
    
    if (!is.null(d1)) {
      polygon(c(d1$x[1], d1$x, d1$x[length(d1$x)]),
              c(0,        d1$y, 0),
              col = adjustcolor(col_h1, alpha.f = 0.20), border = NA)
      lines(d1$x, d1$y, col = col_h1, lwd = prior_lwd)
    }
    
    if (identical(type, "direction")) {
      abline(v = p0, lty = 2, col = "grey40")
    }
    
    legend("topright",
           legend = c(lab0, lab1),
           col    = c(col_h0, col_h1),
           lty    = c(1, 1),
           lwd    = prior_lwd,
           bty    = "n",
           cex    = 0.80)
  }
  
  ## panel 3: design priors
  .draw_prior_panel(
    a0         = x$inputs$da0,
    b0         = x$inputs$db0,
    a1         = x$inputs$da1,
    b1         = x$inputs$db1,
    p0         = x$inputs$p0,
    type       = x$inputs$type,
    main_title = "Design priors"
  )
  
  ## panel 4: analysis priors
  .draw_prior_panel(
    a0         = x$inputs$a0,
    b0         = x$inputs$b0,
    a1         = x$inputs$a1,
    b1         = x$inputs$b1,
    p0         = x$inputs$p0,
    type       = x$inputs$type,
    main_title = "Analysis priors"
  )
  
  invisible(x)
}


#' Convert a one-stage single-arm BF design to a data frame
#'
#' @param x An object of class \code{"singlearm_onestage_bf_design"}.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Currently unused.
#'
#' @return A data frame with the search results.
#' @export
as.data.frame.singlearm_onestage_bf_design <- function(
    x, row.names = NULL, optional = FALSE, ...
) {
  x$search_results
}