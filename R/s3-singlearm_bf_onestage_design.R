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
  cat("Status     :", x$status, "\n")
  
  if (isTRUE(x$feasible)) {
    cat("Selected sample size :", x$design["n"], "\n")
    cat("Evidence threshold k (efficacy)         :", x$inputs$k, "\n")
    if (!is.null(x$inputs$k_ce)) {
      cat("Evidence threshold k_f (futility)     :", x$inputs$k_ce, "\n")
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
    call = object$call,
    feasible = object$feasible,
    status = object$status,
    calibration = object$calibration,
    design = object$design,
    inputs = object$inputs,
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
  cat("Feasible   :", x$feasible, "\n")
  cat("Status     :", x$status, "\n")
  
  if (isTRUE(x$feasible)) {
    oc <- x$operating_characteristics
    cat("\nSelected design\n")
    cat("  n          :", x$design["n"], "\n")
    cat("  k          :", formatC(x$inputs$k, digits = digits, format = "fg"), "\n")
    if (!is.null(x$inputs$k_ce)) {
      cat("  k_ce       :", formatC(x$inputs$k_ce, digits = digits, format = "fg"), "\n")
    }
    
    cat("\nOperating characteristics\n")
    cat("  Bayes power      :", formatC(oc$pfineff, digits = digits, format = "f"), "\n")
    cat("  Bayes type-I     :", formatC(oc$pfineff0, digits = digits, format = "f"), "\n")
    cat("  CE(H0)           :", formatC(oc$pce0_corr, digits = digits, format = "f"), "\n")
    cat("  Freq power       :", formatC(oc$pfineff_freq, digits = digits, format = "f"), "\n")
    cat("  Freq type-I      :", formatC(oc$pfineff_freq0, digits = digits, format = "f"), "\n")
  }
  
  invisible(x)
}

#' Plot a one-stage single-arm BF design
#'
#' @param x An object of class \code{"singlearm_onestage_bf_design"}.
#' @param what Character string; currently one of \code{"all"} or \code{"oc"}.
#' @param legend_pos Position passed to \code{\link[graphics:legend]{legend}}.
#'   Either a keyword such as \code{"topright"} or a numeric vector
#'   \code{c(x, y)}.
#' @param legend_inset Numeric inset passed to \code{legend()} when
#'   \code{legend_pos} is a keyword.
#' @param ... Currently unused.
#'
#' @return Invisibly returns \code{x}.
#' @importFrom graphics abline
#' @export
plot.singlearm_onestage_bf_design <- function(
    x,
    what = c("all", "oc"),
    legend_pos = "right",
    legend_inset = 0,
    ...
) {
  what <- match.arg(what)
  sr <- x$search_results
  
  if (is.null(sr) || nrow(sr) == 0L) {
    stop("No search results available for plotting.")
  }
  
  show_ce <- isTRUE(x$inputs$target_ce_h0 > 0)
  show_freq_power <- !all(is.na(sr$freq_power))
  show_freq_type1 <- !all(is.na(sr$freq_type1))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.0, 1.1))
  par(mar = c(4, 4, 2.5, 0.5))
  
  ## Left panel: OC curves
  y_candidates <- c(
    sr$power,
    sr$type1,
    if (show_ce) sr$ce_h0 else NULL,
    if (show_freq_power) sr$freq_power else NULL,
    if (show_freq_type1) sr$freq_type1 else NULL,
    x$inputs$target_power + x$inputs$power_cushion,
    x$inputs$target_type1,
    if (show_ce) x$inputs$target_ce_h0 else NULL,
    if (!is.na(x$inputs$dp)) x$inputs$target_freq_power + x$inputs$power_cushion else NULL,
    if (!is.na(x$inputs$dp)) x$inputs$target_freq_type1 else NULL
  )
  y_candidates <- y_candidates[is.finite(y_candidates)]
  
  ylim_use <- if (length(y_candidates)) {
    c(max(0, min(y_candidates, na.rm = TRUE) - 0.05),
      min(1, max(y_candidates, na.rm = TRUE) + 0.05))
  } else {
    c(0, 1)
  }
  
  plot(sr$n, sr$power, type = "l", lwd = 2, lty = 1, col = "#D55E00",
       ylim = ylim_use,
       xlab = "Total sample size n",
       ylab = "Operating characteristic",
       main = "Operating characteristics")
  
  lines(sr$n, sr$type1, lwd = 2, lty = 1, col = "#0072B2")
  
  if (show_ce) {
    lines(sr$n, sr$ce_h0, lwd = 2, lty = 1, col = "#CC79A7")
  }
  
  if (show_freq_power) {
    lines(sr$n, sr$freq_power, lwd = 2, lty = 1, col = "#E69F00")
  }
  
  if (show_freq_type1) {
    lines(sr$n, sr$freq_type1, lwd = 2, lty = 1, col = "#56B4E9")
  }
  
  ## Target/reference lines remain dotted for orientation
  abline(h = x$inputs$target_power + x$inputs$power_cushion,
         lty = 3, col = "#D55E00")
  abline(h = x$inputs$target_type1,
         lty = 3, col = "#0072B2")
  
  if (show_ce) {
    abline(h = x$inputs$target_ce_h0,
           lty = 3, col = "#CC79A7")
  }
  
  if (!is.na(x$inputs$dp)) {
    abline(h = x$inputs$target_freq_power + x$inputs$power_cushion,
           lty = 3, col = "#E69F00")
    abline(h = x$inputs$target_freq_type1,
           lty = 3, col = "#56B4E9")
  }
  
  if (!is.na(x$design["n"])) {
    abline(v = x$design["n"], lty = 3, col = "grey40")
  }
  
  legend_labels <- c("Bayes power", "Bayes type-I")
  legend_cols <- c("#D55E00", "#0072B2")
  legend_lty <- c(1, 1)
  
  if (show_ce) {
    legend_labels <- c(legend_labels, "CE(H0)")
    legend_cols <- c(legend_cols, "#CC79A7")
    legend_lty <- c(legend_lty, 1)
  }
  
  if (show_freq_power) {
    legend_labels <- c(legend_labels, "Freq power")
    legend_cols <- c(legend_cols, "#E69F00")
    legend_lty <- c(legend_lty, 1)
  }
  
  if (show_freq_type1) {
    legend_labels <- c(legend_labels, "Freq type-I")
    legend_cols <- c(legend_cols, "#56B4E9")
    legend_lty <- c(legend_lty, 1)
  }
  
  legend_args <- list(
    legend = legend_labels,
    col = legend_cols,
    lty = legend_lty,
    lwd = 2,
    bty = "n",
    cex = 0.9,
    inset = legend_inset
  )
  
  if (is.character(legend_pos) && length(legend_pos) == 1L) {
    do.call(graphics::legend, c(list(x = legend_pos), legend_args))
  } else if (is.numeric(legend_pos) && length(legend_pos) == 2L) {
    do.call(graphics::legend, c(list(x = legend_pos[1], y = legend_pos[2]), legend_args))
  } else {
    stop("'legend_pos' must be a single character keyword or a numeric vector of length 2.")
  }
  
  ## Right panel: text summary
  par(mar = c(4, 0.5, 2.5, 0.5))
  plot.new()
  title("Selected design")
  
  oc <- x$operating_characteristics
  
  test_label <- if (identical(x$inputs$type, "point")) "two-sided" else "directional"
  
  h0_line <- if (identical(x$inputs$type, "point")) {
    paste0("H0: p = ", formatC(x$inputs$p0, digits = 3, format = "fg"))
  } else {
    paste0("H0: p <= ", formatC(x$inputs$p0, digits = 3, format = "fg"))
  }
  
  h1_line <- if (identical(x$inputs$type, "point")) {
    paste0("H1: p != ", formatC(x$inputs$p0, digits = 3, format = "fg"))
  } else {
    paste0("H1: p > ", formatC(x$inputs$p0, digits = 3, format = "fg"))
  }
  
  freq_lines <- if (is.na(x$inputs$dp)) {
    c("Frequentist power calculations carried out",
      "under: not requested")
  } else {
    c("Frequentist power calculations carried out",
      paste0("under: dp = ", formatC(x$inputs$dp, digits = 3, format = "fg")))
  }
  
  text_lines <- c(
    paste0("Calibration: ", x$calibration),
    paste0("Sustain: ", x$inputs$sustain_n, " future n"),
    paste0("Status: ", x$status),
    "",
    paste0("Selected n: ", x$design["n"]),
    paste0("Evidence threshold k (efficacy): ", formatC(x$inputs$k, digits = 3, format = "fg")),
    if (!is.null(x$inputs$k_ce) && show_ce)
      paste0("Evidence threshold k_f (futility): ", formatC(x$inputs$k_ce, digits = 3, format = "fg"))
    else NULL,
    "",
    paste0("Bayes power: ",
           ifelse(is.null(oc), "NA",
                  formatC(oc$pfineff, digits = 3, format = "f"))),
    paste0("Bayes type-I: ",
           ifelse(is.null(oc), "NA",
                  formatC(oc$pfineff0, digits = 3, format = "f"))),
    if (show_ce)
      paste0("CE(H0): ",
             ifelse(is.null(oc) || is.na(oc$pce0_corr), "NA",
                    formatC(oc$pce0_corr, digits = 3, format = "f")))
    else NULL,
    paste0("Freq power: ",
           ifelse(is.null(oc) || is.na(oc$pfineff_freq), "NA",
                  formatC(oc$pfineff_freq, digits = 3, format = "f"))),
    paste0("Freq type-I: ",
           ifelse(is.null(oc) || is.na(oc$pfineff_freq0), "NA",
                  formatC(oc$pfineff_freq0, digits = 3, format = "f"))),
    "",
    paste0("Test: ", test_label),
    h0_line,
    h1_line,
    "",
    freq_lines[1],
    freq_lines[2]
  )
  
  y <- seq(0.95, 0.05, length.out = length(text_lines))
  text(x = 0.02, y = y, labels = text_lines, adj = c(0, 1), cex = 0.95)
  
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