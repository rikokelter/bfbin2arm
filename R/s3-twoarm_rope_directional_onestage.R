# S3 methods for one-stage two-arm directional (NI / superiority / nested)
# ROPE designs.

.format_rope2arm_directional_number <- function(x, digits = 4L) {
  if (is.null(x) || length(x) != 1L || !is.finite(x)) return("NA")
  formatC(x, digits = digits, format = "f")
}

.rope2arm_directional_selected_fit <- function(x) {
  if (is.null(x$selected_fit)) {
    stop("No feasible design was found; there is no selected fit to summarize.",
         call. = FALSE)
  }
  x$selected_fit
}

#' @export
print.twoarm_rope_directional_onestage_design <- function(x, ...) {
  direction_label <- switch(
    x$inputs$direction,
    noninferiority = "One-stage two-arm ROPE non-inferiority design",
    superiority = "One-stage two-arm ROPE superiority design",
    nested = "One-stage two-arm nested ROPE non-inferiority/superiority design"
  )
  
  cat(direction_label, "\n")
  cat(strrep("-", nchar(direction_label)), "\n\n")
  
  cat("Non-inferiority margin delta_ni:", x$inputs$delta_ni, "\n")
  cat("Superiority margin delta_sup:", x$inputs$delta_sup, "\n")
  cat(
    "Thresholds: gamma_ni =", x$inputs$gamma_ni,
    ", gamma_sup =", x$inputs$gamma_sup,
    ", gamma_inf =", x$inputs$gamma_inf, "\n\n"
  )
  
  if (is.null(x$selected)) {
    cat("No feasible per-arm sample size was found in [",
        x$inputs$nmin, ", ", x$inputs$nmax, "].\n", sep = "")
    return(invisible(x))
  }
  
  fit <- .rope2arm_directional_selected_fit(x)
  
  cat("Selected design:\n")
  cat("  nC =", fit$nC, ", nT =", fit$nT, ", N =", fit$N, "\n\n")
  
  if (!is.null(x$inputs$target_ni_power) && is.finite(x$inputs$target_ni_power)) {
    cat(
      "Target non-inferiority (or superiority) power:",
      .format_rope2arm_directional_number(x$inputs$target_ni_power), "\n"
    )
  }
  
  if (!is.null(x$inputs$target_sup_power) && is.finite(x$inputs$target_sup_power)) {
    cat(
      "Target superiority power:",
      .format_rope2arm_directional_number(x$inputs$target_sup_power), "\n"
    )
  }
  
  if (!is.null(x$inputs$target_false_ni) && is.finite(x$inputs$target_false_ni)) {
    cat(
      "Maximum predictive false non-inferiority probability:",
      .format_rope2arm_directional_number(x$inputs$target_false_ni), "\n"
    )
  }
  
  if (!is.null(x$inputs$target_false_sup) && is.finite(x$inputs$target_false_sup)) {
    cat(
      "Maximum predictive false superiority probability:",
      .format_rope2arm_directional_number(x$inputs$target_false_sup), "\n"
    )
  }
  
  cat("\n")
  
  cat(
    "Non-inferiority (or superiority) power, P(NI or sup | G_ni):",
    .format_rope2arm_directional_number(fit$ni_power), "\n"
  )
  cat(
    "Superiority power, P(sup | G_sup):",
    .format_rope2arm_directional_number(fit$sup_power), "\n"
  )
  cat(
    "False non-inferiority probability, P(NI or sup | G_inf):",
    .format_rope2arm_directional_number(fit$false_ni), "\n"
  )
  cat(
    "False superiority probability, P(sup | G_not_sup):",
    .format_rope2arm_directional_number(fit$false_sup), "\n"
  )
  
  if (!is.null(fit$freq_ni_power) && is.finite(fit$freq_ni_power)) {
    cat(
      "\nPointwise frequentist non-inferiority power:",
      .format_rope2arm_directional_number(fit$freq_ni_power), "\n"
    )
  }
  
  if (!is.null(fit$freq_sup_power) && is.finite(fit$freq_sup_power)) {
    cat(
      "Pointwise frequentist superiority power:",
      .format_rope2arm_directional_number(fit$freq_sup_power), "\n"
    )
  }
  
  invisible(x)
}

#' @export
summary.twoarm_rope_directional_onestage_design <- function(object, ...) {
  fit <- .rope2arm_directional_selected_fit(object)
  
  out <- list(
    inputs = object$inputs,
    nC = fit$nC, nT = fit$nT, N = fit$N,
    
    ni_power = fit$ni_power, sup_power = fit$sup_power,
    false_ni = fit$false_ni, false_sup = fit$false_sup,
    
    decision_g_inf = c(
      superiority = fit$p_sup_g_inf,
      noninferiority_only = fit$p_ni_only_g_inf,
      inferiority = fit$p_inf_g_inf,
      inconclusive = fit$p_inc_g_inf
    ),
    decision_g_ni = c(
      superiority = fit$p_sup_g_ni,
      noninferiority_only = fit$p_ni_only_g_ni,
      inferiority = fit$p_inf_g_ni,
      inconclusive = fit$p_inc_g_ni
    ),
    decision_g_sup = c(
      superiority = fit$p_sup_g_sup,
      noninferiority_only = fit$p_ni_only_g_sup,
      inferiority = fit$p_inf_g_sup,
      inconclusive = fit$p_inc_g_sup
    ),
    
    freq_ni_power = fit$freq_ni_power,
    freq_sup_power = fit$freq_sup_power,
    freq_false_ni = fit$freq_false_ni,
    freq_false_sup = fit$freq_false_sup
  )
  
  class(out) <- "summary.twoarm_rope_directional_onestage_design"
  out
}

#' @export
print.summary.twoarm_rope_directional_onestage_design <- function(x, ...) {
  cat("Summary: one-stage two-arm directional ROPE design\n")
  cat("-------------------------------------------------\n\n")
  
  cat("Direction:", x$inputs$direction, "\n")
  cat("delta_ni =", x$inputs$delta_ni, ", delta_sup =", x$inputs$delta_sup, "\n")
  cat("nC =", x$nC, ", nT =", x$nT, ", N =", x$N, "\n\n")
  
  cat("Non-inferiority (or superiority) power:",
      .format_rope2arm_directional_number(x$ni_power), "\n")
  cat("Superiority power:",
      .format_rope2arm_directional_number(x$sup_power), "\n")
  cat("False non-inferiority probability:",
      .format_rope2arm_directional_number(x$false_ni), "\n")
  cat("False superiority probability:",
      .format_rope2arm_directional_number(x$false_sup), "\n\n")
  
  cat("Decision probabilities under G_inf (superiority, NI-only, inferiority, inconclusive):\n")
  cat(" ", paste(vapply(x$decision_g_inf, .format_rope2arm_directional_number, character(1)),
                 collapse = "  "), "\n\n")
  
  cat("Decision probabilities under G_ni:\n")
  cat(" ", paste(vapply(x$decision_g_ni, .format_rope2arm_directional_number, character(1)),
                 collapse = "  "), "\n\n")
  
  cat("Decision probabilities under G_sup:\n")
  cat(" ", paste(vapply(x$decision_g_sup, .format_rope2arm_directional_number, character(1)),
                 collapse = "  "), "\n")
  
  invisible(x)
}

#' @export
as.data.frame.twoarm_rope_directional_onestage_design <- function(x, ...) {
  if (is.null(x$grid)) {
    stop("returngrid = TRUE is required to extract the operating-characteristic grid.",
         call. = FALSE)
  }
  x$grid
}

#' Plot a two-arm directional ROPE design
#'
#' @param x An object of class \code{twoarm_rope_directional_onestage_design}.
#' @param which Character, one of \code{"default"} (six-panel layout),
#'   \code{"oc"}, \code{"summary"}, \code{"analysis_prior"},
#'   \code{"design_prior"}, or \code{"heatmap"}.
#' @param ... Currently unused.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_hline
#'   geom_tile scale_fill_manual labs theme_minimal theme element_text
#'   element_blank annotate coord_cartesian
#' @importFrom patchwork plot_layout wrap_plots plot_annotation
#' @export
plot.twoarm_rope_directional_onestage_design <- function(
    x,
    which = c("default", "oc", "summary", "analysis_prior", "design_prior", "heatmap"),
    ...
) {
  which <- match.arg(which)
  
  panel_oc      <- .plot_directional_oc_panel(x)
  panel_summary <- .plot_directional_summary_panel(x)
  panel_ana     <- .plot_directional_analysis_prior_panel(x)
  panel_design  <- .plot_directional_design_prior_panel(x)
  panel_heat    <- .plot_directional_heatmap_panel(x)
  
  if (which == "oc")             return(panel_oc)
  if (which == "summary")        return(panel_summary)
  if (which == "analysis_prior") return(panel_ana)
  if (which == "design_prior")   return(panel_design)
  if (which == "heatmap")        return(panel_heat)
  
  layout <- c(
    patchwork::area(t = 1, l = 1, b = 1, r = 2),  # OC panel: top row, spans both columns
    patchwork::area(t = 2, l = 1, b = 2, r = 1),  # analysis prior: middle-left
    patchwork::area(t = 2, l = 2, b = 2, r = 2),  # design prior: middle-right
    patchwork::area(t = 3, l = 1, b = 3, r = 1),  # summary: bottom-left
    patchwork::area(t = 3, l = 2, b = 3, r = 2)   # heatmap: bottom-right
  )
  
  (panel_oc + panel_ana + panel_design + panel_summary + panel_heat) +
    patchwork::plot_layout(
      design = layout,
      widths = c(1, 1),
      heights = c(1.2, 1.1, 1.4)
    ) +
    patchwork::plot_annotation(
      title = sprintf(
        "Two-arm directional ROPE design (%s), n = %d per arm",
        x$inputs$direction, x$nstar
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
    )
}

#' @keywords internal
.plot_directional_summary_panel <- function(x) {
  fit <- x$selected_fit
  
  fmt <- function(v, digits = 4L) {
    if (is.null(v) || length(v) == 0L || is.na(v)) return("NA")
    formatC(v, digits = digits, format = "f")
  }
  
  rows <- list(
    c("Direction", x$inputs$direction),
    c("Non-inferiority margin", fmt(x$inputs$delta_ni)),
    c("Superiority margin", fmt(x$inputs$delta_sup)),
    c("Control arm sample size", as.character(fit$nC)),
    c("Treatment arm sample size", as.character(fit$nT)),
    c("Total sample size", as.character(fit$N)),
    c("NI (or superiority) power", fmt(fit$ni_power)),
    c("Superiority power", fmt(fit$sup_power)),
    c("False NI probability", fmt(fit$false_ni)),
    c("False superiority probability", fmt(fit$false_sup)),
    c("P(superiority | non-inferiority)", fmt(fit$p_sup_g_ni)),
    c("P(NI only | non-inferiority)", fmt(fit$p_ni_only_g_ni))
  )
  
  has_freq <- !all(vapply(
    list(fit$freq_ni_power, fit$freq_sup_power, fit$freq_false_ni, fit$freq_false_sup),
    function(v) is.null(v) || is.na(v), logical(1)
  ))
  
  if (has_freq) {
    rows <- c(rows, list(
      c("Freq. NI power", fmt(fit$freq_ni_power)),
      c("Freq. superiority power", fmt(fit$freq_sup_power)),
      c("Freq. false NI (type I error)", fmt(fit$freq_false_ni)),
      c("Freq. false superiority (type I error)", fmt(fit$freq_false_sup))
    ))
  }
  
  df <- data.frame(
    Quantity = vapply(rows, `[`, character(1), 1L),
    Value    = vapply(rows, `[`, character(1), 2L),
    stringsAsFactors = FALSE
  )
  
  tbl_theme <- gridExtra::ttheme_minimal(
    core = list(
      fg_params = list(hjust = 0, x = 0.02, fontsize = 8),
      bg_params = list(fill = c("white", "grey95"), col = NA)
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 8.5),
      bg_params = list(fill = "grey85", col = NA)
    ),
    padding = grid::unit(c(3, 2), "mm")
  )
  
  tbl <- gridExtra::tableGrob(df, rows = NULL, theme = tbl_theme)
  
  # Stretch the table grob to fill the full available viewport, so it
  # scales with the panel rather than overflowing it at a fixed size.
  tbl$widths  <- grid::unit(rep(1 / ncol(tbl), ncol(tbl)), "npc")
  tbl$heights <- grid::unit(rep(1 / nrow(tbl), nrow(tbl)), "npc")
  
  title_grob <- grid::textGrob(
    "Selected design summary",
    gp = grid::gpar(fontsize = 10, fontface = "bold")
  )
  
  combined <- gridExtra::arrangeGrob(
    title_grob, tbl,
    ncol = 1,
    heights = grid::unit(c(0.08, 0.92), "npc")
  )
  
  patchwork::wrap_elements(full = combined)
}

#' @keywords internal
.plot_directional_oc_panel <- function(x) {
  grid <- x$grid
  
  if (is.null(grid) || nrow(grid) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No grid stored\n(returngrid = FALSE)") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    )
  }
  
  oc_long <- data.frame(
    n = rep(grid$n, 2),
    value = c(grid$ni_power, grid$false_ni),
    quantity = rep(c("Non-inferiority power", "False non-inferiority probability"), each = nrow(grid))
  )
  
  add_curve <- function(col, label) {
    if (col %in% names(grid) && any(!is.na(grid[[col]]))) {
      oc_long <<- rbind(
        oc_long,
        data.frame(n = grid$n, value = grid[[col]], quantity = label)
      )
    }
  }
  
  add_curve("sup_power",       "Superiority power")
  add_curve("freq_ni_power",   "Freq. NI power")
  add_curve("freq_sup_power",  "Freq. superiority power")
  add_curve("freq_false_ni",   "Freq. false NI (type I error)")
  add_curve("freq_false_sup",  "Freq. false superiority (type I error)")
  
  targets <- list(
    list(value = x$inputs$target_ni_power,      label = "Target: NI power"),
    list(value = x$inputs$target_sup_power,      label = "Target: superiority power"),
    list(value = x$inputs$target_false_ni,       label = "Target: false NI"),
    list(value = x$inputs$target_false_sup,      label = "Target: false superiority"),
    list(value = x$inputs$target_freq_ni_power,  label = "Target: freq. NI power"),
    list(value = x$inputs$target_freq_sup_power, label = "Target: freq. superiority power"),
    list(value = x$inputs$target_freq_false_ni,  label = "Target: freq. false NI"),
    list(value = x$inputs$target_freq_false_sup, label = "Target: freq. false superiority")
  )
  targets <- Filter(function(t) !is.null(t$value) && !is.na(t$value), targets)
  
  target_df <- if (length(targets) > 0L) {
    data.frame(
      value = vapply(targets, function(t) t$value, numeric(1)),
      label = vapply(targets, function(t) t$label, character(1))
    )
  } else {
    NULL
  }
  
  okabe_ito <- c(
    "Non-inferiority power"                 = "#0072B2",
    "False non-inferiority probability"     = "#D55E00",
    "Superiority power"                     = "#009E73",
    "Freq. NI power"                        = "#56B4E9",
    "Freq. superiority power"               = "#CC79A7",
    "Freq. false NI (type I error)"         = "#E69F00",
    "Freq. false superiority (type I error)" = "#F0E442",
    "Target: NI power"                       = "#0072B2",
    "Target: superiority power"              = "#009E73",
    "Target: false NI"                       = "#D55E00",
    "Target: false superiority"              = "#000000",
    "Target: freq. NI power"                 = "#56B4E9",
    "Target: freq. superiority power"        = "#CC79A7",
    "Target: freq. false NI"                 = "#E69F00",
    "Target: freq. false superiority"        = "#F0E442"
  )
  
  p <- ggplot2::ggplot(oc_long, ggplot2::aes(x = n, y = value, colour = quantity)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_vline(xintercept = x$nstar, linetype = "dashed", colour = "grey40")
  
  if (!is.null(target_df)) {
    p <- p +
      ggplot2::geom_hline(
        data = target_df,
        ggplot2::aes(yintercept = value, colour = label),
        linetype = "dotted", linewidth = 0.5
      )
  }
  
  p +
    ggplot2::scale_colour_manual(values = okabe_ito) +
    ggplot2::labs(
      title = "Operating characteristics vs. per-arm sample size",
      x = "Per-arm sample size n", y = "Probability", colour = NULL
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}


#' @keywords internal
.plot_directional_analysis_prior_panel <- function(x) {
  aC <- unlist(x$inputs$analysis_prior_C)
  aT <- unlist(x$inputs$analysis_prior_T)
  
  if (length(aC) != 2L || length(aT) != 2L || anyNA(c(aC, aT))) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Analysis prior not available") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    )
  }
  
  p_seq <- seq(0, 1, length.out = 400)
  
  df <- data.frame(
    p = rep(p_seq, 2),
    density = c(
      stats::dbeta(p_seq, aC[1L], aC[2L]),
      stats::dbeta(p_seq, aT[1L], aT[2L])
    ),
    arm = rep(c("Control analysis prior", "Treatment analysis prior"), each = length(p_seq))
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = p, y = density, colour = arm)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = "Analysis priors",
      x = "Response probability p", y = "Density", colour = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

#' @keywords internal
#' @keywords internal
#' @keywords internal
.plot_directional_design_prior_panel <- function(x) {
  priors <- list(
    "Non-inferiority" = x$inputs$design_prior_ni,
    "Inferiority"      = x$inputs$design_prior_inf,
    "Superiority"      = x$inputs$design_prior_sup,
    "Not superiority"  = x$inputs$design_prior_not_sup
  )
  priors <- priors[!vapply(priors, is.null, logical(1))]
  
  if (length(priors) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No design priors stored") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    )
  }
  
  density_rows <- list()
  spike_rows <- list()
  
  for (nm in names(priors)) {
    nodes <- .directional_design_nodes(priors[[nm]], n_nodes = 128L)
    rd <- nodes$pT - nodes$pC
    w  <- nodes$weights
    
    if (stats::sd(rd) < 1e-4) {
      spike_rows[[nm]] <- data.frame(
        x = stats::weighted.mean(rd, w),
        scenario = nm
      )
    } else {
      density_rows[[nm]] <- data.frame(
        risk_difference = rd, weight = w, scenario = nm
      )
    }
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = -x$inputs$delta_ni, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = x$inputs$delta_sup, linetype = "dotted")
  
  if (length(density_rows) > 0L) {
    df_density <- do.call(rbind, density_rows)
    p <- p +
      ggplot2::geom_density(
        data = df_density,
        ggplot2::aes(x = risk_difference, weight = weight, colour = scenario, fill = scenario),
        alpha = 0.2
      )
  }
  
  if (length(spike_rows) > 0L) {
    df_spike <- do.call(rbind, spike_rows)
    p <- p +
      ggplot2::geom_vline(
        data = df_spike,
        ggplot2::aes(xintercept = x, colour = scenario),
        linetype = "solid", linewidth = 1
      )
  }
  
  p +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      title = "Design priors",
      x = expression(p[T] - p[C]), y = "Density", colour = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

#' @keywords internal
.plot_directional_heatmap_panel <- function(x) {
  fit <- x$selected_fit
  
  decision <- .directional_decision_matrix_twoarm(
    fit$nC, fit$nT,
    delta_ni = x$inputs$delta_ni, delta_sup = x$inputs$delta_sup,
    gamma_ni = x$inputs$gamma_ni, gamma_sup = x$inputs$gamma_sup, gamma_inf = x$inputs$gamma_inf,
    analysis_prior_C = x$inputs$analysis_prior_C, analysis_prior_T = x$inputs$analysis_prior_T,
    integration = "quantile", quad_nodes = 64L
  )
  
  category <- matrix("Inconclusive", nrow(decision$superiority), ncol(decision$superiority))
  category[decision$noninferiority_only] <- "NI only"
  category[decision$superiority] <- "Superiority"
  category[decision$inferiority] <- "Inferiority"
  
  df <- data.frame(
    yC = rep(0:fit$nC, times = fit$nT + 1L),
    yT = rep(0:fit$nT, each = fit$nC + 1L),
    category = as.vector(category)
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = yT, y = yC, fill = category)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      values = c(
        "Superiority" = "#1b9e77",
        "NI only" = "#7570b3",
        "Inferiority" = "#d95f02",
        "Inconclusive" = "grey85"
      )
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::labs(
      title = "Decision regions",
      x = expression(y[T]),
      y = expression(y[C]),
      fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

.rope2arm_directional_plot_oc <- function(x) {
  grid <- x$grid
  
  if (is.null(grid)) {
    stop("returngrid = TRUE is required to plot operating characteristics.",
         call. = FALSE)
  }
  
  cols <- c(
    ni_power = "#0072B2",
    sup_power = "#009E73",
    false_ni = "#D55E00",
    false_sup = "#CC79A7"
  )
  
  yvals <- c(grid$ni_power, grid$sup_power, grid$false_ni, grid$false_sup)
  ylim <- range(yvals, na.rm = TRUE)
  
  graphics::plot(
    grid$n, grid$ni_power,
    type = "l", lwd = 1.6, col = cols["ni_power"],
    xlab = "Per-arm sample size n", ylab = "Predictive probability",
    ylim = ylim, main = "Operating characteristics"
  )
  
  if (any(is.finite(grid$sup_power))) {
    graphics::lines(grid$n, grid$sup_power, lwd = 1.6, col = cols["sup_power"])
  }
  if (any(is.finite(grid$false_ni))) {
    graphics::lines(grid$n, grid$false_ni, lwd = 1.6, lty = 2, col = cols["false_ni"])
  }
  if (any(is.finite(grid$false_sup))) {
    graphics::lines(grid$n, grid$false_sup, lwd = 1.6, lty = 2, col = cols["false_sup"])
  }
  
  if (!is.null(x$inputs$target_ni_power) && is.finite(x$inputs$target_ni_power)) {
    graphics::abline(h = x$inputs$target_ni_power, lty = 3, col = cols["ni_power"])
  }
  if (!is.null(x$inputs$target_sup_power) && is.finite(x$inputs$target_sup_power)) {
    graphics::abline(h = x$inputs$target_sup_power, lty = 3, col = cols["sup_power"])
  }
  
  graphics::legend(
    "bottomright",
    legend = c("NI (or sup.) power", "Superiority power",
               "False NI probability", "False superiority probability"),
    col = unname(cols), lwd = 1.6, lty = c(1, 1, 2, 2), bty = "n", cex = 0.8
  )
  
  invisible(x)
}

.rope2arm_directional_plot_heatmap <- function(x) {
  fit <- .rope2arm_directional_selected_fit(x)
  
  if (is.null(fit$decision)) {
    stop("returnmatrices = TRUE is required to plot the decision heatmap.",
         call. = FALSE)
  }
  
  decision <- fit$decision
  
  code <- matrix(
    0L,
    nrow(decision$superiority),
    ncol(decision$superiority),
    dimnames = dimnames(decision$superiority)
  )
  
  code[decision$noninferiority_only] <- 1L
  code[decision$superiority] <- 2L
  code[decision$inferiority] <- 3L
  
  cols <- c("grey85", "#56B4E9", "#009E73", "#D55E00")
  
  graphics::image(
    x = decision$yC, y = decision$yT, z = t(code),
    col = cols, breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
    xlab = "yC", ylab = "yT",
    main = "Directional decision regions"
  )
  
  graphics::legend(
    "topleft",
    legend = c("Inconclusive", "Non-inferiority only", "Superiority", "Inferiority"),
    fill = cols, bty = "n", cex = 0.8
  )
  
  invisible(x)
}

.rope2arm_directional_plot_priors <- function(x) {
  fit <- .rope2arm_directional_selected_fit(x)
  
  panels <- list(
    pred_inf = "Design prior G_inf",
    pred_ni = "Design prior G_ni",
    pred_sup = "Design prior G_sup",
    pred_not_sup = "Design prior G_not_sup"
  )
  
  available <- panels[vapply(names(panels), function(nm) !is.null(fit[[nm]]), logical(1))]
  
  if (length(available) == 0L) {
    stop("returnmatrices = TRUE is required to plot design-prior predictive mass.",
         call. = FALSE)
  }
  
  old_par <- graphics::par(mfrow = c(1L, length(available)))
  on.exit(graphics::par(old_par), add = TRUE)
  
  for (nm in names(available)) {
    mat <- fit[[nm]]
    graphics::image(
      x = as.numeric(rownames(mat)), y = as.numeric(colnames(mat)), z = t(mat),
      col = grDevices::hcl.colors(50, "Blues", rev = TRUE),
      xlab = "yC", ylab = "yT", main = available[[nm]]
    )
  }
  
  invisible(x)
}