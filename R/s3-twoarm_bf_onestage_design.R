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

#' Plot a one-stage two-arm Bayes factor design
#'
#' Displays operating characteristics, feasibility, or an overview
#' of design and analysis priors and operating characteristics.
#'
#' @param x An object of class \code{twoarm_onestage_bf_design}.
#' @param type Plot type: \code{"old"} for the combined layout,
#'   \code{"oc"} for operating characteristics, or
#'   \code{"feasibility"} for feasibility.
#' @param h0_design,h0_analysis Optional named numeric vectors
#'   \code{c(a = ..., b = ...)} specifying beta priors for the
#'   common response probability under the equality null.
#'   Defaults are obtained from \code{x$inputs$a_0_d},
#'   \code{x$inputs$b_0_d}, \code{x$inputs$a_0_a}, and
#'   \code{x$inputs$b_0_a}. Explicit values affect the displayed
#'   priors only, not the calibrated design.
#'   Not used for \code{test = "BF+-"}.
#' @param ... Additional arguments; currently unused.
#'
#' @details
#' The combined layout displays a shared mathematical title above
#' both prior rows. For \code{"BF+-"}, the title labels the
#' comparison as H+: p2 > p1 versus H-: p1 >= p2.
#' This label does not change the underlying calculations.
#'
#' A single prior legend appears below the second prior row.
#' Alternative-prior panels display beta basis densities rather
#' than conditional marginal densities induced by directional
#' restrictions on the joint parameter space.
#'
#' Operating-characteristic curves are solid, horizontal target
#' lines are dashed, and a dotted vertical line marks the first
#' jointly feasible design in the stored search results.
#' Individual target-attainment sample sizes appear at the right
#' edge and need not equal the jointly calibrated sample size.
#'
#' @return A \code{ggplot} object for \code{"oc"} or
#'   \code{"feasibility"}, or a \code{patchwork} object
#'   for \code{"old"}.
#'
#' @importFrom rlang .data
#' @export
plot.twoarm_onestage_bf_design <- function(
    x,
    type = c("old", "oc", "feasibility"),
    h0_design = NULL,
    h0_analysis = NULL,
    ...
) {
  type <- match.arg(type)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  
  sr <- x$search_results
  inp <- x$inputs
  
  if (is.null(sr) || nrow(sr) == 0L) {
    stop("No search results available for plotting.", call. = FALSE)
  }
  
  # Standardize column aliases.
  aliases <- list(
    c("n_total", "n"),
    c("type1", "t1e"),
    c("ce_h0", "pceH0")
  )
  
  for (pair in aliases) {
    if (!pair[1L] %in% names(sr) && pair[2L] %in% names(sr)) {
      sr[[pair[1L]]] <- sr[[pair[2L]]]
    }
    if (!pair[2L] %in% names(sr) && pair[1L] %in% names(sr)) {
      sr[[pair[2L]]] <- sr[[pair[1L]]]
    }
  }
  
  required <- c("n_total", "power", "type1", "ce_h0")
  missing_columns <- setdiff(required, names(sr))
  
  if (length(missing_columns)) {
    stop(
      "Missing search-result columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!is.numeric(sr$n_total) || any(!is.finite(sr$n_total))) {
    stop("Sample sizes must be finite numeric values.", call. = FALSE)
  }
  
  sr <- sr[order(sr$n_total), , drop = FALSE]
  ns <- sr$n_total
  
  available <- function(column) {
    column %in% names(sr) && any(is.finite(sr[[column]]))
  }
  
  valid_target <- function(z) {
    is.numeric(z) && length(z) == 1L && is.finite(z)
  }
  
  theme_common <- function() {
    ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "plain",
          size = 14,
          margin = ggplot2::margin(b = 8)
        ),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 10),
        legend.key.width = grid::unit(1.3, "cm"),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(6, 12, 6, 8)
      )
  }
  
  # ------------------------------------------------------------------
  # Feasibility.
  # ------------------------------------------------------------------
  
  if (type == "feasibility") {
    if (!all(c("feasible_pointwise", "feasible") %in% names(sr))) {
      stop("Feasibility columns are unavailable.", call. = FALSE)
    }
    
    feasibility_data <- rbind(
      data.frame(
        n_total = ns,
        status = "Pointwise feasible",
        value = as.numeric(sr$feasible_pointwise)
      ),
      data.frame(
        n_total = ns,
        status = "Sustained feasible",
        value = as.numeric(sr$feasible)
      )
    )
    
    return(
      ggplot2::ggplot(
        feasibility_data,
        ggplot2::aes(
          x = .data$n_total,
          y = .data$value,
          colour = .data$status
        )
      ) +
        ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
        ggplot2::scale_colour_manual(
          values = c(
            "Pointwise feasible" = "#0072B2",
            "Sustained feasible" = "#E69F00"
          )
        ) +
        ggplot2::scale_y_continuous(
          breaks = c(0, 1),
          limits = c(-0.05, 1.05)
        ) +
        ggplot2::labs(
          x = "Total sample size",
          y = "Feasibility",
          colour = NULL
        ) +
        theme_common() +
        ggplot2::theme(legend.position = "top")
    )
  }
  
  # ------------------------------------------------------------------
  # Operating characteristics.
  # ------------------------------------------------------------------
  
  pce_name <- if (identical(inp$test, "BF+-")) {
    "PCE(H-)"
  } else {
    "PCE(H0)"
  }
  
  metric_columns <- c(
    "Bayesian power" = "power",
    "Bayesian type-I error" = "type1",
    "Frequentist power" = "freq_power",
    "Frequentist type-I error" = "freq_type1"
  )
  metric_columns[pce_name] <- "ce_h0"
  metric_columns <- metric_columns[
    vapply(metric_columns, available, logical(1))
  ]
  metric_order <- names(metric_columns)
  
  metric_colours <- c(
    "Bayesian power" = "#000000",
    "Bayesian type-I error" = "#0072B2",
    "Frequentist power" = "#E69F00",
    "Frequentist type-I error" = "#CC79A7"
  )
  metric_colours[pce_name] <- "#009E73"
  
  legend_labels <- function(breaks) {
    as.expression(lapply(breaks, function(z) {
      if (identical(z, "PCE(H0)")) {
        quote(PCE(H[0]))
      } else if (identical(z, "PCE(H-)")) {
        quote(PCE(H["-"]))
      } else {
        z
      }
    }))
  }
  
  oc_data <- do.call(
    rbind,
    lapply(metric_order, function(metric) {
      data.frame(
        n_total = ns,
        metric = metric,
        value = sr[[metric_columns[[metric]]]],
        stringsAsFactors = FALSE
      )
    })
  )
  oc_data$metric <- factor(oc_data$metric, levels = metric_order)
  
  p_oc <- ggplot2::ggplot(
    oc_data,
    ggplot2::aes(
      x = .data$n_total,
      y = .data$value,
      colour = .data$metric
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.95,
      linetype = "solid",
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = metric_colours,
      breaks = metric_order,
      labels = legend_labels
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = 3, byrow = TRUE)
    ) +
    ggplot2::labs(
      title = "Operating characteristic",
      x = "Total sample size",
      y = "Probability",
      colour = NULL
    ) +
    theme_common() +
    ggplot2::theme(legend.position = "top")
  
  if (type == "oc") {
    return(p_oc)
  }
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.", call. = FALSE)
  }
  
  if (
    length(inp$test) != 1L ||
    !inp$test %in% c("BF01", "BF+0", "BF-0", "BF+-")
  ) {
    stop("Unsupported or missing Bayes factor test.", call. = FALSE)
  }
  
  # ------------------------------------------------------------------
  # Target-attainment annotations.
  # ------------------------------------------------------------------
  
  sustain_n <- inp$sustain_n
  if (
    !is.numeric(sustain_n) ||
    length(sustain_n) != 1L ||
    !is.finite(sustain_n)
  ) {
    sustain_n <- 0L
  }
  sustain_n <- as.integer(max(0, sustain_n))
  
  first_sustained <- function(ok) {
    for (i in seq_along(ok)) {
      j <- min(length(ok), i + sustain_n)
      window <- ok[i:j]
      if (all(!is.na(window)) && all(window)) {
        return(i)
      }
    }
    NA_integer_
  }
  
  target_values <- list(
    "Bayesian power" = inp$target_power,
    "Bayesian type-I error" = inp$target_type1,
    "Frequentist power" = inp$target_freq_power,
    "Frequentist type-I error" = inp$target_freq_type1
  )
  target_values[pce_name] <- list(inp$target_ce_h0)
  
  target_rows <- lapply(metric_order, function(metric) {
    target <- target_values[[metric]]
    
    if (!valid_target(target)) {
      return(NULL)
    }
    if (metric == pce_name && target <= 0) {
      return(NULL)
    }
    
    values <- sr[[metric_columns[[metric]]]]
    upper_bound <- grepl("type-I error", metric, fixed = TRUE)
    
    ok <- if (upper_bound) {
      is.finite(values) & values <= target
    } else {
      is.finite(values) & values >= target
    }
    
    idx <- first_sustained(ok)
    
    data.frame(
      metric = metric,
      target = target,
      n = if (is.na(idx)) NA_real_ else ns[idx],
      stringsAsFactors = FALSE
    )
  })
  
  target_rows <- Filter(Negate(is.null), target_rows)
  
  if (length(target_rows)) {
    targets <- do.call(rbind, target_rows)
    
    p_oc <- p_oc +
      ggplot2::geom_hline(
        data = targets,
        ggplot2::aes(
          yintercept = .data$target,
          colour = .data$metric
        ),
        linetype = "dashed",
        linewidth = 0.55,
        alpha = 0.65,
        show.legend = FALSE
      )
    
    targets$label <- vapply(
      targets$n,
      function(z) {
        if (is.finite(z)) {
          sprintf("n = %d", as.integer(z))
        } else {
          "n = NA"
        }
      },
      character(1)
    )
    targets$x <- max(ns) - 0.02 * diff(range(ns))
    targets$y <- if (nrow(targets) == 1L) {
      0.91
    } else {
      seq(0.91, 0.39, length.out = nrow(targets))
    }
    
    p_oc <- p_oc +
      ggplot2::geom_text(
        data = targets,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          label = .data$label,
          colour = .data$metric
        ),
        inherit.aes = FALSE,
        hjust = 1,
        size = 3.8,
        show.legend = FALSE
      )
  }
  
  selected_n <- NA_real_
  
  if ("feasible" %in% names(sr)) {
    feasible_idx <- which(
      !is.na(sr$feasible) & as.logical(sr$feasible)
    )
    if (length(feasible_idx)) {
      selected_n <- ns[feasible_idx[1L]]
    }
  }
  
  if (is.finite(selected_n)) {
    p_oc <- p_oc +
      ggplot2::geom_vline(
        xintercept = selected_n,
        linetype = "dotted",
        colour = "#595959",
        linewidth = 1
      ) +
      ggplot2::annotate(
        "text",
        x = selected_n,
        y = 0.98,
        label = sprintf(
          "Calibrated n = %d",
          as.integer(selected_n)
        ),
        angle = 90,
        hjust = 1,
        vjust = -0.5,
        size = 3.3,
        colour = "#595959"
      )
  }
  
  # ------------------------------------------------------------------
  # Prior data.
  # ------------------------------------------------------------------
  
  xseq <- seq(0.0001, 0.9999, length.out = 1000L)
  
  get_shape <- function(name) {
    value <- inp[[name]]
    if (
      !is.numeric(value) ||
      length(value) != 1L ||
      !is.finite(value) ||
      value <= 0
    ) {
      stop(
        "Missing or invalid prior parameter in x$inputs: ",
        name,
        call. = FALSE
      )
    }
    value
  }
  
  make_basis_data <- function(minus = FALSE) {
    suffix <- if (minus) "_Hminus" else ""
    rows <- list()
    counter <- 0L
    
    for (arm in 1:2) {
      for (role in c("Design", "Analysis")) {
        tag <- if (role == "Design") "d" else "a"
        a <- get_shape(paste0("a_", arm, "_", tag, suffix))
        b <- get_shape(paste0("b_", arm, "_", tag, suffix))
        
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          p = xseq,
          density = stats::dbeta(xseq, a, b),
          prior_type = role,
          parameter = paste0("p[", arm, "]"),
          stringsAsFactors = FALSE
        )
      }
    }
    
    result <- do.call(rbind, rows)
    result$prior_type <- factor(
      result$prior_type,
      levels = c("Design", "Analysis")
    )
    result
  }
  
  prior_theme <- function() {
    theme_common() +
      ggplot2::theme(
        aspect.ratio = NULL,
        legend.position = "bottom",
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(0.1, "cm"),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          size = 12,
          face = "plain"
        ),
        panel.spacing.x = grid::unit(0.7, "cm"),
        axis.title.y = ggplot2::element_text(
          margin = ggplot2::margin(r = 8)
        ),
        plot.margin = ggplot2::margin(3, 12, 3, 8)
      )
  }
  
  make_basis_plot <- function(data, y_label) {
    ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data$p,
        y = .data$density,
        linetype = .data$prior_type
      )
    ) +
      ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
      ggplot2::facet_grid(
        . ~ parameter,
        labeller = ggplot2::label_parsed
      ) +
      ggplot2::scale_linetype_manual(
        values = c("Design" = "solid", "Analysis" = "dotted"),
        breaks = c("Design", "Analysis"),
        labels = c("Design prior", "Analysis prior")
      ) +
      ggplot2::guides(
        linetype = ggplot2::guide_legend(nrow = 1)
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25)
      ) +
      ggplot2::labs(
        x = "Probability",
        y = y_label,
        linetype = NULL
      ) +
      prior_theme()
  }
  
  # ------------------------------------------------------------------
  # Shared prior title.
  # ------------------------------------------------------------------
  
  prior_title <- switch(
    inp$test,
    "BF01" = expression(
      "Design and analysis priors under " *
        H[0] * ": " *
        bgroup("", p[1] == p[2], "") *
        " and " *
        H[1] * ": " *
        bgroup("", p[1] != p[2], "")
    ),
    "BF+0" = expression(
      "Design and analysis priors under " *
        H[0] * ": " *
        bgroup("", p[1] == p[2], "") *
        " and " *
        H["+"] * ": " *
        bgroup("", p[2] > p[1], "")
    ),
    "BF-0" = expression(
      "Design and analysis priors under " *
        H[0] * ": " *
        bgroup("", p[1] == p[2], "") *
        " and " *
        H["-"] * ": " *
        bgroup("", p[1] > p[2], "")
    ),
    "BF+-" = expression(
      "Design and analysis priors under " *
        H["+"] * ": " *
        bgroup("", p[2] > p[1], "") *
        " and " *
        H["-"] * ": " *
        bgroup("", p[1] >= p[2], "")
    )
  )
  
  top_y <- switch(
    inp$test,
    "BF01" = expression("Density under " * H[1]),
    "BF+0" = expression("Density under " * H["+"]),
    "BF-0" = expression("Density under " * H["-"]),
    "BF+-" = expression("Density under " * H["+"])
  )
  
  p_top <- make_basis_plot(
    make_basis_data(minus = identical(inp$test, "BF-0")),
    y_label = top_y
  ) +
    ggplot2::labs(title = prior_title) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "plain",
        size = 14,
        margin = ggplot2::margin(b = 8)
      )
    )
  
  # ------------------------------------------------------------------
  # Second prior row.
  # ------------------------------------------------------------------
  
  if (identical(inp$test, "BF+-")) {
    p_middle <- make_basis_plot(
      make_basis_data(minus = TRUE),
      y_label = expression("Density under " * H["-"])
    )
  } else {
    if (is.null(h0_design)) {
      h0_design <- c(
        a = get_shape("a_0_d"),
        b = get_shape("b_0_d")
      )
    }
    
    if (is.null(h0_analysis)) {
      h0_analysis <- c(
        a = get_shape("a_0_a"),
        b = get_shape("b_0_a")
      )
    }
    
    validate_h0 <- function(z, argument) {
      if (
        !is.numeric(z) ||
        length(z) != 2L ||
        is.null(names(z)) ||
        !all(c("a", "b") %in% names(z)) ||
        any(!is.finite(z)) ||
        any(z <= 0)
      ) {
        stop(
          argument,
          " must be a positive named vector c(a = ..., b = ...).",
          call. = FALSE
        )
      }
      z[c("a", "b")]
    }
    
    h0_design <- validate_h0(h0_design, "h0_design")
    h0_analysis <- validate_h0(h0_analysis, "h0_analysis")
    
    null_data <- rbind(
      data.frame(
        p = xseq,
        density = stats::dbeta(
          xseq,
          h0_design[["a"]],
          h0_design[["b"]]
        ),
        prior_type = "Design prior"
      ),
      data.frame(
        p = xseq,
        density = stats::dbeta(
          xseq,
          h0_analysis[["a"]],
          h0_analysis[["b"]]
        ),
        prior_type = "Analysis prior"
      )
    )
    
    null_data$prior_type <- factor(
      null_data$prior_type,
      levels = c("Design prior", "Analysis prior")
    )
    
    p_middle <- ggplot2::ggplot(
      null_data,
      ggplot2::aes(x = .data$p, y = .data$density)
    ) +
      ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
      ggplot2::facet_grid(. ~ prior_type) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25)
      ) +
      ggplot2::labs(
        x = expression("Common response probability " * p),
        y = expression("Density under " * H[0])
      ) +
      prior_theme()
  }
  
  # Collect the prior legends below both prior rows.
  prior_block <- patchwork::wrap_plots(
    p_top,
    p_middle,
    ncol = 1,
    heights = c(1, 1),
    guides = "collect"
  )
  
  prior_block <- prior_block &
    ggplot2::theme(legend.position = "bottom")
  
  # Keep the OC legend separate and preserve its additional height.
  patchwork::wrap_plots(
    prior_block,
    p_oc,
    ncol = 1,
    heights = c(3.2, 3.0),
    guides = "keep"
  )
}