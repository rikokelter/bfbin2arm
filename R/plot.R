#' Plot an optimal two-stage two-arm Bayes factor design
#'
#' Given the result from \code{optimal_twostage_2arm_bf()}, this function
#' produces a six-panel base R plot showing the design schematic, operating
#' characteristics, and the design and analysis priors under \eqn{H_0} and
#' \eqn{H_1}.
#'
#' @param res A list returned by \code{optimal_twostage_2arm_bf()}, containing
#'   components \code{$design}, \code{$naive_oc}, \code{$occ} and \code{$priors}.
#' @param main Character string with the main title of the plot.
#'
#' @return Invisibly returns \code{NULL}; called for its side effect of
#'   producing a plot.
#'
#' @export
#'
#' @examples
#' res <- optimal_twostage_2arm_bf(
#'   alpha = 0.10, beta = 0.20, k = 1/3, k_f = 3,
#'   n1_min = c(3, 3), n2_max = c(8, 8),
#'   alloc1 = 0.5, alloc2 = 0.5,
#'   power_cushion = 0,
#'   interim_fraction = c(0.5, 0.5),
#'   grid_step = 1,
#'   progress = FALSE,
#'   max_iter = 16,
#'   test = "BF01",
#'   a_0_d = 1, b_0_d = 1,
#'   a_0_a = 1, b_0_a = 1,
#'   a_1_d = 1, b_1_d = 1,
#'   a_2_d = 1, b_2_d = 1,
#'   a_1_a = 1, b_1_a = 1,
#'   a_2_a = 1, b_2_a = 1
#' )
#' if (is.numeric(res$design) && length(res$design) == 4 && !anyNA(res$design)) {
#'   plot_twostage_2arm_bf(res)
#' }
plot_twostage_2arm_bf <- function(
    res,
    main = "Optimal two-stage two-arm Bayes factor design"
) {
  
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
  
  get_occ <- function(x, aliases) {
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
  
  text(5.0, 5.2, bquote(k == .(formatC(k, format = "f", digits = 3)) ~ "," ~
                          k[f] == .(formatC(k_f, format = "f", digits = 3))), cex = 0.95)
  if (!is.na(enh0_val)) {
    text(5.0, 4.1,
         bquote(E[H[0]](N) == .(formatC(enh0_val, format = "f", digits = 2))),
         cex = 0.95)
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
    d1 <- dbeta(p_grid, a_1_d, b_1_d)
    d2 <- dbeta(p_grid, a_2_d, b_2_d)
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