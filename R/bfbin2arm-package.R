#' bfbin2arm: Bayesian Power for Two-Arm Binomial Bayes Factors
#'
#' Functions for Bayesian power and sample size calculation in two-arm binomial
#' trials using Bayes factors for point-null and directional hypotheses.
#'
#' The package provides simulation-free and simulation-based methods to
#' calibrate Bayes factor designs in two-arm phase II trials with binary
#' endpoints, including power and type-I-error calculations and visualization
#' tools.
#'
#' @docType package
#' @name bfbin2arm
#'
#' @keywords internal
#'
#' @importFrom stats integrate dbeta pbeta dbinom density
#' @importFrom utils flush.console globalVariables
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline annotate labs
#'   theme_minimal theme element_text scale_color_manual coord_cartesian
#'   facet_grid scale_linetype_manual guides guide_legend unit element_blank
#'   labeller label_parsed
#' @importFrom dplyr `%>%`
#' @importFrom patchwork plot_layout
#' @importFrom graphics axis layout mtext plot.new points text title
"_PACKAGE"

utils::globalVariables(
  c("n", "power", "t1e", "pceH0", "freq_power",
    "x", "density", "prior_type")
)
