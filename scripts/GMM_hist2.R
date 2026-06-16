#!/usr/bin/env Rscript

# =============================================================================
# plot_gmm_hist.R
#
# Visualizes a Gaussian Mixture Model (GMM) fit for a single CpG site by
# overlaying a density-filled histogram with per-component Gaussian curves,
# crossover split lines, and a legend summarizing component means and mixing
# proportions.
#
# Dependencies: mclust
#
# Usage:
#   source("plot_gmm_hist.R")
#   plot_gmm_hist("cg14531093", spline_intercepts, int_cols)
# =============================================================================


#' Plot a GMM-annotated histogram for a single CpG site
#'
#' Fits a Gaussian Mixture Model (GMM) via \code{Mclust} to the per-subject
#' spline intercepts (BLUPs) for a given CpG, then plots:
#'   - A histogram scaled to density
#'   - Solid colored lines and fills for each model-estimated Gaussian component
#'     (i.e. props[k] * dnorm(x, mean[k], sd[k])); this is what the GMM fitted
#'   - A single dashed black line for the empirical kernel density estimate,
#'     shown as a reference only
#'   - Dotted vertical lines marking the crossover split points between adjacent
#'     components (where adjacent weighted Gaussians intersect)
#'   - A legend reporting each component's mean (mu) and mixing proportion (pi)
#'
#' The best number of components G is selected automatically by BIC over the
#' range specified in \code{G_range}. The plot title reports the selected G.
#'
#' Variance handling: mclust may fit an equal-variance model (one shared
#' sigma^2 across components, model type "E") or an unequal-variance model
#' (per-component sigma^2, model type "V"). Both cases are handled
#' automatically.
#'
#' @param cpg Character string. The CpG identifier to plot (must be present in
#'   the \code{CpG} column of \code{spline_intercepts}).
#'
#' @param spline_intercepts data.table. Wide-format table of per-subject spline
#'   intercepts (BLUPs). Must contain:
#'   \describe{
#'     \item{CpG}{Character column of CpG identifiers (one row per CpG).}
#'     \item{intercept_* columns}{Numeric columns holding per-subject intercept
#'       values; these are identified via \code{int_cols}.}
#'   }
#'
#' @param int_cols Character vector. Names of the intercept columns in
#'   \code{spline_intercepts} to use as the data vector for GMM fitting.
#'   Typically produced by:
#'   \code{grep("^intercept_", names(spline_intercepts), value = TRUE)}
#'
#' @param G_range Integer vector. The candidate numbers of Gaussian components
#'   to evaluate. BIC is used to select the best G within this range.
#'   Default: \code{1:3}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   producing a base R plot.
#'
#' @examples
#' # Plot a single CpG, allowing up to 3 components
#' plot_gmm_hist("cg14531093", spline_intercepts, int_cols)
#'
#' # Force evaluation over 2-3 components only
#' plot_gmm_hist("cg14531093", spline_intercepts, int_cols, G_range = 2:3)

plot_gmm_hist <- function(
    cpg, 
    spline_intercepts, 
    int_cols, 
    G_range = 1:3,
    legend_pos = "topright",
    show_empirical = TRUE
) {
  
  # ---- 1. Extract and clean data for the requested CpG --------------------
  
  x <- as.numeric(spline_intercepts[CpG == cpg, ..int_cols])
  x <- x[is.finite(x)]   # drop NA / Inf / -Inf
  
  # ---- 2. Fit GMM via BIC-selected number of components -------------------
  
  fit    <- Mclust(x, G = G_range, verbose = FALSE)
  best_G <- fit$G   # number of components chosen by BIC
  
  # Sort all component parameters by ascending mean so colors are applied
  # consistently left-to-right across the x-axis
  ord   <- order(fit$parameters$mean)
  means <- fit$parameters$mean[ord]
  props <- fit$parameters$pro[ord]
  
  # Handle equal-variance (scalar sigmasq, model "E") vs
  # unequal-variance (vector sigmasq, model "V") mclust parameterizations
  sigmasq <- fit$parameters$variance$sigmasq
  if (length(sigmasq) == 1) {
    sds <- rep(sqrt(sigmasq), best_G)   # shared SD across all components
  } else {
    sds <- sqrt(sigmasq[ord])           # component-specific SDs
  }
  
  # ---- 3. Find crossover points between adjacent components ---------------
  #
  # The crossover is the x value where two adjacent weighted Gaussian densities
  # are equal: p_k * N(x | mu_k, sd_k) = p_{k+1} * N(x | mu_{k+1}, sd_{k+1})
  # Solved numerically with uniroot() over the interval [mu_k, mu_{k+1}].
  
  find_crossover <- function(m1, s1, p1, m2, s2, p2) {
    f <- function(x) p1 * dnorm(x, m1, s1) - p2 * dnorm(x, m2, s2)
    uniroot(f, interval = c(m1, m2))$root
  }
  
  splits <- numeric(best_G - 1)   # one split per adjacent component pair
  for (k in seq_len(best_G - 1)) {
    splits[k] <- find_crossover(
      means[k],   sds[k],   props[k],
      means[k+1], sds[k+1], props[k+1])
  }
  
  # ---- 4. Define per-component colors (teal, coral, purple) ---------------
  
  fill_palette   <- c("#1D9E75", "#D85A30", "#7F77DD")
  border_palette <- c("#0F6E56", "#993C1D", "#534AB7")
  
  fills   <- adjustcolor(fill_palette[seq_len(best_G)],  alpha.f = 0.3)
  borders <- border_palette[seq_len(best_G)]
  
  # ---- 5. Base plot: histogram on density scale ---------------------------
  
  d      <- density(x)
  breaks <- seq(min(x), max(x), length.out = 31)
  
  hist(x,
       breaks = breaks,
       freq   = FALSE,
       main   = paste0(cpg, "  [G = ", best_G, "]"),
       xlab   = "Intercepts (BLUP)",
       ylab   = "Density",
       col    = adjustcolor("#D3D1C7", alpha.f = 0.55),
       border = "#888780",
       las    = 1)
  
  # ---- 6. Fill and outline each model-estimated Gaussian component ---------
  #
  # Each component curve is the weighted Gaussian: props[k] * dnorm(x, mean, sd)
  # This is what the GMM actually fitted; solid lines = model estimate.
  # Filled polygons are drawn first so outlines render on top.
  
  xseq <- seq(min(x), max(x), length.out = 300)
  
  for (k in seq_len(best_G)) {
    y_comp <- props[k] * dnorm(xseq, means[k], sds[k])
    
    # Fill under the model Gaussian curve
    polygon(c(xseq[1], xseq, xseq[length(xseq)]),
            c(0,       y_comp, 0),
            col    = fills[k],
            border = NA)
    
    # Solid colored line = model-estimated component curve
    lines(xseq, y_comp, col = borders[k], lwd = 1.5)
  }
  
  # ---- 7. Empirical kernel density as dashed reference --------------------
  #
  # Drawn after the fills so it sits on top and is clearly visible.
  # Dashed black line = empirical density (data only, no model).
  # ---- 7. Empirical kernel density as dashed reference --------------------
  if (show_empirical) {
    lines(d, col = "#2C2C2A", lwd = 1.4, lty = 2)
  }
  
  # ---- 8. Mark crossover split points -------------------------------------
  
  abline(v = splits, col = "#444441", lwd = 1.2, lty = 3)
  
  # ---- 9. Legend with component summary statistics ------------------------
  #
  # Colored filled boxes = GMM model components (solid line = model estimate)
  # Dashed line entry    = empirical kernel density (reference only)
  
  comp_labels <- paste0("G", seq_len(best_G),
                        "  \u03bc=", round(means, 2),
                        "  \u03c0=", round(props, 2))
  
  emp_label <- if (show_empirical) "Empirical density" else NULL
  
  legend(legend_pos,
         legend = c(comp_labels, emp_label),
         fill   = c(fills,   if (show_empirical) NA),
         border = c(borders, if (show_empirical) NA),
         lty    = c(rep(NA, best_G), if (show_empirical) 2),
         lwd    = c(rep(NA, best_G), if (show_empirical) 1.4),
         col    = c(rep(NA, best_G), if (show_empirical) "#2C2C2A"),
         seg.len = 1.5,
         bty    = "n",
         cex    = 0.85)
  
  invisible(NULL)
}
