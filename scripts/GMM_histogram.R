

plot_gmm_hist <- function(cpg, spline_intercepts, int_cols, G_range = 1:3) {
  
  x <- as.numeric(spline_intercepts[CpG == cpg, ..int_cols])
  x <- x[is.finite(x)]
  
  fit <- Mclust(x, G = G_range, verbose = FALSE)
  best_G <- fit$G
  
  # extract and sort components by mean
  ord   <- order(fit$parameters$mean)
  means <- fit$parameters$mean[ord]
  props <- fit$parameters$pro[ord]
  
  # handle scalar (E model) vs vector (V model) sigmasq
  sigmasq <- fit$parameters$variance$sigmasq
  if (length(sigmasq) == 1) {
    sds <- rep(sqrt(sigmasq), best_G)
  } else {
    sds <- sqrt(sigmasq[ord])
  }
  
  # find crossover points between adjacent components
  find_crossover <- function(m1, s1, p1, m2, s2, p2) {
    f <- function(x) p1 * dnorm(x, m1, s1) - p2 * dnorm(x, m2, s2)
    uniroot(f, interval = c(m1, m2))$root
  }
  
  splits <- numeric(best_G - 1)
  for (k in seq_len(best_G - 1)) {
    splits[k] <- find_crossover(
      means[k], sds[k], props[k],
      means[k+1], sds[k+1], props[k+1])
  }
  
  # color ramps: teal, coral, purple
  fills   <- c("#1D9E75", "#D85A30", "#7F77DD")
  borders <- c("#0F6E56", "#993C1D", "#534AB7")
  fills   <- adjustcolor(fills[seq_len(best_G)], alpha.f = 0.40)
  borders <- borders[seq_len(best_G)]
  
  d <- density(x)
  breaks <- seq(min(x), max(x), length.out = 31)
  
  hist(x, breaks = breaks, freq = FALSE,
       main = paste0(cpg, "  [G = ", best_G, "]"),
       xlab = "Intercepts (BLUP)", ylab = "Density",
       col  = adjustcolor("#D3D1C7", alpha.f = 0.55),
       border = "#888780", las = 1)
  
  # filled density regions per component
  boundaries <- c(-Inf, splits, Inf)
  for (k in seq_len(best_G)) {
    idx <- d$x >= boundaries[k] & d$x <= boundaries[k+1]
    # extend to exactly hit boundary
    xi <- d$x[idx]; yi <- d$y[idx]
    polygon(c(xi[1], xi, xi[length(xi)]),
            c(0,     yi, 0),
            col = fills[k], border = NA)
    lines(xi, yi, col = borders[k], lwd = 2)
  }
  
  # overlay individual GMM component curves
  xseq <- seq(min(x), max(x), length.out = 300)
  for (k in seq_len(best_G)) {
    ycomp <- props[k] * dnorm(xseq, means[k], sds[k])
    lines(xseq, ycomp, col = borders[k], lwd = 1.2, lty = 2)
  }
  
  # vertical split lines
  abline(v = splits, col = "#444441", lwd = 1.2, lty = 3)
  
  legend("topright",
         legend = paste0("G", seq_len(best_G),
                         "  μ=", round(means, 2),
                         "  π=", round(props, 2)),
         fill   = fills,
         border = borders,
         bty    = "n", cex = 0.85)
}