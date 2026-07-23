##############################################################################
# Compare ns() vs bs() fitted trajectories for a single CpG
# Focus: boundary behavior in sparse late-age region
##############################################################################

library(ggplot2)

##############################################################################
# STEP 1: Build a prediction grid across the observed age range
##############################################################################
# Use a FINE, evenly spaced grid so boundary wiggle (if present) is visible
age_grid <- seq(min(dat$age_yrs), max(dat$age_yrs), length.out = 200)

##############################################################################
# STEP 2: Hold all other covariates at fixed, representative values
##############################################################################
# Continuous covariates -> mean; categorical -> most common level (mode)
# This isolates the age effect for visualization -- exactly what you want
# when comparing the SHAPE of the two spline bases

get_mode <- function(x) names(sort(table(x), decreasing = TRUE))[1]

newdat <- data.frame(
  age_yrs  = age_grid,
  gender   = get_mode(dat$gender),
  cc       = get_mode(dat$cc),
  new_Bcell = mean(dat$new_Bcell),
  new_CD4T  = mean(dat$new_CD4T),
  new_CD8T  = mean(dat$new_CD8T),
  new_Mono  = mean(dat$new_Mono),
  new_NK    = mean(dat$new_NK)
)

##############################################################################
# STEP 3: Get population-level predictions (exclude random effects)
##############################################################################
# re.form = NA -> predict the FIXED-EFFECTS-only trajectory (population
# average), which is what you want for comparing curve SHAPE between
# ns() and bs() -- including random effects would just add maskid-specific
# noise that has nothing to do with the spline basis choice

newdat$pred_ns <- predict(m1, newdata = newdat, re.form = NA)
newdat$pred_bs <- predict(m2, newdata = newdat, re.form = NA)

##############################################################################
# STEP 4 (optional but recommended): Add confidence bands via bootstrap
##############################################################################
# lme4 doesn't give analytic SEs for predict() by default; bootMer() is the
# standard way to get them. This can be slow -- reduce nsim for a quick look,
# increase (e.g., 500-1000) for a figure you'll actually report.

library(lme4)

pred_fun_ns <- function(fit) predict(fit, newdata = newdat, re.form = NA)
pred_fun_bs <- function(fit) predict(fit, newdata = newdat, re.form = NA)

boot_ns <- bootMer(m1, FUN = pred_fun_ns, nsim = 200, use.u = FALSE, seed = 123)
boot_bs <- bootMer(m2, FUN = pred_fun_bs, nsim = 200, use.u = FALSE, seed = 123)

newdat$lwr_ns <- apply(boot_ns$t, 2, quantile, probs = 0.025)
newdat$upr_ns <- apply(boot_ns$t, 2, quantile, probs = 0.975)
newdat$lwr_bs <- apply(boot_bs$t, 2, quantile, probs = 0.025)
newdat$upr_bs <- apply(boot_bs$t, 2, quantile, probs = 0.975)

##############################################################################
# STEP 5: Reshape to long format for a clean overlay plot
##############################################################################

plot_df <- data.frame(
  age_yrs = rep(newdat$age_yrs, 2),
  pred    = c(newdat$pred_ns, newdat$pred_bs),
  lwr     = c(newdat$lwr_ns, newdat$lwr_bs),
  upr     = c(newdat$upr_ns, newdat$upr_bs),
  basis   = rep(c("Natural spline (ns)", "B-spline (bs)"), each = nrow(newdat))
)

##############################################################################
# STEP 6: Plot -- overlay both curves, with a rug showing where data actually
#         exist (this is the key diagnostic for boundary behavior)
##############################################################################

ggplot(plot_df, aes(x = age_yrs, y = pred, color = basis, fill = basis)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_rug(data = dat, aes(x = age_yrs), inherit.aes = FALSE,
           sides = "b", alpha = 0.3) +
  labs(
    title = paste0("ns() vs bs() fitted trajectory: ", probe),
    subtitle = "Rug marks show observed age density -- watch boundary behavior where data are sparse",
    x = "Age (years)", y = "Predicted methylation (fixed effects only)",
    color = "Spline basis", fill = "Spline basis"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

##############################################################################
# NOTES
##############################################################################
# - The rug plot (bottom tick marks) is the key diagnostic here: look at
#   whether ns() and bs() start to visibly diverge in regions where rug
#   marks are sparse (e.g., your oldest ages, ~3.5-6 yrs based on your
#   earlier pred_age_* columns). If bs() swings/curves sharply while ns()
#   stays comparatively linear in that sparse region, that's direct visual
#   confirmation of the boundary-behavior concern.
# - re.form = NA gives the population-average curve. If you instead want to
#   see a specific child's trajectory (including their random intercept),
#   use re.form = NULL and provide that child's maskid in newdat.
# - If bootMer is too slow for exploration, first run steps 1-3 and 5-6
#   without the confidence bands (skip step 4) to get a quick visual;
#   add bootstrap bands only once you're ready to finalize a figure.
# - Repeat this same plot for a handful of other CpGs (ideally ones near
#   your significant/candidate list) to see if the pattern is consistent,
#   rather than judging basis choice off a single probe.