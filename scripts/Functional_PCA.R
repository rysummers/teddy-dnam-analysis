##############################################################################
# Functional PCA (fPCA) for DNA methylation trajectories
# Goal: reduce each CpG's fitted spline trajectory to a small number of
#       fPCA scores that summarize its SHAPE, for downstream clustering
#       (increasing / decreasing / transient patterns)
##############################################################################

# install.packages(c("fda", "refund"))  # run once if needed
library(fda)
library(refund)

##############################################################################
# STEP 1: Reconstruct each CpG's trajectory on a common, evenly-spaced grid
##############################################################################
# You do NOT need the original model objects for this -- your saved spline
# coefficients (fixed_intercept, ns1_est, ns2_est, ns3_est) fully define the
# curve. You DO need the exact knots/Boundary.knots used when the ns() basis
# was originally created -- make sure these are saved/hardcoded somewhere in
# your pipeline metadata, since they're required to reconstruct the basis
# identically.

# --- Example placeholders: replace with your actual saved values ---
ns_basis <- ns(pheno_complete$age_yrs[pheno_complete$rgName %in% sample_ids], 
               df = 3) 

bs_basis <- bs(pheno_complete$age_yrs[pheno_complete$rgName %in% sample_ids], 
   df = 3)

attr(bs_basis, "knots")
attr(bs_basis, "Boundary.knots")
attr(ns_basis, "knots")
attr(ns_basis, "Boundary.knots")

og_knots <- attr(ns_basis, "knots")  
og_boundary <- attr(ns_basis, "Boundary.knots")

# Evenly spaced, fine grid across the observed age range
# (evenly spaced matters here -- unlike your saved pred_age_* columns,
#  which were unevenly spaced and would bias distance/integration below)
fine_ages <- seq(og_boundary[1], og_boundary[2], length.out = 100)

# Rebuild the spline basis ONCE (same basis is reused for every CpG)
ns_basis_fine <- ns(fine_ages, knots = og_knots,
                    Boundary.knots = og_boundary)

# Build a matrix: rows = CpGs, columns = fine_ages (100 evenly spaced points)
# This is the "dense functional data" input to fPCA
traj_matrix <- with(spline_summary_BLUPs,
                    outer(fixed_intercept, rep(1, length(fine_ages))) +
                      outer(ns1_est, ns_basis_fine[, 1]) +
                      outer(ns2_est, ns_basis_fine[, 2]) +
                      outer(ns3_est, ns_basis_fine[, 3]))
rownames(traj_matrix) <- spline_summary_BLUPs$CpG
# traj_matrix: n_CpGs x 100 matrix of predicted methylation (M or beta) values

##############################################################################
# OPTION A: fPCA using the `refund` package (fpca.face) -- fast, recommended
#           for dense, evenly-spaced grids like this one
##############################################################################

fpca_fit <- fpca.face(
  Y = traj_matrix, 
  argvals = fine_ages, 
  #pve = 0.95,
  npc = 4
  )
# pve = 0.95 -> keep enough components to explain 95% of variance in shapes

# Variance explained by each functional PC ("eigenfunction")
fpca_fit$evalues / sum(fpca_fit$evalues)

# The fPCA scores -- one row per CpG, one column per retained component
# These are your new, LOW-DIMENSIONAL, UNCORRELATED shape features
fpca_scores <- as.matrix(fpca_fit$scores)
colnames(fpca_scores) <- paste0("fPC", seq_len(ncol(fpca_scores)))
rownames(fpca_scores) <- rownames(traj_matrix)

head(fpca_scores)

sum(fpca_fit$evalues[1:4]) / sum(fpca_fit$evalues)

##############################################################################
# STEP 2: Interpret the eigenfunctions (what shape does each fPC represent?)
##############################################################################
# This is the functional-data equivalent of looking at PCA "loadings" --
# plot each eigenfunction to see what pattern of variation it captures
# (e.g., fPC1 might capture "overall increasing vs decreasing",
#  fPC2 might capture "early rise then plateau" vs "late rise", etc.)

par(mfrow = c(2, 2))
for (k in 1:ncol(fpca_fit$efunctions)) {
  plot(fine_ages, fpca_fit$efunctions[, k], type = "l", lwd = 2,
       main = paste0("Eigenfunction ", k,
                     " (", round(100 * fpca_fit$evalues[k] /
                                   sum(fpca_fit$evalues), 4), "% var)"),
       xlab = "Age", ylab = "Loading")
  abline(h = 0, lty = 2, col = "gray50")
}
par(mfrow = c(1, 1))

# You can also visualize what "high" vs "low" scores on fPC1 look like,
# by adding/subtracting the mean curve +/- eigenfunction * sd(score):
mean_curve <- fpca_fit$mu
pc1_sd <- sd(fpca_scores[, "fPC1"])

plot(fine_ages, mean_curve, type = "l", lwd = 2, ylim = range(traj_matrix),
     main = "Mean trajectory +/- 2 SD along fPC1", xlab = "Age", ylab = "Methylation")
lines(fine_ages, mean_curve + 2 * pc1_sd * fpca_fit$efunctions[, 1],
      col = "firebrick", lwd = 2, lty = 2)
lines(fine_ages, mean_curve - 2 * pc1_sd * fpca_fit$efunctions[, 1],
      col = "steelblue", lwd = 2, lty = 2)
legend("topleft", c("Mean", "+2SD fPC1", "-2SD fPC1"),
       col = c("black", "firebrick", "steelblue"), lty = c(1, 2, 2), lwd = 2)

##############################################################################
# STEP 3: Cluster CpGs on the fPCA scores (not the raw 100-point trajectories)
##############################################################################

set.seed(123)

# Decide k -- e.g., via elbow method on within-cluster sum of squares
wss <- sapply(2:8, function(k) kmeans(fpca_scores, centers = k, nstart = 25)$tot.withinss)
plot(2:8, wss, type = "b", xlab = "Number of clusters (k)", ylab = "Within-cluster SS")

# Fit k-means with your chosen k (e.g., k = 4)
km_fit <- kmeans(fpca_scores, centers = 4, nstart = 25)

spline_summary_BLUPs$fpca_cluster <- km_fit$cluster[
  match(spline_summary_BLUPs$CpG, rownames(fpca_scores))
]

##############################################################################
# STEP 4: Characterize each cluster's average trajectory shape
##############################################################################
# This is how you'll actually LABEL clusters as "increasing", "decreasing",
# "transient", etc. -- by inspecting the mean shape per cluster

library(ggplot2)

traj_df <- as.data.frame(traj_matrix)
traj_df$CpG <- rownames(traj_matrix)
traj_df$cluster <- spline_summary_BLUPs$fpca_cluster[
  match(traj_df$CpG, spline_summary_BLUPs$CpG)
]

traj_long <- reshape(traj_df, direction = "long",
                     varying = as.character(seq_along(fine_ages)),
                     v.names = "methylation", timevar = "age_idx")
traj_long$age <- fine_ages[traj_long$age_idx]

cluster_means <- aggregate(methylation ~ age + cluster, data = traj_long, FUN = mean)

ggplot(cluster_means, aes(x = age, y = methylation, color = factor(cluster))) +
  geom_line(linewidth = 1.2) +
  labs(color = "Cluster", title = "Mean trajectory shape per fPCA-derived cluster") +
  theme_minimal()

##############################################################################
# ALTERNATIVE: Classic `fda` package workflow (more control, more manual)
##############################################################################
# Useful if you want to fit fPCA on a B-spline basis representation
# directly (rather than the dense-grid approach above), or need finer
# control over smoothing/roughness penalties.

# 1. Represent each CpG's curve as a functional data (fd) object
bspl_basis <- create.bspline.basis(rangeval = range(fine_ages), nbasis = 15)
fd_obj <- smooth.basis(argvals = fine_ages, y = t(traj_matrix), fdParobj = bspl_basis)$fd

# 2. Run functional PCA
pca_fd_result <- pca.fd(fd_obj, nharm = 4)

# 3. Proportion of variance explained by each harmonic (eigenfunction)
pca_fd_result$varprop

# 4. Plot harmonics (equivalent to eigenfunction plots above)
plot(pca_fd_result)

# 5. Scores for clustering
fd_scores <- pca_fd_result$scores
rownames(fd_scores) <- rownames(traj_matrix)

##############################################################################
# NOTES / GOTCHAS
##############################################################################
# - fpca.face (Option A) is generally faster and simpler for a dense, evenly
#   spaced grid like this one -- recommended as your default.
# - The `fda` package (Option B) gives more explicit control (roughness
#   penalties, basis choice) but requires more manual setup.
# - Always cluster on fPCA SCORES, never on the raw 100-point matrix directly
#   -- the scores are uncorrelated and low-dimensional, avoiding the
#   redundancy/curse-of-dimensionality issues from raw dense trajectories.
# - Standardize/scale fPCA scores before k-means if the variance across
#   components differs a lot (kmeans() is sensitive to scale, same as with
#   any Euclidean-distance-based method).







##############################################################################
# STEP 1b: Filter flat CpGs, then row-center (and scale) each CpG's own
#          trajectory -- CRITICAL before running fPCA
##############################################################################
# fpca.face(center = TRUE) only subtracts the CROSS-SECTIONAL mean (the
# average curve across ALL CpGs at each age) -- it does NOT remove each
# CpG's own baseline methylation level. Since CpG baselines naturally range
# from near 0 to near 1 (beta scale), and true age-related change is only a
# few percentage points, leaving raw levels in means the eigen-decomposition
# gets dominated by "which CpGs sit at high vs low absolute methylation" --
# NOT by the shape of change you actually care about. This is almost
# certainly why Eigenfunction 1 captured 99.965% of variance while its
# y-axis range was nearly flat (-0.101 to -0.099).

# 1. Filter out flat/near-flat CpGs FIRST (using your effect-size or LRT
#    filter from earlier) -- this must happen before scaling, since dividing
#    a genuinely flat curve by its near-zero SD blows noise up to Inf/NaN
flat_threshold <- 0.01   # e.g., trajectory range < 1% beta-value change
traj_range <- apply(traj_matrix, 1, function(x) max(x) - min(x))
keep <- traj_range >= flat_threshold
traj_matrix_nonflat <- traj_matrix[keep, ]

# 2. Row-center: remove each CpG's own baseline level, leaving only ITS OWN
#    shape of change relative to itself
traj_matrix_centered <- traj_matrix_nonflat - rowMeans(traj_matrix_nonflat)

# 3. Row-scale (optional but recommended, given how small/variable the raw
#    magnitude of methylation change is across CpGs): standardizes each
#    curve to unit variance, so fPCA captures RELATIVE SHAPE rather than
#    being dominated by whichever CpGs happen to swing the most in absolute
#    terms. Skip this step if you specifically want magnitude-weighted shape.
traj_matrix_scaled <- traj_matrix_centered /
  apply(traj_matrix_centered, 1, sd)

##############################################################################
# OPTION A: fPCA using the `refund` package (fpca.face) -- fast, recommended
#           for dense, evenly-spaced grids like this one
##############################################################################

fpca_fit <- fpca.face(Y = traj_matrix_scaled, argvals = fine_ages, pve = 0.95)
# pve = 0.95 -> keep enough components to explain 95% of variance in shapes
# NOTE: now that baseline level and magnitude are removed, expect variance
# to spread more meaningfully across the first several components, rather
# than collapsing into one -- re-check the scree plot / evalues after this fix

# Variance explained by each functional PC ("eigenfunction")
fpca_fit$evalues / sum(fpca_fit$evalues)

# The fPCA scores -- one row per CpG, one column per retained component
# These are your new, LOW-DIMENSIONAL, UNCORRELATED shape features
# Force to matrix in case only 1 component is retained (drop = TRUE issue)
fpca_scores <- as.matrix(fpca_fit$scores)
rownames(fpca_scores) <- rownames(traj_matrix_scaled)   # note: filtered set
colnames(fpca_scores) <- paste0("fPC", seq_len(ncol(fpca_scores)))

head(fpca_scores)

##############################################################################
# STEP 2: Interpret the eigenfunctions (what shape does each fPC represent?)
##############################################################################
# This is the functional-data equivalent of looking at PCA "loadings" --
# plot each eigenfunction to see what pattern of variation it captures
# (e.g., fPC1 might capture "overall increasing vs decreasing",
#  fPC2 might capture "early rise then plateau" vs "late rise", etc.)

par(mfrow = c(1, 3))
for (k in 1:min(3, ncol(fpca_fit$efunctions))) {
  plot(fine_ages, fpca_fit$efunctions[, k], type = "l", lwd = 2,
       main = paste0("Eigenfunction ", k,
                     " (", round(100 * fpca_fit$evalues[k] /
                                   sum(fpca_fit$evalues), 1), "% var)"),
       xlab = "Age", ylab = "Loading")
  abline(h = 0, lty = 2, col = "gray50")
}
par(mfrow = c(1, 1))

# You can also visualize what "high" vs "low" scores on fPC1 look like,
# by adding/subtracting the mean curve +/- eigenfunction * sd(score):
mean_curve <- fpca_fit$mu
pc1_sd <- sd(fpca_scores[, "fPC1"])

plot(fine_ages, mean_curve, type = "l", lwd = 2, ylim = range(traj_matrix_scaled),
     main = "Mean CENTERED/SCALED trajectory +/- 2 SD along fPC1 (not raw methylation)",
     xlab = "Age", ylab = "Standardized methylation (centered/scaled)")
lines(fine_ages, mean_curve + 2 * pc1_sd * fpca_fit$efunctions[, 1],
      col = "firebrick", lwd = 2, lty = 2)
lines(fine_ages, mean_curve - 2 * pc1_sd * fpca_fit$efunctions[, 1],
      col = "steelblue", lwd = 2, lty = 2)
legend("topleft", c("Mean", "+2SD fPC1", "-2SD fPC1"),
       col = c("black", "firebrick", "steelblue"), lty = c(1, 2, 2), lwd = 2)

##############################################################################
# STEP 3: Cluster CpGs on the fPCA scores (not the raw 100-point trajectories)
##############################################################################

set.seed(123)

# Decide k -- e.g., via elbow method on within-cluster sum of squares
wss <- sapply(2:8, function(k) kmeans(fpca_scores, centers = k, nstart = 25)$tot.withinss)
plot(2:8, wss, type = "b", xlab = "Number of clusters (k)", ylab = "Within-cluster SS")

# Fit k-means with your chosen k (e.g., k = 4)
km_fit <- kmeans(fpca_scores, centers = 4, nstart = 25)

spline_summary_BLUPs$fpca_cluster <- km_fit$cluster[
  match(spline_summary_BLUPs$CpG, rownames(fpca_scores))
]
# Note: CpGs filtered out as "flat" in Step 1b will have NA cluster here --
# that's expected and correct, they were excluded from the fPCA entirely

##############################################################################
# STEP 4: Characterize each cluster's average trajectory shape
##############################################################################
# This is how you'll actually LABEL clusters as "increasing", "decreasing",
# "transient", etc. -- by inspecting the mean shape per cluster
# IMPORTANT: plot using the RAW traj_matrix (actual predicted methylation),
# not traj_matrix_scaled, for this step -- standardized units aren't
# interpretable as "% methylation change"; clustering happened on scaled
# shape, but characterization/reporting should use real methylation values

library(ggplot2)

traj_df <- as.data.frame(traj_matrix[rownames(fpca_scores), ])  # non-flat CpGs only
traj_df$CpG <- rownames(traj_df)
traj_df$cluster <- spline_summary_BLUPs$fpca_cluster[
  match(traj_df$CpG, spline_summary_BLUPs$CpG)
]

traj_long <- reshape(traj_df, direction = "long",
                     varying = as.character(seq_along(fine_ages)),
                     v.names = "methylation", timevar = "age_idx")
traj_long$age <- fine_ages[traj_long$age_idx]

cluster_means <- aggregate(methylation ~ age + cluster, data = traj_long, FUN = mean)

ggplot(cluster_means, aes(x = age, y = methylation, color = factor(cluster))) +
  geom_line(linewidth = 1.2) +
  labs(color = "Cluster", title = "Mean trajectory shape per fPCA-derived cluster") +
  theme_minimal()

##############################################################################
# ALTERNATIVE: Classic `fda` package workflow (more control, more manual)
##############################################################################
# Useful if you want to fit fPCA on a B-spline basis representation
# directly (rather than the dense-grid approach above), or need finer
# control over smoothing/roughness penalties.

# 1. Represent each CpG's curve as a functional data (fd) object
bspl_basis <- create.bspline.basis(rangeval = range(fine_ages), nbasis = 15)
fd_obj <- smooth.basis(argvals = fine_ages, y = t(traj_matrix), fdParobj = bspl_basis)$fd

# 2. Run functional PCA
pca_fd_result <- pca.fd(fd_obj, nharm = 4)

# 3. Proportion of variance explained by each harmonic (eigenfunction)
pca_fd_result$varprop

# 4. Plot harmonics (equivalent to eigenfunction plots above)
plot(pca_fd_result)

# 5. Scores for clustering
fd_scores <- pca_fd_result$scores
rownames(fd_scores) <- rownames(traj_matrix)

##############################################################################
# NOTES / GOTCHAS
##############################################################################
# - fpca.face (Option A) is generally faster and simpler for a dense, evenly
#   spaced grid like this one -- recommended as your default.
# - The `fda` package (Option B) gives more explicit control (roughness
#   penalties, basis choice) but requires more manual setup.
# - Always cluster on fPCA SCORES, never on the raw 100-point matrix directly
#   -- the scores are uncorrelated and low-dimensional, avoiding the
#   redundancy/curse-of-dimensionality issues from raw dense trajectories.
# - Standardize/scale fPCA scores before k-means if the variance across
#   components differs a lot (kmeans() is sensitive to scale, same as with
#   any Euclidean-distance-based method).

## after running fPCA (e.g., via fda::pca.fd or refund::fpca.sc)
mean_curve <- eval.fd(fine_ages, mean.fd(fd_obj))

for (pc in 1:3) {
  harm <- eval.fd(fine_ages, pca_fit$harmonics[pc])
  scale_factor <- sqrt(pca_fit$values[pc]) * 2  # ~2 SD in score units
  
  plot(fine_ages, mean_curve, type = "l", lwd = 2,
       main = paste("PC", pc), ylim = range(mean_curve) + c(-1,1)*max(abs(harm))*scale_factor)
  lines(fine_ages, mean_curve + scale_factor * harm, col = "blue", lty = 2)
  lines(fine_ages, mean_curve - scale_factor * harm, col = "red", lty = 2)
}
