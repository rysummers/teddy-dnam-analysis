## Clustering


res <- as.data.table(lmeR_spline_results)

# keep non-singular fits only
res_ok <- res[fit_status == "ok"]

res[fit_status == "error"]

# prediction columns
traj_cols <- grep("^pred_age_", names(res_ok), value = TRUE)


#Build trajectory matrix

traj_mat <- as.matrix(res_ok[, ..traj_cols]) # .. => select
rownames(traj_mat) <- res_ok$CpG

dim(traj_mat)



#Center and scale by row for PCA

traj_shape <- scale_rows(traj_mat)
dim(traj_shape)


PCA

set.seed(2522)

pca <- prcomp(traj_shape, center = FALSE, scale. = FALSE)

pc_scores <- pca$x
rownames(pc_scores) <- rownames(traj_shape)

dim(pc_scores)
summary(pca)



# first 3 explain 100%
pc_use <- pc_scores[, 1:3, drop = FALSE]




# function to determine optimum number of k
run_kmeans_grid <- function(x, k_values = 3:12, nstart = 20, iter.max = 100) {
  out <- future_lapply(k_values, function(k) {
    message("Running kmeans for k = ", k)
    km <- kmeans(x, centers = k, nstart = nstart, iter.max = iter.max,
                 algorithm = "Hartigan-Wong")
    
    list(
      k = k,
      km = km,
      tot_withinss = km$tot.withinss,
      betweenss = km$betweenss,
      totss = km$totss,
      ratio = km$betweenss / km$totss)
  }, future.seed = TRUE)
  
  names(out) <- paste0("k", k_values)
  out
}




set.seed(2522)
future::plan(future::multisession, workers = 12)

t0 <- Sys.time()

km_grid <- run_kmeans_grid(
  x = pc_use,
  k_values = 3:12,
  nstart = 20,
  iter.max = 100)

t1 <- Sys.time()
message("lmeR runtime: ", 
        round(as.numeric(difftime(t1, t0, units = "mins")), 2), 
        " minutes")



km_summary <- rbindlist(lapply(km_grid, function(x) {
  data.table(
    k = x$k,
    tot_withinss = x$tot_withinss,
    betweenss = x$betweenss,
    totss = x$totss,
    ratio = x$ratio)
}))

print(km_summary)


#elbow plots

plot(km_summary$k, km_summary$tot_withinss, type = "b",
     xlab = "k", ylab = "Total within-cluster SS")

plot(km_summary$k, km_summary$ratio, type = "b",
     xlab = "k", ylab = "BetweenSS / TotalSS")


#Cluster quality checks


cluster_indices <- rbindlist(
  future_lapply(km_keep, function(z) {
    cl_sub <- z$km$cluster[sub_idx]
    
    crit <- clusterCrit::intCriteria(
      traj = as.matrix(pc_sub),
      part = as.integer(cl_sub),
      crit = c("Calinski_Harabasz", "Davies_Bouldin"))
    
    data.table(
      k = z$k,
      Calinski_Harabasz = crit$calinski_harabasz,
      Davies_Bouldin = crit$davies_bouldin)
  }, future.seed = TRUE),
  fill = TRUE)

print(cluster_indices)


#*	higher CH is better
#*	lower DB is better


set.seed(2522)
sub_idx <- sample(seq_len(nrow(pc_use)), size = min(50000, nrow(pc_use)))
pc_sub <- pc_use[sub_idx, , drop = FALSE]
cpg_sub <- rownames(pc_use)[sub_idx]


#Stability function using repeated k-means and Adjusted Rand Index (ARI)


kmeans_stability <- function(x, k, n_reps = 5, nstart = 20) {
  fits <- vector("list", n_reps)
  
  for (i in seq_len(n_reps)) {
    set.seed(100 + i)
    fits[[i]] <- kmeans(x, centers = k, nstart = nstart)
  }
  
  ari_vals <- c()
  for (i in 1:(n_reps - 1)) {
    for (j in (i + 1):n_reps) {
      ari_vals <- c(ari_vals,
                    mclust::adjustedRandIndex(fits[[i]]$cluster, fits[[j]]$cluster))
    }
  }
  
  data.table(
    k = k,
    mean_ARI = mean(ari_vals),
    median_ARI = median(ari_vals),
    min_ARI = min(ari_vals),
    max_ARI = max(ari_vals))
}


#Run stability for candidate k


stab_tbl <- rbindlist(lapply(3:12, function(k) {
  message("Stability check for k = ", k)
  kmeans_stability(pc_sub, k = k, n_reps = 5, nstart = 25)
}))

print(stab_tbl)


want k values with:
  *	decent CH
*	low DB
*	high ARI stability
*	interpretable trajectories

Choose optimum k and run final kmeans

set.seed(2522)
final_k <- 6

final_km <- kmeans(pc_use, centers = final_k, nstart = 50, iter.max = 100)

cluster_dt <- data.table(
  CpG = rownames(pc_use),
  cluster = final_km$cluster)

res_clustered <- merge(
  res_ok,
  cluster_dt,
  by = "CpG",
  all.x = TRUE)



res_clustered %>% 
  group_by(cluster) %>% 
  summarize(counts = n())



Inspect mean trajectory per cluster


age_vals <- as.numeric(sub("^pred_age_", "", traj_cols))


Cluster mean trajectories

cluster_means <- res_clustered[, lapply(.SD, mean, na.rm = TRUE),
                               by = cluster,
                               .SDcols = traj_cols]

print(cluster_means)


Plot cluster means

matplot(
  x = age_vals,
  y = m2beta(t(as.matrix(cluster_means[, ..traj_cols]))),
  type = "l", lty = 1, lwd = 2,
  xlab = "Age (years)",
  ylab = "Predicted methylation trajectory",
  main = "Mean trajectories by cluster")
legend("bottom", horiz = TRUE, cex=0.8,
       legend = paste("Cluster", cluster_means$cluster),
       col = seq_len(nrow(cluster_means)),
       lty = 1, lwd = 2, bty = "n")


clustering appears to be dominated by high and low baseline methylation CpG values instead
of trajectories


shape-only visualization:
  
cluster_means_shape <- copy(cluster_means)
cluster_means_shape[, (traj_cols) := 
                      as.data.table(scale_rows(
                        as.matrix(cluster_means[, ..traj_cols])))]

# sort clusters numerically
cluster_means_shape <- cluster_means_shape[order(cluster)]

matplot(
  x = age_vals,
  y = t(as.matrix(cluster_means_shape[, ..traj_cols])),
  type = "l", lty = 1, lwd = 2,
  xlab = "Age (years)",
  ylab = "Scaled trajectory",
  main = "Mean scaled trajectories by cluster")
legend("top", horiz = TRUE, cex=0.8,
       legend = paste("Cluster", cluster_means_shape$cluster),
       col = seq_len(nrow(cluster_means_shape)),
       lty = 1, lwd = 2, bty = "n")


Cluster stats

res_clustered2 <- merge(res_clustered, cmp[, .(CpG, delta_AIC, delta_BIC)],
                        by = "CpG", all.x = TRUE)

cluster_annotation <- res_clustered2[, .(
  n_cpg = .N,
  mean_aic = mean(aic, na.rm = TRUE),
  mean_delta_AIC_slope = mean(delta_AIC, na.rm = TRUE),
  pct_slope_help_AIC_gt2 = mean(delta_AIC > 2, na.rm = TRUE) * 100,
  mean_var_rand_slope = mean(var.rand.slope, na.rm = TRUE)), by = cluster]

cluster_annotation <- cluster_annotation[order(cluster)]
print(cluster_annotation)





Are clusters:
  *	mostly simple linear-ish
*	strongly nonlinear
*	more subject-heterogeneous


two groups : increasing and flat
simulate the data to test the method
20 each - test current clustering method
Test different variables to clustering to see what works best

Since we are controlling for cell proportions - residual epigenetic developmental 
are patterns not explained by cell mixture

Predicted methylation trajectories were evaluated on a pre-specified age grid 
concentrated in infancy and early childhood, where data support was strongest.
Because observations after age 4 years were sparse, those ages were excluded 
from the clustering grid to avoid unstable tail predictions driving cluster 
assignment. A secondary plotting grid included ages up to 6 years to reflect 
clinically relevant T1D screening windows.
