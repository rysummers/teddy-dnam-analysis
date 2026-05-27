## Gaussian Mixture Model (GMM)


mc_fit <- Mclust(pc_use_test)
summary(mc_fit, parameters = TRUE)

mc_clusters <- mc_fit$classification
mc_fit$G
mc_fit$modelName


BIC <- mclustBIC(pc_use_test)
plot(BIC)



mod1 <- Mclust(pc_use_test, x=BIC)
summary(mod1, parameters = TRUE)

plot(mod1, what = "classification")
