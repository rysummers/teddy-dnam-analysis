## FDApace

```{r}
library(fdapace)

tst_fpca <- traj_scaled %>%
  as.data.frame() %>%
  rownames_to_column("CpG") %>%
  pivot_longer(
    cols = starts_with("pred_age_"),
    names_to = "age_label",
    values_to = "M_value") %>%
  mutate(age = as.numeric(sub("pred_age_", "", age_label)))

T1D <- MakeFPCAInputs(tst_fpca$CpG, tst_fpca$age, tst_fpca$M_value)
fpcaObjT1D <- FPCA(T1D$Ly,T1D$Lt,
                   list(plot=TRUE,dataType = "Sparse"))
```

Estimate path plots: create the fitted sample path plot based on the results from the parse data.
```{r}
CreatePathPlot(fpcaObjT1D)
```


```{r}
SelectK(fpcaObjT1D,criterion='AIC')
predSparse <- predict(fpcaObjT1D, T1D$Ly, T1D$Lt, K=2) 
predSparse$scores
```

```{r}
set.seed(2522)

km <- kmeans(predSparse$predCurves, centers = 2, nstart = 50)
cpgs_to_test2$pred <- km$cluster[
  match(cpgs_to_test2$CpG, rownames(traj_scaled))]
```


```{r}
adjustedRandIndex(cpgs_to_test2$group, cpgs_to_test2$pred)
table(True = cpgs_to_test2$group, Pred = cpgs_to_test2$pred)
res <- data.frame(
  CpG = cpgs_to_test2$CpG,
  group = cpgs_to_test2$group,
  pred = pred)

res
```