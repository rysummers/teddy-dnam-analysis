# install.packages("/Users/ryan_summers/GitHub/Teddy-DNAm-Analysis/R/clv_0.3-2.5.tar.gz", 
#                  repos = NULL, type = "source")
# 
# install.packages("/Users/ryan_summers/GitHub/Teddy-DNAm-Analysis/R/longitudinalData_2.4.7.tar.gz", 
#                  repos = NULL, type = "source")
# 
# install.packages("/Users/ryan_summers/GitHub/Teddy-DNAm-Analysis/R/kml_2.5.0.tar.gz", 
#                  repos = NULL, type = "source")

##### Load Libraries #####
options(rgl.useNULL = TRUE)
library("kml")

##### Load Data #####
data("epipageShort", package = "kml")
head(epipageShort)


##### Impute for Missing Values #####
imputation(as.matrix(epipageShort[, 3:6]))
imputation(as.matrix(epipageShort[, 3:6]), method = "linearInterpol")

## Trajectories of sdq are in columns 3 to 6 of the data.frame. So we can build the object of
## class ‘cld’:
cldSDQ <- cld(epipageShort, timeInData = 3:6)
cldSDQ


## Partitioning the data is then performed using the function kml(). At first, we
## want to “see” the clustering process. Since it might be slow, we will ask only
## for two redrawings with the graph display.
kml(cldSDQ, nbRedraw = 2, toPlot = "both")


## Then we want more rerolling, without the graphical display (that slow down the process):
kml(cldSDQ)

















