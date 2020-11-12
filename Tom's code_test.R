

library(caret)
library(randomForest)
library(ranger)
library(mlr)
library(tibble)
library(raster)
library(sp)
library(rgdal)
library(hexbin)
library(lattice)
library(RColorBrewer)
library(viridis)
library(Metrics)

ksat_df<-read.csv("E:/Ksat_dataset_mapping.csv")

## Unique IDs for 1 degrees by 1 degrees

## Spatial ID:
pol.100km = readOGR("E:/Ksat/tiles_ll_100km_mask.shp")
sp.pnts = ksat_df[,c("longitude_decimal_degrees", "latitude_decimal_degrees")]
ov.ID = sp::over(SpatialPoints(sp.pnts, proj4string = CRS(proj4string(pol.100km))), pol.100km["ID"])
summary(is.na(ov.ID$ID))
ksat_df$ID = ov.ID$ID


source("E:/OpenLandMap/R/saveRDS_functions.R")
source("E:/OpenLandMap/R/LandGIS_functions.R")

## saveRDS_functions.R and LandGIS_functions.R available at https://github.com/Envirometrix/LandGISmaps/tree/477460d1d0099646c508f65e68769b9edf050ce8/functions

## 3D modeling (see Hengl, T., & MacMillan, R. A. (2019). Predictive soil mapping with R. Lulu. com.)

dfs <- hor2xyd(ksat_df, U="hzn_top", L="hzn_bot")


I.vars = make.names(unique(unlist(sapply(c("clm_", "dtm_", "lcv", "veg_","lat","long", "olm_c", "olm_s", "olm_bd", "ID","DEPTH"), function(i){names(dfs)[grep(i, names(dfs))]}))))

t.vars = c("log_ksat")
sel.n <- c(t.vars,I.vars)
sel.r <- complete.cases(dfs[,sel.n])
PTF_temp2 <- dfs[sel.r,sel.n]


unique(PTF_temp2$ID)



library(mlr)
# coordinates needed for the spatial partitioning
coords = PTF_temp2[, c("longitude_decimal_degrees", "latitude_decimal_degrees")]
# select response and predictors to use in the modeling
data = dplyr::select( PTF_temp2, -longitude_decimal_degrees, -latitude_decimal_degrees,-ID,-FID_Fish_n)


# create task
task = makeRegrTask(data = data, target = "log_ksat", coordinates = coords, blocking = as.factor( PTF_temp2$ID))

lrn.rf = mlr::makeLearner("regr.ranger", num.threads = parallel::detectCores(), mtry=6, num.trees=200, importance="impurity")

ctrl = makeFeatSelControlRandom(maxit = 2)

perf_level = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 2)

lrn1 = mlr::makeFeatSelWrapper(lrn.rf, resampling = perf_level, control = ctrl, show.info=TRUE)

result = mlr::resample(learner =lrn1, 
                       task = task,
                       resampling = perf_level,
                       measures = mlr::mse)

hh<-result$pred$data



