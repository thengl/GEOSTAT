## spatial predictions LoaLoa data set
## by: Tom.Hengl@isric.org

library(geostatsp)
library(GSIF)
library(plotKML)
library(rgdal)
library(quantregForest)
library(gstat)
library(randomForest)

## observations:
data(loaloa)
## covariates:
s <- list(elevationLoa, eviLoa, tempLoa, ltLoa)
r <- as(elevationLoa, "SpatialGridDataFrame")
grid2km <- spsample(type="regular", cellsize=2e3, r)
## TH: not all covariates are aligned to the same grid/extent
ov <- list(NULL)
for(i in 1:length(s)){
  if(!proj4string(s[[i]])==proj4string(elevationLoa)){
    s[[i]] <- projectRaster(s[[i]], crs=proj4string(elevationLoa), method="ngb")
  }
  ov[[i]] <- extract(s[[i]], grid2km)
}
grid2km <- SpatialPointsDataFrame(grid2km, data.frame(do.call(cbind, ov)))
gridded(grid2km) = TRUE
grid2km <- as(grid2km, "SpatialPixelsDataFrame")
names(grid2km) <- c("elevationLoa", "eviLoa", "tempLoa", "ltLoa")
## clean up / reformat covs:
grid2km$ltLoa <- as.factor(grid2km$ltLoa)
grid2km$tempLoa <- ifelse(is.nan(grid2km$tempLoa), NA, grid2km$tempLoa)
str(grid2km@data)
plot(raster(grid2km))
points(loaloa, pch="+")

## target variable:
loaloa$pre <- (loaloa$y/loaloa$N)*100
## randomForest-kriging:
m.pre <- fit.gstatModel(loaloa, pre~elevationLoa+eviLoa+tempLoa+ltLoa, grid2km, method="randomForest")
plot(m.pre)
varImpPlot(m.pre@regModel)
rk.pre <- predict(m.pre, grid2km, zmin=0, zmax=100)
plotKML(rk.pre, colour_scale=SAGA_pal[[1]])
