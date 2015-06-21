# title         : fishcamp.R
# purpose       : Pre-processing of elevation data, extraction of landform units and predictive mapping of soil mapping units
# reference     : Hengl, T. 2009. A Practical Guide to Geostatistical Mapping, 2nd Edt. University of Amsterdam, www.lulu.com, 291 p. ISBN: ISBN 978-90-9024981-0
# producer      : Prepared by T. Hengl
# inputs        : study area fishcamp (37.46353 N; 119.6119 W) comprises: (1) LIDAR points, (2) contour lines, (3) SRTM DEM; These can be obtained from [http://geomorphometry.org/content/fishcamp];
# outputs       : DEMs generated from points, stream networks, DEM-parameters, predicted soil mapping units and landforms;
# remarks 1     : local coordinates system used - UTM zone 11N with North American Datum 83;

library(maptools)
library(rgdal)
library(gstat)
library(spatstat)
library(mda) # kappa statistics
library(vcd)
library(plotKML)
library(raster)
library(RSAGA); rsaga.env()

# ------------------------------------------------------------
# Download of the maps:
# ------------------------------------------------------------

## Download the datasets:
if(!file.exists("fishcamp.zip")){
  download.file("http://geomorphometry.org/system/files/fishcamp.zip", destfile=paste(getwd(), "fishcamp.zip", sep="/"))
}
## Extract all data sets:
unzip("fishcamp.zip")
system("open README.txt")

## import LiDAR points:
lidar <- readOGR("lidar.shp", "lidar")
str(lidar@data) # 273,028 points (a very large point dataset);
## GRIDS:
grids25m <- readGDAL("DEMSRTM1.asc")
names(grids25m) <- "DEMSRTM1"
grids5m <- readGDAL("soilmu.asc")
names(grids5m) <- "soilmu"
grids5m$soilmu.c <- as.factor(grids5m$soilmu)
spplot(grids5m[1], col.regions=rainbow(8))
## estimate the pixel size (http://dx.doi.org/10.1016/j.cageo.2005.11.008):
pixelsize <- round(2*sqrt(areaSpatialGrid(grids25m)/length(lidar$Z)), 0)
pixelsize
## set-up the correct CRS:
proj4string(lidar) <- CRS("+init=epsg:26911")
proj4string(lidar)
proj4string(grids25m) <- CRS("+init=epsg:26911")
proj4string(grids5m) <- CRS("+init=epsg:26911")

## coordinates of the center:
grids25m.ll <- spTransform(grids25m, CRS("+proj=longlat +ellps=WGS84"))
grids25m.ll@bbox
clon <- mean(grids25m.ll@bbox[1,])
clat <- mean(grids25m.ll@bbox[2,])

# ------------------------------------------------------------
# DEM generation and variogram modelling:
# ------------------------------------------------------------

## variogram modelling

## This is a very large data set so we need to sub-sample:
lidar.sample <- lidar[runif(length(lidar$Z))<0.05,]
#varmap.plt <- plot(variogram(Z~1, lidar.sample, map=TRUE, cutoff=50*pixelsize, width=pixelsize), col.regions=grey(rev(seq(0,1,0.025))))
## obvious anisotropy
Z.svar <- variogram(Z~1, lidar.sample, alpha=c(45,135)) # cutoff=50*dem.pixelsize
Z.vgm <- fit.variogram(Z.svar, vgm(psill=var(lidar.sample$Z), "Gau", sqrt(areaSpatialGrid(grids25m))/4, nugget=0, anis=c(p=135, s=0.6)))
vgm.plt <- plot(Z.svar, Z.vgm, plot.nu=F, cex=2, pch="+", col="black")
Z.vgm
print(vgm.plt)
## plot the two variograms next to each other:
#print(varmap.plt, split=c(1,1,2,1), more=T)
#print(vgm.plt, split=c(2,1,2,1), more=F)

## Generate an initial DEM:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(TARGET_OUT_GRID="DEM5LIDAR.sgrd", INPUT="lidar.shp", FIELD=0, LINE_TYPE=0, TARGET_USER_SIZE=pixelsize, TARGET_USER_XMIN=grids5m@bbox[1,1]+pixelsize/2, TARGET_USER_XMAX=grids5m@bbox[1,2]-pixelsize/2, TARGET_USER_YMIN=grids5m@bbox[2,1]+pixelsize/2, TARGET_USER_YMAX=grids5m@bbox[2,2]-pixelsize/2), check.module.exists = FALSE)
## There are a quite some articifacts and missing pixels - we need to filter them
## first, run the residual analysis and mask out spikes:
rsaga.geoprocessor(lib="statistics_grid", 1, param=list(GRID="DEM5LIDAR.sgrd", RADIUS=5, DIFF="dif_lidar.sgrd"), check.module.exists = FALSE) # all features >7 pixels remain!
## read back to R and mask out all areas with high difference:
grids5m$DEM5LIDAR <- readGDAL("DEM5LIDAR.sdat")$band1
grids5m$dif <- readGDAL("dif_lidar.sdat")$band1
lim.dif <- quantile(grids5m$dif, c(0.025,0.975), na.rm=TRUE)
lim.dif
grids5m$DEM5LIDARf <- ifelse(grids5m$dif<=lim.dif[[1]]|grids5m$dif>=lim.dif[[2]], NA, grids5m$DEM5LIDAR)
summary(grids5m$DEM5LIDARf)[7]/length(grids5m@data[[1]])  # 15% pixels masked out!
writeGDAL(grids5m["DEM5LIDARf"], "DEM5LIDARf.sdat", "SAGA", mvFlag=-99999)

## Filter the missing values (close gaps):
rsaga.geoprocessor(lib="grid_tools", module=25, param=list(GRID="DEM5LIDARf.sgrd", CLOSED="DEM5LIDARf2.sgrd"))
grids5m$DEM5LIDARf <- readGDAL("DEM5LIDARf2.sdat")$band1

## Generate DEM from contours:
rsaga.geoprocessor(lib="grid_spline", module=1, param=list(TARGET_OUT_GRID="DEM25TPS.sgrd", SHAPES="contours.shp", TARGET_DEFINITION=0, SEARCH_POINTS_MAX=10, TARGET_USER_SIZE=25, TARGET_USER_XMIN=grids25m@bbox[1,1]+25/2, TARGET_USER_XMAX=grids25m@bbox[1,2]-25/2, TARGET_USER_YMIN=grids25m@bbox[2,1]+25/2, TARGET_USER_YMAX=grids25m@bbox[2,2]-25/2), check.module.exists = FALSE)

## resample LiDAR DEM to 25 m:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="DEM5LIDARf2.sgrd", TARGET_OUT_GRID="DEM25LIDAR.sgrd", SCALE_UP_METHOD=6, TARGET_USER_SIZE=25, TARGET_USER_XMIN=grids25m@bbox[1,1]+25/2, TARGET_USER_XMAX=grids25m@bbox[1,2]-25/2, TARGET_USER_YMIN=grids25m@bbox[2,1]+25/2, TARGET_USER_YMAX=grids25m@bbox[2,2]-25/2))
## read back to R:
grids25m$DEM25LIDAR <- readGDAL("DEM25LIDAR.sdat")$band1
grids25m$DEM25TPS <- readGDAL("DEM25TPS.sdat")$band1
## difference between the LiDAR DEM and topo DEM:
plot(grids25m$DEM25LIDAR, grids25m$DEM25TPS, asp=1)
sqrt(sum((grids25m$DEM25LIDAR - grids25m$DEM25TPS)^2)/length(grids25m$DEM25LIDAR))

# ------------------------------------------------------------
# Extraction of Land Surface Parameters:
# ------------------------------------------------------------

## Slope and curvature:
rsaga.geoprocessor(lib="ta_morphometry", module=0, param=list(ELEVATION="DEM5LIDARf2.sgrd", SLOPE="SLP.sgrd", C_TOTA="CRV.sgrd"))
## Topographic Wetness Index:
rsaga.geoprocessor(lib="ta_hydrology", module=15, param=list(DEM="DEM5LIDARf2.sgrd", TWI="TWI.sgrd"))
## MrVBF:
rsaga.geoprocessor(lib="ta_morphometry", module=8, param=list(DEM="DEM5LIDARf2.sgrd", MRVBF="VBF.sgrd"))
## valley depth:
rsaga.geoprocessor(lib="ta_channels", module=7, param=list(ELEVATION="DEM5LIDARf2.sgrd", VALLEY_DEPTH="VDP.sgrd"))
## deviation from the mean value:
rsaga.geoprocessor(lib="statistics_grid", module=1, param=list(GRID="DEM5LIDARf2.sgrd", DEVMEAN="DVM.sgrd", RADIUS=11))
## incoming solar radiation:
rsaga.geoprocessor(lib="ta_lighting", module=2, param=list(GRD_DEM="DEM5LIDARf2.sgrd", GRD_TOTAL="INS.sgrd", LATITUDE=clat, DHOUR=2, PERIOD=2, DDAYS=5))  # time-consuming!

## read maps back to R:
LSP.lst <- c("SLP.sdat", "CRV.sdat", "VBF.sdat", "TWI.sdat", "VDP.sdat", "DVM.sdat", "INS.sdat")
grids5mLSP <- stack(LSP.lst)
grids5mLSP <- as(grids5mLSP, "SpatialGridDataFrame")
proj4string(grids5mLSP) <- proj4string(grids5m)
grids5mLSP$DEM5LIDARf <- grids5m$DEM5LIDARf

# ------------------------------------------------------------
# Unsupervised extraction of landforms:
# ------------------------------------------------------------

## Principal components: 
pc.dem <- prcomp(~SLP+VBF+TWI+VDP+DVM+INS, scale=TRUE, grids5mLSP@data)
biplot(pc.dem, arrow.len=0.1, xlabs=rep(".", length(pc.dem$x[,1])), main="PCA biplot")
## LSPs are relatively uncorrelated;

## Run fuzzy k-means classification
## Determine number of clusters
demdata <- as.data.frame(pc.dem$x)
wss <- (nrow(demdata)-1)*sum(apply(demdata,2,var))
for (i in 2:20) {wss[i] <- sum(kmeans(demdata, centers=i)$withinss)}
plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
## does not converge :(

kmeans.dem <- kmeans(demdata, 9)
grids5m$kmeans.dem <- kmeans.dem$cluster
grids5m$landform <- as.factor(kmeans.dem$cluster)
summary(grids5m$landform)
spplot(grids5m["landform"], col.regions=rainbow(12))
write.asciigrid(grids5m["kmeans.dem"], "landform.asc", na.value=-1)
## plot in GE:
plotKML(grids5m["landform"], colour_scale=rainbow(9))

# ------------------------------------------------------------
# Fitting variograms for different landform classes:
# ------------------------------------------------------------

lidar.sample.ov <- over(y=grids5m["landform"], x=lidar.sample)
lidar.sample$landform <- lidar.sample.ov$landform

landform.no <- length(levels(lidar.sample.ov$landform))
landform.vgm <- as.list(rep(NA, landform.no))
landform.par <- data.frame(landform=as.factor(levels(lidar.sample.ov$landform)), Nug=rep(NA, landform.no), Sill=rep(NA, landform.no), range=rep(NA, landform.no))

landform.vgm <- NULL
for(i in 1:length(levels(lidar.sample.ov$landform))) {
   tmp <- lidar.sample[which(lidar.sample$landform==levels(lidar.sample$landform)[i]),]
   landform.vgm[[i]] <- fit.variogram(variogram(Z~1, tmp, cutoff=50*pixelsize), vgm(psill=var(tmp$Z), "Gau", sqrt(areaSpatialGrid(grids25m))/4, nugget=0))
   landform.par$Nug[i] <- round(landform.vgm[[i]]$psill[1], 1)
   landform.par$Sill[i] <- round(landform.vgm[[i]]$psill[2], 1)
   landform.par$range[i] <- round(landform.vgm[[i]]$range[2], 1)
}
landform.par  

# ------------------------------------------------------------
# Spatial prediction of soil mapping units:
# ------------------------------------------------------------

## convert the raster map to polygon map:
rsaga.esri.to.sgrd(in.grids="soilmu.asc", out.sgrd="soilmu.sgrd", in.path=getwd())
rsaga.geoprocessor(lib="shapes_grid", module=6, param=list(GRID="soilmu.sgrd", POLYGONS="soilmu.shp", CLASS_ALL=1))
soilmu <- readOGR("soilmu.shp", "soilmu")
proj4string(soilmu) <- proj4string(grids5m)
plotKML(soilmu["NAME"])
## convert the polygon to line map:
rsaga.geoprocessor(lib="shapes_lines", module=0, param=list(POLYGONS="soilmu.shp", LINES="soilmu_l.shp"))
## derive the buffer map using the shape file:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(TARGET_OUT_GRID="soilmu_r.sgrd", INPUT="soilmu_l.shp", FIELD=0, LINE_TYPE=0, TARGET_USER_SIZE=pixelsize, TARGET_USER_XMIN=grids5m@bbox[1,1]+pixelsize/2, TARGET_USER_XMAX=grids5m@bbox[1,2]-pixelsize/2, TARGET_USER_YMIN=grids5m@bbox[2,1]+pixelsize/2, TARGET_USER_YMAX=grids5m@bbox[2,2]-pixelsize/2))
## buffer distance:
rsaga.geoprocessor(lib="grid_tools", module=10, param=list(SOURCE="soilmu_r.sgrd", DISTANCE="soilmu_dist.sgrd", ALLOC="tmp.sgrd", BUFFER="tmp.sgrd", DIST=sqrt(areaSpatialGrid(grids25m))/3, IVAL=pixelsize))
## surface specific points (medial axes!):
rsaga.geoprocessor(lib="ta_morphometry", module=3, param=list(ELEVATION="soilmu_dist.sgrd", RESULT="soilmu_medial.sgrd", METHOD=1))
## read to R:
grids5m$soilmu_medial <- readGDAL("soilmu_medial.sdat")$band1
## generate training pixels:
grids5m$weight <- abs(ifelse(grids5m$soilmu_medial>=0, 0, grids5m$soilmu_medial))
dens.weight <- as.im(as.image.SpatialGridDataFrame(grids5m["weight"]))
# image(dens.weight)
training.pix <- rpoint(length(grids5m$weight)/10, f=dens.weight)
# plot(training.pix)
training.pix <- data.frame(x=training.pix$x, y=training.pix$y , no=1:length(training.pix$x))
coordinates(training.pix) <- ~x+y
proj4string(training.pix) <- proj4string(grids5m)
# spplot(grids5m["weight"], col.regions=grey(rev(seq(0,0.95,0.05))), sp.layout=list("sp.points", pch="+", training.pix, col="yellow"))
#writeOGR(training.pix, "training_pix.shp", "training.pix", "ESRI Shapefile")

## Run MLR:
training.pix.ov <- cbind(over(y=grids5m, x=training.pix), over(y=grids5mLSP, x=training.pix))
library(nnet)
## fit the model:
mlr.soilmu <- multinom(soilmu.c~DEM5LIDARf+SLP+VBF+TWI+VDP+DVM+INS, training.pix.ov)
# summary(mlr.soilmu)
grids5m$soilmu.mlr <- predict(mlr.soilmu, grids5mLSP@data)
spplot(grids5m["soilmu.mlr"], col.regions=rainbow(length(levels(grids5m$soilmu.c))))
plotKML(grids5m["soilmu.mlr"])

## the kappa statistics (the whole map!)
sel <- !is.na(grids5m$soilmu.c)
Kappa(confusion(grids5m$soilmu.c[sel], grids5m$soilmu.mlr[sel]))
x <- confusion(grids5m$soilmu.c, grids5m$soilmu.mlr)
#agreementplot(x)

# ------------------------------------------------------------
# Extraction of memberships:
# ------------------------------------------------------------

## mask-out classes with <5 points:
mask.c <- as.integer(attr(summary(training.pix.ov$soilmu.c[summary(training.pix.ov$soilmu.c)<5]), "names"))
mask.c
## fuzzy exponent:
fuzzy.e <- 1.2
tvars <- c("DEM5LIDARf", "SLP", "VBF", "TWI", "VDP", "DVM", "INS")
## extract the class centroids:
class.c <- aggregate(training.pix.ov[,tvars], by=list(training.pix.ov$soilmu.c), FUN="mean")
class.sd <- aggregate(training.pix.ov[,tvars], by=list(training.pix.ov$soilmu.c), FUN="sd")
## derive distances in feature space:
distmaps <- as.list(levels(grids5m$soilmu.c)[mask.c])
tmp <- rep(NA, nrow(grids5m))
for(c in (1:length(levels(grids5m$soilmu.c)))[mask.c]){
  distmaps[[c]] <- data.frame(DEM5LIDARf=tmp, SLP=tmp, VBF=tmp, TWI=tmp, VDP=tmp, DVM=tmp, INS=tmp)
  for(j in tvars){
    distmaps[[c]][j] <- ((grids5mLSP@data[j]-class.c[c,j])/class.sd[c,j])^2
  }
}
## sum up distances per class:
distsum <- data.frame(tmp)
for(c in (1:length(levels(grids5m$soilmu.c)))[mask.c]){
  distsum[paste(c)] <- sqrt(rowSums(distmaps[[c]], na.rm=T, dims=1))
}
distsum[1] <- NULL
str(distsum)
totsum <- rowSums(distsum^(-2/(fuzzy.e-1)), na.rm=T, dims=1)
## derive the fuzzy membership:
for(c in (1:length(levels(grids5m$soilmu.c)))[mask.c]){
  grids5m@data[paste("mu_", c, sep="")] <- (distsum[paste(c)]^(-2/(fuzzy.e-1))/totsum)[,1]
}
spplot(grids5m[c("mu_1","mu_2","mu_3","mu_4","mu_5","mu_6")], at=seq(0,1,0.05), col.regions=grey(rev(seq(0,0.95,0.05))))
## this is different of course from the exercises 3
write.asciigrid(grids5m["mu_4"], "mu_4.asc", na.value=-1)
write.asciigrid(grids5m["mu_5"], "mu_5.asc", na.value=-1)

## end of script!