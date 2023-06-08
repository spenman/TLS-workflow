# This R script accompanies the paper:
# "An instructional workflow for using Terrestrial Laser Scanning (TLS) to 
# quantify vegetation structure for wildlife studies"

library(lidR)
library(TreeLS)
library(dplyr)
library(ggplot2)


# Process scans -----------------------------------------------------------

# Load raw point cloud
rawlas<-readTLSLAS("./data/pointcloud-example.las", select="xyzc")

plot(rawlas)


# Normalise plot ----------------------------------------------------------

# find middle of plot
x<-median(rawlas$X)
y<-median(rawlas$Y)

# Initial clip (40 x 40 m) to reduce processing time, but larger than the plot size to reduce edge artefacts
mylas<-lidR::clip_rectangle(rawlas, (x-20), (y-20), (x+20), (y+20))
plot(mylas)


ntls<-TreeLS::tlsNormalize(mylas, keep_ground=T)

# You can visualise the difference by plotting the normalised point cloud
plot(ntls)

# Crop to area of interest (20 x 20 m)
tls<-lidR::clip_rectangle(ntls, (x-10), (y-10), (x+10), (y+10))

rm(rawlas, mylas)


# Derive tree-related metrics ---------------------------------------------

# Detect and delineate trees
thin<-tlsSample(tls, smp.voxelize(0.005))

g<-plot(thin)

## NB default values for classifying trees have been adjusted to suit this vegetation type
# default values are (min_density = 0.1, min_h = 1, max_h = 3)
map<-treeMap(thin, map.hough(min_density = 0.0125, min_h = 4, max_h = 6), 0)
add_treeMap(g, map, color='yellow', size=1)

tls<-treePoints(tls, map, trp.crop(l=0.75))

# add the TreeID labels to plot
add_treeIDs(g, tls, cex = 1, col='pink')

# classify stem points
tls<-stemPoints(tls, stm.hough())
add_stemPoints(g, tls, color='pink', size=3)

# make the plot's inventory
inv<-tlsInventory(tls, d_method=shapeFit(shape = 'circle', algorithm='ransac', n=20))

# calculate plot size in Ha
plotsize<-(20*20)/(100*100)


# Calculate dbh
inv$dbh<-inv$Radius*2

# remove any stems < 0.1 dbh
inv<-dplyr::filter(inv, dbh>=0.1)

# Calculate Basal area (m2/Ha)
inv$BA<-pi*(inv$Radius^2)

tempBA<-sum(inv$BA)
BasalArea<-tempBA/plotsize

# Calculate stem density (stems/Ha)
plotstems<-nrow(inv)

StemDensity<-plotstems/plotsize

# Measure canopy height
Height<-max(inv$H)


# Density profile of plot -------------------------------------------------

# Voxelise point cloud 
tmptls<-voxelize_points(tls, 0.1)

# get the voxel data and assign it a value of 1
atst<-tmptls@data
atst$bin<-1

# filter out ground & litter (points below 0.1 m)
atst<-filter(atst, atst$Z>0.1)

# plot the density of veg voxels at different heights
j<-ggplot(atst, aes(Z)) + geom_density() + labs(title = "Vegetation density plot")
j


# Use the density profile to inform canopy base height and understorey height decisions
# In this case, the density plot suggests that the understorey shrub layer is 
# around 4 m in height as that is where the curve reaches a local minimum, and
# the canopy looks to be around 11 m where the density of vegetation increases

# Canopy base height
cbht<-11

# Understorey height
uht<-4


# Canopy metrics ----------------------------------------------------------

## Canopy Cover ####

# Filter cloud to only include canopy points (i.e. all points above 11 metres)
ctls<-filter_poi(tls, tls$Z>cbht)

xl2<-round(min(ctls$X), 0)
yb2<-round(min(ctls$Y), 0)

# calculate grid metrics for the canopy point cloud
canras<-grid_metrics(ctls, sum(Z), res = 0.1, start=c(xl2, yb2))

plot(canras, col="black", legend=FALSE, yaxt='n', xaxt='n', main="Canopy cover")

# get values from the raster stack (canras)
ccvr<-canras[[1]]@data@values

# assign 1 to cells with a value >0
ccvr[ccvr>0]<-1

# count non-empty cells
fcell<-sum(ccvr, na.rm=T)

# total cells (20 x 20 m plot resolution = 0.1 m)
tcells<-((20/0.1)+1)^2

# Calculate cover % (full cells/total cells)
cvr<-(fcell/tcells)*100


## Canopy height model and heterogeneity ####
a<-0.5

### Canopy height model
chm <- grid_canopy(tls, res = a, p2r())
col <- height.colors(50)

# Calculate rumple index from the chm (i.e. 3D heterogeneity of the canopy)
canrum<-rumple_index(chm)
canrum<-round(canrum, digits=3)

plot(chm, col = col, main=paste0("Canopy height model: resolution ", a, " m \n Rumple Index ", canrum))


## Canopy volume and density ####

# Voxelise the canopy cloud
vctls<-voxelize_points(ctls, 0.2)

atst<-vctls@data
atst$bin<-1

tmpdat<-atst %>%
  group_by(Z, bin) %>%
  summarise(vox=n())

zrange<-(max(tmpdat$Z)-min(tmpdat$Z))/0.2

# canopy density
candens<-(sum(tmpdat$vox))/(10201*zrange)

# canopy gap volume in cubic metres
cgapvol<-(1-candens)*(20*20*(max(tmpdat$Z)-min(tmpdat$Z)))


# Extract metrics for other strata ----------------------------------------

# In this case the understorey

# Filter out tree stem points
utls<-filter_poi(tls, tls$TreeID==0)

# Filter pointcloud to only include understorey (i.e. Z<4 which was previously identified as the understorey height)
utls<-filter_poi(utls, Z<uht)

# Scan coverage is too sparse at edges of this plot, so understorey was best assessed in a 10 x 10 m plot
utls<-clip_rectangle(utls, (x-5), (y-5), (x+5), (y+5))

# resolution for the understorey 'canopy' height model
a<-0.1

### understorey height model (a canopy height model applied to understorey point cloud)
chm <- grid_canopy(utls, res = a, p2r())
col <- height.colors(50)

# 3D heterogeneity of understorey (i.e. calculate rumple index from the chm)
usRI<-rumple_index(chm)
usRI<-round(usRI, 3)

# plot the understorey height model
plot(chm, col = col, main=paste0("Understorey height model \n Rumple index: ", usRI))

## Understorey cover #####

# filter out ground & litter (points below 0.1 m)
utls<-filter_poi(utls, Z>0.1)

# find corner of utls
xl2<-round(min(utls$X), 0)
yb2<-round(min(utls$Y), 0)

utls<-filter_poi(utls, Z>=0.1)

# calculate grid metrics for the understorey point cloud
uras<-grid_metrics(utls, sum(Z), res = 0.1, start=c(xl2, yb2))

plot(uras, col="black", legend=FALSE, yaxt='n', xaxt='n', main="Understorey cover")

# get values from the raster stack (canras)
ucvr<-uras[[1]]@data@values

# assign 1 to cells with a value >0
ucvr[ucvr>0]<-1

# count non-empty cells
fcell<-sum(ucvr, na.rm=T)

# Calculate cover % (full cells/total cells)
ucvr<-(fcell/10201)*100

## Understorey density ####

# Voxelise point cloud
vtls<-voxelize_points(utls, 0.1)

# get the voxel data and assign it a value of 1
atst<-vtls@data
atst$bin<-1

tmpdat<-atst %>%
  group_by(Z, bin) %>%
  summarise(vox=n())

zrange<-(max(tmpdat$Z)-min(tmpdat$Z))/0.1

# understorey density
udens<-(sum(tmpdat$vox))/(10201*zrange)

# Understorey gap volume in cubic metres
ugapvol<-(1-udens)*(10*10*(max(tmpdat$Z)-min(tmpdat$Z)))


# Foliage height diversity ------------------------------------------------

# This includes the full Z profile of the plot (i.e. all strata) 20 x 20 m plot using 0.2 m voxels

fhdtls<-voxelize_points(tls, 0.2)

# get the voxel data and assign it a value of 1
atst<-fhdtls@data
atst$bin<-1

fhddat<-atst %>%
  group_by(Z, bin) %>%
  summarise(vox=n())

# vegetation density
fhddat$dens<-fhddat$vox/10201

# Foliage height diversity (analagous measure)
fhddat$temp<--(fhddat$dens*(log(fhddat$dens)))

fhd<-sum(fhddat$temp)
