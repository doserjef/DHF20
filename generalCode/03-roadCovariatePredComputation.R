
rm(list = ls())
# 25-covariateGVRegion ----------------------------------------------------

# Purpose: compute the road covariate for all locations bounded by the 
# recording sites

# Spatial Functions
# library(rgdal)
library(sp)
library(raster)
# For computing nearest neighbors
# Data manipulation
library(dplyr)

# Read in data 
all.roads <- read.csv("../data/completeRoadData.csv")

# Convert all.roads back to a spatial data frame
coordinates(all.roads) <- ~long + lat
proj4string(all.roads) <- CRS("+proj=utm +zone=18 +units=m +datum=NAD83")


# Get distinct pixels for all roads ---------------------------------------
# Roads broken up into 10 x 10 m blocks
r <- raster(extent(all.roads), nrow = 4485, ncol = 4485)
res(r)
raster.aadt <- rasterize(all.roads, r, "AADT", fun = max)
raster.speed <- rasterize(all.roads, r, "SPEED", fun = max)
raster.truck <- rasterize(all.roads, r, "PERC_TRUCK", fun = max)

# Convert the rasters back to points for the pixels
aadt.points <- as.data.frame(rasterToPoints(raster.aadt))
colnames(aadt.points) <- c("x", "y", "aadt")
speed.points <- as.data.frame(rasterToPoints(raster.speed))
colnames(speed.points) <- c("x", "y", "speed")
truck.points <- as.data.frame(rasterToPoints(raster.truck))
colnames(truck.points) <- c("x", "y", "truck")
tmp <- inner_join(aadt.points, speed.points, 
                  by = c("x", "y"))
dat.points <- inner_join(tmp, truck.points, by = c("x", "y"))
# Standardize variables to have mean of 0 and sd of 1 
dat.points[, -c(1, 2)] <- scale(dat.points[, -c(1, 2)])
dat.points.mat <- as.matrix(dat.points)

# Get coordinates for entire GV region ------------------------------------
# Only predict for points within the radius of which you have road data
# so you don't have problems on the edges
radius <- 600
tmp <- c(radius, -radius, radius, -radius)
road.extent <- as.vector(extent(all.roads)) + tmp

# For sound, break into 250 x 250 m resolution (175 x 175)
r <- raster(extent(road.extent), nrow = 175, ncol = 175)
res(r)
# A data frame with longitudes and latitudes
gv.area <- as.matrix(rasterToPoints(r))
colnames(gv.area) <- c("x", "y")
str(gv.area)

write.table(gv.area, "../data/predLocations-SOUND", row.names = FALSE, col.names = FALSE, sep = "\t")

# Compute Road Covariate --------------------------------------------------

# Dividing this up into a for loop because otherwise the matrix becomes extremely 
# large

road.covariate <- rep(0, dim(gv.area)[1])
my.factor <- 175
for (i in 1:my.factor) {
  print(paste("We are on part ", i, " of ", my.factor, sep = ""))
  # Determine the specific part of the gv location you're currently working wth
  gv.area.part = gv.area[((i -1)*my.factor + 1) : (i *my.factor), ]
  
  # Compute distances from the current locations to road pixels
  dists <- pointDistance(gv.area.part, dat.points.mat[, c("x", "y")], lonlat = FALSE)
  dists.t <- t(dists)
  
  dat <- cbind(dat.points, dists.t)
  
  dat.scaled <- dat
  dat.scaled[, -c(1:5)] <- scale(dat[, -c(1:5)])
  
  # Will want to change this hard coding
  for (j in 5:dim(dat)[2]) {
    vals <- dat[, j] <= 600
    road.values <- dat$aadt[vals] + dat$speed[vals] + dat$truck[vals] - 
      dat.scaled[, j][vals]
    road.covariate[(i - 1)*my.factor + j - 4] <- sum(road.values)
  }
}

# Scale road covariate
road.covariate <- road.covariate / 100

gv.road.cov <- as.data.frame(cbind(gv.area, road.covariate))

names(gv.road.cov) <- c("x", "y", "z")

# # Read out data to a table
write.table(road.covariate, "../data/predRoadCovariate-SOUND", row.names = FALSE,
            col.names = FALSE, sep = "\t")
