rm(list = ls())
library(rgdal)
library(sp)
library(raster)
library(dplyr)
library(ggplot2)


# Read in all Data --------------------------------------------------------
# Note: all of the commented out code is run once, then outputed a final
# csv file that contains everything so you don't have to write in the 
# large shapefiles each time. 
# # AADT and Percent Truck
# aadt2016 <- readOGR("../data/roadCovariateData/AADT2016/", "TDV_Shapefile_AADT_2016")
# 
# # Pavement Type Data (CURRENTLY NOT ABLE TO READ IN DATA)
# # fgdb <- "../data/roadCovariateData/pavementTypeData/RoadwayInventorySystem2017Pub.gdb"
# # subset(ogrDrivers(), grepl("GDB", name))
# # fc_list <- ogrListLayers(fgdb)
# 
# # Road Speed
# roadSpeed <- readOGR("../data/roadCovariateData/roadSpeed/Streets_shp/", "StreetSegment")
# 
# Recording site locations
rec.locs <- read.csv("../data/locations.csv")
coordinates(rec.locs) <- ~longitude + latitude
proj4string(rec.locs) <- CRS("+proj=longlat +datum=WGS84")
# Convert to UTM NAD83
rec.locs <- spTransform(rec.locs, "+proj=utm +zone=18 +units=m +datum=NAD83")
rec.locs.bounds <- as.vector(extent(rec.locs))
# 
# 
# # Setting Boundary for RC Computation -------------------------------------
# radius <- 600 # meters
# tmp <- c(-radius, radius, -radius, radius)
# rec.locs.bounds <- rec.locs.bounds + tmp
# northing.diff <- abs(rec.locs.bounds[4] - rec.locs.bounds[3])
# easting.diff <- abs(rec.locs.bounds[2] - rec.locs.bounds[1])
# change <- (northing.diff - easting.diff) / 2
# rec.locs.bounds[1] <- rec.locs.bounds[1] - change
# rec.locs.bounds[2] <- rec.locs.bounds[2] + change
# 
# 
# #  Crop the data sets to GV Region ----------------------------------------
# 
# gv.road.data <- crop(aadt2016, extent(rec.locs.bounds))
# gv.road.data.f <- fortify(gv.road.data)
# # Joining condition
# gv.road.data$id <- row.names(gv.road.data)
# # Join all the data to the data frame so you have aadt in data frame
# gv.road.data.full <- left_join(gv.road.data.f, gv.road.data@data)
# # Filter out erroneous data occurring over Conesus Lake
# road.dat.lake <- crop(gv.road.data, extent(276395.8, 277573, 4731115, 4742194))
# bad.ids <- road.dat.lake@data$OBJECTID
# gv.road.data.f <- filter(gv.road.data.full, !(OBJECTID %in% bad.ids[1]))
# write.csv(gv.road.data.f, "../data/gvRoadData.csv", row.names = FALSE)
# 
# roadSpeed.gv <- crop(roadSpeed, extent(rec.locs.bounds))
# roadSpeed.gv.f <- fortify(roadSpeed.gv)
# # Specify the joining condition
# roadSpeed.gv$id <- row.names(roadSpeed.gv)
# # Join data to the data frame so you have speed values in the data frame
# roadSpeed.gv <- left_join(roadSpeed.gv.f, roadSpeed.gv@data)
# # Output as csv so you don't have to to read it in again since it takes foreer
# write.csv(roadSpeed.gv, "../data/roadSpeedData.csv", row.names = FALSE)

# Fill in missing values from both data sets ------------------------------
# AADT
# gv.road.data.f <- arrange(gv.road.data.f, long, lat)
# missing.data <- which(with(gv.road.data.f, AADT == 0))
# actual.data <- which(with(gv.road.data.f, AADT != 0))
# # Set min.dist to the number of entries to start off
# # Compare each missing data point to non missing data points, and find the
# # closest one then set equal to that value
# for (i in missing.data) {
#   min.dist <- dim(gv.road.data.f)[1] - 1
#   for (j in actual.data) {
#     if (abs((i - j)) < min.dist) {
#       min.dist <- abs(i - j)
#       new.value <- j
#     }
#   }
#   gv.road.data.f$AADT[i] <- gv.road.data.f$AADT[new.value]
# }
# 
# all.roads <- left_join(gv.road.data.f, roadSpeed.gv, by = c("long", "lat"))
# 
# speed.data.f.new <- arrange(all.roads, long, lat)
# missing.data <- which(with(all.roads, is.na(SPEED)))
# actual.data <- which(with(all.roads, !is.na(SPEED)))
# for (i in missing.data) {
#   min.dist <- dim(all.roads)[1] - 1
#   for (j in actual.data) {
#     if (abs((i - j)) < min.dist) {
#       min.dist <- abs(i - j)
#       new.value <- j
#     }
#   }
#   all.roads$SPEED[i] <- all.roads$SPEED[new.value]
# }
# 
# 
# # Output so you don't have to do any of the above code --------------------
# 
# write.csv(all.roads, "../data/completeRoadData.csv", row.names = FALSE)

all.roads <- read.csv("../data/completeRoadData.csv")
# Convert all.roads to spatial data frame
coordinates(all.roads) <- ~long + lat
proj4string(all.roads) <- CRS("+proj=utm +zone=18 +units=m +datum=NAD83")

# Get unique recording locations
plot.locations <- as.data.frame(rec.locs)
plot.locations.matrix <- unique(as.matrix(plot.locations[, c("longitude", "latitude")]))

# Convert to raster -------------------------------------------------------

# Roads are broken up into approximately 10 x 10 meter blocks. 
# A car is about 4.5 meters long, so this is a little of 2 times the length of a car
r <- raster(extent(all.roads), nrow = 4485, ncol = 4485)
res(r)
raster.aadt <- rasterize(all.roads, r, "AADT", fun = max)
raster.speed <- rasterize(all.roads, r, "SPEED", fun = max)
raster.truck <- rasterize(all.roads, r, "PERC_TRUCK", fun = max)

# Convert rasters to points -----------------------------------------------
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

# Compute value for each point --------------------------------------------

# For all points in dat.points, compute distance between recording locations
# Add this to dat.points data frame

# Gives distance in meters
dists <- pointDistance(plot.locations.matrix, dat.points.mat[, c("x", "y")],
                       lonlat = FALSE)
dists.t <- t(dists)
# dists.t <- scale(dists.t)
dat <- cbind(dat.points, dists.t)
# Names of forest plots in desired order
names <- c("ADP-C", "ADP-E", "BTP-C", "BTP-E", "CLI-C", "CLI-E", "IFP-C", "IFP-E", 
           "IP-C", "IP-E", "JOR-C", "JOR-E", "MCC-C", "MCC-E", "MCP-C", "MCP-E", 
           "RA-C", "RA-E")
names <- gsub("-", "", names)
dat.names <- sapply(names, paste, "dist", sep = ".")
colnames(dat) <- c("x", "y", "aadt", "speed", "truck", dat.names)
dat.scaled <- dat
dat.scaled[, -c(1:5)] <- scale(dat[, -c(1:5)])

# Compute the road covariate site by site
road.covariate <- rep(0, length(names))
for (j in 6:dim(dat)[2]) {
  vals <- dat[, j] <= 600
  road.values <- dat$aadt[vals] + dat$speed[vals] + dat$truck[vals] - 
    dat.scaled[, j][vals]
  road.covariate[j - 5] <- sum(road.values)
}
# Divide by 100 just to scale values to smaller numbers
final.road.covariate <- rep(rep(road.covariate, each = 29), times = 3) / 100
names.ordered <- rep(rep(names, each = 29), times = 3)
plot.locations$roadCovariate <- final.road.covariate
names(plot.locations) <- c("easting", "northing", "road.covariate")
write.csv(plot.locations, "../data/predictors.csv", row.names = FALSE)
write.table(final.road.covariate, "../beta/data/roadCovariate", row.names = FALSE,
            col.names = FALSE, sep = "\t")
