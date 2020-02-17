rm(list = ls())
library(raster)
library(rgdal)
library(seewave)
library(viridis)

pred.rc <- read.table("../data/predRoadCovariate")
pred.locations <- read.table("../data/predLocations")

pred.locations.km <- pred.locations / 1000
gv.road.cov <- as.data.frame(cbind(pred.locations.km, pred.rc))
names(gv.road.cov) <- c("x", "y", "z")
# gv.road.cov$z[gv.road.cov$z > 4] <- 4



# Recording site locations ------------------------------------------------

rec.locs <- read.csv("../data/locations.csv")
coordinates(rec.locs) <- ~longitude + latitude
proj4string(rec.locs) <- CRS("+proj=longlat +datum=WGS84")
# Convert to UTM NAD83
rec.locs <- spTransform(rec.locs, "+proj=utm +zone=18 +units=m +datum=NAD83")
recording.sites <- as.data.frame(rec.locs)
sites <- unique(recording.sites) / 1000



# Road Data ---------------------------------------------------------------

roads <- readOGR("../data/roadData/", "AADT_2015_tdv")
roads <- spTransform(roads, "+proj=utm +zone=18 +units=km +datum=NAD83")
road.gv <- crop(roads, extent(gv.road.cov))
bounds <- as.vector(extent(gv.road.cov))

roc <- data.frame(lat = 43.1566, long = -77.6088)
coordinates(roc) <- ~long + lat
proj4string(roc) <- CRS("+proj=longlat +datum=WGS84")
roc <- spTransform(roc, "+proj=utm +zone=18 +units=m +datum=NAD83")
roc <- as.data.frame(roc) / 1000
png("../figures/nysRoads.png", width = 480, height = 480)
plot(roads)
rect(xleft = bounds[1], ybottom = bounds[3], xright = bounds[2],
     ytop = bounds[4], border = 'red', lwd = 3)
points(roc$long, roc$lat, las = 1, col = 'blue', pch = 19)
text(roc$long, roc$lat + 50, "Rochester")
arrows(x0 = roc$long, y0 = roc$lat + 35, x1 = roc$long, y1 = roc$lat + 5, xpd = TRUE, len = 0.1, col = 'black')
dev.off()

# Plot --------------------------------------------------------------------

png("../figures/roadCovariateRoads.png", width = 480, height = 480)
curr <- par()$mar
par(mar = c(5.1, 6.1, 4.1, 2.1))
plot.raster <- rasterFromXYZ(gv.road.cov)
magma.palette <- magma(25000)
plot(plot.raster, col = magma.palette, lty = 'blank', xlab = 'Easting (km)', 
     axes = FALSE, box = FALSE, 
     cex.lab = 1.5)
points(sites$longitude, sites$latitude, pch = 19, cex = 0.75)
lines(road.gv, lwd = 0.5, col = 'white')
axis(1, seq(260, 300, by = 10), pos = 4732)
axis(2, seq(4735, 4775, by = 10), las = 1, pos = 256)
title(ylab="Northing (km)", mgp=c(4,1,1), cex.lab = 1.5)
text(x = 304.5, y = 4769, "RC", xpd = TRUE)
dev.off()

