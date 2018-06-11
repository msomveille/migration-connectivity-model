## Load required libraries

library(dggridR)
library(rworldmap)
library(rgeos)
library(igraph)
library(raster)

setwd('/Users/mariussomveille/Desktop/Yale-MPIO/White Storks/Data')


##  Construct an hexagon grid covering the region of interest  ##

hexgrid <- dgconstruct(projection="ISEA", topology="HEXAGON", res=8, metric=T)
hexgrid_center <- dgSEQNUM_to_GEO(hexgrid, 1:65612) # 21872 / 65612
hexgrid_centroids <- cbind(hexgrid_center$lon_deg, hexgrid_center$lat_deg)
hex_sel <- which(hexgrid_centroids[,1] > -29 & hexgrid_centroids[,1] < 59 & hexgrid_centroids[,2] > -40 & hexgrid_centroids[,2] < 65)
hexgrid2 <- dgcellstogrid(hexgrid, hex_sel, frame=F,wrapcells=TRUE)
newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, proj4string(hexgrid2))
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- gUnaryUnion(newmap, id=newmap@data$world)
hexgrid2 <- gIntersection(hexgrid2, newmap, byid=T)
hexgrid2_centroids <- matrix(unlist(lapply(hexgrid2@polygons, function(x) x@labpt)), byrow=T, ncol=2)


##  Load gthe geographical distribution of white storks (downloaded from BirdLife International)  ##

white.storks.distribution <- readOGR("distribution","Ciconia_ciconia_3835_BL", verbose=FALSE)
white.storks.distribution <- spTransform(white.storks.distribution, proj4string(hexgrid2))
nonbreeding.grounds <- gIntersects(white.storks.distribution[which(white.storks.distribution$SEASONAL == 3),], hexgrid2, byid=T)
breeding.grounds <- gIntersects(white.storks.distribution[which(white.storks.distribution$SEASONAL == 2),], hexgrid2, byid=T)
resident.grounds <- gIntersects(white.storks.distribution[which(white.storks.distribution$SEASONAL == 1),], hexgrid2, byid=T)
res2 <- which(nonbreeding.grounds + breeding.grounds + resident.grounds > 1)
nonbreeding.grounds[res2] <- FALSE
breeding.grounds[res2] <- FALSE
resident.grounds[res2] <- TRUE


##  Compute the shortest migration path between every pair of breeding/wintering sites (i.e. seasonally occupied hexagons)  ##

##  Compute transition/conductance matrix for western and eastern flyways  ##

neighbours.west <- gTouches(hexgrid2, byid=T)
neighbours.west[378,379] <- TRUE  # Connect at gibraltar
neighbours.west[379,378] <- TRUE
neighbours.west[4526,which(neighbours.west[4526,] == TRUE)] <- FALSE  # Disconnect eastern flyway
neighbours.west[4527,which(neighbours.west[4527,] == TRUE)] <- FALSE
neighbours.west[4528,which(neighbours.west[4528,] == TRUE)] <- FALSE
neighbours.west[4529,which(neighbours.west[4529,] == TRUE)] <- FALSE
neighbours.west[4530,which(neighbours.west[4530,] == TRUE)] <- FALSE
neighbours.west[4531,which(neighbours.west[4531,] == TRUE)] <- FALSE
neighbours.west[4532,which(neighbours.west[4532,] == TRUE)] <- FALSE
neighbours.west[4533,which(neighbours.west[4533,] == TRUE)] <- FALSE
neighbours.west[4701,which(neighbours.west[4701,] == TRUE)] <- FALSE
neighbours.west[4702,which(neighbours.west[4702,] == TRUE)] <- FALSE
neighbours.west[4703,which(neighbours.west[4703,] == TRUE)] <- FALSE
neighbours.west[4704,which(neighbours.west[4704,] == TRUE)] <- FALSE
neighbours.west[4705,which(neighbours.west[4705,] == TRUE)] <- FALSE
neighbours.west[4706,which(neighbours.west[4706,] == TRUE)] <- FALSE
neighbours.west[4707,which(neighbours.west[4707,] == TRUE)] <- FALSE
Conductance.west <- neighbours.west
Conductance.west[which(Conductance.west == FALSE)] <- 0
Conductance.west[which(Conductance.west == TRUE)] <- 1

neighbours.east <- gTouches(hexgrid2, byid=T)
Conductance.east <- neighbours.east
Conductance.east[which(Conductance.east == FALSE)] <- 0
Conductance.east[which(Conductance.east == TRUE)] <- 1


# Simulate shortest paths between pairs of breeding and non-breeding hexagons
breeding.hex <- which(breeding.grounds==TRUE | resident.grounds==TRUE)
nonbreeding.hex <- which(nonbreeding.grounds==TRUE | resident.grounds==TRUE)
resident.hex <- which(resident.grounds==TRUE)
# Western flyway
g.west <- graph.adjacency(Conductance.west, weighted=T)
distance.matrix.west <- shortest.paths(g.west, v=nonbreeding.hex, to=breeding.hex, mode="out")
# Eastern flyway
g.east <- graph.adjacency(Conductance.east, weighted=T)
distance.matrix.east <- shortest.paths(g.east, v=nonbreeding.hex, to=breeding.hex, mode="out")

write.csv(distance.matrix.west, "distanceMatrix_west.csv", row.names=F, col.names=F)
write.csv(distance.matrix.east, "distanceMatrix_east.csv", row.names=F, col.names=F)



##  Spatial distribution of energy supply derived from NDVI  ##

NDVI_june <- raster("MOD_NDVI_M_2010-06.TIFF")
NDVI_dec <- raster("MOD_NDVI_M_2010-12.TIFF")
NDVI_june.hex <- extract(NDVI_june, hexgrid2[breeding.hex], mean, na.rm=T)
NDVI_dec.hex <- extract(NDVI_dec, hexgrid2[nonbreeding.hex], mean, na.rm=T)
energySupply_summer = NDVI_june.hex * 10 
energySupply_winter = NDVI_dec.hex * 10

# PLOT
plot(hexgrid2, col="dark grey", border="dark grey", bg="light grey")
rbPal <- colorRampPalette(c("yellow", "red"))
datcol <- rbPal(10)[as.numeric(cut(energySupply_summer, breaks= quantile(energySupply_summer, seq(0,1,0.1))))]
plot(hexgrid2[which(breeding.grounds == TRUE | resident.grounds==TRUE)], add=T, col=datcol, border=datcol)
rbPal <- colorRampPalette(c("light blue", "dark blue"))
datcol <- rbPal(10)[as.numeric(cut(energySupply_winter, breaks= quantile(energySupply_winter, seq(0,1,0.1))))]
plot(hexgrid2[which(nonbreeding.grounds == TRUE | resident.grounds==TRUE)], add=T, col=datcol, border=datcol) 

write.csv(energySupply_summer, "energySupply_summer.csv", row.names=F, col.names=F)
write.csv(energySupply_winter, "energySupply_winter.csv", row.names=F, col.names=F)