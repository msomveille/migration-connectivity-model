## Load required libraries

library(dggridR)
library(rworldmap)
library(rgeos)
library(igraph)

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


setwd('/Users/mariussomveille/Desktop/Yale-MPIO/White Storks/migration-connectivity-model')

##  Import model outputs  ##

pathways.west <- read.csv("modelOutput_whiteStorks_west_lambda2.csv", header=F)
pathways.east <- read.csv("modelOutput_whiteStorks_east_lambda2.csv", header=F)

pathways.west <- apply(pathways.west, c(1,2), function(x) ifelse(x<0.1, 0, x))
pathways.east <- apply(pathways.east, c(1,2), function(x) ifelse(x<0.1, 0, x))

pathways = pathways.west + pathways.east


##  Compute the prevalence of western migrants at each breeding site, and plot its spatial distribution  ##

prevalence.west <- apply(pathways.west, 2, sum) / apply(pathways, 2, sum)

plot(hexgrid2, col="dark grey", border="dark grey", bg="light grey")
plot(hexgrid2[which(breeding.grounds == TRUE | resident.grounds==TRUE)][which(prevalence.west > 0.66)], add=T, col="orange", border="orange")
plot(hexgrid2[which(breeding.grounds == TRUE | resident.grounds==TRUE)][which(prevalence.west < 0.33)], add=T, col="blue", border="blue")
plot(hexgrid2[which(breeding.grounds == TRUE | resident.grounds==TRUE)][which(prevalence.west >= 0.33 & prevalence.west <= 0.66)], add=T, col="green", border="green")


# Summary statistics indicating the location of the migratory divide
divide = max(hexgrid2_centroids[which(breeding.grounds == TRUE | resident.grounds==TRUE),][which(prevalence.west > 0.66),1])


hist(apply(pathways, 2, sum))   # distribution of the number of individuals at breeding sites
hist(apply(pathways, 1, sum))   # distribution of the number of individuals at wintering sites



##  Plot the winter redistribution of individuals from a selected breeding site  ##
 
s = 285   #40, 190 , 500
plot(hexgrid2, col="dark grey", border="dark grey", bg="light grey")
rbPal <- colorRampPalette(c("light blue", "dark blue"))
datcol <- rbPal(10)[as.numeric(cut(pathways[,s], breaks=c(0.1,1,2,5,10,20,50,100,200,500,max(pathways))))]
plot(hexgrid2[which(nonbreeding.grounds == TRUE | resident.grounds==TRUE)], add=T, col=datcol, border=datcol)
plot(hexgrid2[which(breeding.grounds == TRUE | resident.grounds==TRUE)][s], add=T, col="red", border="red")




##  Compute transition/conductance matrix for western and eastern flyways in order to plot migratory paths ##

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

breeding.hex <- which(breeding.grounds==TRUE | resident.grounds==TRUE)
nonbreeding.hex <- which(nonbreeding.grounds==TRUE | resident.grounds==TRUE)
resident.hex <- which(resident.grounds==TRUE)
# Western flyway
g.west <- graph.adjacency(Conductance.west, weighted=T)
# Eastern flyway
g.east <- graph.adjacency(Conductance.east, weighted=T)



##  Plot the migratory path to the wintering sites where most individuals go from the selected breeding site  ##

w = which(pathways[,s] == max(pathways[,s]))
s.paths <- get.shortest.paths(g.east, from=nonbreeding.hex[w], to=breeding.hex[s])
s.paths <- match(s.paths$vpath[[1]], V(g.east))
lines(hexgrid2_centroids[s.paths,], col="yellow", lwd=3)
plot(hexgrid2[which(breeding.grounds == TRUE | resident.grounds==TRUE)][s], add=T, col="red", border="red")

