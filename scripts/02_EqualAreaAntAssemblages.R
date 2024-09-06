#### Create ant 1x1 assemblages with equal area projection
## Gabriel Muñoz 
## Nov 2021
#############


# Load data
occ_new<-read.csv("Data/occ_withsitedensity_taxocorrect.csv",stringsAsFactors = F)

# Remove duplicates

dupRec <- duplicated(
  paste0(occ_new$dec_lat,
         occ_new$dec_long,
         occ_new$valid_species_name))


# Extent of observations 
rangeLat <- range(unique(round(occ_new$dec_lat)))
rangeLon <- range(unique(round(occ_new$dec_long)))

# Visualize data

plot(occ_new$dec_long, 
     occ_new$dec_lat, 
     pch = ".", col = "red")

################################
head(AntPoints)
# First set up a projection to the spatial points df. 

crsAn <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
AntPoints <- SpatialPointsDataFrame(data.frame((occ_unique$dec_long), 
                                               (occ_unique$dec_lat)),
                                    data.frame(occ_unique))
AntPoints@proj4string <- crsAn
################ Right projection Equal Area ? ################
crs.laea <- CRS("+proj=laea +lat_0=90 +lon_0=0
               +x_0=0 +y_0=0 +ellps=WGS84 
               +datum=WGS84 +units=m +no_defs")


AntReproj <- spTransform(AntPoints, crs.laea)

library(raster)
r <- raster(ext = extent(-180, 180, -70, 70), res=c(1,1))
values(r) <- 1:ncell(r)
plot(r)
points(AntPoints, pch = ".")
rA <- projectRaster(r, 
                    crs=crs.laea)
plot(rA)
points(AntReproj, pch = ".")

library(rgdal)
p <- rasterToPolygons(r) 
#p <- as(r, "SpatialPolygonsDataFrame")
pA <- spTransform(p, crs.laea)

# plot(pA, add=TRUE)
# 

GridPoly <- as(pA,"SpatialPolygons")
occurrences<- over(AntReproj,GridPoly)
corEA <- coordinates(GridPoly)
AntReproj$ID_ea 

occurrences
dim(PAMat_eqArea)
length(occurrences)
AntReproj$ID_ea <- occurrences
PAMat_eqArea <- table(AntReproj$ID_ea,
                      AntReproj$valid_species_name)

GridCentroids <- data.frame("x" = corEA[,1],"y" = corEA[,2])
GridCentroids <- GridCentroids[as.numeric(rownames(PAMat_eqArea)),]
GridID <- paste0(c(GridCentroids$x),"_", c(GridCentroids$y))
rownames(PAMat_eqArea) <- GridID

GridCentroids <- data.frame("x" = corEA[,1],"y" = corEA[,2])

AntData <- data.frame(AntReproj,
                      GridCentroids[occurrences,])


names(AntData)[c(21,22)] <- c("x_grid", "y_grid")

# remove previously climate data and extract 
AntData$MAP.point <- NULL
AntData$MAT.point <- NULL

#Load climate data
MAT<-raster("Data/wc2.0_bio_10m_01.tif")
TAP<-raster("Data/wc2.0_bio_10m_12.tif")

coords <- SpatialPoints(data.frame(AntData$dec_lat, 
                                   AntData$dec_long))

# Note I extract from the globally scaled environment, not "ant distribution" scaled 
MAT.point<-extract(MAT,coords)
TAP.point<-extract(TAP,coords)

# add climatic data to ant data.frame 
AntData$MAT <- MAT.point
AntData$TAP <- TAP.point












