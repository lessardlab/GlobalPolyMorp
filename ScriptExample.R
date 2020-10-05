DATA PREPARATION
#datamatching = polymorphism database
#lat-longGABI = GABI database with lat/long coordinates

1. Creating data frame with poly_id
specieslist<-read.csv("matched-polymorphism_bg.csv",stringsAsFactors = TRUE)
occ<-read.csv("Lat-Long_Data_GABI.csv",stringsAsFactors = TRUE)

occ$valid_species_name <- factor(occ$valid_species_name, levels = levels(specieslist$taxon_code))
occ$poly_id <- 0

for(i in 1:nrow(occ)){
temp<-subset(specieslist,specieslist$taxon_code==occ$valid_species_name[i])
occ$poly_id[i]<-temp$poly_id[1]
}

occ_match<-occ[,-6] #elevation has too many NAs
occ_match<-na.omit(occ_match)#remove NAs to be able to attribute continent

2. Matching Continents Per Occurrence
library(rworldmap)
library(rgeos)
library(maptools)
library(cleangeo) ## For clgeo_Clean()
library(raster)
library(sp)

coords<-as.data.frame(cbind(occ_match$dec_long,occ_match$dec_lat))
colnames(coords)<-c("Lat","Long")
coords<-na.omit(coords)

###########
# Error starts here:
########

coords<-SpatialPointsDataFrame(coords=coords,data=coords,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
#convert coordinates into spatial object

cont<-rgdal::readOGR("~/Documents/Work/Concordia/TEACHING/Supervision/MSc/Frédérique Larichelière/Masters Project (Frederique)/FINAL SCRIPT/ne_50m_admin_0_countries")
proj4string(coords)<- proj4string(cont)

points<-over(coords,cont)#must remove NAs from dataframe

occ_match["Continent"]<-points$CONTINENT
