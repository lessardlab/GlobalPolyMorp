#### DATA PREPARATION #### 

#datamatching = polymorphism database 
#lat-longGABI = GABI database with lat/long coordinates
### 1. Creating data frame with poly_id ###

specieslist<-read.csv("matched-polymorphism_bg.csv",stringsAsFactors = TRUE, header = T)
occ<-read.csv("Lat-Long_Data_GABI.csv",stringsAsFactors = TRUE, header = T)

# match species names and add to a poly_id column
occ_match <- occ[match(paste0(specieslist$genus,".",specieslist$species),
                       occ$valid_species_name),]

occ_match<-na.omit(occ_match)#remove NAs to be able to attribute continent

### 2. Matching Continents Per Occurrence ###

library(rworldmap)
library(rgeos)
library(maptools)
library(cleangeo)  ## For clgeo_Clean()
library(raster)
library(sp)

coords<-as.data.frame(cbind(occ_match$dec_long,occ_match$dec_lat))
colnames(coords)<-c("Lat","Long")
coords<-na.omit(coords)
head(coords)
coords<-SpatialPointsDataFrame(coords=coords,data=coords,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
#convert coordinates into spatial object

#######
###
# if devtools not installed run `install.packages("devtools")´
devtools::install_github("ropenscilabs/rnaturalearthdata")
# adding a programable interface to the natural earth data 
install.packages("rnaturalearthhires",
                 repos = "http://packages.ropensci.org",
                 type = "source")

library(rnaturalearth)
library(sp)

# add the shapefile with world countries directly from naturalearth see: https://cran.r-project.org/web/packages/rnaturalearth/README.html
cont <- ne_countries()  ## get the shapefile
# transforming the coords object to share the same projection as the country shapfile
coords <- spTransform(coords,CRSobj = cont@proj4string )
points<-over(coords,cont)#must remove NAs from dataframe

occ_match$continent <- points$continent[match(occ_match$country,points$name)]


# write.csv(occ_match,file="sites_withcontinents.csv")


### 3. Create density/site ID dataframe ### 

occ_match<-na.omit(occ_match)

new_id<-paste(ceiling(occ_match$dec_long),ceiling(occ_match$dec_lat),sep=".")#ceiling function rounds down

occ_id<-cbind(occ_match,new_id)

id.table<-as.data.frame(table(new_id))
# add point density as a new column of occ_id
occ_id$Freq <- id.table$Freq[match(occ_id$new_id, id.table$new_id)]


# write.csv(occ_id,file="occ_withsitedensity.csv")
#must add in genus column manually in the file 

### 4. Extract Climate ###

library(sp)
library(raster)
library(vegan)


# programatically download worldclim data
WORLDCLIM <-raster::getData("worldclim",var="bio",res=10 )

MAT<-WORLDCLIM[[1]] # BIO_01 mean annual temperature
MAP<-WORLDCLIM[[12]] # BIO_12 precipitation

### Extract variables ###

coords<-as.data.frame(cbind(occ_id$dec_long,occ_id$dec_lat))

coords<-SpatialPointsDataFrame(coords=coords,data=coords,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

MAT.point<-extract(MAT,coords)
MAP.point<-extract(MAP,coords)

### Climate Data Frame ### 

climate<-as.data.frame(cbind(MAP.point,MAT.point))
# write.csv(climate,file="climate_perpoint.csv")
stand.climate<-scale(climate)

### Combine both Datasets ### 

occ_full<-as.data.frame(cbind(occ_id,stand.climate))

head(occ_full)

## separate species name column and reconstruct genus and species columns
gen_sp <- data.frame(stringr::str_split(occ_full$valid_species_name, pattern  = "\\.", simplify = T))
names(gen_sp) <- c("genus", "species")
# bind into the dataset
occ_full <- cbind(occ_full, gen_sp)

##### From here below I don't understand what the code wants, is a bit criptic and convoluted. 

#### 5. FINDING CENTROID ####

#library(geosphere)
#library(dplyr)
#library(ggplot2)

#occ<-read.csv("occ_withsitedensity.csv",stringsAsFactors = FALSE)

#for(z in 1:2) {

#for(i in 1:length(unique(occ$genus))) {
#x<-subset(occ,occ$genus==unique(occ$genus)[i])
# precip5[i,2*z-1] <- mean(x[which(x[,z + 1] >= quantile(x[,z+1], 0.95)),z + 1], na.rm = TRUE)
# precip5[i,2*z] <- mean(x[which(x[,z + 1] <= quantile(x[,z+1], 0.05)),z + 1], na.rm = TRUE)
#}
#}


###############################
# 6. Species per grid cell avg#
###############################

occ<-read.csv("occ_withsitedensity.csv",stringsAsFactors = TRUE)

trial <- aggregate(occ$MAP.point, list(occ$valid_species_name, occ$new_id), mean)
colnames(trial) <- c("valid_species_name", "new_id", "MAP.point")
trial$MAT.point <- aggregate(occ$MAT.point, list(occ$valid_species_name, occ$new_id), mean)$x
trial$poly_id <- aggregate(occ$poly_id, list(occ$valid_species_name, occ$new_id), mean)$x

write.csv(trial,file="new_occ_df.csv")
