## GM update on Freddie's base script. 
## DIC 2020

#### DATA PREPARATION #### 

#datamatching = polymorphism database 
#lat-longGABI = GABI database with lat/long coordinates
### 1. Creating data frame with poly_id ###

specieslist<-read.csv("Data/matched-polymorphism_bg_JP_taxocorrect.csv")

unique(specieslist$Subfamily)
length(unique(specieslist$taxon_code))

subFa <- reshape::melt(stringr::str_split(specieslist$Subfamily, "\\t"))$value
specieslist$Subfamily <- subFa[!stringr::str_count(subFa) == 0]

occ<-read.csv("Data/Lat-Long_Data_GABI.csv",stringsAsFactors = FALSE)
occ$valid_species_name <- as.factor(occ$valid_species_name)
occ$poly_id <- 0

specieslist$taxon_code <- paste0(specieslist$genus,".",specieslist$species)

#for(i in 1:nrow(occ)){
#  print(paste(i, "of", nrow(occ)))
#  temp<-subset(specieslist,specieslist$taxon_code==occ$valid_species_name[i])
#  occ$poly_id[i]<-temp$poly_id[1]
#  cat("\014")
#}
############
# match species names and add to a poly_id column


# fast alternative  
occ$poly_id <- specieslist$poly_id[match(occ$valid_species_name,paste0(specieslist$genus,".",specieslist$species))]

occ_match<-occ[,-6] #elevation has too many NAs  
occ_match<-na.omit(occ_match)#remove NAs to be able to attribute continent

new_id<-paste(ceiling(occ_match$dec_long),ceiling(occ_match$dec_lat),sep=".")#ceiling function rounds down



#############################################
### 2. Matching Continents Per Occurrence ###
#############################################

library(rworldmap)
library(rgeos)
library(maptools)
library(cleangeo)  ## For clgeo_Clean()
library(raster)
library(sp)

coords<-as.data.frame(cbind(occ_match$dec_long,occ_match$dec_lat))
colnames(coords)<-c("Long", "Lat")
coords<-na.omit(coords)

#convert coordinates into spatial object
coords<-SpatialPoints(coords=coords,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
st_crs(coords)
# http://rgdal.r-forge.r-project.org/articles/PROJ6_GDAL3.html
cont<-rgdal::readOGR("Data/ne_50m_admin_0_countries.shp") 

proj4string(coords)<- proj4string(cont) 


### Download ocean layer
# URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
# fil <- basename(URL)
# if (!file.exists(fil)) download.file(URL, fil)
# fils <- unzip(fil)
# oceans <- rgdal::readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
# stringsAsFactors=FALSE, verbose=FALSE)
# 
# 
# 

oceans <-rgdal::readOGR("Data/ne_110m_ocean.shp", "ne_110m_ocean",stringsAsFactors=FALSE, verbose=FALSE)


points<-over(coords,cont)#must remove NAs from dataframe
occ_match["Continent"]<-points$CONTINENT

# remove points that lay in the ocean
OcenDAta <- sp::over(coords, oceans)
coords <- coords[which(is.na(OcenDAta$scalerank)),]

# remove Hawaii
clean_coords<-data.frame(coordinates(coords))
clean_coords<-clean_coords[!((round(clean_coords$Lat) %in% seq(18, 28, by = 1)) & (round(clean_coords$Long) %in% seq(-178, -154, by = 1))),]

occ_id<-cbind(occ_match,new_id)
occ_id2<-(occ_id[(occ_id$new_id %in%  paste0(round(clean_coords$Long), ".", round(clean_coords$Lat)) ),])
occ_id2<-occ_id2[complete.cases(occ_id2),]


#plot(occ_id2$dec_lat~occ_id2$dec_long, pch = ".", cex = 1)

write.csv(occ_id2,file="Data/occ_withsitedensity.csv")
#must add in genus column manually in the file 





### 3. Create density/site ID dataframe ### 

id.table<-as.data.frame(table(new_id))
occ_id2$point_density<-NA

for(i in 1:nrow(id.table)){
  occ_id2$point_density[which(occ_id2$new_id==id.table$new_id[i])]<-id.table[i,2]
}

#write.csv(occ_id2,file="Data/occ_withsitedensity.csv")
#must add in genus column manually in the file 








###########################
### 4. Extract Climate ###
###########################

library(sp)
library(raster)
library(vegan)

MAT<-raster("Data/wc2.0_bio_10m_01.tif")
MAP<-raster("Data/wc2.0_bio_10m_12.tif")


### Extract variables ###

coords<-as.data.frame(cbind(occ_id2$dec_long,occ_id2$dec_lat))
coords<-SpatialPointsDataFrame(coords=coords,data=coords,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Note I extract from the globally scaled environment, not "ant distribution" scaled 
MAT.point<-extract(scale(MAT),coords)
MAP.point<-extract(scale(MAP),coords)

### Climate Data Frame ### 

climate<-as.data.frame(cbind(MAP.point,MAT.point))
# write.csv(climate,file="climate_perpoint.csv")

### Combine both Datasets ### 

occ_full<-as.data.frame(cbind(occ_id2, climate))

### Subset all Data to Eliminate Seven Seas as Continent ### 

occ<-subset(occ_full,occ_full$Continent!="Seven seas (open ocean)")
occ$Continent<-factor(occ$Continent)
occ$new_id<-factor(occ$new_id)
occ$genus <- specieslist$genus[match(occ$valid_species_name, paste0(specieslist$genus, ".",specieslist$species))]
occ$subFamily <- specieslist$Subfamily[match(occ$valid_species_name,paste0(specieslist$genus, ".", specieslist$species))]

write.csv(occ,file="Data/occ_withsitedensity.csv") # good file to use 

plot(occ$dec_lat ~ occ$dec_long, pch = ".")


############################
# Fix database large polymorphism 
#####

install.packages("AntWeb")
library(AntWeb)


PlyBG <- readxl::read_xlsx("Data/Polymorphism_BG.xlsx", sheet = 1, skip = 3)
PlyBG2 <- readxl::read_xlsx("Data/Polymorphism_BG.xlsx", sheet = 2, skip = 0)

head(PlyBG2)


PlyBG$SubFamily <- specieslist$Subfamily[match(PlyBG$`Row Labels`, specieslist$genus)]
PlyBG$SubFamily  <- ifelse(is.na(PlyBG$SubFamily), FixFam(PlyBG$`Row Labels`), PlyBG$SubFamily)


PlyBG2$SubFamily <- specieslist$Subfamily[match(PlyBG2$genus, specieslist$genus)]
PlyBG2$SubFamily  <- ifelse(is.na(PlyBG2$SubFamily), FixFam(PlyBG2$genus), PlyBG$SubFamily)


write.csv(PlyBG,file = "Data/Polymorphism_BG_TaxCorrect.csv", row.names = F)
write.csv(PlyBG2,file = "Data/Polymorphism_BG_TaxCorrectb.csv", row.names = F)










