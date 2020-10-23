### Functions

### custom function to colorize https://gist.github.com/fgabriel1891
f <- function(x,n=10, pal, rev = F){
  if(rev == F){ 
    rev(RColorBrewer::brewer.pal(n, pal))[cut(x,n)]
  }else{
    (RColorBrewer::brewer.pal(n, pal))[cut(x,n)]
  }
}


DataCleaner <- function(occ, specieslist){
  
  # match species names and add to a poly_id column
  occ$poly_id <- specieslist$poly_id..0.monomorphic..1.polymorphic.[match(occ$valid_species_name,
                                                                          paste0(specieslist$genus,".",specieslist$species))]
  
  
  
  occ_match<-na.omit(occ)#remove NAs to be able to attribute continent
  
  ### 2. Matching Continents Per Occurrence ###
  
  
  coords<-as.data.frame(cbind(occ_match$dec_long,occ_match$dec_lat))
  colnames(coords)<-c("Lat","Long")
  coords<-na.omit(coords)
  head(coords)
  coords<-sp::SpatialPointsDataFrame(coords=coords,data=coords,proj4string = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  #convert coordinates into spatial object
  
  #######
  ###
  # # if devtools not installed run `install.packages("devtools")´
  # devtools::install_github("ropenscilabs/rnaturalearthdata")
  # # adding a programable interface to the natural earth data
  # install.packages("rnaturalearthhires",
  #                  repos = "http://packages.ropensci.org",
  #                  type = "source")
  # devtools::install_github("ropenscilabs/rnaturalearth")
  
  
  # add the shapefile with world countries directly from naturalearth see: https://cran.r-project.org/web/packages/rnaturalearth/README.html
  cont <- ne_countries()  ## get the shapefile of countries
  
  ### Substract data from oceans
  # download ocean file from naturalearth
  URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
  fil <- basename(URL)
  if (!file.exists(fil)) download.file(URL, fil)
  fils <- unzip(fil)
  oceans <- rgdal::readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
                           stringsAsFactors=FALSE, verbose=FALSE)
  
  # plot the global spatial distribution of data
  raster::plot(oceans, col = "blue")
  points(coords, col = "red", pch = ".")
  
  # remove points that lay in the ocean
  OcenDAta <- sp::over(coords, oceans)
  coords <- coords[which(is.na(OcenDAta$scalerank)),]
  
  # check if we remove correctly
  raster::plot(oceans, col = "blue")
  points(coords, col = "red", pch = ".")
  
  off <- over(coords, cont)
  off[c("region_un", "continent"),]
  
  coords$continent <-  off$continent
  coords$region_un <-  off$region_un
  # add unique id 
  occ_match$uniID <-paste0(occ_match$dec_lat,"_",occ_match$dec_long)
  
  # remove points in the ocean from the occ dataset
  occ_match <- occ_match[match(paste0(coords@coords[,2],"_",coords@coords[,1]), occ_match$uniID),]
  
  # add continent as variable
  occ_match$continent <- coords$continent[match(occ_match$uniID,paste0(coords@coords[,2],"_",coords@coords[,1]))]
  # add region as variable
  occ_match$region_un <- coords$region_un[match(occ_match$uniID,paste0(coords@coords[,2],"_",coords@coords[,1]))]
 
  # check if we remove correctly
  raster::plot(oceans, col = "blue")
  points(occ_match$dec_lat~occ_match$dec_long, col = "red", pch = ".")
  
  # visualize data structure
  data.frame(coords)
  print(occ_match)
  
  # write.csv(occ_match,file="sites_withcontinents.csv")
  
  ### 3. Create density/site ID dataframe ### 
  
  
  # remove na
  occ_clean<-na.omit(occ_match)
  # make a id vector
  occ_clean$new_id<-paste(ceiling(occ_clean$dec_long),ceiling(occ_clean$dec_lat),sep=".")#ceiling function rounds down
  id.table<-as.data.frame(table(occ_clean$new_id))
  # add point density as a new column of occ_id
  occ_clean$Freq <- id.table$Freq[match(occ_clean$new_id, id.table$Var1)]
  # scale frequency to 0-1
  occ_clean$scalFreq <- occ_clean$Freq/max(occ_clean$Freq)
  
  return(occ_clean)
}


WorldClimMatch <- function(occ_clean){
  ### 4. Extract Climate ###
  
  
  # programatically download worldclim data
  WORLDCLIM <-raster::getData("worldclim",var="bio",res=10 )
  
  MAT<-WORLDCLIM[[1]] # BIO_01 mean annual temperature
  MAP<-WORLDCLIM[[12]] # BIO_12 precipitation
  
  ### Extract variables ###
  
  coords<-as.data.frame(cbind(occ_clean$dec_long,occ_clean$dec_lat))
  coords<-SpatialPointsDataFrame(coords=coords,data=coords,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  plot(MAT.point)
  MAT.point<-extract(MAT,coords)
  MAP.point<-extract(MAP,coords)
  
  plot(coords, pch = "." ,col = f(scale(extract( MAT,coords)), 10, "Spectral"))
  
  ### Observe in relation to climatic variables
  
  #png(filename = "ClimMaps.png", width = 1000, height = 600, pointsize = 20)
  par(mfrow = c(2,3), mar = c(1,1,1,1))
  plot(MAP)
  plot(MAP)
  points(occ_clean$dec_lat[occ_clean$poly_id == 1]~occ_clean$dec_long[occ_clean$poly_id == 1], 
         pch = ".", cex = 0.2, 
         col = scales::alpha(f(log(occ_clean$scalFreq[occ_clean$poly_id == 1]), 10, "Spectral"),0.5))
  
  plot(MAP)
  points(occ_clean$dec_lat[occ_clean$poly_id == 0]~occ_clean$dec_long[occ_clean$poly_id == 0], 
         pch = ".", cex = 0.2, 
         col = scales::alpha(f(log(occ_clean$scalFreq[occ_clean$poly_id == 1]), 10, "Spectral"),0.5))
  
  #######
  plot(MAT)
  plot(MAT)
  points(occ_clean$dec_lat[occ_clean$poly_id == 1]~occ_clean$dec_long[occ_clean$poly_id == 1], 
         pch = ".", cex = 0.2, 
         col = scales::alpha(f(log(occ_clean$scalFreq[occ_clean$poly_id == 1]), 10, "Spectral"),0.5))
  
  plot(MAT)
  points(occ_clean$dec_lat[occ_clean$poly_id == 0]~occ_clean$dec_long[occ_clean$poly_id == 0], 
         pch = ".", cex = 0.2, 
         col = scales::alpha(f(log(occ_clean$scalFreq[occ_clean$poly_id == 1]), 10, "Spectral"),0.5))
  
  dev.off()
  
  ### Climate Data Frame ### 
  
  climate<-as.data.frame(cbind(MAP.point,MAT.point))
  # write.csv(climate,file="climate_perpoint.csv")
  stand.climate<-scale(climate)
  
  ### Combine both Datasets ### 
  
  occ_full<-as.data.frame(cbind(occ_clean,stand.climate))
  return(occ_full)
}

ModelSelection <- function(occ_full){
  par(mfrow = c(2,2), mar = c(5,5,5,5))
  plot(occ_full$MAP.point~log1p(occ_full$scalFreq), 
       col = ifelse(occ_full$poly_id == 1, "firebrick", "skyblue"),
       pch = ".", cex = 2, 
       xlab = "Frequency (log)",
       ylab = "Mean Annual Precipitation")
  
  
  plot(occ_full$MAP.point,occ_full$MAT.point, 
       col = scales::alpha(ifelse(occ_full$poly_id == 1, "firebrick", "skyblue"),0.5),
       pch = ".", cex = log1p(occ_full$scalFreq),
       xlab = "Mean Annual Temperature",
       ylab = "Mean Annual Precipitation")
  
  
  plot(occ_full$MAP.point[occ_full$poly_id == 0],occ_full$MAT.point[occ_full$poly_id == 0], 
       col = "skyblue",
       pch = ".", cex = log1p(occ_full$scalFreq),
       xlab = "Mean Annual Temperature Poly = 0",
       ylab = "Mean Annual Precipitation Poly = 0")
  
  plot(occ_full$MAP.point[occ_full$poly_id == 1],occ_full$MAT.point[occ_full$poly_id == 1], 
       col = "firebrick",
       pch = ".", cex = log1p(occ_full$scalFreq),
       xlab = "Mean Annual Temperature Poly = 1",
       ylab = "Mean Annual Precipitation Poly = 1")
  
  
  
  
  
  ###############################
  # 2. Global models
  ###############################
  occ_new <- occ_full[-which(is.na(occ_full$MAT.point)),] # remove non matching points 
print("nullmodel")  
  nullmodel<-glm(poly_id~1,data=occ_new,family="binomial"(link=logit))
  print("spamodel1")  
  spamodel1<-glmer(poly_id~1+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("spamodel2")  
  spamodel2<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("spamodel3")
  spamodel3<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("spamodel4")
  spamodel4<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  
  ###############################
  # 3. Global models extra
  ###############################
  
  print("MATmodel")
  MATmodel<-glmer(poly_id~MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("MAPmodel")
  MAPmodel<-glmer(poly_id~MAP.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("MAP2model")
  MAP2model<-glmer(poly_id~I(MAP.point^2)+MAP.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("MAT2model")
  MAT2model<-glmer(poly_id~I(MAT.point^2)+MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  
  print("AIC")
  AIC <- AIC(nullmodel,spamodel1,spamodel2,spamodel3,spamodel4, MATmodel,MAPmodel,MAP2model,MAT2model)
  AIC <- AIC[order(AIC$AIC, decreasing = F),]
  AIC2 <- AIC(nullmodel, spamodel1, spamodel2, spamodel3, spamodel4)
  AIC2 <- AIC2[order(AIC2$AIC, decreasing = F),]
  print("Done")
  myMod <- list(nullmodel,spamodel1,spamodel2,spamodel3,spamodel4, MATmodel,MAPmodel,MAP2model,MAT2model, AIC)
  cat("\014")
  return(myMod)
  
}

plotModels <- function(myModelx){
  a<-plot_model(myModelx,type="pred",terms = "MAT.point [all]")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(family="Times New Roman", face="bold", size=20))+
    xlab("Temperature")+
    ylab("Predicted Probabilities of Worker Caste Polymorphism")+
    ggtitle("")+
    scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),
                       labels=scales::percent_format(accuracy = 1),limits=c(0,0.7));a
  
  b<-plot_model(myModelx,type="pred",terms = "MAP.point [all]")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(family="Times New Roman", face="bold", size=20))+
    xlab("Precipitation")+
    ylab("Predicted Probabilities of Worker Caste Polymorphism")+
    ggtitle("")+
    scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels=scales::percent_format(accuracy = 1), limits = c(0,NA));b
  
  ggarrange(a,b)
}
