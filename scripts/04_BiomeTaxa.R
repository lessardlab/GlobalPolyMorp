#====
# Load ecoregions shapes
ecoRegi <- shapefile("Data/EcoRegions/wwf_terr_ecos.shp")
head(ecoRegi)
ecoRegiEA <- spTransform(ecoRegi, crs.laea)

plot(ecoRegiEA)
# make spatial object
occSP <- SpatialPointsDataFrame(data.frame(occ_new$x_grid,
                                           occ_new$y_grid),
                                data = occ_new, 
                                proj4string =crs.laea)
# crossdata with biomes
biomCros <- over(occSP, ecoRegiEA)
occSP@data <- data.frame(occSP@data,biomCros)
# subset datframe
BiomDat <- occSP@data






PolyPorGen <- function(BiomDat,
                       biome,
                       which = c("propW", "propPol", "propAll")){
  geBiom <- BiomDat[BiomDat$BIOME == biome,]
  geBiom <- geBiom[!duplicated(geBiom$valid_species_name),]
  geBiom <- table(geBiom$genus,geBiom$poly_id)
  
  if(which == "propW"){
    geBiom2 <- geBiom[,2] /(geBiom[,1] + geBiom[,2])
    }
  if(which == "propPol"){
    geBiom2 <-geBiom[,2] 
    #/sum(geBiom[,2])
  }
  if(which == "propAll"){
    geBiom2 <-(geBiom[,1] + geBiom[,2])
    
  }
  geBiom <- sort((geBiom[,1] + geBiom[,2]), T)
  cumsum <- cumsum(geBiom)
  cumsum <- cumsum[cumsum<quantile(cumsum,c(0.3))]
  geBiom2 <-  geBiom2[match(names(cumsum), names(geBiom2))]
  return(  data.frame("names" = names(geBiom2),
                      "polyPor" = geBiom2))
}



agSpOcCounts <- function(BiomDat, biome){
  geBiom <- BiomDat[BiomDat$BIOME == biome,]
  agBio <- aggregate(geBiom$X, by= list(geBiom$valid_species_name,geBiom$genus), length)
  agBio <- aggregate(agBio$x, list(agBio$Group.2), median)
  
  
  
  return(agBio)
  
  
}





unique(BiomDat$BIOME.1)
######### Data for heatmap of species richness per genus and biome
PolyxBiome <- reshape2::melt(lapply(1:13,
                                    function(x) PolyPorGen(BiomDat,x, "propW")))
PolyxBiomeTab <- xtabs(PolyxBiome$value~PolyxBiome$names + PolyxBiome$L1)
PolyxBiomeTab <- PolyxBiomeTab[!rowSums(PolyxBiomeTab)< 0.01,]


PolyxBiome2 <- reshape2::melt(lapply(1:13,
                                    function(x) PolyPorGen(BiomDat,x,"propPol")))
PolyxBiomeTab2 <- xtabs(PolyxBiome2$value~PolyxBiome2$names + PolyxBiome2$L1)
PolyxBiomeTab2 <- PolyxBiomeTab2[!rowSums(PolyxBiomeTab2) < 0.01,]


PolyxBiome3 <- reshape2::melt(lapply(1:13,
                                     function(x) PolyPorGen(BiomDat,x,"propAll")))
PolyxBiomeTab3 <- xtabs(PolyxBiome3$value~PolyxBiome3$names + PolyxBiome3$L1)
PolyxBiomeTab3 <- PolyxBiomeTab3[!rowSums(PolyxBiomeTab3) < 0.01,]

############


###########
# Aggregating by observation counts 
###########


SpOcurxBiome <- reshape2::melt(lapply(1:13,
                                    function(x) agSpOcCounts(BiomDat,x)))
SpOcurxBiome <- xtabs(SpOcurxBiome$value~SpOcurxBiome$Group.1 + SpOcurxBiome$L1)
SpOcurxBiome <- SpOcurxBiome[!rowSums(SpOcurxBiome)< 0.01,]




######### Plot heatmaps 

colnames(PolyxBiomeTab2) <- paste0("BM",colnames(PolyxBiomeTab2))



png("reviewProceedingsB/heatmap.png",
    1500, 1500, 
    pointsize = 10)
pheatmap::pheatmap(log1p(PolyxBiomeTab2),
                   col = cpl(10),angle_col = 45,
                   treeheight_col = 50,
                   cellheight = 20,
                   cellwidth = 40,border_color = NA,
                   fontsize_col = 17,
                   fontsize_row = 17,
                   legend = F,treeheight_row = 300,
                   scale = "none")


dev.off()

colnames(PolyxBiomeTab3) <- paste0("BM",colnames(PolyxBiomeTab3))

png("reviewProceedingsB/hetmapAll.png",
    1500, 2000, 
    pointsize = 20)
pheatmap::pheatmap(log1p(PolyxBiomeTab3),
                   col = cpl(10),angle_col = 45,
                   treeheight_col = 50,
                   cellheight = 13,
                   cellwidth = 40,border_color = NA,
                   fontsize_col = 13,
                   fontsize_row = 10,
                   legend = F,treeheight_row = 300,
                   scale = "none")


dev.off()


colnames(SpOcurxBiome) <- paste0("BM",colnames(SpOcurxBiome))

png("reviewProceedingsB/heatmapOccurrences.png",
    3000, 3000, 
    pointsize = 10)
pheatmap::pheatmap(log1p(SpOcurxBiome),
                   col = cpl(10),angle_col = 45,
                   treeheight_col = 50,
                   cellheight = 9,
                   cellwidth = 40,border_color = NA,
                   fontsize_col = 12,
                   fontsize_row = 10,
                   legend = F,treeheight_row = 300,
                   scale = "none")


dev.off()








png("reviewProceedingsB/heatLegend.png", 500,500, pointsize = 10)
plot(0,0, pch= ".",  cex = 0.1, frame = F, axes = NA)
legend("center",
       cex = 1.2, box.col = "white",ncol = 2,
       legend = c("BM1:Tropical rain forest",
                  "BM2:Tropical dry forest",
                  "BM3:Tropical coniferous forest",
                  "BM4:Temperate mixed forest",
                  "BM5:Temperate conifer forest",
                  "BM6:Boreal forest",
                  "BM7:Tropical savanna",
                  "BM8:Temperate grassland",
                  "BM9:Flooded savanna",
                  "BM10:Páramo",
                  "BM11:Tundra",
                  "BM12:Meditarranean scrubland",
                  "BM13:Desert"))
dev.off()












library(magick)

image_write(image_append(c(image_trim(image_read("reviewProceedingsB/heatmap.png")),
               image_border(
                 image_trim(image_read("reviewProceedingsB/heatLegend.png")),
                 color = "white",geometry = "10x10")),T),
            path = "reviewProceedingsB/richnessBiomeHeatmap.png")
