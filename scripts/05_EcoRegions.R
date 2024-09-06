library(raster)


##### Plot the distribution of the interaction between temperature and precipitation globally 
# Results from a PCA-raster between MAP and MAT

PCAras <- RStoolbox::rasterPCA(EnvStack)

PCAras$map$PC1

#====
png("PCA_locallyother.png", 
    width = 1000, 
    height = 1000,pointsize = 5)

par(mar = c(0,0,0,0), oma = c(0,2,2,2))
raster::plot(-(PCAras$map$PC1),
             box = F,
             legend = F,
             axes = F,
             bty="n",
             ylim = c(-60,90),
             breaks = round(quantile(-PCAras$map$PC1,
                                     probs = c(0,0.15,0.3,0.65,1)), 
                            2), 
             col = rev(RColorBrewer::brewer.pal(5, "Greys")))
dev.off()


png("PCA_locallyDrier.png", 
    width = 1000, 
    height = 1000,pointsize = 5)

par(mar = c(0,0,0,0), oma = c(0,2,2,2))

raster::plot((PCAras$map$PC2),
             box = F,
             bty="n",
             axes = F,
             ylim = c(-60,90),
             legend = F,
             breaks = round(quantile(PCAras$map$PC2,
                                     probs = c(0,0.15,0.3,0.65,1)), 2), 
             col = rev(RColorBrewer::brewer.pal(5, "Greys")))

dev.off()

library(magick)


magick::image_write(
  magick::image_append(
    c(image_trim(image_read("PCA_locallyother.png")),
      image_trim(image_read("PCA_locallyDrier.png"))),
    T),
  "FigsPNAS/PCA_MAP_MAT_World.png")
#====


##### Plots the bioregions map and the aggregated probabilities by biome



# Aggregating probabilities by biome at the grid scale
#====
PredVal <- data.frame(replicate(14,1:length(getValues(pred2x))))

for(i in 1:14){
  predBiome <- mask(pred2x,ecoRegi[ecoRegi$BIOME == i,])
  PredVal[,i] <- getValues(predBiome)
  
  
}

## Map the biome geographical distribution into climatic space. 

# CreateSpatialObject with Coordinates from bioclim
poinT <- SpatialPoints(coordinates(MAP), proj4string = crs(ecoRegi))

# Extract biome from points
PointVal <- data.frame(replicate(14,1:length(getValues(pred2x))))

for(i in 1:14){
  print(i)
  predBiome <- mask(MAP,ecoRegi[ecoRegi$BIOME == i,])
  PointVal[,i] <- getValues(predBiome)
  cat("\014")
  
}
# Give appropiate names and reshape the data.frame
names(PointVal) <- 1:14
MAP_Val <- reshape2::melt(PointVal[,-14])

PointVal2 <- data.frame(replicate(14,1:length(getValues(pred2x))))
# Extract values
for(i in 1:14){
  print(i)
  predBiome <- mask(MAT,ecoRegi[ecoRegi$BIOME == i,])
  PointVal2[,i] <- getValues(predBiome)
  cat("\014")
  
}
# proper names and reshape
names(PointVal2) <- 1:14
MAT_Val <- reshape2::melt(PointVal2[,-14])
MAP_Val$MAT <- MAT_Val$value


# make a nice palette that clearly shows biomes 
colpal <- c("#006400",
            "#FFF5D7",
            "#AAC800",
            "#CAFE8F",
            "#9ED003",
            "#008D02",
            "#F7D600",
            "#FFB432",
            "#cacd01", 
            "#8A9000",
            "#AAAAAA", 
            "brown",
            "yellow"  )

## Plot a map of biome in climatic space
#====
png("reviewProceedingsB/EnvMap.png", 
    1820, 
    1820, 
    pointsize = 20)
par(las= 1, mar = c(8,6,0,0))
plot(log(MAP_Val$value)~ MAP_Val$MAT, 
     pch = ".", 
     col = scales::alpha(colpal[MAP_Val$variable], 0.7),
     frame = F,
     cex.lab = 3,
     cex.axis = 2,
     xlab = "",
     xaxt = "n",
     ylab = "Mean annual precipitation (log)"
)

dev.off()
#====


#====
a <- data.frame(coordinates(MAP),
                extract(Stack,coordinates(MAP)))
myOv <- over(SpatialPoints(coordinates(Stack),
                           proj4string = crs(ecoRegi)),
             ecoRegi)
png("FigsPNAS/PredictedNicheLog.png", 
    2000,
    1500, pointsize = 25)
par(mar = c(0,0,0,0))
scatterplot3d::scatterplot3d(x= log1p(a$MAP.point),
                             y = a$MAT.point,
                             z= getValues(predx),
                             pch = ".", 
                             xlab = "MAP (ln)",
                             ylab = "MAT",
                             zlab = "Polymorphism probability",
                             color = scales::alpha(colpal[myOv$BIOME], 0.7),
                             box = F)

dev.off()
#====
range(MAT)
mean(getValues(MAT), na.rm=T)
range(a$MAT.point*8.45, na.rm = T)
#====

png("FigsPNAS/PredictedNicheNorm.png", 
    1500,
    1500, pointsize = 45)
par(las = 2)
scatterplot3d::scatterplot3d(x= (getValues(MAP)),
                             y = (getValues(MAT)),
                             z= getValues(predx),
                             angle = 30,
                             zlim = c(0,0.5),
                             xlim = c(0,10000),
                             ylim = c(-20,30), 
                             pch = ".",scale.y = 1,
                             y.margin.add=0.4,
                             xlab = "",
                             ylab = "",mar = c(5,5,0,0.5),
                             cex.lab = 1.45, cex.axis = 1.4,
                             zlab = "",
                             color = scales::alpha(colpal[myOv$BIOME], 0.4),
                             box = F)
par(new = T, las = 2)
scatterplot3d::scatterplot3d(x= (getValues(MAP)),
                             y = (getValues(MAT)),
                             z= rep(0,length(getValues(predx))),
                             angle = 30,scale.y = 1,
                             zlim = c(0,0.5),
                             xlim = c(0,10000),
                             ylim = c(-20,30),
                             pch = ".",
                             y.margin.add=0.4,
                             xlab = "",
                             ylab = "",mar = c(5,5,0,0.5),
                             cex.lab = 1.45, cex.axis = 1.4,
                             zlab = "",
                             color = "black",
                             cex = 0.3,
                             #color = scales::alpha(colpal[myOv$BIOME], 0.7),
                             box = F)


#text(x = 6, y = 0.3, "MAP", srt = 35, cex = 1.45)


dev.off()


#====

# Plot the map of the global biomes
#====
png("BiomeMap.png",
    3000, 
    1000, 
    pointsize = 20)
par(mar = c(0,0,0,0))
plot(ecoRegi[ecoRegi$BIOME == 1,], col = "#006400", border = "#006400" ,
     xlim = c(-70,70), 
     ylim = c(-70,90)) # Tropical forest
plot(ecoRegi[ecoRegi$BIOME == 2,], col = "#FFF5D7", border = "#FFF5D7", add = T ) # dry forest
plot(ecoRegi[ecoRegi$BIOME == 3,], col = "#AAC800", border = "#AAC800", add = T ) # Coniferous forest tropi
plot(ecoRegi[ecoRegi$BIOME == 4,], col = "#CAFE8F", border = "#CAFE8F", add = T ) # temperate mixed forest
plot(ecoRegi[ecoRegi$BIOME == 5,], col = "#9ED003", border = "#9ED003", add = T ) # temperate conifer
plot(ecoRegi[ecoRegi$BIOME == 6,], col = "#008D02", border = "#008D02", add = T ) # taiga
plot(ecoRegi[ecoRegi$BIOME == 7,], col = "#F7D600", border = "#F7D600", add = T ) # grasslands and savana tropcal
plot(ecoRegi[ecoRegi$BIOME == 8,], col = "#FFB432", border = "#FFB432", add = T ) # temperate grassland
plot(ecoRegi[ecoRegi$BIOME == 9,], col = "#cacd01", border = "#cacd01", add = T ) # flodded grassland savanna
plot(ecoRegi[ecoRegi$BIOME == 10,], col = "#8A9000", border = "#8A9000", add = T ) #montane grasslands (paramo)
plot(ecoRegi[ecoRegi$BIOME == 11,], col = "#AAAAAA", border = "#AAAAAA", add = T ) # tundra
plot(ecoRegi[ecoRegi$BIOME == 12,], col = "brown", border = "brown", add = T ) # #meditarrenan forest
plot(ecoRegi[ecoRegi$BIOME == 13,], col = "yellow", border = "yellow", add = T ) # dessert
maps::map("world", add = T,interior = F, lwd =0.2, color = "grey20")
plot(occSP, pch= ".", add = T, cex = 1.5)
legend(-170,20,
fill = colpal,
cex = 1.2, box.col = "white",
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
#====

# Make a boxplot aggregating the predicted grids by biome

#====
# subset object to remove mangroves
sampDat <- PredVal[,-14]
#sampDat <- ifelse(is.infinite(sampDat), NA, sampDat)
# reorder appropiatly
ord <- order(apply(sampDat, 2, median, na.rm = T))
length(ord)
sam <- sample(1:nrow(PredVal), 100000)
toPlot <- PredVal[,-14]
colnames(toPlot) <- paste0("BM",1:13)
toPlot <- toPlot[,ord]
toPlot <- reshape2::melt(toPlot)
toPlot <- toPlot[-which(is.na(toPlot$value)),]

png("FigsPNAS/BoxplotBiome.png",
    1500, 
    1500, 
    pointsize = 60)
par( mar = c(6,6,1,2))

boxplot(toPlot$value~toPlot$variable,
        frame = F, 
        horizontal = F,
        ylim = c(0,0.5),
        col = colpal[ord],
        cex.axis = 1.1,
        yaxt ="n",
        xaxt = "n",
        xlab = c("Biome class","Precipitation (mm)", "Temperature ºC"),
        ylab = c("Probability of polymorphism"),
        cex.lab = 1.3)
axis(1, 0:13, c("",as.character(unique(toPlot$variable))), 
     las = 1, cex.axis = 0.8,  tck = 0.01, las =  2)
axis(1, 0:13, c("",as.character(unique(toPlot$variable))), 
     las = 2, cex.axis = 0.8,  tck = -0.01)
axis(2, seq(-0.1,0.5,0.1),seq(-0.1,0.5,0.1), 
     las = 1,
     cex.axis = 1.2,tck = -0.01)
axis(2, seq(-0.1,0.5,0.1),seq(-0.1,0.5,0.1), 
     las = 1, cex.axis = 1.2,tck = 0.01)

dev.off()
#====

# Composite the final biome figure

#====
IMG1 <- image_append(
  c(
    image_border(image_scale(
      image_border(image_read("FigsPNAS/BoxplotBiome.png"),
                 color = "white",
                 geometry = "100x100"), "1500"), 
      color = "white",
      geometry = "150x0"),
    image_border(
      image_trim(image_read("FigsPNAS/PredictedNicheNorm.png")),
      color = "white",
      geometry = "100x100")),
  stack = F)



image_write(
  image_append(c(
    image_scale(image_crop(image_trim(image_read("BiomeMap.png")),
                           geometry = "3000x830"),geometry = "3200"),
    IMG1), stack = T), 
  "FigsPNAS/BiomeMap.png")


#====

# Plot predicted probabilities 

#====
f1 <- focal(x = pred2x,w = matrix(1,5,5), fun = median)


png("FigsPNAS/AllSpecies_YlOrRd.png",
     width = 1200, 
     height = 1000, pointsize = 20)

plot((f1),
      box = F,
      bty="n",
      axes=F,
      ylim = c(-60,70),
      breaks = c(0.06,0.1,0.15,0.20,0.23,0.26,0.35,0.4,0.5,1),
      col = (scales::alpha(RColorBrewer::brewer.pal(9, "YlOrRd"))))
maps::map("world", regions = "Greenland",
          col = "white", add = T, fill = T, border="grey80")

maps::map("world", col = "black", lwd = 0.2,
          add = T, fill = F, interior = F, ylim = c(-60,70))

dev.off()

image_write(image_trim(image_read("FigsPNAS/AllSpecies_YlOrRd.tiff")), 
            "FigsPNAS/AllSpecies_YlOrRd.png")



png("FigsPNAS/LegFig2.png",
    200,200, pointsize = 20)
par(las = 1, mar = c(0,0,0,0))
plot(c(0.06,0.1,0.15,0.20,0.23,0.26,0.35,0.4,0.5)~rep(1,9),
     pch = 15, 
     cex = 5,
     frame = F, xlim = c(0.9,1.1),
     xaxt  = "n", yaxt = "n", xlab = "", ylab = "",
     col = RColorBrewer::brewer.pal(9, "YlOrRd"))
axis(2,pos = 1)
dev.off()


image_write(
  image_composite(
    image_read("FigsPNAS/AllSpecies_YlOrRd.tiff"),
    image_trim(image_read("FigsPNAS/LegFig2.png")),
    offset = "-100-200"), "FigsPNAS/Figure2.png", format = "png")

#====

a <- c()

# Richness per family
# make spatial object of richness per grid distribution by family
#====
subSetFam <-function(occ_new, fam){
  occSP2 <- SpatialPointsDataFrame(data.frame(round(occ_new$dec_long), round(occ_new$dec_lat)),
                                   data = occ_new, proj4string =crs(ecoRegi) )
  corF <- occ_new[occ_new$subFamily == fam & occ_new$poly_id == 1,]
  
  
  
  corF$ID <- paste0(round(corF$dec_long),"_",round(corF$dec_lat))
  corF <- aggregate(corF$valid_species_name, by = list(corF$ID), length)
  
  spcorF <- SpatialPointsDataFrame(data.frame(
    apply(stringr::str_split(corF$Group.1, "_" ,
                             simplify = T), 2, as.numeric)+1/2),
    data.frame(log(corF$x)))
  return(spcorF)
  
  
  
}

spcorF <- subSetFam(occ_new, "Formicinae")
spcorM <- subSetFam(occ_new, "Myrmicinae")
spcorD <- subSetFam(occ_new, "Dolichoderinae")
#====

# plot results
#====
png("RichFormicinae.png",
    1000, 
    500, 
    pointsize = 40)
par(mar = c(0,0,0,0))
maps::map("world", interior = F, border =F,
          col = scales::alpha("black",0.5),
          fill = T, bg = "white")
image(rasterFromXYZ(spcorF),
      col =  rev(RColorBrewer::brewer.pal(9, "PuOr")), 
      frame= F, 
      xaxt = "n", add = T, 
      yaxt = "n")
dev.off()


png("RichMyrmicinae.png",
    1000, 
    500, 
    pointsize = 40)
par(mar = c(0,0,0,0))

maps::map("world", 
          interior = F, border =F, 
          col = scales::alpha("black", 0.5),fill = T, bg = "white")
image(rasterFromXYZ(spcorM),
      col =   rev(RColorBrewer::brewer.pal(9, "PuOr")), 
      frame= F, 
      xaxt = "n", add = T, 
      yaxt = "n")
dev.off()

png("RichDolicho.png",
    1000, 
    500, 
    pointsize = 40)
par(mar = c(0,0,0,0))

maps::map("world", interior = F, 
          border =F, col = scales::alpha("black",0.5) ,
          fill = T, bg = "white")
image(rasterFromXYZ(spcorD),
      col = rev(RColorBrewer::brewer.pal(9, "PuOr")),
      frame= F, 
      xaxt = "n", add = T, 
      yaxt = "n")
dev.off()
#====

# Composite image
#====
image_write(
  image_append(
    c(image_crop(image_trim(image_read("RichFormicinae.png")), "600x250"),
      image_crop(image_trim(image_read("RichMyrmicinae.png")), "600x250"),
      image_crop(image_trim(image_read("RichDolicho.png")), "600x250")), 
    stack = T), "FigsPNAS/RichnessAll.png")
#====

