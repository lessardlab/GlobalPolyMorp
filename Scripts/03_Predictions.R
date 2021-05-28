spamodel4



#### All species 

plot(MAP)
Stack <- raster::stack(scale(MAP), scale(MAT))
Stack <- mask(Stack,
              wrld_simpl[wrld_simpl@data$ISO3 != "ATA",])


names(Stack) <- c("MAP.point", "MAT.point")

predx <- predict(Stack, 
                spamodel3,
                re.form=~0, 
                type = "response")
pred2x <- predx


###### Normalize environment
MAT <- mask(MAT,
            wrld_simpl[wrld_simpl@data$ISO3 != "ATA",])

MAP <- mask(MAP,
            wrld_simpl[wrld_simpl@data$ISO3 != "ATA",])


MATminmax <- range(getValues(MAT), na.rm = T)

Stack2 <- raster::stack(MAP/max(getValues(MAP), na.rm = T), 
                        (MAT-MATminmax[1])/(MATminmax[2]-MATminmax[1]))

raster::res(MAP)
## build a ocean raster
Or <- raster(oceans, res = res(MAP))
Or <- rasterize(oceans, Or)

## artic raster

artic <- raster::raster("Data/RasterCAVM/raster_cavm_v1.tif")

plot(artic )


wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
artic_wgs84 = projectRaster(artic, crs = wgs84, method = "ngb")
hist(artic_wgs84)



artic_wgs84_2 <- extend(artic_wgs84, pred2)
artic_wgs84_2 <- crop(artic_wgs84_2, pred2)
extent(artic_wgs84_2) <- round(extent(artic_wgs84_2))-c(-1,1,0,0)

art <- projectRaster(artic,pred2)

plot(art)

art[is.na(art),] <- 80
art[art<79] <- NA
art[art>79,] <- 1



pred2 <- mask(pred2,art)




png("AllSpecies_RdYlBu.png",width = 700, 
    height = 500, 
    pointsize = 12)
pdf("AllSpecies_YlOrRd.pdf",width = 200, 
    height = 100, 
    pointsize = 200)
tiff("AllSpecies_YlOrRd.tiff",width = 1500, 
     height = 700, pointsize = 20)
raster::plot(Or, 
             box = F,
             bty="n",
             axes = F,
             legend = F,
             ylim = c(-60,70),
             col = "white")
f1 <- focal(x = pred2x,w = matrix(1,5,5), fun = median)

raster::plot(f1,
             box = F,
             bty="n",
             axes=F,
             add = T,
             ylim = c(-60,70),
             breaks = round(quantile(f1),2),
             # col = c("#335c67","#fff3b0","#e09f3e",
             #         "#9e2a2b","#540b0e"),
             col = scales::alpha((RColorBrewer::brewer.pal(4, "YlOrRd"))))
maps::map("world", regions = "Greenland",
          col = "white", add = T, fill = T, border="grey80")

maps::map("world", col = "black", lwd = 0.5,
          add = T, fill = F, interior = F, ylim = c(-60,70))

dev.off()

image_trim(image_read("AllSpecies_YlOrRd.tiff"))


#### By Subfamilies 


unique(occ_new$SubFamily)
myrT <- occ_new[occ_new$SubFamily == "Myrmicinae",]
plot(myrT$poly_id~myrT$MAP.point)



summary(spamodel4.a[[9]])

# Model 4
for(i in c(1:12)){
  
   if(unique(occ_new$SubFamily)[i] %in% c("Dorylinae", "Myrmicinae", "Formicinae", "Dolichoderinae") == F){
     next
   }
  
  i = 2
  if(length(summary(spamodel4.a[[i]])) == 18){
    nam <- paste(c("FigSubFam/SubFamilyProb", 
                   unique(occ_new$SubFamily)[i], 
                   ".pdf"), collapse = "")
    print(nam)
    
    
    Form_pred <- predict(Stack, 
                         spamodel4.a[[i]],
                         re.form=~0, 
                         type = "response")
    Form_pred2 <- Form_pred
    Form_pred2 <- mask(Form_pred2,art)
    
    
    
    pdf(nam,
        width = 200, 
        height = 100, 
        pointsize = 200)
    #par(mfrow=c(1,1), mar = c(4,0,4,4), oma = c(2,2,2,2))
    
    ########
    raster::plot(Or, 
                 box = F,
                 bty="n",
                 axes = F,
                 legend = F,
                 ylim = c(-60,70),
                 col = "white")
    f1 <- focal(x = Form_pred2,w = matrix(1,5,5), fun = median)
    
    raster::plot(f1,
                 box = F,
                 bty="n",
                 axes=F,
                 add = T,
                 ylim = c(-60,70),
                 breaks = seq(0,1,length.out = 9),
                 # col = c("#335c67","#fff3b0","#e09f3e",
                 #         "#9e2a2b","#540b0e"),
                 col = scales::alpha((RColorBrewer::brewer.pal(9, "PuBu"))))
    maps::map("world", regions = "Greenland",
              col = "white", add = T, fill = T, border="grey90")
    
    maps::map("world", col = "black", border = "grey20", lwd = 0.5,
              add = T, fill = F, interior = F, ylim = c(-60,70))
    
    legend("left",   unique(occ_new$SubFamily)[i], bty = "n", text.font = 2)
    dev.off()
    
    cat("\014")
  }
  
}


# Model 3

for(i in c(1,3,)){
  
  if(unique(occ_new$SubFamily)[i] %in% c("Dorylinae", "Myrmicinae", "Formicinae", "Dolichoderinae") == F){
   next
     }
  
  
  if(length(summary(spamodel3.a[[i]])) == 18){
    nam <- paste(c("FigSubFam/SubFamilyProb", 
                   unique(occ_new$SubFamily)[i], 
                   ".pdf"), collapse = "")
    print(nam)
    
    
    Form_pred <- predict(Stack, 
                         spamodel3.a[[i]],
                         re.form=~0, 
                         type = "response")
    Form_pred2 <- Form_pred
    Form_pred2 <- mask(Form_pred2,art)
    
    
    
    pdf(nam,
        width = 200, 
        height = 100, 
        pointsize = 200)
    #par(mfrow=c(1,1), mar = c(4,0,4,4), oma = c(2,2,2,2))
    
    ########
    raster::plot(Or, 
                 box = F,
                 bty="n",
                 axes = F,
                 legend = F,
                 ylim = c(-60,70),
                 col = "white")
    f1 <- focal(x = Form_pred2,w = matrix(1,5,5), fun = median)
    
    raster::plot(f1,
                 box = F,
                 bty="n",
                 axes=F,
                 add = T,
                 ylim = c(-60,70),
                 breaks = seq(0,1,length.out = 9),
                 # col = c("#335c67","#fff3b0","#e09f3e",
                 #         "#9e2a2b","#540b0e"),
                 col = scales::alpha((RColorBrewer::brewer.pal(9, "PuBu"))))
    maps::map("world", regions = "Greenland",
              col = "white", add = T, fill = T, border="grey90")
    
    maps::map("world", col = "black", border = "grey20", lwd = 0.5,
              add = T, fill = F, interior = F, ylim = c(-60,70))
    
    legend("left",   unique(occ_new$SubFamily)[i], bty = "n", text.font = 2)
    dev.off()
    
    cat("\014")
  }
  
}








library(magick)



image_write(
  magick::image_append(
    c(image_border(image_trim(image_read("FigSubFam/PAST/SubfamilyProbDolichoderinae.tiff")),
                   "white", "10x10"),
      image_border(image_trim(image_read("FigSubFam/PAST/SubfamilyProbFormicinae.tiff")), 
                   "white", "10x10"),
      image_border(image_trim(image_read("FigSubFam/PAST/SubfamilyProbMyrmicinae.tiff")),
                   "white","10x10")),
    stack = T),
  "FigSubFam/FullImageSubFamMod3.tiff")










image_write(
  magick::image_append(
      c(image_border(image_trim(image_read("FigSubFam/Mod3a.png")),
                     "white", "10x10"),
        image_border(image_trim(image_read("FigSubFam/Mod3b.png")), 
                       "white", "10x10"),
          image_border(image_trim(image_read("FigSubFam/Mod3c.png")),
                       "white","10x10")),
        stack = T),
  "FigSubFam/FullImageSubFamMod3.tiff")




image_write(
  magick::image_append(
    c(image_border(image_trim(image_read("FigSubFam/Mod4a.png")),
                   "white", "10x10"),
      image_border(image_trim(image_read("FigSubFam/Mod4b.png")), 
                   "white", "10x10"),
      image_border(image_trim(image_read("FigSubFam/Mod4c.png")),
                   "white","10x10")),
    stack = T),
  "FigSubFam/FullImageSubFamMod4.tiff")










Mod3 <- image_read("FigSubFam/FullImageSubFamMod3.tiff")

Mod4 <- image_read("FigSubFam/FullImageSubFamMod4.tiff")

image_append(c(Mod3, Mod4))





img <- image_composite(img, 
                       image_scale(image_read("FigSubFam/Linephitema.png"), 
                                   geometry = "200x200"), 
                       offset = "+30-220")

img <- image_composite(img, 
                       image_transparent(image_rotate(image_scale(image_read("FigSubFam/Formica.png"), 
                                   geometry = "200x200"),30), "white"),
                       offset = "+30-610")

img <- image_composite(img, 
                       image_scale(image_read("FigSubFam/Messor.jpeg"), 
                                   geometry = "200x200"), 
                       offset = "+20-1080")

image_write(img, path = "FullMapImage.tiff", format = "tiff")
img2 <- image_read("FullMapImage.tiff")

pdf("Map.pdf", width = 125, height = 250, pointsize = 20)
par(oma = c(0,0,0,0), mar = c(0,0,0,0))
plot(img2)
dev.off()


PCA <- RStoolbox::rasterPCA(Stack2, spca = T)

world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")


axis <- PCA$map$PC1

PCA$model$loadings

png("PCA_locallyother.png",width = 700, 
    height = 500, 
    pointsize = 10)
pdf("MAPS_high.pdf", 
    width = 175*0.50, 
    height = 125*0.50, 
    pointsize = 12
)
tiff("MAPS_high.tiff", 
     width = 175, units = "cm",
     height = 125,res = 140, pointsize = 5)

par(mfrow = c(1,2) , mar = c(0,0,0,0), oma = c(0,2,2,2))
raster::plot(-(axis),
             box = F,
             legend = F,
             axes = F,
             bty="n",
             ylim = c(-60,70),
             breaks = round(quantile(-axis,
                                     probs = c(0,0.15,0.22,0.4,0.6,0.9,1)), 
                            2), 
             col = rev(RColorBrewer::brewer.pal(6, "Greys")))
# maps::map("world", 
#           add = T,
#           ylim = c(-60,70),
#           interior = F)
# legend(-170,0, 
#        cex = 1.5,
#        legend = rev(c("Colder + Drier", rep("",4),"Hotter + Humid")), 
#        bty = "n",y.intersp =0.7,
#        fill =rev(RColorBrewer::brewer.pal(7, "Greys") ))

# mtext("A", 3, outer = T, cex = 30, adj = 0, line = -7)




# axis2 <- PCA$map$PC2
# 
# png("PCA_locallyDrier.png",width = 700, 
#     height = 500, 
#     pointsize = 10)
# pdf("PCA_locallyDrier.pdf", 
#     width = 175*0.75, 
#     height = 125*0.75, 
#     pointsize = 10
# )
# par(mar = c(0,0,0,0), oma = c(0,0,2,2))

raster::plot((axis2),
             box = F,
             bty="n",
             axes = F,
             ylim = c(-60,70),
             legend = F,
             breaks = round(quantile(axis2, probs = c(0,0.15,0.22,0.4,0.6,0.9,1)), 2), 
             col = rev(RColorBrewer::brewer.pal(7, "Greys")))

# maps::map("world", 
#           add = T,
# #           ylim = c(-60,70),
# #           interior = F)
# legend(-170,0,
#        legend = rev(c("Colder + Humid", rep("",4),"Hotter + Drier")), 
#        bty = "n",y.intersp =0.7, cex = 1.5,
#        fill =rev(RColorBrewer::brewer.pal(7, "Greys") ))
#mtext("B", 3, outer = T, cex = 30, adj = 0, line = -7)

dev.off()

library(magick)
magick::image_write(
  magick::image_append(
    c(image_read("PCA_locallyother.png"),
      image_read("PCA_locallyDrier.png"))),
  "PCA_Axis.World.png")




magick::image_write(
  magick::image_append(
    c(image_read_pdf("PCA_locallyother.pdf"),
      image_read_pdf("PCA_locallyDrier.pdf"))),
  "PCA_Axis.World.pdf")


