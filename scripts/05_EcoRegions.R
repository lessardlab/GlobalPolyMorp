# Load necessary library for raster data manipulation
library(raster)

##### Plot the distribution of the interaction between temperature and precipitation globally 
# Perform PCA (Principal Component Analysis) on a raster stack (EnvStack) containing climate variables (MAP and MAT)
PCAras <- RStoolbox::rasterPCA(EnvStack)

# Extract the first principal component (PC1) from the PCA result
PCAras$map$PC1

#====
# Save the plot of PC1 to a PNG file
png("PCA_locallyother.png", 
    width = 1000, 
    height = 1000, pointsize = 5)

# Set plot margins and outer margins for cleaner visualization
par(mar = c(0,0,0,0), oma = c(0,2,2,2))

# Plot the negative of PC1 without axis or legend, using a grayscale palette and quantile breaks
raster::plot(-(PCAras$map$PC1),
             box = F,
             legend = F,
             axes = F,
             bty="n",
             ylim = c(-60,90),  # Set the latitude limits from -60 to 90
             breaks = round(quantile(-PCAras$map$PC1,  # Use quantile-based breaks for better contrast
                                     probs = c(0, 0.15, 0.3, 0.65, 1)), 
                            2), 
             col = rev(RColorBrewer::brewer.pal(5, "Greys")))  # Use a grayscale color palette
dev.off()

# Save the plot of PC2 to a different PNG file
png("PCA_locallyDrier.png", 
    width = 1000, 
    height = 1000, pointsize = 5)

# Set plot margins and outer margins for cleaner visualization
par(mar = c(0,0,0,0), oma = c(0,2,2,2))

# Plot PC2, similar to PC1, without axis or legend
raster::plot((PCAras$map$PC2),
             box = F,
             bty="n",
             axes = F,
             ylim = c(-60,90),  # Set the latitude limits from -60 to 90
             legend = F,
             breaks = round(quantile(PCAras$map$PC2,  # Use quantile-based breaks for better contrast
                                     probs = c(0, 0.15, 0.3, 0.65, 1)), 2), 
             col = rev(RColorBrewer::brewer.pal(5, "Greys")))  # Use a grayscale color palette
dev.off()

# Load magick library to handle image processing
library(magick)

# Combine the two PCA maps (PC1 and PC2) into a single image and save it as a PNG file
magick::image_write(
  magick::image_append(
    c(image_trim(image_read("PCA_locallyother.png")),  # Trim and read the first PCA image
      image_trim(image_read("PCA_locallyDrier.png"))),  # Trim and read the second PCA image
    T),  # Stack the images horizontally
  "FigsPNAS/PCA_MAP_MAT_World.png")
#====

##### Plots the bioregions map and the aggregated probabilities by biome

# Aggregate probabilities for different biomes on a grid scale
# Create a data frame to store values from 14 biomes (1 column per biome)
PredVal <- data.frame(replicate(14, 1:length(getValues(pred2x))))

# Loop over all 14 biomes, mask the raster data by biome, and extract the values
for(i in 1:14){
  predBiome <- mask(pred2x, ecoRegi[ecoRegi$BIOME == i,])  # Mask the raster by the current biome
  PredVal[,i] <- getValues(predBiome)  # Store the extracted values for the biome
}

## Map the biome's geographical distribution into climatic space

# Create spatial points from the coordinates of the MAP raster and assign the projection system
poinT <- SpatialPoints(coordinates(MAP), proj4string = crs(ecoRegi))

# Create a data frame to store values from MAP for each biome
PointVal <- data.frame(replicate(14, 1:length(getValues(pred2x))))

# Loop over all biomes, mask the MAP raster by biome, and extract the values
for(i in 1:14){
  print(i)  # Print progress for debugging
  predBiome <- mask(MAP, ecoRegi[ecoRegi$BIOME == i,])  # Mask MAP raster by the current biome
  PointVal[,i] <- getValues(predBiome)  # Store the extracted values for the biome
  cat("\014")  # Clear the console output
}

# Assign names to columns and reshape the data frame for MAP values
names(PointVal) <- 1:14
MAP_Val <- reshape2::melt(PointVal[,-14])  # Reshape the data, excluding the 14th column

# Create a data frame to store values from MAT for each biome
PointVal2 <- data.frame(replicate(14, 1:length(getValues(pred2x))))

# Extract values for each biome from the MAT raster
for(i in 1:14){
  print(i)  # Print progress for debugging
  predBiome <- mask(MAT, ecoRegi[ecoRegi$BIOME == i,])  # Mask MAT raster by the current biome
  PointVal2[,i] <- getValues(predBiome)  # Store the extracted values for the biome
  cat("\014")  # Clear the console output
}

# Assign names and reshape the data frame for MAT values
names(PointVal2) <- 1:14
MAT_Val <- reshape2::melt(PointVal2[,-14])  # Reshape the data, excluding the 14th column
MAP_Val$MAT <- MAT_Val$value  # Add MAT values to the MAP data frame

# Create a color palette for the 13 biomes to use in the plots
colpal <- c("#006400",  # Tropical rainforest
            "#FFF5D7",  # Tropical dry forest
            "#AAC800",  # Tropical coniferous forest
            "#CAFE8F",  # Temperate mixed forest
            "#9ED003",  # Temperate conifer forest
            "#008D02",  # Boreal forest (taiga)
            "#F7D600",  # Tropical savanna
            "#FFB432",  # Temperate grassland
            "#cacd01",  # Flooded grassland
            "#8A9000",  # Montane grassland (Páramo)
            "#AAAAAA",  # Tundra
            "brown",    # Mediterranean scrubland
            "yellow")   # Desert


## Plot a map of biome in climatic space
#====
# Create a PNG file to visualize the map of biomes in climatic space
png("reviewProceedingsB/EnvMap.png", 
    1820, 
    1820, 
    pointsize = 20)

# Set plotting parameters for margins and axis labeling
par(las = 1, mar = c(8,6,0,0))

# Plot a scatter plot of log-transformed MAP values against MAT (mean annual temperature)
# Color points by biome, with transparency for better visualization
plot(log(MAP_Val$value) ~ MAP_Val$MAT, 
     pch = ".",  # Use small points
     col = scales::alpha(colpal[MAP_Val$variable], 0.7),  # Apply biome colors with transparency
     frame = F,  # Remove the plot frame
     cex.lab = 3,  # Increase the label size
     cex.axis = 2,  # Increase the axis label size
     xlab = "",  # Omit the x-axis label
     xaxt = "n",  # Remove the x-axis ticks
     ylab = "Mean annual precipitation (log)")  # Label for the y-axis

# Close the PNG device to save the file
dev.off()
#====


#====
# Extract coordinates and raster data from the MAP raster and the environmental stack
a <- data.frame(coordinates(MAP),
                extract(Stack, coordinates(MAP)))

# Overlay the extracted spatial points with the ecoregion shapefile to assign biomes
myOv <- over(SpatialPoints(coordinates(Stack),
                           proj4string = crs(ecoRegi)),
             ecoRegi)

# Save the predicted niche space plot as a PNG file
png("FigsPNAS/PredictedNicheLog.png", 
    2000,
    1500, pointsize = 25)

# Set plotting margins to remove extra whitespace
par(mar = c(0, 0, 0, 0))

# Create a 3D scatter plot of MAP (logged), MAT, and the polymorphism probability
scatterplot3d::scatterplot3d(x = log1p(a$MAP.point),  # Log-transform MAP
                             y = a$MAT.point,  # Use MAT values
                             z = getValues(predx),  # Use predicted polymorphism probabilities
                             pch = ".",  # Small points
                             xlab = "MAP (ln)",  # X-axis label
                             ylab = "MAT",  # Y-axis label
                             zlab = "Polymorphism probability",  # Z-axis label
                             color = scales::alpha(colpal[myOv$BIOME], 0.7),  # Color points by biome
                             box = F)  # Remove the box around the plot

# Close the PNG device to save the file
dev.off()
#====


# Show the range and mean of MAT values for reference
range(MAT)  # Show the range of MAT values
mean(getValues(MAT), na.rm = T)  # Calculate the mean of MAT values, ignoring NAs
range(a$MAT.point * 8.45, na.rm = T)  # Show the range of rescaled MAT values
#====


# Create a 3D scatter plot of MAP, MAT, and predicted polymorphism probability
png("FigsPNAS/PredictedNicheNorm.png", 
    1500, 
    1500, pointsize = 45)

# Set the axis orientation and layout
par(las = 2)

# Plot MAP, MAT, and predicted values with a specified angle and axis limits
scatterplot3d::scatterplot3d(x = (getValues(MAP)),
                             y = (getValues(MAT)),
                             z = getValues(predx),  # Z-axis shows polymorphism probability
                             angle = 30,  # Set the viewing angle
                             zlim = c(0, 0.5),  # Z-axis limits
                             xlim = c(0, 10000),  # X-axis limits for MAP
                             ylim = c(-20, 30),  # Y-axis limits for MAT
                             pch = ".",  # Small points
                             scale.y = 1,  # Scale Y-axis
                             y.margin.add = 0.4,  # Add extra margin for the y-axis
                             xlab = "",  # Omit the x-axis label
                             ylab = "",  # Omit the y-axis label
                             mar = c(5,5,0,0.5),  # Set plot margins
                             cex.lab = 1.45,  # Label size
                             cex.axis = 1.4,  # Axis label size
                             zlab = "",  # Omit the z-axis label
                             color = scales::alpha(colpal[myOv$BIOME], 0.4),  # Color by biome with transparency
                             box = F)  # Remove the box

# Add a second scatterplot with Z-axis set to 0 for baseline comparison
par(new = T, las = 2)
scatterplot3d::scatterplot3d(x = (getValues(MAP)),
                             y = (getValues(MAT)),
                             z = rep(0, length(getValues(predx))),  # Set z-values to 0
                             angle = 30, 
                             scale.y = 1,
                             zlim = c(0, 0.5),
                             xlim = c(0, 10000),
                             ylim = c(-20, 30),
                             pch = ".",
                             y.margin.add = 0.4,
                             xlab = "",
                             ylab = "", 
                             mar = c(5,5,0,0.5),
                             cex.lab = 1.45, cex.axis = 1.4,
                             zlab = "",
                             color = "black",  # Set point color to black
                             cex = 0.3,
                             box = F)

# Close the PNG device to save the file
dev.off()
#====


# Plot the map of the global biomes
#====
# Save the biome map as a PNG file
png("BiomeMap.png",
    3000, 
    1000, 
    pointsize = 20)

# Set margins to remove whitespace
par(mar = c(0, 0, 0, 0))

# Plot each biome, using different colors for each biome region
plot(ecoRegi[ecoRegi$BIOME == 1,], col = "#006400", border = "#006400", 
     xlim = c(-70, 70), 
     ylim = c(-70, 90))  # Tropical rainforest
plot(ecoRegi[ecoRegi$BIOME == 2,], col = "#FFF5D7", border = "#FFF5D7", add = T)  # Dry forest
plot(ecoRegi[ecoRegi$BIOME == 3,], col = "#AAC800", border = "#AAC800", add = T)  # Tropical coniferous forest
plot(ecoRegi[ecoRegi$BIOME == 4,], col = "#CAFE8F", border = "#CAFE8F", add = T)  # Temperate mixed forest
plot(ecoRegi[ecoRegi$BIOME == 5,], col = "#9ED003", border = "#9ED003", add = T)  # Temperate coniferous forest
plot(ecoRegi[ecoRegi$BIOME == 6,], col = "#008D02", border = "#008D02", add = T)  # Boreal forest
plot(ecoRegi[ecoRegi$BIOME == 7,], col = "#F7D600", border = "#F7D600", add = T)  # Tropical savanna
plot(ecoRegi[ecoRegi$BIOME == 8,], col = "#FFB432", border = "#FFB432", add = T)  # Temperate grassland
plot(ecoRegi[ecoRegi$BIOME == 9,], col = "#cacd01", border = "#cacd01", add = T)  # Flooded savanna
plot(ecoRegi[ecoRegi$BIOME == 10,], col = "#8A9000", border = "#8A9000", add = T)  # Montane grasslands (Páramo)
plot(ecoRegi[ecoRegi$BIOME == 11,], col = "#AAAAAA", border = "#AAAAAA", add = T)  # Tundra
plot(ecoRegi[ecoRegi$BIOME == 12,], col = "brown", border = "brown", add = T)  # Mediterranean scrubland
plot(ecoRegi[ecoRegi$BIOME == 13,], col = "yellow", border = "yellow", add = T)  # Desert

# Overlay a world map and add biome points from 'occSP'
maps::map("world", add = T, interior = F, lwd = 0.2, color = "grey20")
plot(occSP, pch = ".", add = T, cex = 1.5)

# Add a legend for the biome classes
legend(-170, 20, fill = colpal, cex = 1.2, box.col = "white", 
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

# Close the PNG device to save the file
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
#' Subset Family and Calculate Species Richness by Grid
#'
#' This function subsets the input occurrence data by a given ant subfamily and then calculates species richness per grid cell.
#' It returns a `SpatialPointsDataFrame` object with spatial coordinates and the logarithm of species richness for the selected subfamily.
#'
#' @param occ_new A `data.frame` containing the occurrence data. This data frame should have columns for longitude (`dec_long`), latitude (`dec_lat`), `subFamily`, `poly_id`, and `valid_species_name`.
#' @param fam A `character` string specifying the subfamily to filter by (e.g., "Formicinae", "Myrmicinae", "Dolichoderinae").
#'
#' @return A `SpatialPointsDataFrame` object where each point represents a grid cell, and the data contains the log-transformed species richness for the selected subfamily.
#'
#' @details The function:
#' - Filters the occurrence data to retain only rows where the `subFamily` matches the input `fam` and `poly_id` equals 1.
#' - Aggregates the number of valid species names by grid cell, using the rounded longitude and latitude values.
#' - Constructs a `SpatialPointsDataFrame` with coordinates and the log-transformed species richness for the grid cells.
#'
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom stringr str_split
#' @importFrom stats aggregate
#' @importFrom raster crs
#' @examples
#' # Example usage:
#' # Assuming `occ_data` is a data frame with the appropriate structure
#' result <- subSetFam(occ_new = occ_data, fam = "Formicinae")
#' 
#' @export
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

