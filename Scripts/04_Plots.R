### GPM Figures 
## Gabriel Muñoz _ dic 2020. 


## Please note that functions to produce plots feed on the "occ" dataset, you need to have that one loaded first. 

# source("GM_Scripts/1_PrepData.R") # run this line if you haven't gone through the 1. Script. 

#source("4_Functions.R")

occ<-read.csv("Data/occ_withsitedensity.csv",stringsAsFactors = F)
occ <- occ[complete.cases(occ),]

#https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html

library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()

# make the global proportion plot 
pdf("GlobalMapProp.pdf", width = 6, height = 4, pointsize = 10)
par(mar = c(2,2,2,2))
makeWorldPlot(occ, "prop", colbin = 11, "RdYlBu", plot = T, rev = T)
dev.off()

pdf("GlobalMapDensity.pdf", width = 6, height = 4, pointsize = 15)
par(mar = c(2,2,2,2))
makeWorldPlot(occ, "density", colbin = 9, "YlOrRd", plot = T, rev = F)
dev.off()

pdf("GlobalMapPoly.pdf", width = 6, height = 4, pointsize = 15)
par(mar = c(2,2,2,2))
makeWorldPlot(occ, "poly", colbin = 11,"RdYlGn" , plot = T, rev = F)
dev.off()

pdf("GlobalMapRichness.pdf", width = 6, height = 4, pointsize = 15)
par(mar = c(2,2,2,2))
makeWorldPlot(occ, "richness", colbin = 9, "Greens", plot = T, rev = F)
dev.off()

pdf("GlobalMapEffort.pdf", width = 6, height = 4, pointsize = 15)
par(mar = c(2,2,2,2))
makeWorldPlot(occ, "effort", colbin = 11, "RdGy", plot = T, rev = F)
dev.off()


# Make the partial plots

## Myrmicinae
png("FiguresPaper/GlobalMapRatio_Myrm.png", width = 900, height = 450, pointsize = 10)
par(mar = c(2,2,2,2))
makeWorldPlot(droplevels(occ[occ$genus == "Camponotus",]), "prop", "Spectral", colbin = 9, plot = T, rev = T)
dev.off()


png("FiguresPaper/GlobalMapDensity_Myr.png", width = 900, height = 450, pointsize = 10)
par(mar = c(2,2,2,2))
makeWorldPlot(droplevels(occ[occ$subFamily == "Myrmicinae",]), "density", "YlOrRd", plot = T, rev = F)
dev.off()

## Formicinae
png("FiguresPaper/GlobalMapRatio_Form.png", width = 900, height = 450, pointsize = 10)
par(mar = c(2,2,2,2))
makeWorldPlot(droplevels(occ[occ$subFamily == "Formicinae",]), "ratio", "Spectral", plot = T, rev = T)
dev.off()


png("FiguresPaper/GlobalMapDensity_Form.png", width = 900, height = 450, pointsize = 10)
par(mar = c(2,2,2,2))
makeWorldPlot(droplevels(occ[occ$subFamily == "Formicinae",]), "density", "YlOrRd", plot = T, rev = F)
dev.off()

###############




##############################


## Plot the global distribution of ant polymorphism in climatic space 
# dataframe with all species
mt <- makeWorldPlot(occ, "prop", plot = F)
mt <- data.frame(mt,makeWorldPlot(occ, "density", plot = F)[c("dens", "PenDens")])
mt$Temp <- MAT1$x[match(rownames(mt), stringr::str_replace( MAT1$Group.1, "_", "."))]
mt$Prec <- MAP1$x[match(rownames(mt), stringr::str_replace( MAP1$Group.1, "_", "."))]

mt$MAT_round <- round(mt$MAT.point, 2)
mt$MAP_round <- round(mt$MAP.point,2)

agg_mt <- aggregate(mt$prop, list(mt$MAT_round, mt$MAP_round), median)
agg_mt_dens <- aggregate(mt$dens, list(mt$MAT_round, mt$MAP_round), median)
agg_mt_Pdens <- aggregate(mt$PenDens, list(mt$MAT_round, mt$MAP_round), median)
###############

par(las = 1, oma = c(2,2,2,2))
### proportion poly
plot(MAP1$x~MAT1$x, cex = 4, 
     xaxt = "n",
     yaxt = "n",
     col = "grey94", frame = F,
     xlab = "", 
     ylab = "",
     pch = ".", 
     ylim = c(-1,6))
points(agg_mt$Group.2~ agg_mt$Group.1, cex = 4, 
       col = scales::alpha(f(round(agg_mt$x), 11, "RdYlBu", F),1),
       #col = ifelse(mt$rat >1, "yellow", "red"),
       pch = ".")

abline(h = 0, v = 0, col = "grey50", lty = 2)
axis(1,seq(-2,1.5, 0.5), labels = c(-3.82-12.07*4, -3.82-12.07*3,-3.82-12.07*2,-3.82-12.07,-3.82,-3.82+12.07,-3.82+(12.07*2),-3.82+(12.07*3)))
axis(2,seq(-1,6, 1), labels = c(0,550.05,550.05+652.48,550.05+(652.48*2),550.05+(652.48*3),550.05+(652.48*4),550.05+(652.48*5),550.05+(652.48*6)))


### density 
plot(MAP1$x~MAT1$x, cex = 5, 
     xaxt = "n",
     yaxt = "n",
     col = "grey90", frame = F,
     xlab = "", 
     ylab = "",
     pch = ".", 
     ylim = c(-1,6))
points(agg_mt_Pdens$Group.2[agg_mt$x>0]~ agg_mt_Pdens$Group.1[agg_mt$x>0], cex = 5, 
       col = scales::alpha(f(log1p(agg_mt_Pdens$x[agg_mt_Pdens$x>0]),9, "YlOrRd", T),1),
       #col = ifelse(mt$rat >1, "yellow", "red"),
       pch = ".")
abline(h = 0, v = 0, col = "grey50", lty = 2)
axis(1,seq(-2,1.5, 0.5), labels = c(-3.82-12.07*4, -3.82-12.07*3,-3.82-12.07*2,-3.82-12.07,-3.82,-3.82+12.07,-3.82+(12.07*2),-3.82+(12.07*3)))
axis(2,seq(-1,6, 1), labels = c(0,550.05,550.05+652.48,550.05+(652.48*2),550.05+(652.48*3),550.05+(652.48*4),550.05+(652.48*5),550.05+(652.48*6)))



################
## Myrmicinae
#################

# dataframe with all species
mt_myr <- makeWorldPlot(occ[occ$subFamily == "Myrmicinae",], "ratio", "Spectral", plot = F, rev = T)
mt_myr <- data.frame(mt_myr,makeWorldPlot(occ[occ$subFamily == "Myrmicinae",], plotType = "density", plot = F, rev = T)[c("dens", "PenDens")])
mt_myr$Temp <- MAT1$x[match(rownames(mt_myr), stringr::str_replace( MAT1$Group.1, "_", "."))]
mt_myr$Prec <- MAP1$x[match(rownames(mt_myr), stringr::str_replace( MAP1$Group.1, "_", "."))]

mt_myr$MAT_round <- round(mt_myr$MAT.point, 2)
mt_myr$MAP_round <- round(mt_myr$MAP.point, 2)

agg_mt_myr <- aggregate(mt_myr$rat, list(mt_myr$MAT_round, mt_myr$MAP_round), median)
agg_mt_dens_myr <- aggregate(mt_myr$dens, list(mt_myr$MAT_round, mt_myr$MAP_round), median)
agg_mt_Pdens_myr <- aggregate(mt_myr$PenDens, list(mt_myr$MAT_round, mt_myr$MAP_round), median)
###############

par(las = 1, oma = c(2,2,2,2))
### ratio 
plot(MAP1$x~MAT1$x, cex = 4, 
     xaxt = "n",
     yaxt = "n",
     col = "grey90", frame = F,
     xlab = "", 
     ylab = "",
     pch = ".", 
     ylim = c(-1,6))
points(agg_mt_myr$Group.2~ agg_mt_myr$Group.1, cex = 4, 
       col = scales::alpha(f(log1p(round(agg_mt_myr$x)), 6, "Spectral", F),0.5),
       #col = ifelse(mt$rat >1, "yellow", "red"),
       pch = ".")

abline(h = 0, v = 0, col = "grey50", lty = 2)
axis(1,seq(-2,1.5, 0.5), labels = c(-3.82-12.07*4, -3.82-12.07*3,-3.82-12.07*2,-3.82-12.07,-3.82,-3.82+12.07,-3.82+(12.07*2),-3.82+(12.07*3)))
axis(2,seq(-1,6, 1), labels = c(0,550.05,550.05+652.48,550.05+(652.48*2),550.05+(652.48*3),550.05+(652.48*4),550.05+(652.48*5),550.05+(652.48*6)))


### density 
plot(MAP1$x~MAT1$x, cex = 4, 
     xaxt = "n",
     yaxt = "n",
     col = "grey90", frame = F,
     xlab = "", 
     ylab = "",
     pch = ".", 
     ylim = c(-1,6))
points(agg_mt_Pdens_myr$Group.2[agg_mt$x>0]~ agg_mt_Pdens_myr$Group.1[agg_mt$x>0], cex = 3, 
       col = scales::alpha(f((log1p(agg_mt_Pdens_myr$x[agg_mt$x>0])), 6, "YlOrRd", T),0.6),
       #col = ifelse(mt$rat >1, "yellow", "red"),
       pch = ".")
abline(h = 0, v = 0, col = "grey50", lty = 2)
axis(1,seq(-2,1.5, 0.5), labels = c(-3.82-12.07*4, -3.82-12.07*3,-3.82-12.07*2,-3.82-12.07,-3.82,-3.82+12.07,-3.82+(12.07*2),-3.82+(12.07*3)))
axis(2,seq(-1,6, 1), labels = c(0,550.05,550.05+652.48,550.05+(652.48*2),550.05+(652.48*3),550.05+(652.48*4),550.05+(652.48*5),550.05+(652.48*6)))


################
## Formycinae
#################


# dataframe with all species
mt_for <- makeWorldPlot(occ[occ$subFamily == "Formicinae",], "ratio", "Spectral", plot = F, rev = T)
mt_for <- data.frame(mt_for,makeWorldPlot(occ[occ$subFamily == "Formicinae",], plotType = "density", plot = F, rev = T)[c("dens", "PenDens")])
mt_for$Temp <- MAT1$x[match(rownames(mt_for), stringr::str_replace( MAT1$Group.1, "_", "."))]
mt_for$Prec <- MAP1$x[match(rownames(mt_for), stringr::str_replace( MAP1$Group.1, "_", "."))]

mt_for$MAT_round <- round(mt_for$MAT.point, 2)
mt_for$MAP_round <- round(mt_for$MAP.point, 2)

agg_mt_for <- aggregate(mt_for$rat, list(mt_for$MAT_round, mt_for$MAP_round), median)
agg_mt_dens_for <- aggregate(mt_for$dens, list(mt_for$MAT_round, mt_for$MAP_round), median)
agg_mt_Pdens_for <- aggregate(mt_for$PenDens, list(mt_for$MAT_round, mt_for$MAP_round), median)
###############

par(las = 1, oma = c(2,2,2,2))
### ratio 
plot(MAP1$x~MAT1$x, cex = 4, 
     xaxt = "n",
     yaxt = "n",
     col = "grey90", frame = F,
     xlab = "", 
     ylab = "",
     pch = ".", 
     ylim = c(-1,6))
points(agg_mt_for$Group.2~ agg_mt_for$Group.1, cex = 4, 
       col = scales::alpha(f(round(agg_mt_for$x), 6, "Spectral", F),0.5),
       #col = ifelse(mt$rat >1, "yellow", "red"),
       pch = ".")

abline(h = 0, v = 0, col = "grey50", lty = 2)
axis(1,seq(-2,1.5, 0.5), labels = c(-3.82-12.07*4, -3.82-12.07*3,-3.82-12.07*2,-3.82-12.07,-3.82,-3.82+12.07,-3.82+(12.07*2),-3.82+(12.07*3)))
axis(2,seq(-1,6, 1), labels = c(0,550.05,550.05+652.48,550.05+(652.48*2),550.05+(652.48*3),550.05+(652.48*4),550.05+(652.48*5),550.05+(652.48*6)))


### density 
plot(MAP1$x~MAT1$x, cex = 4, 
     xaxt = "n",
     yaxt = "n",
     col = "grey90", frame = F,
     xlab = "", 
     ylab = "",
     pch = ".", 
     ylim = c(-1,6))
points(agg_mt_Pdens_for$Group.2[agg_mt$x>0]~ agg_mt_Pdens_for$Group.1[agg_mt$x>0], cex = 4, 
       col = scales::alpha(f((log1p(agg_mt_Pdens_for$x[agg_mt$x>0])), 6, "YlOrRd", T),0.5),
       #col = ifelse(mt$rat >1, "yellow", "red"),
       pch = ".")
abline(h = 0, v = 0, col = "grey50", lty = 2)
axis(1,seq(-2,1.5, 0.5), labels = c(-3.82-12.07*4, -3.82-12.07*3,-3.82-12.07*2,-3.82-12.07,-3.82,-3.82+12.07,-3.82+(12.07*2),-3.82+(12.07*3)))
axis(2,seq(-1,6, 1), labels = c(0,550.05,550.05+652.48,550.05+(652.48*2),550.05+(652.48*3),550.05+(652.48*4),550.05+(652.48*5),550.05+(652.48*6)))

#####################


head(mt)


par(mfrow = c(1,3), mar = c(4,4,2,2))
plot(log(mt$poly1),log(mt$poly0), 
     col = f(mt$rat, 6, "Spectral"),
     xlim = c(0,7), ylim = c(0,7),
     pch  =15, main = "All species", 
     cex.lab = 1.5,
     ylab = "Monomorphic richness (Log)", 
     xlab = "Polymorphic richness (Log)",
     frame = F)
abline(a=0,b=1)


plot(log(mt_myr$poly1),log(mt_myr$poly0), 
     col = f(mt_myr$rat, 6, "Spectral"),
     xlim = c(0,7), ylim = c(0,7),
     cex.lab = 1.5,
     pch  =15, main = "Myrmicinae",
     ylab = "Monomorphic richness (Log)", 
     xlab = "Polymorphic richness (Log)",
     frame = F)
abline(a=0,b=1)


plot(log(mt_for$poly1),log(mt_for$poly0), 
     col = f(mt_for$rat, 6, "Spectral"),
     xlim = c(0,7), ylim = c(0,7),
     cex.lab = 1.5,
     pch  =15, main = "Formicinae",
     ylab = "Monomorphic richness (Log)", 
     xlab = "Polymorphic richness (Log)",
     frame = F)
abline(a=0,b=1)


######## Richness vs environment tradeoffs

library(RStoolbox)


EnvStack <- raster::stack(MAP, MAT)

PCAras <- RStoolbox::rasterPCA(EnvStack)

dev.off()

##########
# ALL SPECIES
##########


png("FiguresPaper/EnvMapsALL.png", width = 1000, height = 1000, pointsize = 10)
par(mfrow = c(2,2), mar = c(1,1,1,1), oma = c(1,1,1,1))
plot(log(PCAras$map[[2]]+50), box=F, axes = F,
     legend = F, main = "Polymorphic species")
points(X2~X1, 
       data = mt[mt$rat>0,], pch = ".", 
       cex = log1p(PenDens)*2)
plot(log(PCAras$map[[2]]+50), box=F, axes = F, legend = F, main = "Monomorphic species")
points(X2~X1, 
       data = mt[mt$rat<0,], pch = ".", 
       cex = log1p(PenDens)*2)

plot(-log(PCAras$map[[1]]+1000), box=F, axes = F, legend = F, main = "Polymorphic species")
points(X2~X1, 
       data = mt[mt$rat>0,], pch = ".", 
       cex = log1p(PenDens)*2)

plot(-log(PCAras$map[[1]]+1000), box=F, axes = F, legend = F, main = "Monomorphic species")
points(X2~X1, 
       data = mt[mt$rat<0,], pch = ".", 
       cex = log1p(PenDens)*2)
dev.off()

##########
# Myrmicinae
##########

png("FiguresPaper/EnvMapsMyr.png", width = 1000, height = 1000, pointsize = 10)
par(mfrow = c(2,2), mar = c(1,1,1,1), oma = c(1,1,1,1))
plot(log(PCAras$map[[2]]+50), box=F, axes = F,
     legend = F, main = "Polymorphic species")
points(X2~X1, 
       data = mt_myr[mt_myr$rat>0,], pch = ".", 
       cex = log1p(PenDens)*2)
plot(log(PCAras$map[[2]]+50), box=F, axes = F, legend = F, main = "Monomorphic species")
points(X2~X1, 
       data = mt_myr[mt_myr$rat<0,], pch = ".", 
       cex = log1p(PenDens)*2)

plot(-log(PCAras$map[[1]]+1000), box=F, axes = F, legend = F, main = "Polymorphic species")
points(X2~X1, 
       data = mt_myr[mt_myr$rat>0,], pch = ".", 
       cex = log1p(PenDens)*2)

plot(-log(PCAras$map[[1]]+1000), box=F, axes = F, legend = F, main = "Monomorphic species")
points(X2~X1, 
       data = mt_myr[mt_myr$rat<0,], pch = ".", 
       cex = log1p(PenDens)*2)
dev.off()



##########
# Formycinae
##########

png("FiguresPaper/EnvMapsFor.png", width = 1000, height = 1000, pointsize = 10)
par(mfrow = c(2,2), mar = c(1,1,1,1), oma = c(1,1,1,1))
plot(log(PCAras$map[[2]]+50), box=F, axes = F,
     legend = F, main = "Polymorphic species")
points(X2~X1, 
       data = mt_for[mt_for$rat>0,], pch = ".", 
       cex = log1p(PenDens)*2)
plot(log(PCAras$map[[2]]+50), box=F, axes = F, legend = F, main = "Monomorphic species")
points(X2~X1, 
       data = mt_for[mt_for$rat<0,], pch = ".", 
       cex = log1p(PenDens)*2)

plot(-log(PCAras$map[[1]]+1000), box=F, axes = F, legend = F, main = "Polymorphic species")
points(X2~X1, 
       data = mt_for[mt_for$rat>0,], pch = ".", 
       cex = log1p(PenDens)*2)

plot(-log(PCAras$map[[1]]+1000), box=F, axes = F, legend = F, main = "Monomorphic species")
points(X2~X1, 
       data = mt_for[mt_for$rat<0,], pch = ".", 
       cex = log1p(PenDens)*2)
dev.off()

