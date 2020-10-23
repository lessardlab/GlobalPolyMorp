### Condensed script to replicate the data preparation and data fitting with different datasets. 

library(raster)
library(sp)
library(rworldmap)
library(rgeos)
library(maptools)
library(rnaturalearth)
library(vegan)
library(scales)
library(adespatial)
library(nlme)
library(lme4)
library(ggplot2)


# load functions 
source("Functions.R")
# load data
occ2 <- read.csv("Lat-Long_Data_GABI.csv",stringsAsFactors = TRUE, header = T)
occ<-read.csv("NewData/Lat-Long_Data_GABI.csv",stringsAsFactors = TRUE, header = T)
specieslist<-read.csv("matched-polymorphism_bg.csv",stringsAsFactors = TRUE, header = T)

# Replicate data prepartion steps
occ_clean <- DataCleaner(occ, specieslist)
occ_clean2 <- DataCleaner(occ2, specieslist)

# Replicate climate matching steps
env_occ_match <- WorldClimMatch(occ_clean)
env_occ_match2 <- WorldClimMatch(occ_clean2)

# Replicate model fitting steps
MyModels1 <- ModelSelection(env_occ_match)
MyModels2 <- ModelSelection(env_occ_match2)

# See AIC differences
MyModels1[[length(MyModels1)]]
MyModels2[[length(MyModels1)]]

# Plot the fitted models 
plotModels(MyModels1[[5]])
plotModels(MyModels2[[5]])
