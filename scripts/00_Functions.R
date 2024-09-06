### Functions 


### custom function to colorize https://gist.github.com/fgabriel1891
f <- function(x,n=10, pal, rev = F){
  if(rev == F){ 
    rev(RColorBrewer::brewer.pal(n, pal))[cut(x,n)]
  }else{
    (RColorBrewer::brewer.pal(n, pal))[cut(x,n)]
  }
}

### Function to plot a mixed effect model output with ggplot

#' This function creates two plots from the given model, displaying the predicted probabilities 
#' of worker caste polymorphism as a function of temperature (MAT01) and precipitation (TAP01).
#' It uses `plot_model` from the `sjPlot` package and customizes the plots with the `ggplot2` 
#' theme and labels. The function returns a combined figure with the two plots arranged side by side.
#'
#' @param myModelx A fitted model object. This is the model from which predicted probabilities 
#' will be extracted and plotted.
#'
#' @return A combined ggplot object with two panels: 
#' - The first panel shows predicted probabilities as a function of temperature (MAT01).
#' - The second panel shows predicted probabilities as a function of precipitation (TAP01).
#' 
#' @details
#' This function produces two plots using the `plot_model` function from the `sjPlot` package, 
#' with predictions based on two covariates: MAT01 (temperature) and TAP01 (precipitation). 
#' The output plots show the predicted probabilities for worker caste polymorphism with custom 
#' formatting using `ggplot2`. 
#' 
#' @examples
#' # Assuming you have a model `my_model` fitted:
#' plotModels(my_model)
#'
#' @importFrom sjPlot plot_model
#' @importFrom ggplot2 theme element_blank element_line element_text xlab ylab ggtitle scale_y_continuous
#' @importFrom scales percent_format
#' @importFrom ggpubr ggarrange
#' 
plotModels <- function(myModelx){
  a<-plot_model(myModelx,type="pred",terms = "MAT01 [all]")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(family="Times New Roman", face="bold", size=20))+
    xlab("Temperature")+
    ylab("Predicted Probabilities of Worker Caste Polymorphism")+
    ggtitle("")+
    scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),
                       labels=scales::percent_format(accuracy = 1),limits=c(0,0.7));a
  
  b<-plot_model(myModelx,type="pred",terms = "TAP01 [all]")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text=element_text(family="Times New Roman", face="bold", size=20))+
    xlab("Precipitation")+
    ylab("Predicted Probabilities of Worker Caste Polymorphism")+
    ggtitle("")+
    scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels=scales::percent_format(accuracy = 1), limits = c(0,NA));b
  
  ggarrange(a,b)
}


## Function to estimate the number of species per grid 
#' Get Polymorphism Richness
#'
#' This function calculates the species richness for polymorphic or non-polymorphic species
#' based on a provided occurrence dataset. It returns the total number of unique species for 
#' each identifier (new_id) where polymorphism is either present or absent.
#'
#' @param occ A data frame containing occurrence records. It must include columns `new_id` 
#' (a factor representing unique identifiers) and `valid_species_name` (the valid species name).
#' @param poly A numeric vector indicating the polymorphism status to filter by. Default is 
#' `c(1,0)` where 1 represents polymorphic species and 0 represents non-polymorphic species.
#'
#' @return A numeric vector representing the species richness for each `new_id` based on 
#' polymorphism status. The values are the row sums of species occurrences (with duplicates removed).
#'
#' @details
#' This function subsets the occurrence data based on the provided polymorphism status (`poly`).
#' It then creates a table of `new_id` against `valid_species_name`, where each cell indicates 
#' whether a species is present for a given identifier. Duplicate species entries are counted only once.
#' The function returns the sum of unique species for each identifier (row).
#'
#' @examples
#' # Assuming `occ_data` is a data frame with the required columns:
#' poly_richness <- getPolyRich(occ_data, poly = 1)
#'
#' @export

getPolyRich <- function(occ, poly = c(1,0)){
  occ$new_id <- as.factor(occ$new_id)
  dat <- (occ[occ$poly_id == poly,])
  myTab <- table(dat$new_id,dat$valid_species_name)
  myTab[myTab>1]<-1
  return(rowSums(myTab))
  
  
  
  
}



# Function to make a plot of the global distribution of polymorphism in ants 

#' Generate World Map Based on Polymorphism, Density, Richness, or Effort
#'
#' This function creates a world map based on a variety of metrics, including the proportion of polymorphic species, 
#' density of occurrences, total richness, or sampling effort. It outputs a rasterized map which can be smoothed and 
#' colored according to the user's preferences. The function supports multiple plot types and a variety of color palettes.
#'
#' @param occ A data frame containing occurrence data. It should include columns for unique identifiers (`new_id`), 
#' climate variables (e.g., `TAP01`, `MAT01`), point density, and other relevant fields for polymorphism and richness.
#' @param plotType A character string specifying the type of plot to generate. Options are:
#' - `"prop"`: Creates a map showing the proportion of polymorphic species to total species richness.
#' - `"density"`: Creates a map based on the density of occurrences.
#' - `"poly"`: Creates a map based on the density of polymorphic species.
#' - `"richness"`: Creates a map based on species richness (total species).
#' - `"effort"`: Creates a map based on the sampling effort (density of records).
#' @param colbin Integer. Number of color bins to use for the map's color palette.
#' @param palName A character string specifying the name of the color palette to use (from `RColorBrewer`).
#' @param plot Logical. If `TRUE`, the function will plot the map; otherwise, it will only return the computed table.
#' @param rev Logical. If `TRUE`, the color palette will be reversed.
#'
#' @return A data frame (`myTab`) containing processed values for the chosen plot type (polymorphism proportion, density, richness, or effort)
#' and other related columns like climate variables (`TAP01`, `MAT01`) and geographic coordinates. If `plot` is `TRUE`, the function also 
#' generates a raster map and displays it.
#'
#' @details
#' The function processes the occurrence data and generates a table (`myTab`) with the relevant values based on the chosen `plotType`. 
#' Depending on the plot type, it calculates the proportion of polymorphic species, density of records, species richness, or sampling effort. 
#' It then rasterizes the values and, if requested, smooths the raster using a median filter with the `focal` function. 
#' The final raster map can be customized with color palettes from `RColorBrewer`, and smoothed for visual clarity.
#'
#' @examples
#' # Assuming `occ_data` is your data frame of occurrences:
#' result <- makeWorldPlot(occ_data, plotType = "richness", colbin = 10, palName = "Blues", rev = TRUE)
#'
#' @importFrom raster rasterFromXYZ disaggregate focal mask
#' @importFrom reshape2 melt
#' @importFrom stringr str_split
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales alpha
#' @importFrom maps map
#' @importFrom ggplot2 theme
#'
#' @export

makeWorldPlot <- function(occ, plotType = c("prop", "density","poly", "richness", "effort"), colbin = 3, palName = "Spectral", plot= T, rev = T){
  
  myTab <- data.frame("poly1" = getPolyRich(occ, 1),"poly0" = getPolyRich(occ, 0))
  
  myTab$poly0[myTab$poly0 == 0] <- 0.1
  myTab$poly1[myTab$poly1 == 0] <- 0.1
  
  if(plotType == "prop"){
    myTab <- data.frame(myTab,
                        apply(stringr::str_split(rownames(myTab), "\\.", simplify = T),2, as.numeric),
                        occ[match(rownames(myTab), occ$new_id),][c("TAP01","MAT01")])
    myTab$dens <- occ$point_density[match(paste0(myTab$X1,".",myTab$X2),occ$new_id)]
    myTab$Rich <- myTab$poly0+myTab$poly1
    myTab$prop <- log(myTab$poly1 / myTab$Rich)
    
    myTab <- data.frame(myTab,
                        apply(stringr::str_split(rownames(myTab), "\\.", simplify = T),2, as.numeric),
                        occ[match(rownames(myTab), occ$new_id),][c("TAP01","MAT01")])
    print("creating ratiomap object")
    myOBj <- reshape2::melt(xtabs(prop~X1+X2, myTab))
    myOBj
    
  }
  if(plotType == "density"){
    
    myTab <- data.frame(myTab,
                        apply(stringr::str_split(rownames(myTab), "\\.", simplify = T),2, as.numeric),
                        occ[match(rownames(myTab), occ$new_id),][c("TAP01","MAT01")])
    
    myTab$dens <- occ$point_density[match(paste0(myTab$X1,".",myTab$X2),occ$new_id)]
    myTab$Rich <- myTab$poly0+myTab$poly1
    # create a polymorphic density per total richness penalized by sampling
    #myTab$PenDens <- (myTab$poly1/myTab$Rich)*(1/log1p(myTab$dens))
    myTab$PenDens <- log((myTab$poly1)*(1/log1p(myTab$dens)))
    
    print("creating map object")
    myOBj <- reshape2::melt(xtabs(PenDens~X1+X2, myTab))
    
  }
  
  if(plotType == "poly"){
    
    myTab <- data.frame(myTab,
                        apply(stringr::str_split(rownames(myTab), "\\.", simplify = T),2, as.numeric),
                        occ[match(rownames(myTab), occ$new_id),][c("TAP01","MAT01")])
    
    # create a polymorphic density per total richness penalized by sampling
    #myTab$PenDens <- (myTab$poly1/myTab$Rich)*(1/log1p(myTab$dens))
    myTab$PenDens <- log(myTab$poly1)
    
    print("creating polymap object")
    myOBj <- reshape2::melt(xtabs(PenDens~X1+X2, myTab))
    
  }
  
  if(plotType == "richness"){
    
    myTab <- data.frame(myTab,
                        apply(stringr::str_split(rownames(myTab), "\\.", simplify = T),2, as.numeric),
                        occ[match(rownames(myTab), occ$new_id),][c("TAP01","MAT01")])
    
    myTab$Rich <- myTab$poly0+myTab$poly1
    # create a polymorphic density per total richness penalized by sampling
    #myTab$PenDens <- (myTab$poly1/myTab$Rich)*(1/log1p(myTab$dens))
    myTab$PenDens <- log(myTab$Rich)
    
    print("creating richnessmap object")
    myOBj <- reshape2::melt(xtabs(PenDens~X1+X2, myTab))
    
  }
  
  if(plotType == "effort"){
    
    myTab <- data.frame(myTab,
                        apply(stringr::str_split(rownames(myTab), "\\.", simplify = T),2, as.numeric),
                        occ[match(rownames(myTab), occ$new_id),][c("TAP01","MAT01")])
    
    myTab$dens <- occ$point_density[match(paste0(myTab$X1,".",myTab$X2),occ$new_id)]
    # create a polymorphic density per total richness penalized by sampling
    #myTab$PenDens <- (myTab$poly1/myTab$Rich)*(1/log1p(myTab$dens))
    myTab$PenDens <- log(myTab$dens)
    
    print("creating effortmap object")
    myOBj <- reshape2::melt(xtabs(PenDens~X1+X2, myTab))
    
  }
  if(plot == T){
    
    
    
    myRas1 <- expand.grid(-179:180, -55:71)
    myRas1$id <- paste(myRas1$Var1, myRas1$Var2, sep ="_")
    myOBj$id <- paste(myOBj$X1, myOBj$X2, sep ="_")
    myRas1$val <- myOBj$value[match(myRas1$id, myOBj$id)]
    myRas1$id <- NULL
    
    myRas1$val[myRas1$val == 0] <- NA
  
    
    
    myRas <- raster::rasterFromXYZ(myRas1, 
                                   crs = "+proj=longlat +datum=WGS84 +no_defs")
    
 
    print("smoothing raster, a little patience here")
    r2 <- disaggregate(myRas,5, method="bilinear")
    r4 <- focal(r2, w= matrix(1,5,5), median)
    
    pal <- RColorBrewer::brewer.pal(n = colbin, name = palName)
    
    if(rev == T){ 
      pal <- rev(pal)
    }
    
    
    plot(r4, legend = F, 
         ylim = c(-66.5,66.5)
         ,axes = F, box = F, 
         col = pal, colNA = scales::alpha("white",1))
    #maps::map("world", add = T, interior = F, border = NA, fill = T, col = scales::alpha("black",0.7))
    maps::map("world", add = T, interior = F, border = NA, fill = T, col = "grey94")
    
    r4[is.na(r4)] <- NA
    mr <- mask(r4, r4)


    plot(mr, legend = F,
         ylim = c(-66.5,66.5)
         ,axes = F, box = F, 
         col = pal, add = TRUE)
    #plot(oceans, ylim = c(-70,70),add = T, col = scales::alpha("white", 0.3), lwd = 0.7)
    #abline(h= seq(-50, 50, 25), v = seq(-150, 150, 25), col = "#4575b4", lwd = 0.2)
    print("done")
    cat("\014")
    
    
  }
  
  
  return(myTab)
}


#' Fit Multiple Generalized Linear Mixed Models (GLMMs)
#'
#' This function fits a series of generalized linear mixed models (GLMMs) to predict the presence 
#' or absence of polymorphism (`poly_id`) in a dataset of occurrences. The models vary in complexity, 
#' with different combinations of predictors and interaction terms, including temperature (`MAT01`) 
#' and precipitation (`TAP01`), with random effects accounted for by `new_id`.
#'
#' @param occ_new A data frame containing occurrence data. It must include columns `poly_id` (binary response), 
#' `new_id` (random effect), and environmental predictors such as `TAP01` (precipitation) and `MAT01` (temperature).
#'
#' @return A list of fitted model objects, including:
#' - `nullmodel`: A null model with no predictors.
#' - `spamodel1`: A random intercept model with no predictors.
#' - `spamodel2`: A model with `TAP01` and `MAT01` as predictors.
#' - `spamodel3`: A model with `TAP01`, `TAP01^2`, and `MAT01`.
#' - `spamodel4`: A model with the interaction between `TAP01` and `MAT01`.
#' - `MATmodel`: A model with `MAT01` as the sole predictor.
#' - `MAPmodel`: A model with `TAP01` as the sole predictor.
#' - `MAP2model`: A model with `TAP01` and `TAP01^2`.
#' - `MAT2model`: A model with `MAT01` and `MAT01^2`.
#'
#' @details
#' This function fits several GLMMs using the `glmer` function from the `lme4` package. Each model 
#' includes a random intercept for `new_id`, and the function explores various combinations of predictors, 
#' including linear and quadratic terms for temperature and precipitation. The models are fitted using 
#' a binomial family with a logit link function. The optimization is controlled using the "Nelder_Mead" method.
#'
#' @examples
#' # Assuming `occ_data` is your data frame of occurrences:
#' models <- fitModels(occ_data)
#'
#' @importFrom lme4 glmer glmerControl
#' @importFrom stats binomial
#'
#' @export

fitModels <- function(occ_new){
  
  print("nullmodel")
  nullmodel<-glm(poly_id~1,data=occ_new,family="binomial"(link=logit))
  print("mod1")
  spamodel1<-glmer(poly_id~1+(1|new_id),
                   data=occ_new,
                   family="binomial"(link=logit),
                   glmerControl(optimizer = "Nelder_Mead"))
  print("mod2")
  spamodel2<-glmer(poly_id~TAP01+MAT01+
                     (1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("mod3")
  spamodel3<-glmer(poly_id~TAP01+I(TAP01^2)+MAT01+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("mod4")
  spamodel4<-glmer(poly_id~TAP01*MAT01+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("MAT_mod1")
  MATmodel<-glmer(poly_id~MAT01+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("MAT_mod2")
  MAPmodel<-glmer(poly_id~TAP01+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("MAT_mod3")
  MAP2model<-glmer(poly_id~I(TAP01^2)+TAP01+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("MAT_mod4")
  MAT2model<-glmer(poly_id~I(MAT01^2)+MAT01+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
  print("done")
  return(list(nullmodel, spamodel1, spamodel2, spamodel3, spamodel4,MATmodel ,MAPmodel,MAP2model,MAT2model))
  cat("\014")
}


#' Fit Multiple Generalized Linear Models (GLMs) Without Random Effects
#'
#' This function fits a series of generalized linear models (GLMs) to predict the presence or absence of polymorphism (`poly_id`) 
#' in a dataset of occurrences. The models vary in complexity, with different combinations of predictors and interaction terms, 
#' including temperature (`MAT01`) and precipitation (`TAP01`). Unlike the `fitModels` function, no random effects are included.
#'
#' @param occ_new A data frame containing occurrence data. It must include columns `poly_id` (binary response), 
#' `TAP01` (precipitation), and `MAT01` (temperature).
#'
#' @return A list of fitted model objects, including:
#' - `nullmodel`: A null model with no predictors.
#' - `spamodel1`: A model with `TAP01` and `MAT01` as predictors.
#' - `spamodel2`: A model with `TAP01`, `TAP01^2`, and `MAT01`.
#' - `spamodel3`: A model with the interaction between `TAP01` and `MAT01`.
#' - `MATmodel`: A model with `MAT01` as the sole predictor.
#' - `MAPmodel`: A model with `TAP01` as the sole predictor.
#' - `MAP2model`: A model with `TAP01` and `TAP01^2`.
#' - `MAT2model`: A model with `MAT01` and `MAT01^2`.
#'
#' @details
#' This function fits several GLMs using the `glm` function from the base R package. The models explore various combinations 
#' of predictors, including linear and quadratic terms for temperature and precipitation. All models are fitted using a binomial 
#' family with a logit link function. The function is designed for cases where random effects are not necessary.
#'
#' @examples
#' # Assuming `occ_data` is your data frame of occurrences:
#' models <- fitModelsNoRan(occ_data)
#'
#' @importFrom stats glm binomial
#'
#' @export

fitModelsNoRan <- function(occ_new){
  
  print("nullmodel")
  nullmodel<-glm(poly_id~1,
                 data=occ_new,
                 family="binomial"(link=logit))

  print("mod1")
  spamodel1<-glm(poly_id~TAP01+MAT01,
                   data=occ_new,family="binomial"(link=logit))
  print("mod2")
  spamodel2<-glm(poly_id~TAP01+I(TAP01^2)+MAT01,
                   data=occ_new,family="binomial"(link=logit))
  print("mod3")
  spamodel3<-glm(poly_id~TAP01*MAT01,
                   data=occ_new,family="binomial"(link=logit))
  print("MAT_mod1")
  MATmodel<-glm(poly_id~MAT01,
                  data=occ_new,family="binomial"(link=logit))
  print("MAT_mod2")
  MAPmodel<-glm(poly_id~TAP01,
                  data=occ_new,family="binomial"(link=logit))
  print("MAT_mod3")
  MAP2model<-glm(poly_id~I(TAP01^2)+TAP01,
                   data=occ_new,family="binomial"(link=logit))
  print("MAT_mod4")
  MAT2model<-glm(poly_id~I(MAT01^2)+MAT01,
                   data=occ_new,family="binomial"(link=logit))
  print("done")
  return(list(nullmodel, spamodel1, spamodel2, spamodel3,MATmodel ,MAPmodel,MAP2model,MAT2model))
  cat("\014")
}


#' Fit Multiple Generalized Linear Models (GLMs) Including Species Richness Without Random Effects
#'
#' This function fits a series of generalized linear models (GLMs) to predict the presence or absence of polymorphism (`poly_id`) 
#' in a dataset of occurrences. The models include species richness (`SpRich`) as a predictor, along with environmental variables 
#' like temperature (`MAT01`) and precipitation (`TAP01`). No random effects are included in these models.
#'
#' @param occ_new A data frame containing occurrence data. It must include columns `poly_id` (binary response), 
#' `TAP01` (precipitation), `MAT01` (temperature), and `SpRich` (species richness).
#'
#' @return A list of fitted model objects, including:
#' - `nullmodel`: A null model with no predictors.
#' - `spamodel1`: A model with species richness (`SpRich`) as the sole predictor.
#' - `spamodel2`: A model with `TAP01`, `MAT01`, and `SpRich` as predictors.
#' - `spamodel3`: A model with `TAP01`, `TAP01^2`, `MAT01`, and `SpRich`.
#' - `spamodel4`: A model with the interaction between `TAP01` and `MAT01`, along with `SpRich`.
#' - `MATmodel`: A model with `MAT01` and `SpRich` as predictors.
#' - `MAPmodel`: A model with `TAP01` and `SpRich` as predictors.
#' - `MAP2model`: A model with `TAP01`, `TAP01^2`, and `SpRich`.
#' - `MAT2model`: A model with `MAT01`, `MAT01^2`, and `SpRich`.
#'
#' @details
#' This function fits several GLMs using the `glm` function from the base R package. Each model includes species richness (`SpRich`) 
#' as a predictor, along with various combinations of environmental predictors such as temperature and precipitation. All models 
#' are fitted using a binomial family with a logit link function. The function is useful when modeling species richness as an additional 
#' explanatory variable without the need for random effects.
#'
#' @examples
#' # Assuming `occ_data` is your data frame of occurrences:
#' models <- fitModelsNoRanplusRich(occ_data)
#'
#' @importFrom stats glm binomial
#'
#' @export

fitModelsNoRanplusRich <- function(occ_new){
  
  print("nullmodel")
  nullmodel<-glm(poly_id~1,
                 data=occ_new,
                 family="binomial"(link=logit))
  spamodel1<-glm(poly_id~ SpRich,
                 data=occ_new,family="binomial"(link=logit))
  print("mod2")
  spamodel2<-glm(poly_id~TAP01+MAT01+ SpRich,
                 data=occ_new,family="binomial"(link=logit))
  print("mod3")
  spamodel3<-glm(poly_id~TAP01+I(TAP01^2)+MAT01+ SpRich,
                 data=occ_new,family="binomial"(link=logit))
  print("mod4")
  spamodel4<-glm(poly_id~TAP01*MAT01+ SpRich,
                 data=occ_new,family="binomial"(link=logit))
  print("MAT_mod1")
  MATmodel<-glm(poly_id~MAT01+ SpRich,
                data=occ_new,family="binomial"(link=logit))
  print("MAT_mod2")
  MAPmodel<-glm(poly_id~TAP01+ SpRich,
                data=occ_new,family="binomial"(link=logit))
  print("MAT_mod3")
  MAP2model<-glm(poly_id~I(TAP01^2)+TAP01+ SpRich,
                 data=occ_new,family="binomial"(link=logit))
  print("MAT_mod4")
  MAT2model<-glm(poly_id~I(MAT01^2)+MAT01+ SpRich,
                 data=occ_new,family="binomial"(link=logit))
  print("done")
  return(list(nullmodel, spamodel1, spamodel2, spamodel3,MATmodel ,MAPmodel,MAP2model,MAT2model))
  cat("\014")
}



