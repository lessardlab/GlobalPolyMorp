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

getPolyRich <- function(occ, poly = c(1,0)){
  occ$new_id <- as.factor(occ$new_id)
  dat <- (occ[occ$poly_id == poly,])
  myTab <- table(dat$new_id,dat$valid_species_name)
  myTab[myTab>1]<-1
  return(rowSums(myTab))
  
  
  
  
}



# Function to make a plot of the global distribution of polymorphism in ants 


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



