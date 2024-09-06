######################
## Fit linear mixed effect models with the effect of ricness of grid as covariate 
## Gabriel Muñoz 
## NOV 2021

library(vegan)
library(adespatial)
library(sp)
library(nlme)
library(lme4)
library(AICcmodavg)
library(sjPlot)


# add richness per grid
SpRich <- rowSums(table(occ_new$new_id,occ_new$valid_species_name))
occ_new$SpRich <- SpRich[match(occ_new$new_id,names(SpRich))]
# Logarith of richness to account for power-laws in richness growth 
occ_new$SpRichLog <- log(occ_new$SpRich)


###############################
# 2. Global models
###############################
names(occ_new)
nullmodel<-glm(poly_id~1,
               data=occ_new,
               family="binomial"(link=logit))
spcmodel1<-glm(poly_id~ SpRich,
               data=occ_new,
               family="binomial"(link=logit))
spcmodel2<-glm(poly_id~TAP01+MAT01+ SpRich,
               data=occ_new,
               family="binomial"(link=logit))
spcmodel3<-glm(poly_id~TAP01+I(TAP01^2)+MAT01+ SpRich,
               data=occ_new,
               family="binomial"(link=logit))
spcmodel4<-glm(poly_id~TAP01*MAT01+ SpRich,
               data=occ_new,family="binomial"(link=logit))



###############################
# 3. Global models extra
###############################

cMATmodel<-glm(poly_id~MAT01+ SpRich,
               data=occ_new,family="binomial"(link=logit))
cMAPmodel<-glm(poly_id~TAP01+ SpRich,
               data=occ_new,family="binomial"(link=logit))
cMAP2model<-glm(poly_id~I(TAP01^2)+TAP01+ SpRich,
                data=occ_new,family="binomial"(link=logit))
cMAT2model<-glm(poly_id~I(MAT01^2)+MAT01+ SpRich,
                data=occ_new,family="binomial"(link=logit))


# Test the difference between convergent and non-covergent models If trying distict optimizers 
# MAPmodel1 <- some alternative model with optimizer
# aa.OK <- list(MAPmodel1,MAPmodel)
# (lliks <- sort(sapply(aa.OK,logLik)))

AIC <- AIC( nullmodel, spcmodel1,
            spcmodel2, spcmodel3,
            cMATmodel,
            cMAPmodel,cMAP2model,
            cMAT2model)
AIC[order(AIC$AIC),]

## Compare both best fit global models
tab_model(cMAT2model, spcmodel3, df.method = "wald")

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'Sans',  
          axis.textcolor.x = "black",
          geom.label.color = "black",
          axis.textcolor.y = "black", #To change the font type
          axis.title.size = 2.0,  #To change axis title size
          axis.textsize.x = 1.3,
          axis.textcolor = "black" ,#To change x axis text size
          axis.textsize.y = 1.3)  #To change y axis text size



sjPlot::plot_models( spcmodel1,
                    spcmodel2, 
                    spcmodel3,
                    cMATmodel,
                    cMAPmodel,
                    cMAP2model,
                    cMAT2model,
                    axis.labels  = c("MAT^2", "TAP^2", "MAT", "TAP", "SpRich"),
                    m.labels= c("~ SpRich",
                               "~ TAP + MAT + SpRich",
                               "~TAP + TAP^2 + MAT + SpRich",
                               "~MAT + SpRich",
                               "~TAP01 + SpRich",
                               "~TAP^2 + TAP+ SpRich",
                              "~MAT^2 + MAT + SpRich"
                               ),
                    show.p = T,
                    p.shape = T ,
                    grid = F)










summary(cMAT2model)
sink("reviewProceedingsB/MAT2_modelNORandom+Rich.txt")
print(summary(spamodel4))
sink()

#install.packages("AICcmodavg")
library(AICcmodavg)

model.names <- c("nullmodel","spcmodel1" , "spcmodel2", "spcmodel3","spcmodel4",
                 "cMATmodel", "cMAPmodel", "cMAP2model", "cMAT2model")
models<-list( nullmodel,spcmodel1 , spcmodel2, spcmodel3, spcmodel4,
              cMATmodel,cMAPmodel,cMAP2model,cMAT2model)
aictab(cand.set = models, modnames = model.names)
aict<-aictab(cand.set = models, modnames = model.names)
write.csv(aict, file="reviewProceedingsB/aict_globalModelsNoRandom+rich.csv")
tab_model(cMAT2model, method = "wald")
best.model<-summary(cMAT2model)



###############################
# 2.1 Global models testing the effect of Pheidole and Pheidole + Camponotus
###############################

cmod_noP <- fitModelsNoRanplusRich(occ_new[occ_new$genus != "Pheidole",])
cmod_noC <- fitModelsNoRanplusRich(occ_new[occ_new$genus != "Camponotus",])
cmod_noP_noC <- fitModelsNoRanplusRich(occ_new[!occ_new$genus %in% c("Pheidole","Camponotus"),])


caic_mod_noC <- aictab(cand.set = cmod_noC[1:8], modnames = model.names[1:8])
caic_mod_noP <- aictab(cand.set = cmod_noP[1:8], modnames = model.names[1:8])
caic_mod_noP_noC <- aictab(cand.set = cmod_noP_noC[1:8], modnames = model.names[1:8])


write.csv(caic_mod_noC, file="reviewProceedingsB/aict_globalModelsNoRam_noCamponotusplusRich.csv")
write.csv(caic_mod_noP, file="reviewProceedingsB/aict_globalModelsNoRam_noPheidoleplusRich.csv")
write.csv(caic_mod_noP_noC, file="reviewProceedingsB/aict_globalModelsNoRam_noPheidole_noCamponotusplusRich.csv")


## get the output of the models without Camponotus  


sink("reviewProceedingsB/sp4_model_noCplusRich.txt")
print(summary(mod_noP[[5]]))
sink()
## get the output of the models without Pheidole
summary(aic_mod_noP)
sink("reviewProceedingsB/sp4_model_noPplusRich.txt")
print(summary(aic_mod_noP[[5]]))
sink()

## get the output of the models without camponotus and pheidole

summary(aic_mod_noP_noC)
sink("reviewProceedingsB/sp4_model_noP_noCplusRich.txt")
print(summary(aic_mod_noP_noC[[5]]))
sink()

summary(aic_mod_noP_noC[[5]])


###############################
# 3.1 Global models extra with the effect of subfamily 
###############################

########
# 3.1.1. 

occ_new$TAP01 <- round(occ_new$TAP01,2)
occ_new$MAT01 <- round(occ_new$MAT01,2)

#alternative best fit models

cspamodel3.a <- environment(list())

for(i in 1:length(unique(occ_new$SubFamily))){
  
  sp <- unique(occ_new$SubFamily)[i]
  
  if(dim(table(occ_new$poly_id[occ_new$SubFamily == sp]))>1){
    
    print(paste("doing", sp))
    
    tryCatch({bspamodel3.a[[i]] <- glm(poly_id~TAP01 + I(TAP01^2)+MAT01,
                                       data=droplevels(occ_new[occ_new$SubFamily ==sp,]),
                                       family="binomial"(link=logit))},
             error = function(e){
               bspamodel3.a[[i]] <- c("NPV")
             })
    
    print("done")
    
    
  }
  else{
    print(paste("Non poly var", sp))
    bspamodel3.a[[i]] <- "NPV"
    
  }
  
  
  if(i == length(unique(occ_new$SubFamily))){
    cat("\014")
  }
}



spamodel4.a <- list()
for(i in 1:length(unique(occ_new$SubFamily))){
  sp <- unique(occ_new$SubFamily)[i]
  
  if(dim(table(occ_new$poly_id[occ_new$SubFamily == sp]))>1){
    print(paste("doing", sp))
    spamodel4.a[[i]] <- glmer(poly_id~TAP01*MAT01+ (1|new_id),
                              data=droplevels(occ_new[occ_new$SubFamily ==sp,]),
                              family="binomial"(link=logit),
                              glmerControl(optimizer = "Nelder_Mead"))
    
    print("done")
    
    
  }
  else{
    print(paste("Non poly var", sp))
    spamodel4.a[[i]] <- "NPV"
    
  }
  
  
  if(i == length(unique(occ_new$SubFamily))){
    cat("\014")
  }
}

spamodel2.a <- list()
for(i in 1:length(unique(occ_new$SubFamily))){
  sp <- unique(occ_new$SubFamily)[i]
  
  if(dim(table(occ_new$poly_id[occ_new$SubFamily == sp]))>1){
    print(paste("doing", sp))
    spamodel2.a[[i]] <- glmer(poly_id~TAP01+MAT01+(1|new_id),
                              data=droplevels(occ_new[occ_new$SubFamily ==sp,]),
                              family="binomial"(link=logit),
                              glmerControl(optimizer = "Nelder_Mead"))
    
    print("done")
    
    
  }
  else{
    print(paste("Non poly var", sp))
    spamodel2.a[[i]] <- "NPV"
    
  }
  
  
  if(i == length(unique(occ_new$SubFamily))){
    cat("\014")
  }
}

MATmodel <- list()
for(i in 1:length(unique(occ_new$SubFamily))){
  sp <- unique(occ_new$SubFamily)[i]
  
  if(dim(table(occ_new$poly_id[occ_new$SubFamily == sp]))>1){
    print(paste("doing", sp))
    MATmodel[[i]] <- glmer(poly_id~MAT01+(1|new_id),
                           data=droplevels(occ_new[occ_new$SubFamily ==sp,]),
                           family="binomial"(link=logit),
                           glmerControl(optimizer = "Nelder_Mead"))
    
    print("done")
    
    
  }
  else{
    print(paste("Non poly var", sp))
    MATmodel[[i]] <- "NPV"
    
  }
  
  
  if(i == length(unique(occ_new$SubFamily))){
    cat("\014")
  }
}



# calculate partial AIC's
#####

aic3_4 <- list()
for(i in 1:12){
  if(is.character(spamodel3.a[[i]]) == F){
    aic3_4[[i]] <- AIC( spamodel3.a[[i]], 
                        spamodel4.a[[i]])
    
  }else{
    aic3_4[[i]] <- "NPV"
    
  }
}




## Compare both best fit global models
#tab_model(spamodel2, spamodel3, spamodel4)
# Observe the output of partial models 
names(aic3_4) <- unique(occ_new$SubFamily)

aic3_4[[4]][order(aic3_4[[4]]$AIC),]

#which subfamilies have poly 1 and 0?
sapply(unique(occ_new$SubFamily), function(x) dim(table(occ_new$poly_id[occ_new$SubFamily == x]))>1)


subSet <- which(apply(sapply(1:length(unique(occ_new$SubFamily)), function(x)
  c(is.character(spamodel2.a[[x]]),
    is.character(spamodel3.a[[x]]),
    is.character(spamodel3.a[[x]]),
    is.character(MATmodel[[x]]))),
  2,mean) == 0)



aic2_3_4 <- list()
for(i in subSet){
  
  aic2_3_4[[i]] <- aictab(cand.set = c(spamodel2.a[[i]],
                                       spamodel3.a[[i]],
                                       spamodel4.a[[i]],
                                       MATmodel[[i]]),
                          modnames = c("model2", "model3", 
                                       "model4", "modelMAT"))
}

names(aic2_3_4) <- unique(occ_new$SubFamily)[1:max(subSet)]


aic_subF <- data.frame(dplyr::bind_rows(aic2_3_4),
                       "names"= c(sapply(names(aic2_3_4)[subSet], function(x) rep(x,4))))

write.csv(aic_subF, file="reviewProceedingsB/aict_SubFamplusRich.csv")

# 
# 
# 
# 
# #########
# ###############################
# # 3.1 Global models extra with the effect of genus 
# ###############################
# 
# 
# unique(occ_new$genus)
# 
# names(occ_new)
# spamodel3.b <- list()
# for(i in 1:length(unique(occ_new$genus))){
#   sp <- unique(occ_new$genus)[i]
#   con <- data.frame("x" = as.numeric(names(table(occ_new$poly_id[occ_new$genus == sp]))),
#                     "b" = as.numeric(names(table(occ_new$poly_id[occ_new$genus == sp])))==1)
#   
#   if( dim(con)[1] == 0){
#     print(paste("Non poly var", sp))
#     spamodel3.b[[i]] <- "NPV"
#   }
#   
#   if(length(c(1) %in% con$x == T & (con$b[con$x == 1] == 1)) != 0){
#     if((c(1) %in% con$x) == T & (con$b[con$x == 1] == 1)){
#       print(paste("doing", sp))
#       spamodel3.b[[i]] <- tryCatch({
#         glmer(poly_id~TAP01+I(TAP01^2)+MAT01+(1|new_id),
#               data=droplevels(occ_new[occ_new$genus ==sp,]),
#               family="binomial"(link=logit),
#               glmerControl(optimizer = "Nelder_Mead"))
#         
#       },error = function(e){
#         return("Error in model")
#         print("Error in model")
#       })
#       
#     }
#     
#     
#     print("done")
#     
#     
#   }
#   else{
#     print(paste("Non poly var", sp))
#     spamodel3.b[[i]] <- "NPV"
#     
#   }
#   
#   
#   if(i == length(unique(occ_new$genus))){
#     cat("\014")
#   }
# }
# names(spamodel3.b) <- unique(occ_new$genus)
# 
# spamodel4.b <- list()
# for(i in 1:length(unique(occ_new$genus))){
#   sp <- unique(occ_new$genus)[i]
#   con <- data.frame("x" = as.numeric(names(table(occ_new$poly_id[occ_new$genus == sp]))),
#                     "b" = as.numeric(names(table(occ_new$poly_id[occ_new$genus == sp])))==1)
#   
#   if( dim(con)[1] == 0){
#     print(paste("Non poly var", sp))
#     spamodel4.b[[i]] <- "NPV"
#   }
#   
#   if(length(c(1) %in% con$x == T & (con$b[con$x == 1] == 1)) != 0){
#     if((c(1) %in% con$x) == T & (con$b[con$x == 1] == 1)){
#       print(paste("doing", sp))
#       spamodel4.b[[i]] <- tryCatch({
#         glmer(poly_id~TAP01*MAT01+(1|new_id),
#               data=droplevels(occ_new[occ_new$genus ==sp,]),
#               family="binomial"(link=logit),
#               glmerControl(optimizer = "Nelder_Mead"))
#         
#       },error = function(e){
#         return("Error in model")
#         print("Error in model")
#       })
#       
#     }
#     
#     
#     print("done")
#     
#     
#   }
#   else{
#     print(paste("Non poly var", sp))
#     spamodel4.b[[i]] <- "NPV"
#     
#   }
#   
#   
#   if(i == length(unique(occ_new$genus))){
#     cat("\014")
#   }
# }
# 
# names(spamodel4.b) <- unique(occ_new$genus)
# 
# spamodel2.b <- list()
# for(i in 1:length(unique(occ_new$genus))){
#   sp <- unique(occ_new$genus)[i]
#   con <- data.frame("x" = as.numeric(names(table(occ_new$poly_id[occ_new$genus == sp]))),
#                     "b" = as.numeric(names(table(occ_new$poly_id[occ_new$genus == sp])))==1)
#   
#   if( dim(con)[1] == 0){
#     print(paste("Non poly var", sp))
#     spamodel2.b[[i]] <- "NPV"
#   }
#   
#   if(length(c(1) %in% con$x == T & (con$b[con$x == 1] == 1)) != 0){
#     if((c(1) %in% con$x) == T & (con$b[con$x == 1] == 1)){
#       print(paste("doing", sp))
#       spamodel2.b[[i]] <- tryCatch({
#         glmer(poly_id~TAP01+MAT01+(1|new_id),
#               data=droplevels(occ_new[occ_new$genus ==sp,]),
#               family="binomial"(link=logit),
#               glmerControl(optimizer = "Nelder_Mead"))
#         
#       },error = function(e){
#         return("Error in model")
#         print("Error in model")
#       })
#       
#     }
#     
#     
#     print("done")
#     
#     
#   }
#   else{
#     print(paste("Non poly var", sp))
#     spamodel2.b[[i]] <- "NPV"
#     
#   }
#   
#   
#   if(i == length(unique(occ_new$genus))){
#     cat("\014")
#   }
# }
# names(spamodel2.b) <- unique(occ_new$genus)
# 
# MATmodel.b <- list()
# for(i in 1:length(unique(occ_new$genus))){
#   sp <- unique(occ_new$genus)[i]
#   con <- data.frame("x" = as.numeric(names(table(occ_new$poly_id[occ_new$genus == sp]))),
#                     "b" = as.numeric(names(table(occ_new$poly_id[occ_new$genus == sp])))==1)
#   
#   if( dim(con)[1] == 0){
#     print(paste("Non poly var", sp))
#     MATmodel.b[[i]] <- "NPV"
#   }
#   
#   if(length(c(1) %in% con$x == T & (con$b[con$x == 1] == 1)) != 0){
#     if((c(1) %in% con$x) == T & (con$b[con$x == 1] == 1)){
#       print(paste("doing", sp))
#       MATmodel.b[[i]] <- tryCatch({
#         glmer(poly_id~MAT01+(1|new_id),
#               data=droplevels(occ_new[occ_new$genus ==sp,]),
#               family="binomial"(link=logit),
#               glmerControl(optimizer = "Nelder_Mead"))
#         
#       },error = function(e){
#         return("Error in model")
#         print("Error in model")
#       })
#       
#     }
#     
#     
#     print("done")
#     
#     
#   }
#   else{
#     print(paste("Non poly var", sp))
#     MATmodel.b[[i]] <- "NPV"
#     
#   }
#   
#   
#   if(i == length(unique(occ_new$genus))){
#     cat("\014")
#   }
# }
# names(MATmodel.b) <- unique(occ_new$genus)
# 
# 
# # calculate partial AIC's for genus
# #####
# 
# 
# subSet <- which(apply(sapply(1:length(unique(occ_new$genus)), function(x)
#   c(is.character(spamodel2.b[[x]]),
#     is.character(spamodel3.b[[x]]),
#     is.character(spamodel2.b[[x]]),
#     is.character(MATmodel.b[[x]])
#     )),
#   2,mean) == 0)
# 
# 
# 
# aic2_3_4 <- list()
# for(i in subSet){
#   
#   aic2_3_4[[i]] <- aictab(cand.set = c(spamodel2.b[[i]],spamodel3.b[[i]],spamodel4.b[[i]], MATmodel.b[[i]]),
#                           modnames = c("model2", "model3", "model4", "MATModel"))
# }
# 
# names(aic2_3_4) <- unique(occ_new$genus)[1:max(subSet)]
# 
# 
# aic_gen <- data.frame(dplyr::bind_rows(aic2_3_4),
#                       "names"= c(sapply(names(aic2_3_4)[subSet], function(x) rep(x,4))))
# 
# 
# 
# aic_gen$SubFamily <- occ$subFamily[match(aic_gen$names, occ$genus)]
# 
# write.csv(aic_gen, file="aict_genLev.csv")
# 
# 
# ###############################
# # 4. Global model plots
# ###############################
# 
# library(ggplot2)
# library(sjPlot)
# library(ggpubr)
# 
# a<-plot_model(spamodel4,type="pred",terms = "MAT01")+
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text=element_text(family="Helvetica", face="bold", size=20))+
#   xlab("Temperature")+
#   ylab("Worker Caste Polymorphism")+
#   ggtitle("")+
#   scale_y_continuous(breaks = c(0,0.5,1),
#                      labels=scales::percent_format(accuracy = 1),limits=c(0,1));a
# 
# b<-plot_model(spamodel4,type="pred",terms = "TAP01")+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text=element_text(family="Helvetica", face="bold", size=20))+
#   xlab("Precipitation")+
#   ylab("")+
#   ggtitle("")+
#   scale_y_continuous(breaks = c(0,0.5, 1),
#                      labels=scales::percent_format(accuracy = 1),
#                      limits = c(0,1));b
# par(oma=c(2,2,2,2))
# ggarrange(a,b)
