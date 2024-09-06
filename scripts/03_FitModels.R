######################
## Fit linear mixed effect models 
## Gabriel Muñoz 
## DIC 2021
######################

source("reviewProceedingsB/02_EqualAreaAntAssemblages.R")
source("reviewProceedingsB/00_Functions.R")
###############################
# Package dependencies
###############################

library(vegan)
library(adespatial)
library(sp)
library(nlme)
library(lme4)
library(AICcmodavg)
library(sjPlot)

occ_new <- AntData
# make sure there is no NAs
occ_new <- occ_new[complete.cases(occ_new),]
str(occ_new)

dev.off()
plot(occ_new$x_grid,
     occ_new$y_grid,
     pch = ".")

names(occ_new)

####
# Transform climate data so it varies on the same range 
#########
occ_new$MAT01 <- vegan::decostand(occ_new$MAT, "range")
occ_new$TAPlog01 <- vegan::decostand(log1p(occ_new$TAP), "range")
occ_new$TAP01 <- vegan::decostand((occ_new$TAP), "range")

occ_new$ID_ea


###############################
# Global models
###############################
names(occ_new)

nullmodel<-glm(poly_id~1,
               data=occ_new,
               family="binomial"(link=logit))

spamodel1<-glmer(poly_id~1+(1|ID_ea),
                 data=occ_new,
                 family="binomial"(link=logit),
                 glmerControl(optimizer = "Nelder_Mead"))

spamodel2<-glmer(poly_id~TAP01+MAT01+(1|ID_ea),
                 data=occ_new,family="binomial"(link=logit),
                 glmerControl(optimizer = "Nelder_Mead"))

spamodel3<-glmer(poly_id~TAP01+I(TAP01^2)+MAT01+(1|ID_ea),
                 data=occ_new,family="binomial"(link=logit),
                 glmerControl(optimizer = "Nelder_Mead"))

spamodel4<-glmer(poly_id~TAP01*MAT01+(1|ID_ea),
                 data=occ_new,family="binomial"(link=logit),
                 glmerControl(optimizer = "Nelder_Mead"))



###############################
# 3. Global models extra
###############################

MATmodel<-glmer(poly_id~MAT01+(1|ID_ea),
                data=occ_new,family="binomial"(link=logit),
                glmerControl(optimizer = "Nelder_Mead"))

MAPmodel<-glmer(poly_id~TAP01+(1|ID_ea),
                data=occ_new,family="binomial"(link=logit),
                glmerControl(optimizer = "Nelder_Mead"))

MAP2model<-glmer(poly_id~I(TAP01^2)+TAP01+(1|ID_ea),
                 data=occ_new,family="binomial"(link=logit),
                 glmerControl(optimizer = "Nelder_Mead"))

MAT2model<-glmer(poly_id~I(MAT01^2)+MAT01+(1|ID_ea),
                 data=occ_new,family="binomial"(link=logit),
                 glmerControl(optimizer = "Nelder_Mead"))




model.names <- c("spamodel1" ,
                 "spamodel2",
                 "spamodel3", 
                 "spamodel4",
                 "MATmodel", 
                 "MAPmodel", 
                 "MAP2model",
                 "MAT2model")

models<-list( spamodel1 , spamodel2, 
              spamodel3, spamodel4, 
              MATmodel,MAPmodel, 
              MAP2model, MAT2model)

aictab(cand.set = models, modnames = model.names)

aict<-aictab(cand.set = models, modnames = model.names)

write.csv(aict, file="reviewProceedingsB/aict_globalModels.csv")
best.model<-summary(spamodel4)
tab_model(spamodel4)

###############################
# 2.1 Global models testing the effect of Pheidole and Pheidole + Camponotus
###############################

mod_noP <- fitModels(occ_new[occ_new$genus != "Pheidole",])
mod_noC <- fitModels(occ_new[occ_new$genus != "Camponotus",])
mod_noP_noC <- fitModels(occ_new[!occ_new$genus %in% c("Pheidole","Camponotus"),])


aic_mod_noC <- aictab(cand.set = mod_noC[2:8], modnames = model.names[2:8])
aic_mod_noP <- aictab(cand.set = mod_noP[2:8], modnames = model.names[2:8])
aic_mod_noP_noC <- aictab(cand.set = mod_noP_noC[2:8], modnames = model.names[2:8])


write.csv(aic_mod_noC, file="reviewProceedingsB/aict_globalModels_noCamponotus.csv")
write.csv(aic_mod_noP, file="reviewProceedingsB/aict_globalModels_noPheidole.csv")
write.csv(aic_mod_noP_noC, file="reviewProceedingsB/aict_globalModels_noPheidole_noCamponotus.csv")


## get the output of the models without Camponotus  
sink("reviewProceedingsB/MATmodelnoC.txt")
print(summary(mod_noP[[5]]))
sink()
## get the output of the models without Pheidole
summary(aic_mod_noP)
sink("reviewProceedingsB/MATmodelnoP.txt")
print(summary(aic_mod_noP[[5]]))
sink()

## get the output of the models without camponotus and pheidole

summary(aic_mod_noP_noC)
sink("reviewProceedingsB/MATmodelnoCnoP.txt")
print(summary(aic_mod_noP_noC[[5]]))
sink()

###############################
# 3.1 Global models extra with the effect of subfamily 
###############################

########
# 3.1.1. 


occ_new$MAT01 <- round(occ_new$MAT01,2)
occ_new$TAP01 <- round(occ_new$TAP01,2)

#alternative best fit models

spamodel3.a <- environment(list())

for(i in 1:length(unique(occ_new$SubFamily))){
  
  sp <- unique(occ_new$SubFamily)[i]
  
  if(dim(table(occ_new$poly_id[occ_new$SubFamily == sp]))>1){
    
    print(paste("doing", sp))
    
    tryCatch({spamodel3.a[[i]] <- glmer(poly_id~TAP01 + I(TAP01^2)+MAT01+ (1|ID_ea),
                              data=droplevels(occ_new[occ_new$SubFamily ==sp,]),
                              family="binomial"(link=logit),
                              glmerControl(optimizer = "Nelder_Mead"))},
            error = function(e){
              spamodel3.a[[i]] <- c("NPV")
             })
    
    print("done")
    
    
  }
  else{
    print(paste("Non poly var", sp))
    spamodel3.a[[i]] <- "NPV"
    
  }
  
  
  if(i == length(unique(occ_new$SubFamily))){
    cat("\014")
  }
}

spamodel3.a[[which(sapply(1:12, function(x) is.null(spamodel3.a[[x]]))==T)]] <- "NPV"

spamodel4.a <- list()
for(i in 1:length(unique(occ_new$SubFamily))){
  sp <- unique(occ_new$SubFamily)[i]
  
  if(dim(table(occ_new$poly_id[occ_new$SubFamily == sp]))>1){
    print(paste("doing", sp))
    spamodel4.a[[i]] <- glmer(poly_id~TAP01*MAT01+ (1|ID_ea),
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
    spamodel2.a[[i]] <- glmer(poly_id~TAP01+MAT01+(1|ID_ea),
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
    MATmodel[[i]] <- glmer(poly_id~MAT01+(1|ID_ea),
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
for(i in 1:11){
  print(i)
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

write.csv(aic_subF, file="reviewProceedingsB/aict_SubFam.csv")
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
#         glmer(poly_id~MAP.point+I(MAP.point^2)+MAT01+(1|ID_ea),
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
#         glmer(poly_id~MAP.point*MAT01+(1|ID_ea),
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
#         glmer(poly_id~MAP.point+MAT01+(1|ID_ea),
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
#         glmer(poly_id~MAT01+(1|ID_ea),
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
# b<-plot_model(spamodel4,type="pred",terms = "MAP.point")+
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
