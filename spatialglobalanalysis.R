#### PROJECT MODELS #### 

library(vegan)
library(adespatial)
library(sp)
library(nlme)
library(lme4)
library(MASS)
library(MuMIn)
library(adephylo)

### 1. LOAD DATA AND PREPARE ### 

occ<-read.csv("new_occ_df.csv",stringsAsFactors = FALSE) #to load from the saved file 
#poly_id (response variable) has 1 or 0 for polymorphism
#continent  has the continent categorical value
#new id is site ID
#point_density is sampling density 
#MAT.point = temperature data that is already scaled
#MAP.point = precipitation data that is already scaled

occ<-trial
### 2. SPATIAL GLOBAL MODEL ### 

nullmodel<-glm(poly_id~1,data=occ,family="binomial"(link=logit))
spamodel1<-glmer(poly_id~1+(1|new_id),data=occ,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=occ,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=occ,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=occ,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

spaMAT<-glmer(poly_id~MAT.point+(1|new_id),data=occ,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAP<-glmer(poly_id~MAP.point+(1|new_id),data=occ,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAP2<-glmer(poly_id~MAP.point+I(MAP.point^2)+(1|new_id),data=occ,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAT2<-glmer(poly_id~MAT.point+I(MAT.point^2)+(1|new_id),data=occ,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))


### 3. PLOTTING DATA AND MODELS ###

library(ggplot2)
library(sjPlot)
library(ggpubr)

hist(s$MAP.point, xlab="Precipitation", main="Precipitation")
hist(s$MAT.point,xlab="Temperature",main="Temperature")

a<-plot_model(spamodel4,type="pred",terms="MAP.point")+
  xlab("Precipitation")+
  ylab("Predicted Probability of Polymorphism")+
  theme_bw()+
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black"))+
  ggtitle("")

b<-plot_model(spamodel4,type="pred",terms="MAT.point")+
  xlab("Temperature")+
  ylab("Predicted Probability of Polymporphism")+
  theme_bw()+
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black"))+
  ggtitle("")

c<-plot_model(spamodel4,type="pred",terms="MAP.point")+
  xlab("Precipitation")+
  ylab("Predicted Probability of Polymorphism")+
  theme_bw()+
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black"))+
  ggtitle("")


ggarrange(plotlist=list(a,b))

### 4. Pheidole and Campontus Models ###
library(lme4)
occ2<-read.csv("new_occ_df_genus.csv",stringsAsFactors = TRUE)

## Camponotus Removed Models ##
nocampo<-occ2[which(occ2$genus != "Camponotus"),]
summary(nocampo)

nullmodelcamp<-glm(poly_id~1,data=nocampo,family="binomial"(link=logit))
spamodel1camp<-glmer(poly_id~1+(1|new_id),data=nocampo,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2camp<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=nocampo,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3camp<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=nocampo,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4camp<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=nocampo,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

spaMATcamp<-glmer(poly_id~MAT.point+(1|new_id),data=nocampo,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAPcamp<-glmer(poly_id~MAP.point+(1|new_id),data=nocampo,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAP2camp<-glmer(poly_id~MAP.point+I(MAP.point^2)+(1|new_id),data=nocampo,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAT2camp<-glmer(poly_id~MAT.point+I(MAT.point^2)+(1|new_id),data=nocampo,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

## Pheidole Removed Models ##
nopheidole<-occ2[which(occ2$genus != "Pheidole"),]
summary(nopheidole)

nullmodelphei<-glm(poly_id~1,data=nopheidole,family="binomial"(link=logit))
spamodel1phei<-glmer(poly_id~1+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2phei<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3phei<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4phei<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

spaMATphei<-glmer(poly_id~MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAPphei<-glmer(poly_id~MAP.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAP2phei<-glmer(poly_id~MAP.point+I(MAP.point^2)+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAT2phei<-glmer(poly_id~MAT.point+I(MAT.point^2)+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

## Pheidole and Camp Removed Models ##

nopc<-occ2[which(occ2$genus != c("Pheidole","Camponotus")),]
summary(nopc)

nullmodelpc<-glm(poly_id~1,data=nopc,family="binomial"(link=logit))
spamodel1pc<-glmer(poly_id~1+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2pc<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3pc<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4pc<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

spaMATpc<-glmer(poly_id~MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAPpc<-glmer(poly_id~MAP.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAP2pc<-glmer(poly_id~MAP.point+I(MAP.point^2)+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spaMAT2pc<-glmer(poly_id~MAT.point+I(MAT.point^2)+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

