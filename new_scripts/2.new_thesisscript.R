##### THESIS MODELS #####

###############################
# 1. Load data 
###############################

library(vegan)
library(adespatial)
library(sp)
library(nlme)
library(lme4)

occ_new<-read.csv("occ_withsitedensity.csv",stringsAsFactors = F)

###############################
# 2. Global models
###############################

nullmodel<-glm(poly_id~1,data=occ_new,family="binomial"(link=logit))
spamodel1<-glmer(poly_id~1+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

###############################
# 3. Global models extra
###############################

MATmodel<-glmer(poly_id~MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAPmodel<-glmer(poly_id~MAP.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAP2model<-glmer(poly_id~I(MAP.point^2)+MAP.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAT2model<-glmer(poly_id~I(MAT.point^2)+MAT.point+(1|new_id),data=occ_new,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

###############################
# 4. Global model plots
###############################

library(ggplot2)
library(sjPlot)
library(ggpubr)

a<-plot_model(spamodel3,type="pred",terms = "MAT.point")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Temperature")+
  ylab("Predicted Probabilities of Worker Caste Polymorphism")+
  ggtitle("")+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=scales::percent_format(accuracy = 1),limits=c(0,0.7));a

b<-plot_model(spamodel3,type="pred",terms = "MAP.point")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Precipitation")+
  ylab("Predicted Probabilities of Worker Caste Polymorphism")+
  ggtitle("")+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels=scales::percent_format(accuracy = 1), limits = c(0,NA));b

ggarrange(a,b)

###############################
# 5. Comparison with Blanchard & Moreau (2017)
###############################

occ_mor<-read.csv("new_occ_df_moreau.csv",stringsAsFactors = F)

nullmodel2<-glm(poly_id_blanchard~1,data=occ_mor,family="binomial"(link=logit))
spamodel1_m<-glmer(poly_id_blanchard~1+(1|new_id),data=occ_mor,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2_m<-glmer(poly_id_blanchard~MAP.point+MAT.point+(1|new_id),data=occ_mor,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3_m<-glmer(poly_id_blanchard~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=occ_mor,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4_m<-glmer(poly_id_blanchard~MAP.point*MAT.point+(1|new_id),data=occ_mor,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

MATmodel_m<-glmer(poly_id_blanchard~MAT.point+(1|new_id),data=occ_mor,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAPmodel_m<-glmer(poly_id_blanchard~MAP.point+(1|new_id),data=occ_mor,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAP2model_m<-glmer(poly_id_blanchard~I(MAP.point^2)+MAP.point+(1|new_id),data=occ_mor,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAT2model_m<-glmer(poly_id_blanchard~I(MAT.point^2)+MAT.point+(1|new_id),data=occ_mor,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

###############################
# 6. test plots
###############################

library(ggplot2)

fit<-fitted(spamodel3)
fit<-na.omit(fit)
resid<-residuals(spamodel3)
resid<-na.omit(resid)

pred<-predict(spamodel3)
pred<-na.omit(pred)

MAP<-na.omit(occ_new$MAP.point)
MAT<-na.omit(occ_new$MAT.point)


mydf<-ggpredict(spamodel3, terms = "MAT.point")

ggplot(mydf, aes(x, predicted)) + 
  geom_line()+
  geom_point(aes(x=MAP,y=fit))

x<-ggplot()+
  geom_point(aes(x=MAT,y=pred))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Temperature")+
  ylab("Predicted Values of Polymorphism")+
  ggtitle("")+
  stat_smooth();x

###############################
# 7. Models without Pheidole
###############################

occ_genus<-read.csv("occ_withsitedensity.csv",stringsAsFactors = F)

nopheidole<-subset(occ_genus,occ_genus$valid_species_name!="Pheidole")
summary(nopheidole) # eliminated 700 entries

nullmodel_nop<-glm(poly_id~1,data=nopheidole,family="binomial"(link=logit))
spamodel1_nop<-glmer(poly_id~1+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2_nop<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3_nop<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4_nop<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

MATmodel_nop<-glmer(poly_id~MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAPmodel_nop<-glmer(poly_id~MAP.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAP2model_nop<-glmer(poly_id~I(MAP.point^2)+MAP.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAT2model_nop<-glmer(poly_id~I(MAT.point^2)+MAT.point+(1|new_id),data=nopheidole,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

library(ggplot2)
library(sjPlot)
library(ggpubr)
library(ggeffects)

a_1<-plot_model(spamodel3_nop,type="pred",terms = "MAT.point")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Temperature")+
  ylab("Predicted Probabilities of Worker Caste Polymorphism")+
  ggtitle("")+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=scales::percent_format(accuracy = 1),limits=c(0,0.7));a

b_1<-plot_model(spamodel3_nop,type="pred",terms = "MAP.point")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Precipitation")+
  ylab("Predicted Probabilities of Worker Caste Polymorphism")+
  ggtitle("")+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels=scales::percent_format(accuracy = 1), limits = c(0,NA));b_1

x_1<-ggarrange(a_1,b_1)

###############################
# 8. Models without Camponotus
###############################

occ_genus<-read.csv("occ_withsitedensity.csv",stringsAsFactors = F)
nocamp<-subset(occ_genus,occ_genus$valid_species_name!="Camponotus")

nullmodel_noc<-glm(poly_id~1,data=nocamp,family="binomial"(link=logit))
spamodel1_noc<-glmer(poly_id~1+(1|new_id),data=nocamp,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2_noc<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=nocamp,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3_noc<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=nocamp,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4_noc<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=nocamp,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

MATmodel_noc<-glmer(poly_id~MAT.point+(1|new_id),data=nocamp,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAPmodel_noc<-glmer(poly_id~MAP.point+(1|new_id),data=nocamp,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAP2model_noc<-glmer(poly_id~I(MAP.point^2)+MAP.point+(1|new_id),data=nocamp,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAT2model_noc<-glmer(poly_id~I(MAT.point^2)+MAT.point+(1|new_id),data=nocamp,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

a_2<-plot_model(spamodel3_noc,type="pred",terms = "MAT.point")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Temperature")+
  ylab("Predicted Probabilities of Worker Caste Polymorphism")+
  ggtitle("")+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=scales::percent_format(accuracy = 1),limits=c(0,0.7));a

b_2<-plot_model(spamodel3_noc,type="pred",terms = "MAP.point")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Precipitation")+
  ylab("Predicted Probabilities of Worker Caste Polymorphism")+
  ggtitle("")+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels=scales::percent_format(accuracy = 1), limits = c(0,NA));b_1

x_2<-ggarrange(a_2,b_2)

###############################
# 9. Models without Pheidole and Camponotus
###############################

occ_genus<-read.csv("occ_withsitedensity.csv",stringsAsFactors = F)
nopc<-subset(nopheidole,nopheidole$valid_species_name!="Camponotus")

nullmodel_nopc<-glm(poly_id~1,data=nopc,family="binomial"(link=logit))
spamodel1_nopc<-glmer(poly_id~1+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2_nopc<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3_nopc<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4_nopc<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

MATmodel_nopc<-glmer(poly_id~MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAPmodel_nopc<-glmer(poly_id~MAP.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAP2model_nopc<-glmer(poly_id~I(MAP.point^2)+MAP.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAT2model_nopc<-glmer(poly_id~I(MAT.point^2)+MAT.point+(1|new_id),data=nopc,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

a_3<-plot_model(spamodel3_nopc,type="pred",terms = "MAT.point")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Temperature")+
  ylab("Predicted Probabilities of Worker Caste Polymorphism")+
  ggtitle("")+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=scales::percent_format(accuracy = 1),limits=c(0,0.7));a

b_3<-plot_model(spamodel3_nopc,type="pred",terms = "MAP.point")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text=element_text(family="Times New Roman", face="bold", size=20))+
  xlab("Precipitation")+
  ylab("Predicted Probabilities of Worker Caste Polymorphism")+
  ggtitle("")+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels=scales::percent_format(accuracy = 1), limits = c(0,NA));b_1

x_3<-ggarrange(a_3,b_3)

###############################
# 10. Models without Monophasic
###############################

#species removed or changed: Cataglyphis.floricola,Crematogaster.coriaria,Erromyrma.latinodis
#Monomorium.centrale,Monomorium.insolescens,Monomorium.rubriceps,Monomorium.rufonigrum,Monomorium.whitei,Paraponera.clavata
#Trichomyrmex.destructor,Trichomyrmex.emeryi,Trichomyrmex.mayri,Trichomyrmex.oscaris,Solenopsis.invicta

minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Cataglyphis.floricola")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Crematogaster.coriaria")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Erromyrma.latinodis")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Monomorium.centrale")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Monomorium.insolescens")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Monomorium.rubriceps")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Monomorium.rufonigrum")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Monomorium.whitei")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Paraponera.clavata")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Trichomyrmex.destructor")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Trichomyrmex.emeryi")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Trichomyrmex.mayri")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Trichomyrmex.oscaris")
minusmonophasic<-subset(occ_new,occ_new$valid_species_name!="Solenopsis.invicta")

nullmodel_nomono<-glm(poly_id~1,data=minusmonophasic,family="binomial"(link=logit))
spamodel1_nomono<-glmer(poly_id~1+(1|new_id),data=minusmonophasic,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2_nomono<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=minusmonophasic,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3_nomono<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=minusmonophasic,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4_nomono<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=minusmonophasic,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

MATmodel_nomono<-glmer(poly_id~MAT.point+(1|new_id),data=minusmonophasic,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAPmodel_nopmono<-glmer(poly_id~MAP.point+(1|new_id),data=minusmonophasic,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAP2model_nomono<-glmer(poly_id~I(MAP.point^2)+MAP.point+(1|new_id),data=minusmonophasic,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAT2model_nomono<-glmer(poly_id~I(MAT.point^2)+MAT.point+(1|new_id),data=minusmonophasic,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

## Changing classification ## 

vector<-c("Cataglyphis.floricola","Crematogaster.coriaria","Erromyrma.latinodis","Monomorium.centrale","Monomorium.insolescens","Monomorium.rubriceps","Monomorium.rufonigrum","Monomorium.whitei","Paraponera.clavata","Trichomyrmex.destructor","Trichomyrmex.emeryi","Trichomyrmex.mayri","Trichomyrmex.oscaris","Solenopsis.invicta")
data<-occ_new
data$poly_id[which(data$valid_species_name %in% vector)] <- 0

nullmodel_mono<-glm(poly_id~1,data=data,family="binomial"(link=logit))
spamodel1_mono<-glmer(poly_id~1+(1|new_id),data=data,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel2_mono<-glmer(poly_id~MAP.point+MAT.point+(1|new_id),data=data,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel3_mono<-glmer(poly_id~MAP.point+I(MAP.point^2)+MAT.point+(1|new_id),data=data,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
spamodel4_mono<-glmer(poly_id~MAP.point*MAT.point+(1|new_id),data=data,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))

MATmodel_mono<-glmer(poly_id~MAT.point+(1|new_id),data=data,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAPmodel_pmono<-glmer(poly_id~MAP.point+(1|new_id),data=data,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAP2model_mono<-glmer(poly_id~I(MAP.point^2)+MAP.point+(1|new_id),data=data,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
MAT2model_mono<-glmer(poly_id~I(MAT.point^2)+MAT.point+(1|new_id),data=data,family="binomial"(link=logit),glmerControl(optimizer = "Nelder_Mead"))
