setwd("~/Desktop/NASSTEC SESIL 2A/Obj 2 SF NASSTEC Germination Niche/Obj 2 data/germination analysis files")
# Install and load the necessary packages
install.packages("effects")
install.packages("binom")
library(effects)
library ("binom")
# Read the dataset
proportionsPEG <- read.table("proportionsPEG.txt", header = T)
# Go through the data, seed population by seed population. Species names are abbreviated with the first two letters each of the genus and specific epithet. Anthemis cotula = ANCO.

#________Anthemis cotula
# Subset data for the seed population of interest
ANCO<-subset(proportionsPEG,Species == "ANCO")
# Define the generalized linear model (full model)
ANCOm1<-glm(cbind(ANCO$Germinated,ANCO$Germinable - ANCO$Germinated) ~  MPa, data = ANCO, family=binomial) #we set up the germinated and not germinated as the response variable

summary(ANCOm1) 
data.frame(Effect(c("MPa"),ANCOm1)$x, 
           Est.=plogis(Effect(c("MPa"),ANCOm1)$fit), lower=plogis(Effect(c("MPa"),ANCOm1)$lower), 
           upper = plogis(Effect(c("MPa"),ANCOm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#Try to group some MPas.

ANCO0.1.2.3newMPa<-ANCO$MPa
levels(ANCO0.1.2.3newMPa)
levels(ANCO0.1.2.3newMPa)[c(1,2,3,4)]<-"0.1.2.3"

ANCOm2<-glm(cbind(ANCO$Germinated,ANCO$Germinable-ANCO$Germinated)~ ANCO0.1.2.3newMPa,data=ANCO,family=binomial)
summary(ANCOm2)
anova(ANCOm1,ANCOm2, test="Chi")
data.frame(Effect(c("ANCO0.1.2.3newMPa"),ANCOm2)$x, 
           Est.=plogis(Effect(c("ANCO0.1.2.3newMPa"),ANCOm2)$fit), lower=plogis(Effect(c("ANCO0.1.2.3newMPa"),ANCOm2)$lower), 
           upper = plogis(Effect(c("ANCO0.1.2.3newMPa"),ANCOm2)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#keep the grouping and
#see about grouping .8 with the other group. 
ANCO0.1.2.3.8newMPa<-ANCO$MPa
levels(ANCO0.1.2.3.8newMPa)
levels(ANCO0.1.2.3.8newMPa)[c(1,2,3,4,7)]<-"0.1.2.3.8"
levels(ANCO0.1.2.3.8newMPa)
ANCOm3<-glm(cbind(ANCO$Germinated,ANCO$Germinable-ANCO$Germinated)~ ANCO0.1.2.3.8newMPa,data=ANCO,family=binomial)
summary(ANCOm3)
anova(ANCOm2,ANCOm3, test="Chi")
library(effects)
data.frame(Effect(c("ANCO0.1.2.3.8newMPa"),ANCOm3)$x, 
           Est.=plogis(Effect(c("ANCO0.1.2.3.8newMPa"),ANCOm3)$fit), lower=plogis(Effect(c("ANCO0.1.2.3.8newMPa"),ANCOm3)$lower), 
           upper = plogis(Effect(c("ANCO0.1.2.3.8newMPa"),ANCOm3)$upper))

# See about grouping any further MPas.

#see about grouping 0, .1, .2, .3, .6, .8  

ANCO0.1.2.3.6.8newMPa<-ANCO$MPa
levels(ANCO0.1.2.3.6.8newMPa)
levels(ANCO0.1.2.3.6.8newMPa)[c(1,2,3,4,6,7)]<-"0.1.2.3.6.8"
levels(ANCO0.1.2.3.6.8newMPa)
ANCOm4<-glm(cbind(ANCO$Germinated,ANCO$Germinable-ANCO$Germinated)~ ANCO0.1.2.3.6.8newMPa,data=ANCO,family=binomial)
summary(ANCOm4)
anova(ANCOm3,ANCOm4, test="Chi")
library(effects)
data.frame(Effect(c("ANCO0.1.2.3.6.8newMPa"),ANCOm4)$x, 
           Est.=plogis(Effect(c("ANCO0.1.2.3.6.8newMPa"),ANCOm4)$fit), lower=plogis(Effect(c("ANCO0.1.2.3.6.8newMPa"),ANCOm4)$lower), 
           upper = plogis(Effect(c("ANCO0.1.2.3.6.8newMPa"),ANCOm4)$upper))

#looking into the Pr values and the way I have grouped the MPas.  Let's try one more model.
#ANCOm5 groups 0.1.2.3.8 and 1

ANCO0.1.2.3.6.81newMPa<-ANCO$MPa
levels(ANCO0.1.2.3.6.81newMPa)
levels(ANCO0.1.2.3.6.81newMPa)[c(1,2,3,4,6,7,8)]<-"0.1.2.3.6.81"
levels(ANCO0.1.2.3.6.81newMPa)
ANCOm5<-glm(cbind(ANCO$Germinated,ANCO$Germinable-ANCO$Germinated)~ ANCO0.1.2.3.6.81newMPa,data=ANCO,family=binomial)
summary(ANCOm5)
#compare it to m4
anova(ANCOm4,ANCOm5, test="Chi")
library(effects)
data.frame(Effect(c("ANCO0.1.2.3.6.81newMPa"),ANCOm5)$x, 
           Est.=plogis(Effect(c("ANCO0.1.2.3.6.81newMPa"),ANCOm5)$fit), lower=plogis(Effect(c("ANCO0.1.2.3.6.81newMPa"),ANCOm5)$lower), 
           upper = plogis(Effect(c("ANCO0.1.2.3.6.81newMPa"),ANCOm5)$upper))


#________Anthyllis vulneraria

# Subset data for the seed population of interest
ANVU<-subset(proportionsPEG,Species == "ANVU")
# Define the generalized linear model (full model)
ANVUm1<-glm(cbind(ANVU$Germinated,ANVU$Germinable - ANVU$Germinated) ~  MPa, data = ANVU, family=binomial) #we set up the germinated and not germinated as the response variable
library(effects)
summary(ANVUm1) 
data.frame(Effect(c("MPa"),ANVUm1)$x, 
           Est.=plogis(Effect(c("MPa"),ANVUm1)$fit), lower=plogis(Effect(c("MPa"),ANVUm1)$lower), 
           upper = plogis(Effect(c("MPa"),ANVUm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#Try to group some MPas.

ANVU0.1.2.3newMPa<-ANVU$MPa
levels(ANVU0.1.2.3newMPa)
levels(ANVU0.1.2.3newMPa)[c(1,2,3,4)]<-"0.1.2.3"

ANVUm2<-glm(cbind(ANVU$Germinated,ANVU$Germinable-ANVU$Germinated)~ ANVU0.1.2.3newMPa,data=ANVU,family=binomial)
summary(ANVUm2)
anova(ANVUm1,ANVUm2, test="Chi")
data.frame(Effect(c("ANVU0.1.2.3newMPa"),ANVUm2)$x, 
           Est.=plogis(Effect(c("ANVU0.1.2.3newMPa"),ANVUm2)$fit), lower=plogis(Effect(c("ANVU0.1.2.3newMPa"),ANVUm2)$lower), 
           upper = plogis(Effect(c("ANVU0.1.2.3newMPa"),ANVUm2)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#keep the grouping and
#see about grouping .6 with the other group. 
ANVU0.1.2.3.6newMPa<-ANVU$MPa
levels(ANVU0.1.2.3.6newMPa)
levels(ANVU0.1.2.3.6newMPa)[c(1,2,3,4,6)]<-"0.1.2.3.6"
levels(ANVU0.1.2.3.6newMPa)
ANVUm3<-glm(cbind(ANVU$Germinated,ANVU$Germinable-ANVU$Germinated)~ ANVU0.1.2.3.6newMPa,data=ANVU,family=binomial)
summary(ANVUm3)
anova(ANVUm2,ANVUm3, test="Chi")
library(effects)
data.frame(Effect(c("ANVU0.1.2.3.6newMPa"),ANVUm3)$x, 
           Est.=plogis(Effect(c("ANVU0.1.2.3.6newMPa"),ANVUm3)$fit), lower=plogis(Effect(c("ANVU0.1.2.3.6newMPa"),ANVUm3)$lower), 
           upper = plogis(Effect(c("ANVU0.1.2.3.6newMPa"),ANVUm3)$upper))

# See about grouping any further MPas.

#see about grouping 0, .1, .2, .3, .6 and grouping .4 and .8

ANVU0.1.2.3.6and.4.8newMPa<-ANVU$MPa
levels(ANVU0.1.2.3.6and.4.8newMPa)
levels(ANVU0.1.2.3.6and.4.8newMPa)[c(1,2,3,4,6)]<-"0.1.2.3.6"
levels(ANVU0.1.2.3.6and.4.8newMPa)
levels(ANVU0.1.2.3.6and.4.8newMPa)[c(2,3)]<-".4.8"
levels(ANVU0.1.2.3.6and.4.8newMPa)
ANVUm4<-glm(cbind(ANVU$Germinated,ANVU$Germinable-ANVU$Germinated)~ ANVU0.1.2.3.6and.4.8newMPa,data=ANVU,family=binomial)
summary(ANVUm4)
anova(ANVUm3,ANVUm4, test="Chi")
library(effects)
data.frame(Effect(c("ANVU0.1.2.3.6and.4.8newMPa"),ANVUm4)$x, 
           Est.=plogis(Effect(c("ANVU0.1.2.3.6and.4.8newMPa"),ANVUm4)$fit), lower=plogis(Effect(c("ANVU0.1.2.3.6and.4.8newMPa"),ANVUm4)$lower), 
           upper = plogis(Effect(c("ANVU0.1.2.3.6and.4.8newMPa"),ANVUm4)$upper))

#________Cleonia lusitanica


# Subset data for the seed population of interest
CLLU<-subset(proportionsPEG,Species == "CLLU")
# Define the generalized linear model (full model)
CLLUm1<-glm(cbind(CLLU$Germinated,CLLU$Germinable - CLLU$Germinated) ~  MPa, data = CLLU, family=binomial) #we set up the germinated and not germinated as the response variable
summary(CLLUm1) 
data.frame(Effect(c("MPa"),CLLUm1)$x, 
           Est.=plogis(Effect(c("MPa"),CLLUm1)$fit), lower=plogis(Effect(c("MPa"),CLLUm1)$lower), 
           upper = plogis(Effect(c("MPa"),CLLUm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#group 0.1.2.3.4
CLLU0.1.2.3.4newMPa<-CLLU$MPa
levels(CLLU0.1.2.3.4newMPa)
levels(CLLU0.1.2.3.4newMPa)[c(1,2,3,4,5)]<-"0.1.2.3.4"
levels(CLLU0.1.2.3.4newMPa)
CLLUm2<-glm(cbind(CLLU$Germinated,CLLU$Germinable-CLLU$Germinated)~ CLLU0.1.2.3.4newMPa,data=CLLU,family=binomial)
summary(CLLUm2)
anova(CLLUm1,CLLUm2, test="Chi")
library(effects)
data.frame(Effect(c("CLLU0.1.2.3.4newMPa"),CLLUm2)$x, 
           Est.=plogis(Effect(c("CLLU0.1.2.3.4newMPa"),CLLUm2)$fit), lower=plogis(Effect(c("CLLU0.1.2.3.4newMPa"),CLLUm2)$lower), 
           upper = plogis(Effect(c("CLLU0.1.2.3.4newMPa"),CLLUm2)$upper))

# m2 is better.  Finally, see if we can also group -0.6



CLLU0.1.2.3.4.6newMPa<-CLLU$MPa
levels(CLLU0.1.2.3.4.6newMPa)
levels(CLLU0.1.2.3.4.6newMPa)[c(1,2,3,4,5,6)]<-"0.1.2.3.4.6"
levels(CLLU0.1.2.3.4.6newMPa)
CLLUm3<-glm(cbind(CLLU$Germinated,CLLU$Germinable-CLLU$Germinated)~ CLLU0.1.2.3.4.6newMPa,data=CLLU,family=binomial)
summary(CLLUm3)
anova(CLLUm2,CLLUm3, test="Chi")
library(effects)
data.frame(Effect(c("CLLU0.1.2.3.4.6newMPa"),CLLUm3)$x, 
           Est.=plogis(Effect(c("CLLU0.1.2.3.4.6newMPa"),CLLUm3)$fit), lower=plogis(Effect(c("CLLU0.1.2.3.4.6newMPa"),CLLUm3)$lower), 
           upper = plogis(Effect(c("CLLU0.1.2.3.4.6newMPa"),CLLUm3)$upper))
#m3 isn't better than m2.  test m4 with 0.1.2.3.4 and .6.8 and 1




CLLU0.1.2.3.4and.6.8newMPa<-CLLU$MPa
levels(CLLU0.1.2.3.4and.6.8newMPa)
levels(CLLU0.1.2.3.4and.6.8newMPa)[c(1,2,3,4,5)]<-"0.1.2.3.4"
levels(CLLU0.1.2.3.4and.6.8newMPa)
levels(CLLU0.1.2.3.4and.6.8newMPa)[c(2,3)]<-".6.8"
levels(CLLU0.1.2.3.4and.6.8newMPa)
CLLUm4<-glm(cbind(CLLU$Germinated,CLLU$Germinable-CLLU$Germinated)~ CLLU0.1.2.3.4and.6.8newMPa,data=CLLU,family=binomial)
summary(CLLUm4)
anova(CLLUm2,CLLUm4, test="Chi")
library(effects)
data.frame(Effect(c("CLLU0.1.2.3.4and.6.8newMPa"),CLLUm4)$x, 
           Est.=plogis(Effect(c("CLLU0.1.2.3.4and.6.8newMPa"),CLLUm4)$fit), lower=plogis(Effect(c("CLLU0.1.2.3.4and.6.8newMPa"),CLLUm4)$lower), 
           upper = plogis(Effect(c("CLLU0.1.2.3.4and.6.8newMPa"),CLLUm4)$upper))


#________Echium plantagineum

# Subset data for the seed population of interest
ECPL<-subset(proportionsPEG,Species == "ECPL")
# Define the generalized linear model (full model)
ECPLm1<-glm(cbind(ECPL$Germinated,ECPL$Germinable - ECPL$Germinated) ~  MPa, data = ECPL, family=binomial) #we set up the germinated and not germinated as the response variable
summary(ECPLm1) 
data.frame(Effect(c("MPa"),ECPLm1)$x, 
           Est.=plogis(Effect(c("MPa"),ECPLm1)$fit), lower=plogis(Effect(c("MPa"),ECPLm1)$lower), 
           upper = plogis(Effect(c("MPa"),ECPLm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
# see if it is better to group 0.1.2.3.4.6

ECPL0.1.2.3.4.6newMPa<-ECPL$MPa
levels(ECPL0.1.2.3.4.6newMPa)
levels(ECPL0.1.2.3.4.6newMPa)[c(1,2,3,4,5,6)]<-"0.1.2.3.4.6"
levels(ECPL0.1.2.3.4.6newMPa)
ECPLm2<-glm(cbind(ECPL$Germinated,ECPL$Germinable-ECPL$Germinated)~ ECPL0.1.2.3.4.6newMPa,data=ECPL,family=binomial)
summary(ECPLm2)
anova(ECPLm1,ECPLm2, test="Chi")
library(effects)
data.frame(Effect(c("ECPL0.1.2.3.4.6newMPa"),ECPLm2)$x, 
           Est.=plogis(Effect(c("ECPL0.1.2.3.4.6newMPa"),ECPLm2)$fit), lower=plogis(Effect(c("ECPL0.1.2.3.4.6newMPa"),ECPLm2)$lower), 
           upper = plogis(Effect(c("ECPL0.1.2.3.4.6newMPa"),ECPLm2)$upper))
#keep m2



#________Moricandia moricandioides

# Subset data for the seed population of interest
MOMO<-subset(proportionsPEG,Species == "MOMO")
# Define the generalized linear model (full model)
MOMOm1<-glm(cbind(MOMO$Germinated,MOMO$Germinable - MOMO$Germinated) ~  MPa, data = MOMO, family=binomial) #we set up the germinated and not germinated as the response variable
summary(MOMOm1)
data.frame(Effect(c("MPa"),MOMOm1)$x, 
           Est.=plogis(Effect(c("MPa"),MOMOm1)$fit), lower=plogis(Effect(c("MPa"),MOMOm1)$lower), 
           upper = plogis(Effect(c("MPa"),MOMOm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
#Looks like to try and group .1.2.4.6.8

MOMO.1.2.3.4.6.8newMPa<-MOMO$MPa
levels(MOMO.1.2.3.4.6.8newMPa)
levels(MOMO.1.2.3.4.6.8newMPa)[c(2,3,4,5,6,7)]<-".1.2.3.4.6.8"
levels(MOMO.1.2.3.4.6.8newMPa)
MOMOm2<-glm(cbind(MOMO$Germinated,MOMO$Germinable-MOMO$Germinated)~ MOMO.1.2.3.4.6.8newMPa,data=MOMO,family=binomial)
summary(MOMOm2)
anova(MOMOm1,MOMOm2, test="Chi")
library(effects)
data.frame(Effect(c("MOMO.1.2.3.4.6.8newMPa"),MOMOm2)$x, 
           Est.=plogis(Effect(c("MOMO.1.2.3.4.6.8newMPa"),MOMOm2)$fit), lower=plogis(Effect(c("MOMO.1.2.3.4.6.8newMPa"),MOMOm2)$lower), 
           upper = plogis(Effect(c("MOMO.1.2.3.4.6.8newMPa"),MOMOm2)$upper))

#keep MOMOm2



#________Nigella damascena
# Subset data for the seed population of interest
NIDA<-subset(proportionsPEG,Species == "NIDA")
# Define the generalized linear model (full model)
NIDAm1<-glm(cbind(NIDA$Germinated,NIDA$Germinable - NIDA$Germinated) ~  MPa, data = NIDA, family=binomial) #we set up the germinated and not germinated as the response variable
summary(NIDAm1) 
data.frame(Effect(c("MPa"),NIDAm1)$x, 
           Est.=plogis(Effect(c("MPa"),NIDAm1)$fit), lower=plogis(Effect(c("MPa"),NIDAm1)$lower), 
           upper = plogis(Effect(c("MPa"),NIDAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#group all treatments against control

NIDA.1.2.3.4.6.81newMPa<-NIDA$MPa
levels(NIDA.1.2.3.4.6.81newMPa)
levels(NIDA.1.2.3.4.6.81newMPa)[c(2,3,4,5,6,7,8)]<-".1.2.3.4.6.81"
levels(NIDA.1.2.3.4.6.81newMPa)
NIDAm2<-glm(cbind(NIDA$Germinated,NIDA$Germinable-NIDA$Germinated)~ NIDA.1.2.3.4.6.81newMPa,data=NIDA,family=binomial)
summary(NIDAm2)
anova(NIDAm1,NIDAm2, test="Chi")
library(effects)
data.frame(Effect(c("NIDA.1.2.3.4.6.81newMPa"),NIDAm2)$x, 
           Est.=plogis(Effect(c("NIDA.1.2.3.4.6.81newMPa"),NIDAm2)$fit), lower=plogis(Effect(c("NIDA.1.2.3.4.6.81newMPa"),NIDAm2)$lower), 
           upper = plogis(Effect(c("NIDA.1.2.3.4.6.81newMPa"),NIDAm2)$upper))

#NIDAm2 is the better model.




#________Scabiosa atropurpurea
# Subset data for the seed population of interest
SCAT<-subset(proportionsPEG,Species == "SCAT")
# Define the generalized linear model (full model)
SCATm1<-glm(cbind(SCAT$Germinated,SCAT$Germinable - SCAT$Germinated) ~  MPa, data = SCAT, family=binomial) #we set up the germinated and not germinated as the response variable
summary(SCATm1) 
data.frame(Effect(c("MPa"),SCATm1)$x, 
           Est.=plogis(Effect(c("MPa"),SCATm1)$fit), lower=plogis(Effect(c("MPa"),SCATm1)$lower), 
           upper = plogis(Effect(c("MPa"),SCATm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#group .2.3

SCAT.2.3newMPa<-SCAT$MPa
levels(SCAT.2.3newMPa)
levels(SCAT.2.3newMPa)[c(3,4)]<-".2.3"
levels(SCAT.2.3newMPa)
SCATm2<-glm(cbind(SCAT$Germinated,SCAT$Germinable-SCAT$Germinated)~ SCAT.2.3newMPa,data=SCAT,family=binomial)
summary(SCATm2)
anova(SCATm1,SCATm2, test="Chi")
library(effects)
data.frame(Effect(c("SCAT.2.3newMPa"),SCATm2)$x, 
           Est.=plogis(Effect(c("SCAT.2.3newMPa"),SCATm2)$fit), lower=plogis(Effect(c("SCAT.2.3newMPa"),SCATm2)$lower), 
           upper = plogis(Effect(c("SCAT.2.3newMPa"),SCATm2)$upper))

#SCATm2 is the better model.

#SCATm3 see if I can also combine 0 and .1

SCAT0.1and.2.3newMPa<-SCAT$MPa
levels(SCAT0.1and.2.3newMPa)
levels(SCAT0.1and.2.3newMPa)[c(3,4)]<-".2.3"
levels(SCAT0.1and.2.3newMPa)
levels(SCAT0.1and.2.3newMPa)[c(1,2)]<-"0.1"
levels(SCAT0.1and.2.3newMPa)
SCATm3<-glm(cbind(SCAT$Germinated,SCAT$Germinable-SCAT$Germinated)~ SCAT0.1and.2.3newMPa,data=SCAT,family=binomial)
summary(SCATm3)
anova(SCATm2,SCATm3, test="Chi")
library(effects)
data.frame(Effect(c("SCAT0.1and.2.3newMPa"),SCATm3)$x, 
           Est.=plogis(Effect(c("SCAT0.1and.2.3newMPa"),SCATm3)$fit), lower=plogis(Effect(c("SCAT0.1and.2.3newMPa"),SCATm3)$lower), 
           upper = plogis(Effect(c("SCAT0.1and.2.3newMPa"),SCATm3)$upper))
#m4 see if 0.1.2.3 can be a group

SCAT0.1.2.3newMPa<-SCAT$MPa
levels(SCAT0.1.2.3newMPa)
levels(SCAT0.1.2.3newMPa)[c(1,2,3,4)]<-"0.1.2.3"
levels(SCAT0.1.2.3newMPa)
SCATm4<-glm(cbind(SCAT$Germinated,SCAT$Germinable-SCAT$Germinated)~ SCAT0.1.2.3newMPa,data=SCAT,family=binomial)
summary(SCATm4)
anova(SCATm1,SCATm4, test="Chi")
library(effects)
data.frame(Effect(c("SCAT0.1.2.3newMPa"),SCATm4)$x, 
           Est.=plogis(Effect(c("SCAT0.1.2.3newMPa"),SCATm4)$fit), lower=plogis(Effect(c("SCAT0.1.2.3newMPa"),SCATm4)$lower), 
           upper = plogis(Effect(c("SCAT0.1.2.3newMPa"),SCATm4)$upper))
#m5 0.1.2.3 and .4.6 and .8 and 1

SCAT0.1.2.3and.4.6newMPa<-SCAT$MPa
levels(SCAT0.1.2.3and.4.6newMPa)
levels(SCAT0.1.2.3and.4.6newMPa)[c(1,2,3,4)]<-"0.1.2.3"
levels(SCAT0.1.2.3and.4.6newMPa)
levels(SCAT0.1.2.3and.4.6newMPa)[c(2,3)]<-".4.6"
levels(SCAT0.1.2.3and.4.6newMPa)
SCATm5<-glm(cbind(SCAT$Germinated,SCAT$Germinable-SCAT$Germinated)~ SCAT0.1.2.3and.4.6newMPa,data=SCAT,family=binomial)
summary(SCATm5)
anova(SCATm1,SCATm5, test="Chi")
library(effects)
data.frame(Effect(c("SCAT0.1.2.3and.4.6newMPa"),SCATm5)$x, 
           Est.=plogis(Effect(c("SCAT0.1.2.3and.4.6newMPa"),SCATm5)$fit), lower=plogis(Effect(c("SCAT0.1.2.3and.4.6newMPa"),SCATm5)$lower), 
           upper = plogis(Effect(c("SCAT0.1.2.3and.4.6newMPa"),SCATm5)$upper))
#m4 is better than m5
#m6 0.1.2.3 and .4 and .6.8 and 1


SCAT0.1.2.3and.4and.6.8newMPa<-SCAT$MPa
levels(SCAT0.1.2.3and.4and.6.8newMPa)
levels(SCAT0.1.2.3and.4and.6.8newMPa)[c(1,2,3,4)]<-"0.1.2.3"
levels(SCAT0.1.2.3and.4and.6.8newMPa)
levels(SCAT0.1.2.3and.4and.6.8newMPa)[c(3,4)]<-".6.8"
levels(SCAT0.1.2.3and.4and.6.8newMPa)
SCATm6<-glm(cbind(SCAT$Germinated,SCAT$Germinable-SCAT$Germinated)~ SCAT0.1.2.3and.4and.6.8newMPa,data=SCAT,family=binomial)
summary(SCATm6)
anova(SCATm1,SCATm6, test="Chi")
library(effects)
data.frame(Effect(c("SCAT0.1.2.3and.4and.6.8newMPa"),SCATm6)$x, 
           Est.=plogis(Effect(c("SCAT0.1.2.3and.4and.6.8newMPa"),SCATm6)$fit), lower=plogis(Effect(c("SCAT0.1.2.3and.4and.6.8newMPa"),SCATm6)$lower), 
           upper = plogis(Effect(c("SCAT0.1.2.3and.4and.6.8newMPa"),SCATm6)$upper))
#Keep m6

#________Stachys arvensis

# Subset data for the seed population of interest
STAR<-subset(proportionsPEG,Species == "STAR")
# Define the generalized linear model (full model)
STARm1<-glm(cbind(STAR$Germinated,STAR$Germinable - STAR$Germinated) ~  MPa, data = STAR, family=binomial) #we set up the germinated and not germinated as the response variable
summary(STARm1) 
data.frame(Effect(c("MPa"),STARm1)$x, 
           Est.=plogis(Effect(c("MPa"),STARm1)$fit), lower=plogis(Effect(c("MPa"),STARm1)$lower), 
           upper = plogis(Effect(c("MPa"),STARm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#group 0.1

STAR0.1newMPa<-STAR$MPa
levels(STAR0.1newMPa)
levels(STAR0.1newMPa)[c(1,2)]<-"0.1"
levels(STAR0.1newMPa)
STARm2<-glm(cbind(STAR$Germinated,STAR$Germinable-STAR$Germinated)~ STAR0.1newMPa,data=STAR,family=binomial)
summary(STARm2)
anova(STARm1,STARm2, test="Chi")
library(effects)
data.frame(Effect(c("STAR0.1newMPa"),STARm2)$x, 
           Est.=plogis(Effect(c("STAR0.1newMPa"),STARm2)$fit), lower=plogis(Effect(c("STAR0.1newMPa"),STARm2)$lower), 
           upper = plogis(Effect(c("STAR0.1newMPa"),STARm2)$upper))

#STARm2 is the better model.
#STARm3 add the group of .4.8

STAR0.1and.4.8newMPa<-STAR$MPa
levels(STAR0.1and.4.8newMPa)
levels(STAR0.1and.4.8newMPa)[c(1,2)]<-"0.1"
levels(STAR0.1and.4.8newMPa)
levels(STAR0.1and.4.8newMPa)[c(4,6)]<-".4.8"
levels(STAR0.1and.4.8newMPa)

STARm3<-glm(cbind(STAR$Germinated,STAR$Germinable-STAR$Germinated)~ STAR0.1and.4.8newMPa,data=STAR,family=binomial)
summary(STARm3)
anova(STARm2,STARm3, test="Chi")
library(effects)
data.frame(Effect(c("STAR0.1and.4.8newMPa"),STARm3)$x, 
           Est.=plogis(Effect(c("STAR0.1and.4.8newMPa"),STARm3)$fit), lower=plogis(Effect(c("STAR0.1and.4.8newMPa"),STARm3)$lower), 
           upper = plogis(Effect(c("STAR0.1and.4.8newMPa"),STARm3)$upper))

#STARm4 group everything


STARnullMPa<-STAR$MPa
levels(STARnullMPa)
levels(STARnullMPa)[c(1,2,3,4,5,6,7)]<-"null"
levels(STARnullMPa)


STARm4<-glm(cbind(STAR$Germinated,STAR$Germinable-STAR$Germinated)~ STARnullMPa,data=STAR,family=binomial)
summary(STARm4)
anova(STARm3,STARm4, test="Chi")
library(effects)
data.frame(Effect(c("STARnullMPa"),STARm4)$x, 
           Est.=plogis(Effect(c("STARnullMPa"),STARm4)$fit), lower=plogis(Effect(c("STARnullMPa"),STARm4)$lower), 
           upper = plogis(Effect(c("STARnullMPa"),STARm4)$upper))

#nope STARm4 isn't possible

#STARm5 group 0.1 and .2.4.6.81

STAR0.1and.2.4.6.81MPa<-STAR$MPa
levels(STAR0.1and.2.4.6.81MPa)
levels(STAR0.1and.2.4.6.81MPa)[c(1,2)]<-"0.1"
levels(STAR0.1and.2.4.6.81MPa)
levels(STAR0.1and.2.4.6.81MPa)[c(2,4,5,6,7)]<-".2.4.6.81"
levels(STAR0.1and.2.4.6.81MPa)


STARm5<-glm(cbind(STAR$Germinated,STAR$Germinable-STAR$Germinated)~ STAR0.1and.2.4.6.81MPa,data=STAR,family=binomial)
summary(STARm5)
anova(STARm3,STARm5, test="Chi")
library(effects)
data.frame(Effect(c("STAR0.1and.2.4.6.81MPa"),STARm5)$x, 
           Est.=plogis(Effect(c("STAR0.1and.2.4.6.81MPa"),STARm5)$fit), lower=plogis(Effect(c("STAR0.1and.2.4.6.81MPa"),STARm5)$lower), 
           upper = plogis(Effect(c("STAR0.1and.2.4.6.81MPa"),STARm5)$upper))

#Keep STARm5 as the best model


#________Tolpis barbata
# Subset data for the seed population of interest
TOBA<-subset(proportionsPEG,Species == "TOBA")
# Define the generalized linear model (full model)
TOBAm1<-glm(cbind(TOBA$Germinated,TOBA$Germinable - TOBA$Germinated) ~  MPa, data = TOBA, family=binomial) #we set up the germinated and not germinated as the response variable
summary(TOBAm1) 
data.frame(Effect(c("MPa"),TOBAm1)$x, 
           Est.=plogis(Effect(c("MPa"),TOBAm1)$fit), lower=plogis(Effect(c("MPa"),TOBAm1)$lower), 
           upper = plogis(Effect(c("MPa"),TOBAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
#TOBAm2 group 0.1.2.3

TOBA0.1.2.3newMPa<-TOBA$MPa
levels(TOBA0.1.2.3newMPa)
levels(TOBA0.1.2.3newMPa)[c(1,2,3,4)]<-"0.1.2.3"
levels(TOBA0.1.2.3newMPa)
TOBAm2<-glm(cbind(TOBA$Germinated,TOBA$Germinable-TOBA$Germinated)~ TOBA0.1.2.3newMPa,data=TOBA,family=binomial)
summary(TOBAm2)
anova(TOBAm1,TOBAm2, test="Chi")
library(effects)
data.frame(Effect(c("TOBA0.1.2.3newMPa"),TOBAm2)$x, 
           Est.=plogis(Effect(c("TOBA0.1.2.3newMPa"),TOBAm2)$fit), lower=plogis(Effect(c("TOBA0.1.2.3newMPa"),TOBAm2)$lower), 
           upper = plogis(Effect(c("TOBA0.1.2.3newMPa"),TOBAm2)$upper))
#keep TOBAm2

#________Tordylium maximum

# Subset data for the seed population of interest
TOMA<-subset(proportionsPEG,Species == "TOMA")
# Define the generalized linear model (full model)
TOMAm1<-glm(cbind(TOMA$Germinated,TOMA$Germinable - TOMA$Germinated) ~  MPa, data = TOMA, family=binomial) #we set up the germinated and not germinated as the response variable
summary(TOMAm1) 
data.frame(Effect(c("MPa"),TOMAm1)$x, 
           Est.=plogis(Effect(c("MPa"),TOMAm1)$fit), lower=plogis(Effect(c("MPa"),TOMAm1)$lower), 
           upper = plogis(Effect(c("MPa"),TOMAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
#TOMAm2 group 0.1

TOMA0.1newMPa<-TOMA$MPa
levels(TOMA0.1newMPa)
levels(TOMA0.1newMPa)[c(1,2)]<-"0.1"
levels(TOMA0.1newMPa)
TOMAm2<-glm(cbind(TOMA$Germinated,TOMA$Germinable-TOMA$Germinated)~ TOMA0.1newMPa,data=TOMA,family=binomial)
summary(TOMAm2)
anova(TOMAm1,TOMAm2, test="Chi")
library(effects)
data.frame(Effect(c("TOMA0.1newMPa"),TOMAm2)$x, 
           Est.=plogis(Effect(c("TOMA0.1newMPa"),TOMAm2)$fit), lower=plogis(Effect(c("TOMA0.1newMPa"),TOMAm2)$lower), 
           upper = plogis(Effect(c("TOMA0.1newMPa"),TOMAm2)$upper))
#keep TOMAm2

#with TOMAm3 group 0.1 and .2.3
TOMA0.1and2.3newMPa<-TOMA$MPa
levels(TOMA0.1and2.3newMPa)
levels(TOMA0.1and2.3newMPa)[c(1,2)]<-"0.1"
levels(TOMA0.1and2.3newMPa)
levels(TOMA0.1and2.3newMPa)[c(2,3)]<-"2.3"
levels(TOMA0.1and2.3newMPa)
TOMAm3<-glm(cbind(TOMA$Germinated,TOMA$Germinable-TOMA$Germinated)~ TOMA0.1and2.3newMPa,data=TOMA,family=binomial)
summary(TOMAm3)
anova(TOMAm2,TOMAm3, test="Chi")
library(effects)
data.frame(Effect(c("TOMA0.1and2.3newMPa"),TOMAm3)$x, 
           Est.=plogis(Effect(c("TOMA0.1and2.3newMPa"),TOMAm3)$fit), lower=plogis(Effect(c("TOMA0.1and2.3newMPa"),TOMAm3)$lower), 
           upper = plogis(Effect(c("TOMA0.1and2.3newMPa"),TOMAm3)$upper))
#TOMAm3 is better.

#TOMAm4 groups 0.1 and .2.3.4

TOMA0.1and2.3.4newMPa<-TOMA$MPa
levels(TOMA0.1and2.3.4newMPa)
levels(TOMA0.1and2.3.4newMPa)[c(1,2)]<-"0.1"
levels(TOMA0.1and2.3.4newMPa)
levels(TOMA0.1and2.3.4newMPa)[c(2,3,4)]<-"2.3.4"
levels(TOMA0.1and2.3.4newMPa)
TOMAm4<-glm(cbind(TOMA$Germinated,TOMA$Germinable-TOMA$Germinated)~ TOMA0.1and2.3.4newMPa,data=TOMA,family=binomial)
summary(TOMAm4)
anova(TOMAm2,TOMAm4, test="Chi")
library(effects)
data.frame(Effect(c("TOMA0.1and2.3.4newMPa"),TOMAm4)$x, 
           Est.=plogis(Effect(c("TOMA0.1and2.3.4newMPa"),TOMAm4)$fit), lower=plogis(Effect(c("TOMA0.1and2.3.4newMPa"),TOMAm4)$lower), 
           upper = plogis(Effect(c("TOMA0.1and2.3.4newMPa"),TOMAm4)$upper))

#keep TOMAm4



