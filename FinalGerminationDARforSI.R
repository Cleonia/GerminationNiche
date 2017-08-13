setwd("~/Desktop/NASSTEC SESIL 2A/Obj 2 SF NASSTEC Germination Niche/Obj 2 data/germination analysis files")
# Install and load the necessary packages
install.packages("effects")
install.packages("binom")
library(effects)
library(binom)
# Read the dataset
proportionsDAR <- read.table("proportionsDAR.txt", header = T)
proportionsDAR
# Go through the data, seed population by seed population. Species names are abbreviated with the first two letters each of the genus and specific epithet. Anthemis cotula = ANCO.

#________Anthemis cotula
# Subset the data for the seed population of interest
ANCO<-subset(proportionsDAR,Species == "ANCO")
# Define the generalized linear model (full model)
ANCOm1<-glm(cbind(ANCO$Germinated,ANCO$Germinable - ANCO$Germinated) ~  DAR, data = ANCO, family=binomial) 
summary(ANCOm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),ANCOm1)$x, 
           Est.=plogis(Effect(c("DAR"),ANCOm1)$fit), lower=plogis(Effect(c("DAR"),ANCOm1)$lower), 
           upper = plogis(Effect(c("DAR"),ANCOm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.


#________Anthyllis vulneraria
# Subset the data for the seed population of interest
ANVU<-subset(proportionsDAR,Species == "ANVU")
# Define the generalized linear model (full model)
ANVUm1<-glm(cbind(ANVU$Germinated,ANVU$Germinable - ANVU$Germinated) ~  DAR, data = ANVU, family=binomial) 
summary(ANVUm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),ANVUm1)$x, 
           Est.=plogis(Effect(c("DAR"),ANVUm1)$fit), lower=plogis(Effect(c("DAR"),ANVUm1)$lower), 
           upper = plogis(Effect(c("DAR"),ANVUm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.


binom.confint(sum(ANVU$Germinated), sum(ANVU$Germinable), methods="logit") #generate the summary of estimates all the same for all the reps and treatments, because there was no significant effect of T or Oscillation.


#________Cleonia lusitanica

# Subset data for the seed population of interest
CLLU<-subset(proportionsDAR,Species == "CLLU")
# Define the generalized linear model (full model)
CLLUm1<-glm(cbind(CLLU$Germinated,CLLU$Germinable - CLLU$Germinated) ~  DAR, data = CLLU, family=binomial) 
summary(CLLUm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),CLLUm1)$x, 
           Est.=plogis(Effect(c("DAR"),CLLUm1)$fit), lower=plogis(Effect(c("DAR"),CLLUm1)$lower), 
           upper = plogis(Effect(c("DAR"),CLLUm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.


#________Echium plantagineum

# Subset data for the seed population of interest
ECPL<-subset(proportionsDAR,Species == "ECPL")
# Define the generalized linear model (full model)
ECPLm1<-glm(cbind(ECPL$Germinated,ECPL$Germinable - ECPL$Germinated) ~  DAR, data = ECPL, family=binomial) 
summary(ECPLm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),ECPLm1)$x, 
           Est.=plogis(Effect(c("DAR"),ECPLm1)$fit), lower=plogis(Effect(c("DAR"),ECPLm1)$lower), 
           upper = plogis(Effect(c("DAR"),ECPLm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.


#________Moricandia moricandioides

# Subset data for the seed population of interest

MOMO<-subset(proportionsDAR,Species == "MOMOINC")
MOMOm1<-glm(cbind(MOMO$Germinated,MOMO$Germinable - MOMO$Germinated) ~  DAR, data = MOMO, family=binomial) 
summary(MOMOm1) 


data.frame(Effect(c("DAR"),MOMOm1)$x, 
           Est.=plogis(Effect(c("DAR"),MOMOm1)$fit), lower=plogis(Effect(c("DAR"),MOMOm1)$lower), 
           upper = plogis(Effect(c("DAR"),MOMOm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#________Nigella damascena

# Subset data for the seed population of interest
NIDA<-subset(proportionsDAR,Species == "NIDA")
# Define the generalized linear model (full model)
NIDAm1<-glm(cbind(NIDA$Germinated,NIDA$Germinable - NIDA$Germinated) ~  DAR, data = NIDA, family=binomial) 
summary(NIDAm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),NIDAm1)$x, 
           Est.=plogis(Effect(c("DAR"),NIDAm1)$fit), lower=plogis(Effect(c("DAR"),NIDAm1)$lower), 
           upper = plogis(Effect(c("DAR"),NIDAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#________Scabiosa atropurpurea

# Subset data for the seed population of interest
SCAT<-subset(proportionsDAR,Species == "SCAT")
# Define the generalized linear model (full model)
SCATm1<-glm(cbind(SCAT$Germinated,SCAT$Germinable - SCAT$Germinated) ~  DAR, data = SCAT, family=binomial) 
summary(SCATm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),SCATm1)$x, 
           Est.=plogis(Effect(c("DAR"),SCATm1)$fit), lower=plogis(Effect(c("DAR"),SCATm1)$lower), 
           upper = plogis(Effect(c("DAR"),SCATm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.


#________Stachys arvensis

# Subset data for the seed population of interest
STAR<-subset(proportionsDAR,Species == "STAR")
# Define the generalized linear model (full model)
STARm1<-glm(cbind(STAR$Germinated,STAR$Germinable - STAR$Germinated) ~  DAR, data = STAR, family=binomial) 
summary(STARm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),STARm1)$x, 
           Est.=plogis(Effect(c("DAR"),STARm1)$fit), lower=plogis(Effect(c("DAR"),STARm1)$lower), 
           upper = plogis(Effect(c("DAR"),STARm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#________Tolpis barbata

# Subset data for the seed population of interest
TOBA<-subset(proportionsDAR,Species == "TOBA")
# Define the generalized linear model (full model)
TOBAm1<-glm(cbind(TOBA$Germinated,TOBA$Germinable - TOBA$Germinated) ~  DAR, data = TOBA, family=binomial) 
summary(TOBAm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),TOBAm1)$x, 
           Est.=plogis(Effect(c("DAR"),TOBAm1)$fit), lower=plogis(Effect(c("DAR"),TOBAm1)$lower), 
           upper = plogis(Effect(c("DAR"),TOBAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.



binom.confint(sum(TOBA$Germinated), sum(TOBA$Germinable), methods="logit") #generate the summary of estimates all the same for all the reps and treatments, because there was no significant effect of T or Oscillation.



#________Tordylium maximum

# Subset data for the seed population of interest
TOMA<-subset(proportionsDAR,Species == "TOMA")
# Define the generalized linear model (full model)
TOMAm1<-glm(cbind(TOMA$Germinated,TOMA$Germinable - TOMA$Germinated) ~  DAR, data = TOMA, family=binomial) 
summary(TOMAm1) 

# Generate the estimates for final germination proportion based on glm. 
data.frame(Effect(c("DAR"),TOMAm1)$x, 
           Est.=plogis(Effect(c("DAR"),TOMAm1)$fit), lower=plogis(Effect(c("DAR"),TOMAm1)$lower), 
           upper = plogis(Effect(c("DAR"),TOMAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.







