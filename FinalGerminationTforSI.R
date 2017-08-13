setwd("~/Desktop/20170424 desktop stuff/NASSTEC SESIL 2A/Obj 2 SF NASSTEC Germination Niche/Obj 2 data/germination analysis files")
# Load the libraries
library(car)
library(effects)
library(MASS)
library(plyr)
library(binom)

# Read the table
proportions <-read.table("proportionsT.txt", header=T)

# Go through the data, seed population by seed population. Species names are abbreviated with the first two letters each of the genus and specific epithet e.g. Anthemis cotula = ANCO.

#________Anthemis cotula

# Subset data for the seed population of interest
ANCO<-subset(proportions,Species == "ANCO")
# Define the generalized linear model (full model)
ANCOm1<-glm(cbind(ANCO$Germinated,ANCO$Germinable - ANCO$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = ANCO, family=binomial)  
# Model simplification, updating the previous model by stepwise removal of interactions and/or explanatory variables
ANCOm2 <-update(ANCOm1,~.-Temperature : Oscillation) 
anova(ANCOm1,ANCOm2, test = "Chi")
# Pr > 0.05 so reject the fuller model and continue with model simplification.
ANCOm3 <-update(ANCOm2,~.-Temperature) 
anova(ANCOm2,ANCOm3, test = "Chi")
# Pr > 0.05 so reject the fuller model and continue with model simplification.
ANCOm4 <-update(ANCOm3,~.-Oscillation) 
anova(ANCOm3,ANCOm4, test = "Chi")
# Pr > 0.05 so reject the fuller model and keep m4 (in this case, the null model.)

# View parameters of the best, simplest model and interpret the output: the intercept is 10C alternate and since it is categorical data, the Estimate is saying the difference (+ or - depending on the sign) from 10C alternating.
summary(ANCOm4) 

# Generate the estimates for a model with no effect of treatment (all the same for all the reps and treatments, because there was no significant effect of Temperature or Oscillation.)

binom.confint(sum(ANCO$Germinated), sum(ANCO$Germinable), methods="logit") 


#________Anthyllis vulneraria
# Subset data for the seed population of interest
ANVU<-subset(proportions,Species == "ANVU")
# Define the generalized linear model (full model)
ANVUm1<-glm(cbind(ANVU$Germinated,ANVU$Germinable - ANVU$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = ANVU, family=binomial) #we set up the germinated and not germinated as the response variable
summary(ANVUm1) 
data.frame(Effect(c("Temperature","Oscillation"),ANVUm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),ANVUm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),ANVUm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),ANVUm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
# Model simplification, updating the previous model by stepwise removal of interactions and/or explanatory variables
ANVUm2 <-update(ANVUm1,~.-Temperature : Oscillation) #Model simplification
anova(ANVUm1,ANVUm2, test = "Chi")
# Pr < 0.05 so ANVUm1 is better. Keep it.
# Based on the similar values for the estimates, next set up a model that groups 10C and 15C and see if it improves the model.
ANVU10.15Temperature<-ANVU$Temperature
levels(ANVU10.15Temperature)
levels(ANVU10.15Temperature)[c(1,2)]<-"10and15C"

ANVUm5<-glm(cbind(ANVU$Germinated,ANVU$Germinable-ANVU$Germinated)~ANVU10.15Temperature+Oscillation
            +ANVU10.15Temperature:Oscillation,data=ANVU,family=binomial)
summary(ANVUm5)
anova(ANVUm1,ANVUm5, test="Chi")
# Pr > 0.05 it so keep m5 (simpler model) and try additional groupings. 

#Next try m6, which has 5,10and15C merged and compare it to m5.

ANVU5.10.15Temperature<-ANVU10.15Temperature
levels(ANVU5.10.15Temperature)
levels(ANVU5.10.15Temperature)[c(1,4)]<-"5and10and15C"
levels(ANVU5.10.15Temperature)

ANVUm6<-glm(cbind(ANVU$Germinated,ANVU$Germinable-ANVU$Germinated)~ANVU5.10.15Temperature+Oscillation
            +ANVU5.10.15Temperature:Oscillation,data=ANVU,family=binomial)
summary(ANVUm6)
anova(ANVUm5,ANVUm6, test="Chi")

# Pr > 0.05 so m6 was better than m5. Next test if grouping 20C and 25C improve the model.

ANVU20.25Temperature<-ANVU$Temperature
levels(ANVU20.25Temperature)
levels(ANVU20.25Temperature)[c(1,2,5)]<-"5and10and15C"
levels(ANVU20.25Temperature)[c(2,3)]<-"20and25C"
ANVUm7<-glm(cbind(ANVU$Germinated,ANVU$Germinable-ANVU$Germinated)~ANVU20.25Temperature+Oscillation
            +ANVU20.25Temperature:Oscillation,data=ANVU,family=binomial)
summary(ANVUm7)
anova(ANVUm6,ANVUm7, test="Chi")
# Pr < 0.05 so m7 (grouping 20C and 25C) was not an improvement. m6 is the best model. 

# Generate the final germination estimates under the best model.
data.frame(Effect(c("ANVU5.10.15Temperature","Oscillation"),ANVUm6)$x, 
           Est.=plogis(Effect(c("ANVU5.10.15Temperature","Oscillation"), ANVUm6)$fit), 
           lower=plogis(Effect(c("ANVU5.10.15Temperature","Oscillation"), ANVUm6)$lower), 
           upper = plogis(Effect(c("ANVU5.10.15Temperature","Oscillation"), ANVUm6)$upper))

#________Cleonia lusitanica

# Subset data for the seed population of interest
CLLU<-subset(proportions,Species == "CLLU")
# Define the generalized linear model (full model)
CLLUm1<-glm(cbind(CLLU$Germinated,CLLU$Germinable - CLLU$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = CLLU, family=binomial) #we set up the germinated and not germinated as the response variable
summary(CLLUm1) 
data.frame(Effect(c("Temperature","Oscillation"),CLLUm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),CLLUm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),CLLUm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),CLLUm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
# Model simplification, updating the previous model by stepwise removal of interactions and/or explanatory variables
CLLUm2 <-update(CLLUm1,~.-Temperature : Oscillation) #Model simplification
anova(CLLUm1,CLLUm2, test = "Chi")
# Pr < 0.05. Cannot further simplify model. Keep the full model.
# Try grouping 10C, 15C, and 20C together.
CLLU10.15.20Temperature<-CLLU$Temperature
levels(CLLU10.15.20Temperature)
levels(CLLU10.15.20Temperature)[c(1,2,3)]<-"10,15and20"

CLLUm5<-glm(cbind(CLLU$Germinated,CLLU$Germinable-CLLU$Germinated)~CLLU10.15.20Temperature+Oscillation
            +CLLU10.15.20Temperature:Oscillation,data=CLLU,family=binomial)
summary(CLLUm5)
anova(CLLUm1,CLLUm5, test="Chi")
# Pr < 0.05 so grouping does not improve the model. Keep the full model, CLLUm1.

#________Echium plantagineum

# Subset data for the seed population of interest
ECPL<-subset(proportions,Species == "ECPL")
# Define the generalized linear model (full model)
ECPLm1<-glm(cbind(ECPL$Germinated,ECPL$Germinable - ECPL$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = ECPL, family=binomial) #we set up the germinated and not germinated as the response variable
summary(ECPLm1) 
data.frame(Effect(c("Temperature","Oscillation"),ECPLm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),ECPLm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),ECPLm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),ECPLm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
# Model simplification, updating the previous model by stepwise removal of interactions and/or explanatory variables
ECPLm2 <-update(ECPLm1,~.-Temperature : Oscillation) #Model simplification
anova(ECPLm1,ECPLm2, test = "Chi")
# Pr < 0.05 so reject the simplification and keep the full model, m1.
# Try to combine the temperatures 20C and 25C and see if that improves the model.

ECPL20.25CTemperature<-ECPL$Temperature
levels(ECPL20.25CTemperature)
levels(ECPL20.25CTemperature)[c(3,4)]<-"20and25C"
ECPLm5<-glm(cbind(ECPL$Germinated,ECPL$Germinable-ECPL$Germinated)~ ECPL20.25CTemperature +Oscillation
            + ECPL20.25CTemperature:Oscillation,data=ECPL,family=binomial)
summary(ECPLm5)
anova(ECPLm1,ECPLm5, test="Chi")
#Pr > 0.05 so keep this grouping. 
# Next see if 15C can also be combined. 

ECPL15.20.25CTemperature<-ECPL20.25CTemperature
levels(ECPL15.20.25CTemperature)
levels(ECPL15.20.25CTemperature)[c(2,3)]<-"15and20and25C"
levels(ECPL15.20.25CTemperature)
ECPLm6<-glm(cbind(ECPL$Germinated,ECPL$Germinable-ECPL$Germinated)~ ECPL15.20.25CTemperature +Oscillation
            + ECPL15.20.25CTemperature:Oscillation,data=ECPL,family=binomial)
summary(ECPLm6)
anova(ECPLm5,ECPLm6, test="Chi")

#Pr < 0.05 so reject the grouping. m5 is the best. 



data.frame(Effect(c("ECPL20.25CTemperature","Oscillation"),ECPLm5)$x, 
           Est.=plogis(Effect(c("ECPL20.25CTemperature","Oscillation"), ECPLm5)$fit), 
           lower=plogis(Effect(c("ECPL20.25CTemperature","Oscillation"), ECPLm5)$lower), 
           upper = plogis(Effect(c("ECPL20.25CTemperature","Oscillation"), ECPLm5)$upper))

#________Moricandia moricandioides

# Subset data for the seed population of interest
MOMO<-subset(proportions,Species == "MOMO")
# Define the generalized linear model (full model)
MOMOm1<-glm(cbind(MOMO$Germinated,MOMO$Germinable - MOMO$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = MOMO, family=binomial) #we set up the germinated and not germinated as the response variable
summary(MOMOm1) 
data.frame(Effect(c("Temperature","Oscillation"),MOMOm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),MOMOm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),MOMOm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),MOMOm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
# Model simplification, updating the previous model by stepwise removal of interactions and/or explanatory variables
MOMOm2 <-update(MOMOm1,~.-Temperature : Oscillation) #Model simplification
anova(MOMOm1,MOMOm2, test = "Chi")
# Pr < 0.05 so reject simplification. 
# Test grouping temperatures. 
MOMO15.20Temperature<-MOMO$Temperature
levels(MOMO15.20Temperature)
levels(MOMO15.20Temperature)[c(2,3)]<-"15and20"

MOMOm5<-glm(cbind(MOMO$Germinated,MOMO$Germinable-MOMO$Germinated)~MOMO15.20Temperature+Oscillation
            +MOMO15.20Temperature:Oscillation,data=MOMO,family=binomial)
summary(MOMOm5)
anova(MOMOm1,MOMOm5, test="Chi")
# Pr < 0.05 so reject the grouping. m1 is the best model. 


#________Nigella damascena

# Subset data for the seed population of interest
NIDA<-subset(proportions,Species == "NIDA")
# Define the generalized linear model (full model)
NIDAm1<-glm(cbind(NIDA$Germinated,NIDA$Germinable - NIDA$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = NIDA, family=binomial) #we set up the germinated and not germinated as the response variable
summary(NIDAm1) 
data.frame(Effect(c("Temperature","Oscillation"),NIDAm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),NIDAm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),NIDAm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),NIDAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

NIDAm2 <-update(NIDAm1,~.-Temperature : Oscillation) #Model simplification
anova(NIDAm1,NIDAm2, test = "Chi")
# Pr < 0.05 reject the simplification and keep the full model. 

#________Scabiosa atropurpurea

# Subset data for the seed population of interest
SCAT<-subset(proportions,Species == "SCAT")
# Define the generalized linear model (full model)
SCATm1<-glm(cbind(SCAT$Germinated,SCAT$Germinable - SCAT$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = SCAT, family=binomial) #we set up the germinated and not germinated as the response variable
summary(SCATm1) 
data.frame(Effect(c("Temperature","Oscillation"),SCATm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),SCATm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),SCATm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),SCATm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.
# Model simplification, updating the previous model by stepwise removal of interactions and/or explanatory variables
SCATm2 <-update(SCATm1,~.-Temperature : Oscillation) #Model simplification
anova(SCATm1,SCATm2, test = "Chi")
# m1 is better. Next try grouping temperatures.

SCAT20.25Temperature<-SCAT$Temperature
levels(SCAT20.25Temperature)
levels(SCAT20.25Temperature)[c(3,4)]<-"20and25C"
SCATm7<-glm(cbind(SCAT$Germinated,SCAT$Germinable-SCAT$Germinated)~SCAT20.25Temperature+Oscillation
            +SCAT20.25Temperature:Oscillation,data=SCAT,family=binomial)
summary(SCATm7)
anova(SCATm1,SCATm7, test="Chi")

#grouping 20and25 is a better model.  Now compare with grouping 10,20and 25.
levels(SCAT20.25Temperature)
levels(SCAT20.25Temperature)[c(1,3)]<-"10and20and25C"
levels(SCAT20.25Temperature)
SCATm8<-glm(cbind(SCAT$Germinated,SCAT$Germinable-SCAT$Germinated)~SCAT20.25Temperature+Oscillation
            +SCAT20.25Temperature:Oscillation,data=SCAT,family=binomial)
summary(SCATm8)
anova(SCATm7,SCATm8, test="Chi")
# Pr < 0.05 reject further grouping. m7 is the best model.
# Estimate the final germination proportion under the best model.
data.frame(Effect(c("SCAT20.25Temperature","Oscillation"),SCATm7)$x, 
           Est.=plogis(Effect(c("SCAT20.25Temperature","Oscillation"), SCATm7)$fit), 
           lower=plogis(Effect(c("SCAT20.25Temperature","Oscillation"), SCATm7)$lower), 
           upper = plogis(Effect(c("SCAT20.25Temperature","Oscillation"), SCATm7)$upper))


#________Stachys arvensis

# Subset data for the seed population of interest
STAR<-subset(proportions,Species == "STAR")
# Define the generalized linear model (full model)
STARm1<-glm(cbind(STAR$Germinated,STAR$Germinable - STAR$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = STAR, family=binomial) #we set up the germinated and not germinated as the response variable
summary(STARm1) 
data.frame(Effect(c("Temperature","Oscillation"),STARm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),STARm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),STARm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),STARm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

STARm2 <-update(STARm1,~.-Temperature : Oscillation) #Model simplification
anova(STARm1,STARm2, test = "Chi")
# Pr < 0.05 reject the simplification and keep the full model. 

#________Tolpis barbata

# Subset data for the seed population of interest
TOBA<-subset(proportions,Species == "TOBA")
# Define the generalized linear model (full model)
TOBAm1<-glm(cbind(TOBA$Germinated,TOBA$Germinable - TOBA$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = TOBA, family=binomial) #we set up the germinated and not germinated as the response variable
summary(TOBAm1) 
data.frame(Effect(c("Temperature","Oscillation"),TOBAm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),TOBAm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),TOBAm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),TOBAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.


# Model simplification, updating the previous model by stepwise removal of interactions and/or explanatory variables
TOBAm2 <-update(TOBAm1,~.-Temperature : Oscillation) #Model simplification
anova(TOBAm1,TOBAm2, test = "Chi")
#Keep m1.
#now see about grouping temps 15 and 20
TOBA15.20Temperature<-TOBA$Temperature
levels(TOBA15.20Temperature)
levels(TOBA15.20Temperature)[c(2,3)]<-"15and20C"
levels(TOBA15.20Temperature)
TOBAm3<-glm(cbind(TOBA$Germinated,TOBA$Germinable-TOBA$Germinated)~TOBA15.20Temperature+Oscillation
            +TOBA15.20Temperature:Oscillation,data=TOBA,family=binomial)
summary(TOBAm3)
anova(TOBAm1,TOBAm3, test="Chi")
#Pr > 0.05 so Keep m3.  Next check if 10 and 25 can be combined.

levels(TOBA15.20Temperature)
levels(TOBA15.20Temperature)[c(1,3)]<-"10and25C"
levels(TOBA15.20Temperature)
TOBAm4<-glm(cbind(TOBA$Germinated,TOBA$Germinable-TOBA$Germinated)~TOBA15.20Temperature+Oscillation
            +TOBA15.20Temperature:Oscillation,data=TOBA,family=binomial)
summary(TOBAm4)
anova(TOBAm3,TOBAm4, test="Chi")
#TOBAm3 is better.

data.frame(Effect(c("TOBA15.20Temperature","Oscillation"),TOBAm3)$x, Est.=plogis(Effect(c("TOBA15.20Temperature","Oscillation"),TOBAm3)$fit), lower=plogis(Effect(c("TOBA15.20Temperature","Oscillation"),TOBAm3)$lower), upper = plogis(Effect(c("TOBA15.20Temperature","Oscillation"),TOBAm3)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

#________Tordylium maximum

# Subset data for the seed population of interest
TOMA<-subset(proportions,Species == "TOMA")
# Define the generalized linear model (full model)
TOMAm1<-glm(cbind(TOMA$Germinated,TOMA$Germinable - TOMA$Germinated) ~  Oscillation + Temperature + Oscillation : Temperature, data = TOMA, family=binomial) #we set up the germinated and not germinated as the response variable
summary(TOMAm1) 
data.frame(Effect(c("Temperature","Oscillation"),TOMAm1)$x, Est.=plogis(Effect(c("Temperature","Oscillation"),TOMAm1)$fit), lower=plogis(Effect(c("Temperature","Oscillation"),TOMAm1)$lower), upper = plogis(Effect(c("Temperature","Oscillation"),TOMAm1)$upper)) #the glm works in "odds" as units and the plogis function converts this back to proportion.

TOMAm2 <-update(TOMAm1,~.-Temperature : Oscillation) #Model simplification
anova(TOMAm1,TOMAm2, test = "Chi")
# Pr < 0.05 reject the simplification and keep the full model. 




