#----------------Longitudinal linear mixed-effects models (LMM)---------------#

# The first part of the script computes LMM analyzing the influence of disease
# duration on visual acuity and VEP latency over time. To account for 
# confounding influences, age and/or sex were subsequently added to the model 
# if the model fit could be improved. The results are presented in the result 
# section of the manuscript.

# The second part of the script computes LMM analyzing the annualized rate of  
# thickness change (globalpRNFL, macular RNFL, GCIPL, INL) for the different
# periods of the disease duration. To account for confounding influences, age 
# and/or sex were subsequently added to the model if the model fit could be improved. 
# The results are presented in Supplementary Table 4. 
# Usual runtime ~= 5 minutes 

#---------------------------------R packages---------------------------------#  

library(lme4)
library(lmerTest) #needed for p-values
library(haven)
library(performance)
library(car)


#-----------------------Data loading & data management-----------------------#

OCT_data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)

OCT_data$Latenz[OCT_data$Latenz == 777] <- NA

#Correct diesease durations, if VEP dates do not match OCT dates
library(lubridate)
OCT_data$Disease_Duration_Latency <- (as_date(OCT_data$VEP) - as_date(OCT_data$OCTDatum))/365 + OCT_data$Disease_Duration_FINAL
OCT_data$Age_Latency <- (as_date(OCT_data$VEP) - as_date(OCT_data$OCTDatum))/365 + OCT_data$Age


Data_longitudinal <- subset(OCT_data, OCT_included==1)
Data_longitudinal$Verlaufsform <- as.factor(Data_longitudinal$Verlaufsform)
Data_longitudinal$Geschlecht <- as.factor(Data_longitudinal$Geschlecht)
Data_Longitudinal_RRMS <- subset(Data_longitudinal, Verlaufsform==1)
Data_Longitudinal_RRMS$Geschlecht <- as.factor(Data_Longitudinal_RRMS$Geschlecht)
Data_Longitudinal_PPMS <- subset(Data_longitudinal, Verlaufsform==2)
Data_Longitudinal_PPMS$Geschlecht <- as.factor(Data_Longitudinal_PPMS$Geschlecht)
Data_Longitudinal_SPMS <- subset(Data_longitudinal, Verlaufsform==3)
Data_Longitudinal_SPMS$Geschlecht <- as.factor(Data_Longitudinal_SPMS$Geschlecht)
Data_LMM_HC <- subset(OCT_data, OCT_included==1 & Verlaufsform==0)
Data_LMM_RRMS <- subset(OCT_data, OCT_included==1 & Verlaufsform==1)
Data_LMM_PPMS <- subset(OCT_data, OCT_included==1 & Verlaufsform==2)
Data_LMM_SPMS <- subset(OCT_data, OCT_included==1 & Verlaufsform==3)

#Verlaufsform=1 --> RRMS
#Verlaufsform=2 --> PPMS
#Verlaufsform=3 --> SPMS


#--------------------------------First part-----------------------------------#
#--------------------Visus HC (high contrast visual acuity)-------------------#  
#---RRMS---#
LMM_DD_predicts_RR <- lmer(VisusHC~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID) + (1|ID/ID_Eye), 
                     data=Data_Longitudinal_RRMS)
summary(LMM_DD_predicts_RR)
isSingular(LMM_DD_predicts_RR)
#singular fit

#Random effects:
# Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.020380 0.14276 
#ID        (Intercept)            0.014253 0.11939 
#ID.1      Disease_Duration_FINAL 0.000000 0.00000    

#RE of DiseaseDuration (Disease_DUration/ID) explains no variance --> reduce complexity of the model by
#eliminating this component

LMM_DD_predicts_RR1 <- lmer(VisusHC~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                      data=Data_Longitudinal_RRMS)
summary(LMM_DD_predicts_RR1)
#model fit without singularity
#best model fit - presented in paper

#to be sure that this is a better model compared to a model with a random slope, 
#check alternative model as well and compare both models

LMM_DD_predicts_RR2 <- lmer(VisusHC~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                            data=Data_Longitudinal_RRMS)
summary(LMM_DD_predicts_RR2)
anova(LMM_DD_predicts_RR2, LMM_DD_predicts_RR1)
#random intercept more important than random slope

#add age
LMM_DD_predicts_RR1_Age <- lmer(VisusHC~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                          data=Data_Longitudinal_RRMS)
summary(LMM_DD_predicts_RR1_Age)
anova(LMM_DD_predicts_RR1_Age, LMM_DD_predicts_RR1)
#age doesn't improve model fit

#add sex
LMM_DD_predicts_RR1_sex <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                          data=Data_Longitudinal_RRMS)
summary(LMM_DD_predicts_RR1_sex)
anova(LMM_DD_predicts_RR1_sex, LMM_DD_predicts_RR1)
#sex doesn't improve model fit

#add age+sex
LMM_DD_predicts_RR_both <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Longitudinal_RRMS)
summary(LMM_DD_predicts_RR_both)
anova(LMM_DD_predicts_RR_both, LMM_DD_predicts_RR1)
#adding both age and sex doesn't improve model fit 


#---PPMS---#
LMM_DD_predicts_PP <- lmer(VisusHC~Disease_Duration_FINAL +(0+Disease_Duration_FINAL|ID) + (1|ID/ID_Eye), 
                           data=Data_Longitudinal_PPMS)
summary(LMM_DD_predicts_PP)
isSingular(LMM_DD_predicts_PP)
#no singular fit
#best model fit - presented in paper

#Random effects:
#Groups    Name                   Variance  Std.Dev.
#ID_Eye.ID (Intercept)            6.150e-03 0.07842 
#ID        (Intercept)            1.862e-02 0.13647 
#ID.1      Disease_Duration_FINAL 3.697e-05 0.00608 
#Residual                         8.263e-03 0.09090 

#RE of DiseaseDuration (Disease_DUration/ID) explains very little variance, but converges

#add age
LMM_DD_predicts_PP_age <- lmer(VisusHC~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID) +(1|ID/ID_Eye), 
                         data=Data_Longitudinal_PPMS)
summary(LMM_DD_predicts_PP_age)
anova(LMM_DD_predicts_PP_age, LMM_DD_predicts_PP)
#age doesn't improve model fit

#add sex
LMM_DD_predicts_PP_sex <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht+ (0+Disease_Duration_FINAL|ID) + (1|ID/ID_Eye), 
                         data=Data_Longitudinal_PPMS)
summary(LMM_DD_predicts_PP_sex)
anova(LMM_DD_predicts_PP_sex, LMM_DD_predicts_PP)
#sex doesn't improve model fit

#add sex+age
LMM_DD_predicts_PP_both <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID)+(1|ID/ID_Eye), 
                          data=Data_Longitudinal_PPMS)
summary(LMM_DD_predicts_PP_both)
anova(LMM_DD_predicts_PP_both, LMM_DD_predicts_PP)
#adding both sex and age doesn't improve model fit 


#---SPMS---#
LMM_DD_predicts_SP <- lmer(VisusHC~Disease_Duration_FINAL + +(0+Disease_Duration_FINAL|ID) + (1|ID/ID_Eye), 
                           data=Data_Longitudinal_SPMS)
summary(LMM_DD_predicts_SP)
#In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#Model failed to converge with max|grad| = 0.00392908 (tol = 0.002, component 1)

#Random effects:
#Groups    Name                   Variance  Std.Dev.
#ID_Eye.ID (Intercept)            1.353e-02 0.116301
#ID        (Intercept)            2.283e-02 0.151103
#ID.1      Disease_Duration_FINAL 3.782e-06 0.001945
#Residual                         5.256e-03 0.072496

#RE of DiseaseDuration (Disease_Duration/ID) explains little variance --> reduce complexity of the model by
#eliminating this component

LMM_DD_predicts_SP2 <- lmer(VisusHC~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                     data=Data_Longitudinal_SPMS)
summary(LMM_DD_predicts_SP2)
#best model fit - presented in paper

#to be sure that this is a better model compared to a model with a random slope, 
#check alternative model as well and compare both models
LMM_DD_predicts_SP3 <- lmer(VisusHC~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                      data=Data_Longitudinal_SPMS)
summary(LMM_DD_predicts_SP3)
isSingular(LMM_DD_predicts_SP3)
anova(LMM_DD_predicts_SP3, LMM_DD_predicts_SP2)
#random slope more important than random intercept

#add age
LMM_DD_predicts_SP_age <- lmer(VisusHC~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                         data=Data_Longitudinal_SPMS)
summary(LMM_DD_predicts_SP_age)
anova(LMM_DD_predicts_SP_age, LMM_DD_predicts_SP2)
#age doesn't improve model fit

#add sex
LMM_DD_predicts_SP_sex <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                         data=Data_Longitudinal_SPMS)
summary(LMM_DD_predicts_SP_sex)
anova(LMM_DD_predicts_SP_sex, LMM_DD_predicts_SP2)
#sex doesn't improve model fi

#add both
LMM_DD_predicts_SP_both <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                          data=Data_Longitudinal_SPMS)
summary(LMM_DD_predicts_SP_both)
anova(LMM_DD_predicts_SP_both, LMM_DD_predicts_SP2)
#adding both sex and age doesn't improve model fit
#best fitting model= LMM_DD_predicts_SP2


#--------------------LATENCY-------------------#  
#---RRMS---#
LMM_Latenz_DD_RR <- lmer(Latenz~Disease_Duration_Latency + (0+Disease_Duration_Latency|ID) + (1|ID/ID_Eye), 
                         data=Data_Longitudinal_RRMS)
summary(LMM_Latenz_DD_RR)
#best model fit - presented in paper

#add age
LMM_Latenz_DD_RR_age <- lmer(Latenz~Disease_Duration_Latency + Age_Latency + (0+Disease_Duration_Latency|ID) + (1|ID/ID_Eye), 
                             data=Data_Longitudinal_RRMS)
summary(LMM_Latenz_DD_RR_age)
anova(LMM_Latenz_DD_RR_age, LMM_Latenz_DD_RR)
#age doesn't improve model fit

#add sex
LMM_Latenz_DD_RR_sex <- lmer(Latenz~Disease_Duration_Latency + Geschlecht + (0+Disease_Duration_Latency|ID) + (1|ID/ID_Eye), 
                             data=Data_Longitudinal_RRMS)
summary(LMM_Latenz_DD_RR_sex)
anova(LMM_Latenz_DD_RR_sex, LMM_Latenz_DD_RR)
#sex doesn't improve model fit

#add age+sex
LMM_Latenz_DD_RR_both <- lmer(Latenz~Disease_Duration_Latency + Geschlecht + Age_Latency + (0+Disease_Duration_Latency|ID) + (1|ID/ID_Eye), 
                              data=Data_Longitudinal_RRMS)
summary(LMM_Latenz_DD_RR_both)
anova(LMM_Latenz_DD_RR_both, LMM_Latenz_DD_RR)
#adding both sex and age doesn't improve model fit 

#---PPMS---#
LMM_Latenz_DD_PP <- lmer(Latenz~Disease_Duration_Latency + (0+Disease_Duration_Latency|ID)+ (1|ID/ID_Eye), 
                         data=Data_Longitudinal_PPMS)
summary(LMM_Latenz_DD_PP)
isSingular(LMM_Latenz_DD_PP)
#singular fit, no random effect for DD*ID necessary

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             58.46    7.646  
#ID        (Intercept)            100.26   10.013  
#ID.1      Disease_Duration_FINAL   0.00    0.000  
#Residual                          64.66    8.041 

#RE of DiseaseDuration (Disease_Duration/ID) explains no variance --> reduce complexity of the model by
#eliminating this component

LMM_Latenz_DD_PP2 <- lmer(Latenz~Disease_Duration_Latency + (1|ID/ID_Eye), 
                          data=Data_Longitudinal_PPMS)
summary(LMM_Latenz_DD_PP2)
#best model fit - presented in paper

#to be sure that this is a better model compared to a model with a random slope, 
#check alternative model as well and compare both models
LMM_Latenz_DD_PP3 <- lmer(Latenz~Disease_Duration_Latency + (0+Disease_Duration_Latency|ID), 
                          data=Data_Longitudinal_PPMS)
summary(LMM_Latenz_DD_PP3)
anova(LMM_Latenz_DD_PP3, LMM_Latenz_DD_PP2)
#P2 much better fit than P3

#add age
LMM_Latenz_DD_PP_age <- lmer(Latenz~Disease_Duration_Latency + Age_Latency + (1|ID/ID_Eye), 
                             data=Data_Longitudinal_PPMS)
summary(LMM_Latenz_DD_PP_age)
anova(LMM_Latenz_DD_PP_age, LMM_Latenz_DD_PP2)
#age doesn't improve model fit

#add sex
LMM_Latenz_DD_PP_sex <- lmer(Latenz~Disease_Duration_Latency + Geschlecht + (1|ID/ID_Eye), 
                             data=Data_Longitudinal_PPMS)
summary(LMM_Latenz_DD_PP_sex)
anova(LMM_Latenz_DD_PP_sex, LMM_Latenz_DD_PP2)
#sex doesn't improve model fit

#add age+sex
LMM_Latenz_DD_PP_both <- lmer(Latenz~Disease_Duration_Latency + Geschlecht + Age_Latency + (1|ID/ID_Eye), 
                              data=Data_Longitudinal_PPMS)
summary(LMM_Latenz_DD_PP_both)
anova(LMM_Latenz_DD_PP_both, LMM_Latenz_DD_PP2)
#adding both age + sex doesn't improve model fit

#---SPMS---#
LMM_Latenz_DD_SP <- lmer(Latenz~Disease_Duration_Latency + (0+Disease_Duration_Latency|ID)+ (1|ID/ID_Eye), 
                         data=Data_Longitudinal_SPMS)
summary(LMM_Latenz_DD_SP)
isSingular(LMM_Latenz_DD_SP)
#singular fit, no random effect for DD*ID necessary

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            150.8    12.28   
#ID        (Intercept)            277.1    16.65   
#ID.1      Disease_Duration_FINAL   0.0     0.00   
#Residual                         116.6    10.80  

#RE of DiseaseDuration (Disease_Duration/ID) explains no variance --> reduce complexity of the model by
#eliminating this component

LMM_Latenz_DD_SP2 <- lmer(Latenz~Disease_Duration_Latency + (1|ID/ID_Eye), 
                          data=Data_Longitudinal_SPMS)
summary(LMM_Latenz_DD_SP2)
#best model fit - presented in paper

#to be sure that this is a better model compared to a model with a random slope, 
#check alternative model as well and compare both models
LMM_Latenz_DD_SP3 <- lmer(Latenz~Disease_Duration_Latency + (0+Disease_Duration_Latency|ID), 
                          data=Data_Longitudinal_SPMS)
summary(LMM_Latenz_DD_SP3)
anova(LMM_Latenz_DD_SP3, LMM_Latenz_DD_SP2)
#P2 much better fit than P3

#add age
LMM_Latenz_DD_SP_age <- lmer(Latenz~Disease_Duration_Latency + Age_Latency + (1|ID/ID_Eye), 
                             data=Data_Longitudinal_SPMS)
summary(LMM_Latenz_DD_SP_age)
anova(LMM_Latenz_DD_SP_age, LMM_Latenz_DD_SP2)
#age doesn't improve model fit

#add sex
LMM_Latenz_DD_SP_sex <- lmer(Latenz~Disease_Duration_Latency + Geschlecht + (1|ID/ID_Eye), 
                             data=Data_Longitudinal_SPMS)
summary(LMM_Latenz_DD_SP_sex)
anova(LMM_Latenz_DD_SP_sex, LMM_Latenz_DD_SP2)
#sex doesn't improve model fit

#add age+sex
LMM_Latenz_DD_SP_both <- lmer(Latenz~Disease_Duration_Latency + Geschlecht + Age_Latency + (1|ID/ID_Eye), 
                              data=Data_Longitudinal_SPMS)
summary(LMM_Latenz_DD_SP_both)
anova(LMM_Latenz_DD_SP_both, LMM_Latenz_DD_SP2)
#adding both age + sex doesn't improve model fit


#--------------------------------Second part----------------------------------#
#-----------------------------------pRNFL-------------------------------------#  
#----HC----#

#-FIRST (ONLY) INTERVAL 0-3,5 YEARS "DISEASE DURATION" (=Time since baseline)-# 
Data_LMM_HC$Verlaufsform<-as.factor(Data_LMM_HC$Verlaufsform)
data_HC_Cat1 <- subset(Data_LMM_HC, DiseaseDurationPPMSCategories_LMM==0)
data_HC_Cat1$Verlaufsform<-as.factor(data_HC_Cat1$Verlaufsform)

LMM_globalpRNFL_HC_Cat1 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_HC_Cat1)
summary(LMM_globalpRNFL_HC_Cat1)
confint(LMM_globalpRNFL_HC_Cat1)
#best model fit - presented in paper

#Add Covariates
#a: Age
LMM_globalpRNFL_HC_Cat1age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_HC_Cat1)
summary(LMM_globalpRNFL_HC_Cat1age)
anova(LMM_globalpRNFL_HC_Cat1age, LMM_globalpRNFL_HC_Cat1)
#adding age doesn't sign. improve model fit

#b: Sex
data_HC_Cat1$Geschlecht<-as.factor(data_HC_Cat1$Geschlecht)
LMM_globalpRNFL_HC_Cat1sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_HC_Cat1)
summary(LMM_globalpRNFL_HC_Cat1sex)
anova(LMM_globalpRNFL_HC_Cat1sex, LMM_globalpRNFL_HC_Cat1)
#sex no sign. effect

#c: age+sex
LMM_globalpRNFL_HC_Cat1both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                    data=data_HC_Cat1)
summary(LMM_globalpRNFL_HC_Cat1both)
anova(LMM_globalpRNFL_HC_Cat1both, LMM_globalpRNFL_HC_Cat1)
#adding both sex and age doesn't improve model fit

#---RRMS---#
#---------------FIRST INTERVAL 0-3,5 YEARS DISEASE DURATION------------------#
data_RRMS_Cat1 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==0)
data_RRMS_Cat1$Geschlecht<-as.factor(data_RRMS_Cat1$Geschlecht)

LMM_globalpRNFL_RR_Cat1 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_RRMS_Cat1)
summary(LMM_globalpRNFL_RR_Cat1)
confint(LMM_globalpRNFL_RR_Cat1)
#best model fit - presented in paper

#add age as covariate
LMM_globalpRNFL_RR_Cat1age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_RRMS_Cat1)
summary(LMM_globalpRNFL_RR_Cat1age)
anova(LMM_globalpRNFL_RR_Cat1, LMM_globalpRNFL_RR_Cat1age)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_RR_Cat1sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_RRMS_Cat1)
summary(LMM_globalpRNFL_RR_Cat1sex)
anova(LMM_globalpRNFL_RR_Cat1, LMM_globalpRNFL_RR_Cat1sex)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_RR_Cat1both <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_RRMS_Cat1)
summary(LMM_globalpRNFL_RR_Cat1both)
anova(LMM_globalpRNFL_RR_Cat1, LMM_globalpRNFL_RR_Cat1both)
#adding both sex and age doesn't improve model fit 

#---------------SECOND INTERVAL 3,6-5,5 YEARS DISEASE DURATION---------------# 
data_RRMS_Cat2 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==1)
data_RRMS_Cat2$Geschlecht<-as.factor(data_RRMS_Cat2$Geschlecht)

LMM_globalpRNFL_RR_Cat2 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_RRMS_Cat2)
summary(LMM_globalpRNFL_RR_Cat2)
confint(LMM_globalpRNFL_RR_Cat2)
#best model fit - presented in paper

#add age
LMM_globalpRNFL_RR_Cat2age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age +(1|ID/ID_Eye)  + (0+Disease_Duration_FINAL|ID), 
                                   data=data_RRMS_Cat2)
summary(LMM_globalpRNFL_RR_Cat2age)
anova(LMM_globalpRNFL_RR_Cat2age,LMM_globalpRNFL_RR_Cat2)
#age doesn't improve model fit

#sex
LMM_globalpRNFL_RR_Cat2sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye)  + (0+Disease_Duration_FINAL|ID), 
                                   data=data_RRMS_Cat2)
summary(LMM_globalpRNFL_RR_Cat2sex)
anova(LMM_globalpRNFL_RR_Cat2sex,LMM_globalpRNFL_RR_Cat2)
#age doesn't improve model fit

#both
LMM_globalpRNFL_RR_Cat2both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye)  + (0+Disease_Duration_FINAL|ID), 
                                    data=data_RRMS_Cat2)
summary(LMM_globalpRNFL_RR_Cat2both)
anova(LMM_globalpRNFL_RR_Cat2both,LMM_globalpRNFL_RR_Cat2)
#adding both sex and age doesn't improve model fit

#---------------Third INTERVAL 5,6-7,5  YEARS DISEASE DURATION---------------#
data_RRMS_Cat3 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==2)
data_RRMS_Cat3$Geschlecht<-as.factor(data_RRMS_Cat3$Geschlecht)

LMM_globalpRNFL_RR_Cat3 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_RRMS_Cat3)
summary(LMM_globalpRNFL_RR_Cat3)
confint(LMM_globalpRNFL_RR_Cat3)
#best model fit - presented in paper

#add age
LMM_globalpRNFL_RR_Cat3age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                   data=data_RRMS_Cat3)
summary(LMM_globalpRNFL_RR_Cat3age)
anova(LMM_globalpRNFL_RR_Cat3, LMM_globalpRNFL_RR_Cat3age)
#age doesn't improve model fit

#add sex
LMM_globalpRNFL_RR_Cat3sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                   data=data_RRMS_Cat3)
summary(LMM_globalpRNFL_RR_Cat3sex)
anova(LMM_globalpRNFL_RR_Cat3sex, LMM_globalpRNFL_RR_Cat3)

#add both
LMM_globalpRNFL_RR_Cat3both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                    data=data_RRMS_Cat3)
summary(LMM_globalpRNFL_RR_Cat3both)
anova(LMM_globalpRNFL_RR_Cat3both, LMM_globalpRNFL_RR_Cat3)
#adding both sex and age doesn't improve model fit

#--------------FOURTH INTERVAL 7,6-10,5 YEARS DISEASE DURATION----------------# 
data_RRMS_Cat4 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==3)
data_RRMS_Cat4$Geschlecht<-as.factor(data_RRMS_Cat4$Geschlecht)

LMM_globalpRNFL_RRMS_Cat4 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_RRMS_Cat4)
summary(LMM_globalpRNFL_RRMS_Cat4)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_globalpRNFL_RRMS_Cat4)
#singular fit

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            3.206e+01 5.662e+00
#ID        (Intercept)            8.237e+01 9.076e+00
#ID.1      Disease_Duration_FINAL 6.536e-09 8.084e-05
#Residual                         2.572e+00 1.604e+00

#Variance of RE of disease duration very low ---> Try without RE of DD

LMM_globalpRNFL_RRMS_Cat4_2 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                                    data=data_RRMS_Cat4)
summary(LMM_globalpRNFL_RRMS_Cat4_2)
confint(LMM_globalpRNFL_RRMS_Cat4_2)
#model without RE converges
#best model fit - presented in paper

#a:add age
LMM_globalpRNFL_RRMS_Cat4age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age +(1|ID/ID_Eye), 
                                     data=data_RRMS_Cat4)
summary(LMM_globalpRNFL_RRMS_Cat4age)
anova(LMM_globalpRNFL_RRMS_Cat4age, LMM_globalpRNFL_RRMS_Cat4_2)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_RRMS_Cat4sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye), 
                                     data=data_RRMS_Cat4)
summary(LMM_globalpRNFL_RRMS_Cat4sex)
anova(LMM_globalpRNFL_RRMS_Cat4sex, LMM_globalpRNFL_RRMS_Cat4_2)
#adding sex doesn't improve model fit

#b: both
LMM_globalpRNFL_RRMS_Cat4sexage <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                                        data=data_RRMS_Cat4)
summary(LMM_globalpRNFL_RRMS_Cat4sexage)
anova(LMM_globalpRNFL_RRMS_Cat4_2, LMM_globalpRNFL_RRMS_Cat4sexage)
#adding both sex and age doesn't improve model fit

#---------------FIFTH INTERVAL 10,6-13,5 YEARS DISEASE DURATION--------------# 
data_RRMS_Cat5 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==4)
data_RRMS_Cat5$Geschlecht<-as.factor(data_RRMS_Cat5$Geschlecht)

LMM_globalpRNFL_RR_Cat5 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_RRMS_Cat5)
summary(LMM_globalpRNFL_RR_Cat5)
confint(LMM_globalpRNFL_RR_Cat5)
#best model fit - presented in paper

#a: add age
LMM_globalpRNFL_RR_Cat5age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_RRMS_Cat5)
summary(LMM_globalpRNFL_RR_Cat5age)
anova(LMM_globalpRNFL_RR_Cat5, LMM_globalpRNFL_RR_Cat5age)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_RR_Cat5sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_RRMS_Cat5)
summary(LMM_globalpRNFL_RR_Cat5sex)
anova(LMM_globalpRNFL_RR_Cat5, LMM_globalpRNFL_RR_Cat5sex)
#adding sex doesn't improve model fit

#b: both
LMM_globalpRNFL_RR_Cat5sexage <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=data_RRMS_Cat5)
summary(LMM_globalpRNFL_RR_Cat5sexage)
anova(LMM_globalpRNFL_RR_Cat5sexage, LMM_globalpRNFL_RR_Cat5)
#adding both sex and age doesn't improve model fit

#-------------SIXTH INTERVAL 13,6-16,5 YEARS DISEASE DURATION-------------# 
data_RRMS_SPCat1 <- subset(Data_LMM_RRMS, DiseaseDurationSPMSCategories_LMM==1)
data_RRMS_SPCat1$Geschlecht<-as.factor(data_RRMS_SPCat1$Geschlecht)

LMM_globalpRNFL_RR_SPCat1 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_RRMS_SPCat1)
summary(LMM_globalpRNFL_RR_SPCat1)
confint(LMM_globalpRNFL_RR_SPCat1)
#best model fit - presented in paper

#add age as covariate
LMM_globalpRNFL_RR_SPCat1age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_RRMS_SPCat1)
summary(LMM_globalpRNFL_RR_SPCat1age)
anova(LMM_globalpRNFL_RR_SPCat1, LMM_globalpRNFL_RR_SPCat1age)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_RR_SPCat1_sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=data_RRMS_SPCat1)
summary(LMM_globalpRNFL_RR_SPCat1_sex)
anova(LMM_globalpRNFL_RR_SPCat1, LMM_globalpRNFL_RR_SPCat1_sex)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_RR_SPCat1_both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                       data=data_RRMS_SPCat1)
summary(LMM_globalpRNFL_RR_SPCat1_both)
anova(LMM_globalpRNFL_RR_SPCat1, LMM_globalpRNFL_RR_SPCat1_both)
#adding both sex and age doesn't improve model fit


######################################################################################### 
###########################               PPMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 0-3,5 YEARS DISEASE DURATION      #################
data_PPMS_Cat1 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==0)
data_PPMS_Cat1$Geschlecht<-as.factor(data_PPMS_Cat1$Geschlecht)

LMM_globalpRNFL_PP_Cat1 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_PPMS_Cat1)
summary(LMM_globalpRNFL_PP_Cat1)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_globalpRNFL_PP_Cat1)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             5.476   2.340   
#ID        (Intercept)            95.303   9.762   
#ID.1      Disease_Duration_FINAL  0.000   0.000   
#Residual                          1.274   1.129 

#RE of Disease duration doesn't explain any variance --> reduce complexity of the model by eliminating this RE

LMM_globalpRNFL_PP_Cat1_2 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                                  data=data_PPMS_Cat1)
summary(LMM_globalpRNFL_PP_Cat1_2)
confint(LMM_globalpRNFL_PP_Cat1_2)
#best model fit - presented in paper

#a: add age
LMM_globalpRNFL_PP_Cat1age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                                   data=data_PPMS_Cat1)
summary(LMM_globalpRNFL_PP_Cat1age)
anova(LMM_globalpRNFL_PP_Cat1age, LMM_globalpRNFL_PP_Cat1_2)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_PPCat1sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                               data=data_PPMS_Cat1)
summary(LMM_globalpRNFL_PPCat1sex)
anova(LMM_globalpRNFL_PPCat1sex, LMM_globalpRNFL_PP_Cat1_2)
#adding sex doesn't improve model fit 

#c: both
LMM_globalpRNFL_PPCat1both <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye), 
                                data=data_PPMS_Cat1)
summary(LMM_globalpRNFL_PPCat1both)
anova(LMM_globalpRNFL_PPCat1both, LMM_globalpRNFL_PP_Cat1_2)
#adding age and sex doesn't improve model fit 

##################      SECOND INTERVAL 3,6-5,5 YEARS DISEASE DURATION   ################
data_PPMS_Cat2 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==1)
data_PPMS_Cat2$Geschlecht<-as.factor(data_PPMS_Cat2$Geschlecht)

LMM_globalpRNFL_PP_Cat2 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_PPMS_Cat2)
summary(LMM_globalpRNFL_PP_Cat2)
#boundary (singular) fit: see help('isSingular')

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)              4.440   2.107  
#ID        (Intercept)            127.389  11.287  
#ID.1      Disease_Duration_FINAL   0.000   0.000  
#Residual                           1.726   1.314  

#RE of Disease duration doesn't explain any variance --> reduce complexity of the model by eliminating this RE

LMM_globalpRNFL_PP_Cat2_2 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                                  data=data_PPMS_Cat2)
summary(LMM_globalpRNFL_PP_Cat2_2)
confint(LMM_globalpRNFL_PP_Cat2_2)
#best model fit - presented in paper

##a: add age
LMM_globalpRNFL_PP_Cat2age <- lmer(globalpRNFL~Disease_Duration_FINAL+ Age +(1|ID/ID_Eye), 
                                   data=data_PPMS_Cat2)
summary(LMM_globalpRNFL_PP_Cat2age)
anova(LMM_globalpRNFL_PP_Cat2age, LMM_globalpRNFL_PP_Cat2_2)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_PP_Cat2sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                               data=data_PPMS_Cat2)
summary(LMM_globalpRNFL_PP_Cat2sex)
anova(LMM_globalpRNFL_PP_Cat2sex, LMM_globalpRNFL_PP_Cat2_2)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_PP_Cat2both <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye), 
                                data=data_PPMS_Cat2)
summary(LMM_globalpRNFL_PP_Cat2both)
anova(LMM_globalpRNFL_PP_Cat2both, LMM_globalpRNFL_PP_Cat2_2)
#adding both sex and age doesn't improve model fit 

##################     THIRD INTERVAL 5,6-7,5 YEARS DISEASE DURATION    ################# 
data_PPMS_Cat3 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==2)
data_PPMS_Cat3$Geschlecht<-as.factor(data_PPMS_Cat3$Geschlecht)

LMM_globalpRNFL_PP_Cat3 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_PPMS_Cat3)
summary(LMM_globalpRNFL_PP_Cat3)
confint(LMM_globalpRNFL_PP_Cat3)
#best model fit - presented in paper

#add age
LMM_globalpRNFL_PP_Cat3age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age +  (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_PPMS_Cat3)
summary(LMM_globalpRNFL_PP_Cat3age)
anova(LMM_globalpRNFL_PP_Cat3age, LMM_globalpRNFL_PP_Cat3)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_PP_Cat3sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht +  (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_PPMS_Cat3)
summary(LMM_globalpRNFL_PP_Cat3sex)
anova(LMM_globalpRNFL_PP_Cat3sex, LMM_globalpRNFL_PP_Cat3)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_PP_Cat3both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                    data=data_PPMS_Cat3)
summary(LMM_globalpRNFL_PP_Cat3both)
anova(LMM_globalpRNFL_PP_Cat3both, LMM_globalpRNFL_PP_Cat3)
#adding both sex and age doesn't improve model fit

##################      FOURTH INTERVAL 7,6-10,5 YEARS DISEASE DURATION   ###############
data_PPMS_Cat4 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==3)
data_PPMS_Cat4$Geschlecht<-as.factor(data_PPMS_Cat4$Geschlecht)

LMM_globalpRNFL_PP_Cat4 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_PPMS_Cat4)
summary(LMM_globalpRNFL_PP_Cat4)
confint(LMM_globalpRNFL_PP_Cat4)
#best model fit - presented in paper

#add age
LMM_globalpRNFL_PP_Cat4age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_PPMS_Cat4)
summary(LMM_globalpRNFL_PP_Cat4age)
anova(LMM_globalpRNFL_PP_Cat4age, LMM_globalpRNFL_PP_Cat4)
#Warning message:
#In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#Model failed to converge with max|grad| = 0.00215875 (tol = 0.002, component 1)
#doesn't improve model fit but resuls in an overly complex model with non-convergence

#add sex
LMM_globalpRNFL_PP_Cat4sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_PPMS_Cat4)
summary(LMM_globalpRNFL_PP_Cat4sex)
anova(LMM_globalpRNFL_PP_Cat4sex, LMM_globalpRNFL_PP_Cat4)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_PP_Cat4both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                    data=data_PPMS_Cat4)
summary(LMM_globalpRNFL_PP_Cat4both)
anova(LMM_globalpRNFL_PP_Cat4both, LMM_globalpRNFL_PP_Cat4)
#adding both sex and age doesn't improve model fit

#################      FIFTH INTERVAL 10,6-13,5 YEARS DISEASE DURATION   ################ 
data_PPMS_Cat5 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==4)
data_PPMS_Cat5$Geschlecht<-as.factor(data_PPMS_Cat5$Geschlecht)

LMM_globalpRNFL_PP_Cat5 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_PPMS_Cat5)
summary(LMM_globalpRNFL_PP_Cat5)
confint(LMM_globalpRNFL_PP_Cat5)

#add age
LMM_globalpRNFL_PP_Cat5age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_PPMS_Cat5)
summary(LMM_globalpRNFL_PP_Cat5age)
anova(LMM_globalpRNFL_PP_Cat5age, LMM_globalpRNFL_PP_Cat5)
#sign. better fit including age
confint (LMM_globalpRNFL_PP_Cat5age)
#best model fit - presented in paper

#add sex
LMM_globalpRNFL_PP_Cat5sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                   data=data_PPMS_Cat5)
summary(LMM_globalpRNFL_PP_Cat5sex)
anova(LMM_globalpRNFL_PP_Cat5sex, LMM_globalpRNFL_PP_Cat5)
#adding sex doesn't improve model fit 

#adding both compared to age only
LMM_globalpRNFL_PP_Cat5both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                    data=data_PPMS_Cat5)
summary(LMM_globalpRNFL_PP_Cat5both)
anova(LMM_globalpRNFL_PP_Cat5both, LMM_globalpRNFL_PP_Cat5age)
#adding sex in addition to age doesn't sign. improve model fit. 

##################     SIXTH INTERVAL 13,6-16,5 YEARS DISEASE DURATION   ################ 
data_PPMS_SPCat1 <- subset(Data_LMM_PPMS, DiseaseDurationSPMSCategories_LMM==1)
data_PPMS_SPCat1$Geschlecht<-as.factor(data_PPMS_SPCat1$Geschlecht)

LMM_globalpRNFL_PP_SPCat1 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_PPMS_SPCat1)
summary(LMM_globalpRNFL_PP_SPCat1)
confint(LMM_globalpRNFL_PP_SPCat1)
#best model fit - presented in paper

#add age
LMM_globalpRNFL_PP_SPCat1age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_PPMS_SPCat1)
summary(LMM_globalpRNFL_PP_SPCat1age)
anova(LMM_globalpRNFL_PP_SPCat1age, LMM_globalpRNFL_PP_SPCat1)
#adding age doesn't sign. improve model fit

#add sex
LMM_globalpRNFL_PP_SPCat1sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_PPMS_SPCat1)
summary(LMM_globalpRNFL_PP_SPCat1sex)
anova(LMM_globalpRNFL_PP_SPCat1sex, LMM_globalpRNFL_PP_SPCat1)
#adding sex doesn't improve model fit

#add both
LMM_globalpRNFL_PP_SPCat1both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_PPMS_SPCat1)
summary(LMM_globalpRNFL_PP_SPCat1both)
anova(LMM_globalpRNFL_PP_SPCat1both, LMM_globalpRNFL_PP_SPCat1)
#adding both sex and age doesn't improve model fit

#################      SEVENTH INTERVAL 16,6-20,5 YEARS DISEASE DURATION   ##############
data_PPMS_SPCat2 <- subset(Data_LMM_PPMS, DiseaseDurationSPMSCategories_LMM==2)
data_PPMS_SPCat2$Geschlecht<-as.factor(data_PPMS_SPCat2$Geschlecht)

LMM_globalpRNFL_PP_SPCat2 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_PPMS_SPCat2)
summary(LMM_globalpRNFL_PP_SPCat2)

#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_globalpRNFL_PP_SPCat2)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            14.2741  3.7781  
#ID        (Intercept)             0.0000  0.0000  
#ID.1      Disease_Duration_FINAL  0.1861  0.4314  
#Residual                          1.6294  1.2765 

#Random Intercept per subject doesn't explain any variance --> reduce complexity of the model by eliminating this RE

LMM_globalpRNFL_PP_SPCat2_2 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                    data=data_PPMS_SPCat2)
summary(LMM_globalpRNFL_PP_SPCat2_2)
confint(LMM_globalpRNFL_PP_SPCat2_2)
#best model fit - presented in paper

#a: age
LMM_globalpRNFL_PP_SPCat2age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_PPMS_SPCat2)
summary(LMM_globalpRNFL_PP_SPCat2age)
anova(LMM_globalpRNFL_PP_SPCat2age, LMM_globalpRNFL_PP_SPCat2_2)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_PP_SPCat2sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_PPMS_SPCat2)
summary(LMM_globalpRNFL_PP_SPCat2sex)
anova(LMM_globalpRNFL_PP_SPCat2_2, LMM_globalpRNFL_PP_SPCat2sex)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_PP_SPCat2both <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=data_PPMS_SPCat2)
summary(LMM_globalpRNFL_PP_SPCat2both)
anova(LMM_globalpRNFL_PP_SPCat2_2, LMM_globalpRNFL_PP_SPCat2both)
#adding both sex and age doesn't improve model fit

######################################################################################### 
###########################               SPMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 3,5-12,5 YEARS DISEASE DURATION      ##############

data_SPMS_SPCat1 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==0)
data_SPMS_SPCat1$Geschlecht<-as.factor(data_SPMS_SPCat1$Geschlecht)

LMM_globalpRNFL_SP_SPCat1 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_SPMS_SPCat1)
summary(LMM_globalpRNFL_SP_SPCat1)
confint(LMM_globalpRNFL_SP_SPCat1)

#a: age
LMM_globalpRNFL_SP_SPCat1age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye)  + (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat1)
summary(LMM_globalpRNFL_SP_SPCat1age)
anova(LMM_globalpRNFL_SP_SPCat1age, LMM_globalpRNFL_SP_SPCat1)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_SP_SPCat1sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye)  + (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat1)
summary(LMM_globalpRNFL_SP_SPCat1sex)
anova(LMM_globalpRNFL_SP_SPCat1sex, LMM_globalpRNFL_SP_SPCat1)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_SP_SPCat1both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                      data=data_SPMS_SPCat1)
#boundary (singular) fit: see help('isSingular')

summary(LMM_globalpRNFL_SP_SPCat1both)
anova(LMM_globalpRNFL_SP_SPCat1both, LMM_globalpRNFL_SP_SPCat1)
isSingular(LMM_globalpRNFL_SP_SPCat1both)
# sign. better fit (but singular fit)

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            1.675e+01 4.093e+00
#ID        (Intercept)            1.483e+02 1.218e+01
#ID.1      Disease_Duration_FINAL 3.068e-08 1.751e-04
#Residual                         2.788e+00 1.670e+00

#RE of disease duration explains almost no variance
#try model with both age and sex without second random effect
LMM_globalpRNFL_SP_SPCat1both1 <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                                       data=data_SPMS_SPCat1)
summary(LMM_globalpRNFL_SP_SPCat1both1)
anova(LMM_globalpRNFL_SP_SPCat1both1, LMM_globalpRNFL_SP_SPCat1)
isSingular(LMM_globalpRNFL_SP_SPCat1both1)
#not singular anymore and sign. better fit than model without both covariates but RE of disease duration 
confint(LMM_globalpRNFL_SP_SPCat1both1)
#best model fit - presented in paper

##################     SECOND INTERVAL 12,6-16,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat2 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==1)
data_SPMS_SPCat2$Geschlecht<-as.factor(data_SPMS_SPCat2$Geschlecht)

LMM_globalpRNFL_SP_SPCat2 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_SPMS_SPCat2)
summary(LMM_globalpRNFL_SP_SPCat2)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_globalpRNFL_SP_SPCat2)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             34.348   5.861  
#ID        (Intercept)            132.744  11.521  
#ID.1      Disease_Duration_FINAL   0.000   0.000  
#Residual                           3.142   1.772 

#RE of disease duration doesn't explain any variance --> eliminate this unneccessary RE to reduce model complexity

LMM_globalpRNFL_SP_SPCat2_new <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                                      data=data_SPMS_SPCat2)
summary(LMM_globalpRNFL_SP_SPCat2_new)
confint(LMM_globalpRNFL_SP_SPCat2_new)
#best model fit - presented in paper

#a: age
LMM_globalpRNFL_SP_SPCat2age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                                     data=data_SPMS_SPCat2)
summary(LMM_globalpRNFL_SP_SPCat2age)
anova(LMM_globalpRNFL_SP_SPCat2age, LMM_globalpRNFL_SP_SPCat2_new)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_SP_SPCat2sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat2)
summary(LMM_globalpRNFL_SP_SPCat2sex)
anova(LMM_globalpRNFL_SP_SPCat2sex, LMM_globalpRNFL_SP_SPCat2_new)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_SP_SPCat2both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht +  Age + (1|ID/ID_Eye), 
                                data=data_SPMS_SPCat2)
summary(LMM_globalpRNFL_SP_SPCat2both)
anova(LMM_globalpRNFL_SP_SPCat2both, LMM_globalpRNFL_SP_SPCat2_new)
#adding both sex and age doesn't improve model fit

##################     THIRD INTERVAL 16,6-20,5 YEARS DISEASE DURATION       ############
data_SPMS_SPCat3 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==2)
data_SPMS_SPCat3$Geschlecht<-as.factor(data_SPMS_SPCat3$Geschlecht)

LMM_globalpRNFL_SP_SPCat3 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_SPMS_SPCat3)
summary(LMM_globalpRNFL_SP_SPCat3)
confint(LMM_globalpRNFL_SP_SPCat3)
#best model fit - presented in paper

#a: age
LMM_globalpRNFL_SP_SPCat3age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat3)
summary(LMM_globalpRNFL_SP_SPCat3age)
anova(LMM_globalpRNFL_SP_SPCat3age, LMM_globalpRNFL_SP_SPCat3)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_SP_SPCat3sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat3)
summary(LMM_globalpRNFL_SP_SPCat3sex)
anova(LMM_globalpRNFL_SP_SPCat3sex, LMM_globalpRNFL_SP_SPCat3)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_SP_SPCat3both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                      data=data_SPMS_SPCat3)
summary(LMM_globalpRNFL_SP_SPCat3both)
anova(LMM_globalpRNFL_SP_SPCat3both, LMM_globalpRNFL_SP_SPCat3)
#adding both sex and age doesn't improve model fit

##################     FOURTH INTERVAL 20,6-25,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat4 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==3)
data_SPMS_SPCat4$Geschlecht<-as.factor(data_SPMS_SPCat4$Geschlecht)

LMM_globalpRNFL_SP_SPCat4 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_SPMS_SPCat4)
summary(LMM_globalpRNFL_SP_SPCat4)

#Warning messages:
# In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge with max|grad| = 0.00213345 (tol = 0.002, component 1)
# Model failed to converge with 1 negative eigenvalue: -7.5e-03 
                
#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            2.146e+01 4.6326204
#ID        (Intercept)            4.712e-07 0.0006865
#ID.1      Disease_Duration_FINAL 2.411e-01 0.4910437
#Residual                         5.026e+00 2.2418642

#Random Intercept of ID explains almost no variance --> reduce model to less complex model

LMM_globalpRNFL_SP_SPCat4_new <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                      data=data_SPMS_SPCat4)
summary(LMM_globalpRNFL_SP_SPCat4_new)
confint(LMM_globalpRNFL_SP_SPCat4_new)
isSingular(LMM_globalpRNFL_SP_SPCat4_new)
#convergence and no singularity
#best model fit - presented in paper

#a: age
LMM_globalpRNFL_SP_SPCat4age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age  + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat4)
summary(LMM_globalpRNFL_SP_SPCat4age)
anova(LMM_globalpRNFL_SP_SPCat4age, LMM_globalpRNFL_SP_SPCat4_new)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_SP_SPCat4sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht  + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat4)
summary(LMM_globalpRNFL_SP_SPCat4sex)
anova(LMM_globalpRNFL_SP_SPCat4sex, LMM_globalpRNFL_SP_SPCat4_new)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_SP_SPCat4both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                      data=data_SPMS_SPCat4)
summary(LMM_globalpRNFL_SP_SPCat4both)
anova(LMM_globalpRNFL_SP_SPCat4both, LMM_globalpRNFL_SP_SPCat4_new)
#adding both sex and age doesn't improve model fit

##################     FIFTH INTERVAL 25,6-30,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat5 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==4)
data_SPMS_SPCat5$Geschlecht<-as.factor(data_SPMS_SPCat5$Geschlecht)

LMM_globalpRNFL_SP_SPCat5 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_SPMS_SPCat5)
summary(LMM_globalpRNFL_SP_SPCat5)
confint(LMM_globalpRNFL_SP_SPCat5)

#a: age
LMM_globalpRNFL_SP_SPCat5age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat5)
summary(LMM_globalpRNFL_SP_SPCat5age)
#convergence error
#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            29.6238  5.4428  
#ID        (Intercept)             5.8610  2.4209  
#ID.1      Disease_Duration_FINAL  0.2479  0.4979  
#Residual                          2.7699  1.6643

#no RE appears to explain too little variance --> adding age results in an overly complex model

#b: sex
LMM_globalpRNFL_SP_SPCat5sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat5)
summary(LMM_globalpRNFL_SP_SPCat5sex)
anova(LMM_globalpRNFL_SP_SPCat5sex, LMM_globalpRNFL_SP_SPCat5)
#adding sex doesn't improve model fit

#c: both
LMM_globalpRNFL_SP_SPCat5both <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=data_SPMS_SPCat5)
summary(LMM_globalpRNFL_SP_SPCat5both)
#convergence error

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            2.883e+01 5.3697016
#ID        (Intercept)            1.145e-07 0.0003384
#ID.1      Disease_Duration_FINAL 1.962e-01 0.4428919
#Residual                         2.628e+00 1.6211456

#Random Intercept per subject explains almost no variance --> eliminate this RE to reduce model complexity

LMM_globalpRNFL_SP_SPCat5both_1 <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=data_SPMS_SPCat5)
summary(LMM_globalpRNFL_SP_SPCat5both_1)
anova(LMM_globalpRNFL_SP_SPCat5both_1, LMM_globalpRNFL_SP_SPCat5)
#no convergence error anymore and much better fit than basic model without covariates and additional random intercept for ID
confint(LMM_globalpRNFL_SP_SPCat5both_1)
isSingular(LMM_globalpRNFL_SP_SPCat5both_1)
#no singular fit
#best model fit - presented in paper

##################        SIXTH INTERVAL 30,6+ YEARS DISEASE DURATION        ############
data_SPMS_SPCat6 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==5)
data_SPMS_SPCat6$Geschlecht<-as.factor(data_SPMS_SPCat6$Geschlecht)

LMM_globalpRNFL_SP_SPCat6 <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                  data=data_SPMS_SPCat6)
summary(LMM_globalpRNFL_SP_SPCat6)
#Warning messages:
#1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#unable to evaluate scaled gradient
#2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
#3: Model failed to converge with 1 negative eigenvalue: -2.9e-01

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            1.578e+02 12.563316
#ID        (Intercept)            2.826e-06  0.001681
#ID.1      Disease_Duration_FINAL 6.275e-03  0.079217
#Residual                         2.592e+00  1.609894

#Random Intercept per subject explains almost no variance --> eliminate this RE to reduce model complexity

LMM_globalpRNFL_SP_SPCat6_new <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=data_SPMS_SPCat6)
summary(LMM_globalpRNFL_SP_SPCat6_new)
confint(LMM_globalpRNFL_SP_SPCat6_new)
#best model fit - presented in paper

#a: age
LMM_globalpRNFL_SP_SPCat6age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat6)
summary(LMM_globalpRNFL_SP_SPCat6age)
anova(LMM_globalpRNFL_SP_SPCat6age, LMM_globalpRNFL_SP_SPCat6_new)
#adding age doesn't improve model fit

#b: sex
LMM_globalpRNFL_SP_SPCat6sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=data_SPMS_SPCat6)
summary(LMM_globalpRNFL_SP_SPCat6sex)
##optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_globalpRNFL_SP_SPCat6sex)
#true
#Random effects:
#Groups   Name                   Variance Std.Dev.
#ID_Eye   (Intercept)            162.970  12.766  
#ID       Disease_Duration_FINAL   0.000   0.000  
#Residual                          2.605   1.614

#RE disease duration doesn't explain any variance but sex also doesn't seem to add value to the model
#check, if elimination of RE of disease duration improves model fit if sex is included
LMM_globalpRNFL_SP_SPCat6sex_1 <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID_Eye), 
                                     data=data_SPMS_SPCat6)
summary(LMM_globalpRNFL_SP_SPCat6sex_1)
#Fixed effects:
#Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)            80.35710    4.33945 71.70785  18.518   <2e-16 ***
#Disease_Duration_FINAL -0.13630    0.09181 43.56921  -1.484    0.145    
#Geschlecht1             6.64368    4.87042 31.22278   1.364    0.182

summary(LMM_globalpRNFL_SP_SPCat6_new)
#Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)            82.53812    4.07301 41.58525  20.265   <2e-16 ***
#Disease_Duration_FINAL -0.14275    0.09771  6.39667  -1.461    0.191 

#results of both models are very similar
#theoretically, the model including a RE of disease duration is more appropriate --> stick to this one 

#c: both
LMM_globalpRNFL_SP_SPCat6both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age+ (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=data_SPMS_SPCat6)
summary(LMM_globalpRNFL_SP_SPCat6both)
anova(LMM_globalpRNFL_SP_SPCat6both, LMM_globalpRNFL_SP_SPCat6_new)
#adding both sex and age doesn't improve model fit 

#########################################################################################              
#########################################################################################              
#                                        mRNFL
#########################################################################################              
######################################################################################### 

######################################################################################### 
###########################                HC           #################################
######################################################################################### 

###    FIRST (ONLY) INTERVAL 0-3,5 YEARS "DISEASE DURATION" (=Time since baseline)    ### 
data_HC_Cat1 <- subset(Data_LMM_HC, DiseaseDurationPPMSCategories_LMM==0)
data_HC_Cat1$Verlaufsform<-as.factor(data_HC_Cat1$Verlaufsform)

LMM_RNFL_HC_Cat1 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_HC_Cat1)
summary(LMM_RNFL_HC_Cat1)
confint(LMM_RNFL_HC_Cat1)
#best model fit - presented in paper

#Add Covariates
#a: Age
LMM_RNFL_HC_Cat1age <- lmer(RNFL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_HC_Cat1)
summary(LMM_RNFL_HC_Cat1age)
anova(LMM_RNFL_HC_Cat1age, LMM_RNFL_HC_Cat1)
#no sign. effect

#b: Sex
data_HC_Cat1$Geschlecht<-as.factor(data_HC_Cat1$Geschlecht)
LMM_RNFL_HC_Cat1sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_HC_Cat1)
summary(LMM_RNFL_HC_Cat1sex)
anova(LMM_RNFL_HC_Cat1sex, LMM_RNFL_HC_Cat1)
#no sign. effect

#c: Both
LMM_RNFL_HC_Cat1agesex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_HC_Cat1)
summary(LMM_RNFL_HC_Cat1agesex)
anova(LMM_RNFL_HC_Cat1agesex, LMM_RNFL_HC_Cat1)
#no sign. effect

######################################################################################### 
###########################               RRMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 0-3,5 YEARS DISEASE DURATION      #################
data_RRMS_Cat1 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==0)
data_RRMS_Cat1$Geschlecht<-as.factor(data_RRMS_Cat1$Geschlecht)

LMM_RNFL_RR_Cat1 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat1)
summary(LMM_RNFL_RR_Cat1)
model_performance(LMM_RNFL_RR_Cat1)
confint(LMM_RNFL_RR_Cat1)
#best model fit - presented in paper

#add age as covariate
LMM_RNFL_RR_Cat1age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat1)
summary(LMM_RNFL_RR_Cat1age)
anova(LMM_RNFL_RR_Cat1, LMM_RNFL_RR_Cat1age)
#age doesn't improve model fit

#b: sex
LMM_RNFL_RR1sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_RRMS_Cat1)
summary(LMM_RNFL_RR1sex)
anova(LMM_RNFL_RR1sex, LMM_RNFL_RR_Cat1)
#adding sex doesn't improve model fit

#c: both
LMM_RNFL_RR1both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat1)
summary(LMM_RNFL_RR1both)
anova(LMM_RNFL_RR1both, LMM_RNFL_RR_Cat1)
#adding both sex and age doesn't improve model fit

##################      SECOND INTERVAL 3,6-5,5 YEARS DISEASE DURATION   ################
data_RRMS_Cat2 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==1)
data_RRMS_Cat2$Geschlecht<-as.factor(data_RRMS_Cat2$Geschlecht)

LMM_RNFL_RR_Cat2 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat2)
summary(LMM_RNFL_RR_Cat2)
#convergence error
model_performance(LMM_RNFL_RR_Cat2)

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            10.34940 3.2170  
#ID        (Intercept)            12.59023 3.5483  
#ID.1      Disease_Duration_FINAL  0.02451 0.1566  
#Residual                          0.67489 0.8215  

#Try without RE of DD as this RE explains least variance
LMM_RNFL_RR_Cat2A <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye),
                          data=data_RRMS_Cat2)
summary(LMM_RNFL_RR_Cat2A)
confint(LMM_RNFL_RR_Cat2A)
#model converges
#best model fit - presented in paper

#add covariates
#age
LMM_RNFL_RR_Cat2age <- lmer(RNFL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye), 
                            data=data_RRMS_Cat2)
summary(LMM_RNFL_RR_Cat2age)
anova(LMM_RNFL_RR_Cat2age,LMM_RNFL_RR_Cat2A)
#age doesn't improve model fit

#add sex
LMM_RNFL_RR_Cat2sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye), 
                            data=data_RRMS_Cat2)
summary(LMM_RNFL_RR_Cat2sex)
anova(LMM_RNFL_RR_Cat2sex,LMM_RNFL_RR_Cat2A)
#adding sex doesn't improve model fit

#add both
LMM_RNFL_RR_Cat2both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht+ Age +(1|ID/ID_Eye), 
                             data=data_RRMS_Cat2)
summary(LMM_RNFL_RR_Cat2both)
anova(LMM_RNFL_RR_Cat2both,LMM_RNFL_RR_Cat2A)
#adding both sex and age doesn't improve model fit

##################     THIRD INTERVAL 5,6-7,5 YEARS DISEASE DURATION    ################# 
data_RRMS_Cat3 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==2)
data_RRMS_Cat3$Geschlecht<-as.factor(data_RRMS_Cat3$Geschlecht)

LMM_RNFL_RR_Cat3 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat3)
summary(LMM_RNFL_RR_Cat3)
plot(LMM_RNFL_RR_Cat3)
isSingular(LMM_RNFL_RR_Cat3)
#no singular fit

#add age  
LMM_RNFL_RR_Cat3age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat3)
summary(LMM_RNFL_RR_Cat3age)
confint(LMM_RNFL_RR_Cat3age)
anova(LMM_RNFL_RR_Cat3, LMM_RNFL_RR_Cat3age)
#better fit than without age
#best model fit - presented in paper

#add sex
LMM_RNFL_RR_Cat3sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat3)
summary(LMM_RNFL_RR_Cat3sex)
anova(LMM_RNFL_RR_Cat3sex, LMM_RNFL_RR_Cat3)
#sex doesn't improve model fit

#sex+age
LMM_RNFL_RR_Cat3both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age+ (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat3)
summary(LMM_RNFL_RR_Cat3both)
anova(LMM_RNFL_RR_Cat3both, LMM_RNFL_RR_Cat3age)
#sex doesn't further improve model compared to model including age only

##################      FOURTH INTERVAL 7,6-10,5 YEARS DISEASE DURATION   ###############
data_RRMS_Cat4 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==3)
data_RRMS_Cat4$Geschlecht<-as.factor(data_RRMS_Cat4$Geschlecht)

LMM_RNFL_RRMS_Cat4 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_RRMS_Cat4)
summary(LMM_RNFL_RRMS_Cat4)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_RRMS_Cat4)
#true 

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             3.913   1.978   
#ID        (Intercept)            11.095   3.331   
#ID.1      Disease_Duration_FINAL  0.000   0.000   
#Residual                          1.605   1.267 

#no variance explained by random effect of disease duration --> Try without this RE
LMM_RNFL_RRMS_Cat4_2 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                             data=data_RRMS_Cat4)
summary(LMM_RNFL_RRMS_Cat4_2)
confint(LMM_RNFL_RRMS_Cat4_2)
#model converges
#best model fit - presented in paper

#a:age
LMM_RNFL_RRMS_Cat4age <- lmer(RNFL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye), 
                              data=data_RRMS_Cat4)
summary(LMM_RNFL_RRMS_Cat4age)
anova(LMM_RNFL_RRMS_Cat4age, LMM_RNFL_RRMS_Cat4_2)
#model including age doesn't improve model fit

#b: sex
LMM_RNFL_RRMS_Cat4sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye), 
                              data=data_RRMS_Cat4)
summary(LMM_RNFL_RRMS_Cat4sex)
anova(LMM_RNFL_RRMS_Cat4sex, LMM_RNFL_RRMS_Cat4_2)
#model including sex doesn't improve model fit

#b: both
LMM_RNFL_RRMS_Cat4sexage <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                                 data=data_RRMS_Cat4)
summary(LMM_RNFL_RRMS_Cat4sexage)
anova(LMM_RNFL_RRMS_Cat4_2, LMM_RNFL_RRMS_Cat4sexage)
#including both sex and age doesn't improve model fit

#################      FIFTH INTERVAL 10,6-13,5 YEARS DISEASE DURATION   ################ 
data_RRMS_Cat5 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==4)
data_RRMS_Cat5$Geschlecht<-as.factor(data_RRMS_Cat5$Geschlecht)

LMM_RNFL_RR_Cat5 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat5)
summary(LMM_RNFL_RR_Cat5)
#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_RR_Cat5)
#true

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            5.6538   2.3778  
#ID        (Intercept)            0.0000   0.0000  
#ID.1      Disease_Duration_FINAL 0.1644   0.4054  
#Residual                         0.7102   0.8427  

#Random intercept of ID doesn't explain any variance
#Try without RE of ID --> only Eyespecific effect but no ID-specific
LMM_RNFL_RR_Cat5_2 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_RRMS_Cat5)
summary(LMM_RNFL_RR_Cat5_2)
confint(LMM_RNFL_RR_Cat5_2)
#best model fit - presented in paper

#a: age
LMM_RNFL_RR_Cat5age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat5)
summary(LMM_RNFL_RR_Cat5age)
anova(LMM_RNFL_RR_Cat5_2, LMM_RNFL_RR_Cat5age)
#adding age doesn't improve model fit

#b: sex
LMM_RNFL_RR_Cat5sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat5)
summary(LMM_RNFL_RR_Cat5sex)
anova(LMM_RNFL_RR_Cat5_2, LMM_RNFL_RR_Cat5sex)
#adding sex doesn't improve model fit

#c: both
LMM_RNFL_RR_Cat5both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat5)
summary(LMM_RNFL_RR_Cat5both)
anova(LMM_RNFL_RR_Cat5both, LMM_RNFL_RR_Cat5_2)
#adding both sex and age doesn't improve model fit

##################     SIXTH INTERVAL 13,6-16,5 YEARS DISEASE DURATION   ################
data_RRMS_Cat6 <- subset(Data_LMM_RRMS, DiseaseDurationSPMSCategories_LMM==1)
data_RRMS_Cat6$Geschlecht<-as.factor(data_RRMS_Cat6$Geschlecht)

LMM_RNFL_RR_Cat6 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat6)
summary(LMM_RNFL_RR_Cat6)
confint(LMM_RNFL_RR_Cat6)
#best model fit - presented in paper

#a: age
LMM_RNFL_RR_Cat6age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat6)
summary(LMM_RNFL_RR_Cat6age)
anova(LMM_RNFL_RR_Cat6age, LMM_RNFL_RR_Cat6)
#age doesn't improve model

#b: sex
LMM_RNFL_RR_Cat6sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat6)
summary(LMM_RNFL_RR_Cat6sex)
anova(LMM_RNFL_RR_Cat6sex, LMM_RNFL_RR_Cat6)
#sex doesn't improve model

#c: both
LMM_RNFL_RR_Cat6both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat6)
summary(LMM_RNFL_RR_Cat6both)
anova(LMM_RNFL_RR_Cat6both, LMM_RNFL_RR_Cat6)
#adding both sex and age doesn't improve model

######################################################################################### 
###########################               PPMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 0-3,5 YEARS DISEASE DURATION      #################
data_PPMS_Cat1 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==0)
data_PPMS_Cat1$Geschlecht<-as.factor(data_PPMS_Cat1$Geschlecht)

LMM_RNFL_PP_Cat1 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat1)
summary(LMM_RNFL_PP_Cat1)
isSingular(LMM_RNFL_PP_Cat1)
confint(LMM_RNFL_PP_Cat1)
#best model fit - presented in paper

#a: age
LMM_RNFL_PP_Cat1age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_Cat1)
summary(LMM_RNFL_PP_Cat1age)
anova(LMM_RNFL_PP_Cat1age, LMM_RNFL_PP_Cat1)
#age doesn't improve model fit

#b: sex
LMM_RNFL_RR1sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat1)
summary(LMM_RNFL_RR1sex)
anova(LMM_RNFL_PP_Cat1, LMM_RNFL_RR1sex)
#sex doesn't improve model fit

#c: both
LMM_RNFL_RR1both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat1)
summary(LMM_RNFL_RR1both)
anova(LMM_RNFL_PP_Cat1, LMM_RNFL_RR1both)
#adding both sex and age doesn't improve model fit

##################      SECOND INTERVAL 3,6-5,5 YEARS DISEASE DURATION   ################
data_PPMS_Cat2 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==1)
data_PPMS_Cat2$Geschlecht<-as.factor(data_PPMS_Cat2$Geschlecht)

LMM_RNFL_PP_Cat2 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat2)
summary(LMM_RNFL_PP_Cat2)
#boundary (singular) fit: see help('isSingular')

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            1.008e+00 1.004e+00
#ID        (Intercept)            2.198e+01 4.689e+00
#ID.1      Disease_Duration_FINAL 1.546e-09 3.932e-05
#Residual                         1.048e+00 1.023e+00

#RE of disease duration doesn't explain a significant proportion of variance
#try model without this (apparently unnecessary RE)
LMM_RNFL_PP_Cat2_2 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                           data=data_PPMS_Cat2)
summary(LMM_RNFL_PP_Cat2_2)
confint(LMM_RNFL_PP_Cat2_2)
#best model fit - presented in paper

#a: age
LMM_RNFL_PP_Cat2age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat2)
summary(LMM_RNFL_PP_Cat2age)
anova(LMM_RNFL_PP_Cat2age, LMM_RNFL_PP_Cat2_2)
#adding age doesn't improve model fit

#b: sex
LMM_RNFL_RR2sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                        data=data_PPMS_Cat2)
summary(LMM_RNFL_RR2sex)
anova(LMM_RNFL_PP_Cat2_2,LMM_RNFL_RR2sex)
#adding sex doesn't improve model fit

#c: both
LMM_RNFL_RR2both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                         data=data_PPMS_Cat2)
summary(LMM_RNFL_RR2both)
anova(LMM_RNFL_PP_Cat2_2,LMM_RNFL_RR2both)
#adding both sex and age doesn't improve model fit

##################     THIRD INTERVAL 5,6-7,5 YEARS DISEASE DURATION    ################# 
data_PPMS_Cat3 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==2)
data_PPMS_Cat3$Geschlecht<-as.factor(data_PPMS_Cat3$Geschlecht)

LMM_RNFL_PP_Cat3 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat3)
summary(LMM_RNFL_PP_Cat3)
#boundary (singular) fit: see help('isSingular')

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             1.7968  1.3404  
#ID        (Intercept)            18.5337  4.3051  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.4074  0.6382  

#RE of disease duration doesn't explain a significant proportion of variance
#try model without this (apparently unnecessary RE)

LMM_RNFL_PP_Cat3_2 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                           data=data_PPMS_Cat3)
summary(LMM_RNFL_PP_Cat3_2)
confint(LMM_RNFL_PP_Cat3_2)
anova(LMM_RNFL_PP_Cat3_2, LMM_RNFL_PP_Cat3)

#a: age
LMM_RNFL_PP_Cat3age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat3)
summary(LMM_RNFL_PP_Cat3age)
anova(LMM_RNFL_PP_Cat3age, LMM_RNFL_PP_Cat3_2)
#age doesn't improve model fit

#b: sex
LMM_RNFL_PP_Cat3sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat3)
summary(LMM_RNFL_PP_Cat3sex)
anova(LMM_RNFL_PP_Cat3sex, LMM_RNFL_PP_Cat3_2)
#better fit than without sex
confint(LMM_RNFL_PP_Cat3sex)
#best model fit - presented in paper

#c: both
LMM_RNFL_PP_Cat3sexage <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                               data=data_PPMS_Cat3)
summary(LMM_RNFL_PP_Cat3sexage)
anova(LMM_RNFL_PP_Cat3sexage, LMM_RNFL_PP_Cat3sex)
#adding age in addition to sex doesn't further improve model fit

##################      FOURTH INTERVAL 7,6-10,5 YEARS DISEASE DURATION   ###############
data_PPMS_Cat4 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==3)
data_PPMS_Cat4$Geschlecht<-as.factor(data_PPMS_Cat4$Geschlecht)

LMM_RNFL_PP_Cat4 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat4)
summary(LMM_RNFL_PP_Cat4)
confint(LMM_RNFL_PP_Cat4)
#best model fit - presented in paper

#a: age
LMM_RNFL_PP_Cat4age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_Cat4)
summary(LMM_RNFL_PP_Cat4age)
anova(LMM_RNFL_PP_Cat4age, LMM_RNFL_PP_Cat4)
#adding age doesn't improve model fit

#b: sex
LMM_RNFL_PP_Cat4sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_Cat4)
summary(LMM_RNFL_PP_Cat4sex)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_PP_Cat4sex)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             1.469   1.212   
#ID        (Intercept)            22.706   4.765   
#ID.1      Disease_Duration_FINAL  0.000   0.000   
#Residual                          2.167   1.472  

#RE disease duration doesn't explain any variance
#check, if elimination of RE of disease duration improves model fit if sex is included

LMM_RNFL_PP_Cat4sex_1 <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat4)
summary(LMM_RNFL_PP_Cat4sex_1)
#Fixed effects:
#Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)             29.5230     2.6701 57.3679  11.057 7.63e-16 ***
#Disease_Duration_FINAL   0.2894     0.2662 38.2750   1.087    0.284    
#Geschlecht1             -2.6043     1.5871 40.9837  -1.641    0.108    

anova(LMM_RNFL_PP_Cat4, LMM_RNFL_PP_Cat4sex_1)

summary(LMM_RNFL_PP_Cat4)
#Fixed effects:
#Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)             27.9532     2.4877 44.4405   11.24 1.42e-14 ***
#Disease_Duration_FINAL   0.2751     0.2672 37.6861    1.03     0.31 

#results of both models are very similar 
#theoretically, model including RE is more appropriate --> stick with this one and 
#don't include sex as a covariate

#################      FIFTH INTERVAL 10,6-13,5 YEARS DISEASE DURATION   ################ 
data_PPMS_Cat5 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==4)
data_PPMS_Cat5$Geschlecht<-as.factor(data_PPMS_Cat5$Geschlecht)

LMM_RNFL_PP_Cat5 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat5)
summary(LMM_RNFL_PP_Cat5)
confint(LMM_RNFL_PP_Cat5)
#best model fit - presented in paper

#a: age
LMM_RNFL_PP_Cat5age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_Cat5)
summary(LMM_RNFL_PP_Cat5age)
anova(LMM_RNFL_PP_Cat5age, LMM_RNFL_PP_Cat5)
#age doesn't sign. improve model fit

#b: sex
LMM_RNFL_PP_Cat5sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_Cat5)
summary(LMM_RNFL_PP_Cat5sex)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_PP_Cat5sex)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             2.1327  1.4604  
#ID        (Intercept)            19.3473  4.3986  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.7693  0.8771 

#RE disease duration doesn't explain any variance
#check, if elimination of RE of disease duration improves model fit if sex is included

LMM_RNFL_PP_Cat5sex_1 <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat5)
summary(LMM_RNFL_PP_Cat5sex_1)
#Fixed effects:
#Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)             37.3969     2.1131 57.4557  17.698  < 2e-16 ***
#Disease_Duration_FINAL  -0.5288     0.1496 33.1869  -3.535  0.00123 ** 
#Geschlecht1             -2.2573     1.5316 35.2488  -1.474  0.14939  

summary(LMM_RNFL_PP_Cat5)
#Fixed effects:
#Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)             36.0172     1.9227 45.1701  18.732  < 2e-16 ***
#Disease_Duration_FINAL  -0.5249     0.1522 25.6678  -3.449  0.00196 ** 

#results of both models are very similar 
#theoretically, model including RE is more appropriate --> stick with this one and 
#don't include sex as a covariate

##################     SIXTH INTERVAL 13,6-16,5 YEARS DISEASE DURATION   ################
data_PPMS_Cat6 <- subset(Data_LMM_PPMS, DiseaseDurationSPMSCategories_LMM==1)
data_PPMS_Cat6$Geschlecht<-as.factor(data_PPMS_Cat6$Geschlecht)

LMM_RNFL_PP_Cat6 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat6)
summary(LMM_RNFL_PP_Cat6)
isSingular(LMM_RNFL_PP_Cat6)
#TRUE
#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             3.1330  1.7700  
#ID        (Intercept)            26.1311  5.1119  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.3245  0.5697 

#RE disease duration doesn't explain any variance --> eliminate this 
LMM_RNFL_PP_Cat6_new <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                             data=data_PPMS_Cat6)
summary(LMM_RNFL_PP_Cat6_new)
confint(LMM_RNFL_PP_Cat6_new)

#a: age
LMM_RNFL_PP_Cat6age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat6)
summary(LMM_RNFL_PP_Cat6age)
anova(LMM_RNFL_PP_Cat6age, LMM_RNFL_PP_Cat6_new)
#age doesn't improve model fit

#b: sex
LMM_RNFL_PP_Cat6sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat6)
summary(LMM_RNFL_PP_Cat6sex)
anova(LMM_RNFL_PP_Cat6sex, LMM_RNFL_PP_Cat6_new)
#better fit 
confint(LMM_RNFL_PP_Cat6sex)
#best model fit - presented in paper

#c: both
LMM_RNFL_PP_Cat6both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                             data=data_PPMS_Cat6)
summary(LMM_RNFL_PP_Cat6both)
anova(LMM_RNFL_PP_Cat6both, LMM_RNFL_PP_Cat6sex)
#adding age in addition to sex doesn't further improve model fit

#################      SEVENTH INTERVAL 16,6-20,5 YEARS DISEASE DURATION   ##############
data_PPMS_Cat7 <- subset(Data_LMM_PPMS, DiseaseDurationSPMSCategories_LMM==2)
data_PPMS_Cat7$Geschlecht<-as.factor(data_PPMS_Cat7$Geschlecht)

LMM_RNFL_PP_Cat7 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat7)
summary(LMM_RNFL_PP_Cat7)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_PP_Cat7)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             3.5243  1.8773  
#ID        (Intercept)            31.9254  5.6503  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.4094  0.6399  

#RE of DD doesn't explain any variance --> try model without it
LMM_RNFL_PP_Cat7_new <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                             data=data_PPMS_Cat7)
summary(LMM_RNFL_PP_Cat7_new)
confint(LMM_RNFL_PP_Cat7_new)
#best model fit - presented in paper

#a: age
LMM_RNFL_PP_Cat7age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat7)
summary(LMM_RNFL_PP_Cat7age)
anova(LMM_RNFL_PP_Cat7age, LMM_RNFL_PP_Cat7_new)
#adding age doesn't improve model fit

#b: sex
LMM_RNFL_PP_Cat7sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat7)
summary(LMM_RNFL_PP_Cat7sex)
anova(LMM_RNFL_PP_Cat7sex, LMM_RNFL_PP_Cat7_new)
#adding sex doesn't improve model fit

#b: both
LMM_RNFL_PP_Cat7both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                             data=data_PPMS_Cat7)
summary(LMM_RNFL_PP_Cat7both)
anova(LMM_RNFL_PP_Cat7_new, LMM_RNFL_PP_Cat7both)
#adding both sex and age doesn't improve model fit

######################################################################################### 
###########################               SPMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 3,5-12,5 YEARS DISEASE DURATION      ##############
data_SPMS_SPCat1 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==0)
data_SPMS_SPCat1$Geschlecht<-as.factor(data_SPMS_SPCat1$Geschlecht)

LMM_RNFL_SP_SPCat1 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_SPMS_SPCat1)
summary(LMM_RNFL_SP_SPCat1)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_SP_SPCat1)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             4.5940  2.1434  
#ID        (Intercept)            23.1397  4.8104  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.6904  0.8309

#RE of DD doesn't explain any variance --> try model without it

LMM_RNFL_SP_SPCat1_new <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat1)
summary(LMM_RNFL_SP_SPCat1_new)
confint(LMM_RNFL_SP_SPCat1_new)
#best model fit - presented in paper

#a: age
LMM_RNFL_SP_SPCat1age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat1)
summary(LMM_RNFL_SP_SPCat1age)
anova(LMM_RNFL_SP_SPCat1age, LMM_RNFL_SP_SPCat1_new)
#age doesn't improve model fit

#b: sex
LMM_RNFL_SP_SPCat1sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat1)
summary(LMM_RNFL_SP_SPCat1sex)
anova(LMM_RNFL_SP_SPCat1sex, LMM_RNFL_SP_SPCat1_new)
#sex doesn't improve model fit

#c: both
LMM_RNFL_SP_SPCat1both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat1)
summary(LMM_RNFL_SP_SPCat1both)
anova(LMM_RNFL_SP_SPCat1both, LMM_RNFL_SP_SPCat1_new)
#adding both sex and age doesn't improve model fit

##################     SECOND INTERVAL 12,6-16,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat2 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==1)
data_SPMS_SPCat2$Geschlecht<-as.factor(data_SPMS_SPCat2$Geschlecht)

LMM_RNFL_SP_SPCat2 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_SPMS_SPCat2)
summary(LMM_RNFL_SP_SPCat2)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_SP_SPCat2)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             5.5132  2.3480  
#ID        (Intercept)            18.2691  4.2742  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.4626  0.6801 

#RE of DD doesn't explain any variance --> try model without it

LMM_RNFL_SP_SPCat2_new <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat2)
summary(LMM_RNFL_SP_SPCat2_new)
confint(LMM_RNFL_SP_SPCat2_new)
#best model fit - presented in paper

#a: age
LMM_RNFL_SP_SPCat2age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat2)
summary(LMM_RNFL_SP_SPCat2age)
anova(LMM_RNFL_SP_SPCat2age, LMM_RNFL_SP_SPCat2_new)
#adding age doesn't improve model

#b: sex
LMM_RNFL_RR1sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                        data=data_SPMS_SPCat2)
summary(LMM_RNFL_RR1sex)
anova(LMM_RNFL_RR1sex, LMM_RNFL_SP_SPCat2_new)
#adding sex doesn't improve model fit

#c: both
LMM_RNFL_RR1both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                         data=data_SPMS_SPCat2)
summary(LMM_RNFL_RR1both)
anova(LMM_RNFL_RR1both,LMM_RNFL_SP_SPCat2_new)
#adding both sex and age doesn't improve model fit

##################     THIRD INTERVAL 16,6-20,5 YEARS DISEASE DURATION       ############
data_SPMS_SPCat3 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==2)
data_SPMS_SPCat3$Geschlecht<-as.factor(data_SPMS_SPCat3$Geschlecht)

LMM_RNFL_SP_SPCat3 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_SPMS_SPCat3)
summary(LMM_RNFL_SP_SPCat3)
isSingular(LMM_RNFL_SP_SPCat3)
confint(LMM_RNFL_SP_SPCat3)
#best model fit - presented in paper

#a: age
LMM_RNFL_SP_SPCat3age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                              data=data_SPMS_SPCat3)
summary(LMM_RNFL_SP_SPCat3age)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_SP_SPCat3age)

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            5.74085  2.3960  
#ID        (Intercept)            0.00000  0.0000  
#ID.1      Disease_Duration_FINAL 0.04028  0.2007  
#Residual                         0.66671  0.8165  

#RE of ID Intercept explains 0,0 variance 
#check, if elimination of RE of disease duration improves model fit if age is included

LMM_RNFL_SP_SPCat3age_1 <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                              data=data_SPMS_SPCat3)
summary(LMM_RNFL_SP_SPCat3age_1)
anova(LMM_RNFL_SP_SPCat3age_1, LMM_RNFL_SP_SPCat3)
#doesn't improve model fit --> stick to model without age and with RE of ID

#b: sex
LMM_RNFL_SP_SPCat3sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                              data=data_SPMS_SPCat3)
summary(LMM_RNFL_SP_SPCat3sex)
anova(LMM_RNFL_SP_SPCat3sex, LMM_RNFL_SP_SPCat3)
#sex doesn't improve model fit

##################     FOURTH INTERVAL 20,6-25,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat4 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==3)
data_SPMS_SPCat4$Geschlecht<-as.factor(data_SPMS_SPCat4$Geschlecht)

LMM_RNFL_SP_SPCat4 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_SPMS_SPCat4)
summary(LMM_RNFL_SP_SPCat4)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_RNFL_SP_SPCat4)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             3.2154  1.793   
#ID        (Intercept)            11.7587  3.429   
#ID.1      Disease_Duration_FINAL  0.0000  0.000   
#Residual                          0.4928  0.702 

#RE of disease duration explains 0,0 variance --> try model without it
LMM_RNFL_SP_SPCat4_new <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat4)
summary(LMM_RNFL_SP_SPCat4_new)
confint(LMM_RNFL_SP_SPCat4_new)
#best model fit - presented in paper

#a: age
LMM_RNFL_SP_SPCat4age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat4)
summary(LMM_RNFL_SP_SPCat4age)
anova(LMM_RNFL_SP_SPCat4age, LMM_RNFL_SP_SPCat4_new)
#age doesn't improve model fit

#b: sex
LMM_RNFL_SP_SPCat4sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat4)
summary(LMM_RNFL_SP_SPCat4sex)
anova(LMM_RNFL_SP_SPCat4sex, LMM_RNFL_SP_SPCat4_new)
#sex doesn't improve model fit

#c: both
LMM_RNFL_SP_SPCat4both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat4)
summary(LMM_RNFL_SP_SPCat4both)
anova(LMM_RNFL_SP_SPCat4sex, LMM_RNFL_SP_SPCat4_new)
#adding both age and sex doesn't improve model fit

##################     FIFTH INTERVAL 25,6-30,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat5 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==4)
data_SPMS_SPCat5$Geschlecht<-as.factor(data_SPMS_SPCat5$Geschlecht)

LMM_RNFL_SP_SPCat5 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_SPMS_SPCat5)
summary(LMM_RNFL_SP_SPCat5)
#Warning message:
#In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#Model failed to converge with max|grad| = 0.00207535 (tol = 0.002, component 1)

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             3.97644 1.9941  
#ID        (Intercept)            12.24498 3.4993  
#ID.1      Disease_Duration_FINAL  0.01359 0.1166  
#Residual                          0.65440 0.8090 

#RE of disease duration doesn't explain a lot of variance --> try model without it

LMM_RNFL_SP_SPCat5_new <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat5)
summary(LMM_RNFL_SP_SPCat5_new)
confint(LMM_RNFL_SP_SPCat5_new)

#a: age
LMM_RNFL_SP_SPCat5age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat5)
summary(LMM_RNFL_SP_SPCat5age)
anova(LMM_RNFL_SP_SPCat5age, LMM_RNFL_SP_SPCat5_new)
#age doesn't improve model fit

#b: sex
LMM_RNFL_SP_SPCat5sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat5)
summary(LMM_RNFL_SP_SPCat5sex)
anova(LMM_RNFL_SP_SPCat5sex, LMM_RNFL_SP_SPCat5_new)
#sex doesn't improve model fit

#c: both
LMM_RNFL_SP_SPCat5both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat5)
summary(LMM_RNFL_SP_SPCat5both)
anova(LMM_RNFL_SP_SPCat5both, LMM_RNFL_SP_SPCat5_new)
#model sign. improved
isSingular(LMM_RNFL_SP_SPCat5both)
#no singular fit
confint(LMM_RNFL_SP_SPCat5both)
#best model fit - presented in paper

##################        SIXTH INTERVAL 30,6+ YEARS DISEASE DURATION        ############
data_SPMS_SPCat6 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==5)
data_SPMS_SPCat6$Geschlecht<-as.factor(data_SPMS_SPCat6$Geschlecht)

LMM_RNFL_SP_SPCat6 <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_SPMS_SPCat6)
summary(LMM_RNFL_SP_SPCat6)
confint(LMM_RNFL_SP_SPCat6)

#There were 50 or more warnings (use warnings() to see the first 50)
warnings()
#1: In zetafun(np, ns) : slightly lower deviances (diff=-7.31575e-11) detected
#2: In nextpar(mat, cc, i, delta, lowcut, upcut) :
#  Last two rows have identical or NA .zeta values: using minstep
#3: In zetafun(np, ns) : slightly lower deviances (diff=-9.65201e-11) detected
#4: In nextpar(mat, cc, i, delta, lowcut, upcut) :
#  Last two rows have identical or NA .zeta values: using minstep

isSingular(LMM_RNFL_SP_SPCat6)
#model not singular

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            1.586e+01 3.983e+00
#ID        (Intercept)            5.742e-09 7.577e-05
#ID.1      Disease_Duration_FINAL 4.364e-03 6.606e-02
#Residual                         5.145e-01 7.173e-01

#RE explain litte variance; especially random intercept of ID doesn't seem to add any value --> eliminate
LMM_RNFL_SP_SPCat6_new <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_SPMS_SPCat6)
summary(LMM_RNFL_SP_SPCat6_new)
confint(LMM_RNFL_SP_SPCat6_new)
#best model fit - presented in paper

 #a: age
LMM_RNFL_SP_SPCat6age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID),  
                              data=data_SPMS_SPCat6)
summary(LMM_RNFL_SP_SPCat6age)
anova(LMM_RNFL_SP_SPCat6age, LMM_RNFL_SP_SPCat6_new)
#age doesn't improve model fit

#b: sex
LMM_RNFL_SP_SPCat6sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                              data=data_SPMS_SPCat6)
summary(LMM_RNFL_SP_SPCat6sex)
anova(LMM_RNFL_SP_SPCat6sex, LMM_RNFL_SP_SPCat6_new)
#sex doesn't improve model fit

#c: both
LMM_RNFL_SP_SPCat6both <- lmer(RNFL_new~Disease_Duration_FINAL + Age + Geschlecht + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                               data=data_SPMS_SPCat6)
summary(LMM_RNFL_SP_SPCat6both)
anova(LMM_RNFL_SP_SPCat6both, LMM_RNFL_SP_SPCat6_new)
#adding both age and sex doesn't improve model fit

#########################################################################################              
#########################################################################################              
#                                        GCIPL
#########################################################################################              
######################################################################################### 

######################################################################################### 
###########################                HC           #################################
######################################################################################### 

###    FIRST (ONLY) INTERVAL 0-3,5 YEARS "DISEASE DURATION" (=Time since baseline)    ### 
data_HC_Cat1 <- subset(Data_LMM_HC, DiseaseDurationPPMSCategories_LMM==0)
data_HC_Cat1$Verlaufsform<-as.factor(data_HC_Cat1$Verlaufsform)

LMM_GCIPL_HC_Cat1_2 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_HC_Cat1)
summary(LMM_GCIPL_HC_Cat1_2)
confint(LMM_GCIPL_HC_Cat1_2)

#Add Covariates
#a: Age
LMM_GCIPL_HC_Cat1age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_HC_Cat1)
summary(LMM_GCIPL_HC_Cat1age)
anova(LMM_GCIPL_HC_Cat1age, LMM_GCIPL_HC_Cat1_2)
#model including age sign. better than model without age
#after controlling for age. time since baseline doesnt have a sign. effect
confint(LMM_GCIPL_HC_Cat1age)
#best model fit - presented in paper

#b: Sex
data_HC_Cat1$Geschlecht<-as.factor(data_HC_Cat1$Geschlecht)
LMM_GCIPL_HC_Cat1sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_HC_Cat1)
summary(LMM_GCIPL_HC_Cat1sex)
anova(LMM_GCIPL_HC_Cat1sex, LMM_GCIPL_HC_Cat1_2)
#sex doesn't improve the model fit

#c: Both
LMM_GCIPL_HC_Cat1agesex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_HC_Cat1)
summary(LMM_GCIPL_HC_Cat1agesex)
anova(LMM_GCIPL_HC_Cat1agesex, LMM_GCIPL_HC_Cat1age)
#sex not needed in addition to age

######################################################################################### 
###########################               RRMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 0-3,5 YEARS DISEASE DURATION      #################
data_RRMS_Cat1 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==0)
data_RRMS_Cat1$Geschlecht<-as.factor(data_RRMS_Cat1$Geschlecht)

LMM_GCIPL_RR_Cat1 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat1)
summary(LMM_GCIPL_RR_Cat1)
confint(LMM_GCIPL_RR_Cat1)
#best model fit - presented in paper

#add age as covariate
LMM_GCIPL_RR_Cat1age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat1)
summary(LMM_GCIPL_RR_Cat1age)
anova(LMM_GCIPL_RR_Cat1, LMM_GCIPL_RR_Cat1age)
#age doesn't improve model fit

#b: sex
LMM_GCIPL_RR1sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_RRMS_Cat1)
summary(LMM_GCIPL_RR1sex)
anova(LMM_GCIPL_RR_Cat1, LMM_GCIPL_RR1sex)
#sex doesn't improve model fit

#c: both
LMM_GCIPL_RR1both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat1)
summary(LMM_GCIPL_RR1both)
anova(LMM_GCIPL_RR1both, LMM_GCIPL_RR_Cat1)
#adding both sex and age doesn't improve model fit

##################      SECOND INTERVAL 3,6-5,5 YEARS DISEASE DURATION   ################
data_RRMS_Cat2 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==1)
data_RRMS_Cat2$Geschlecht<-as.factor(data_RRMS_Cat2$Geschlecht)

LMM_GCIPL_RR_Cat2 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_RRMS_Cat2)
summary(LMM_GCIPL_RR_Cat2)
confint(LMM_GCIPL_RR_Cat2)
#best model fit - presented in paper

#age
LMM_GCIPL_RR_Cat2age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat2)
summary(LMM_GCIPL_RR_Cat2age)
anova(LMM_GCIPL_RR_Cat2age,LMM_GCIPL_RR_Cat2)
#adding age doesn't improve model fit

#add sex
LMM_GCIPL_RR_Cat2sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat2)
summary(LMM_GCIPL_RR_Cat2sex)
anova(LMM_GCIPL_RR_Cat2sex,LMM_GCIPL_RR_Cat2)
#adding sex doesn't improve model fit

#add both
LMM_GCIPL_RR_Cat2both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht+ Age +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat2)
summary(LMM_GCIPL_RR_Cat2both)
anova(LMM_GCIPL_RR_Cat2both,LMM_GCIPL_RR_Cat2)
#adding both sex and age doesn't improve model fit

##################     THIRD INTERVAL 5,6-7,5 YEARS DISEASE DURATION    ################# 
data_RRMS_Cat3 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==2)
data_RRMS_Cat3$Geschlecht<-as.factor(data_RRMS_Cat3$Geschlecht)

LMM_GCIPL_RR_Cat3 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_RRMS_Cat3)
summary(LMM_GCIPL_RR_Cat3)
#optimizer (nloptwrap) convergence code: 0 (OK)
#unable to evaluate scaled gradient
#Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            1.701e+01 4.1246757
#ID        (Intercept)            5.946e-07 0.0007711
#ID.1      Disease_Duration_FINAL 8.800e-01 0.9380871
#Residual                         4.533e-01 0.6732437

#Random intercept of ID explains ~0 variance --> eliminate this unnecessary RE
#Try without RE of ID  
LMM_GCIPL_RR_Cat3A <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_RRMS_Cat3)
summary(LMM_GCIPL_RR_Cat3A)
confint(LMM_GCIPL_RR_Cat3A)
#best model fit - presented in paper

#add covariates to best fitting model
#a. age
LMM_GCIPL_RR_Cat3ageA <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                              data=data_RRMS_Cat3)
summary(LMM_GCIPL_RR_Cat3ageA)
anova(LMM_GCIPL_RR_Cat3A, LMM_GCIPL_RR_Cat3ageA)
#age doesn't improve model fit

#sex
LMM_GCIPL_RR_Cat3sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat3)
summary(LMM_GCIPL_RR_Cat3sex)
anova(LMM_GCIPL_RR_Cat3sex, LMM_GCIPL_RR_Cat3A)
#adding sex doesn't improve model fit

#add sex and age to best fitting model 
LMM_GCIPL_RR_Cat3both <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + Geschlecht + (1|ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat3)
summary(LMM_GCIPL_RR_Cat3both)
anova(LMM_GCIPL_RR_Cat3both, LMM_GCIPL_RR_Cat3A)
#adding both sex and age doesn't improve model fit

##################      FOURTH INTERVAL 7,6-10,5 YEARS DISEASE DURATION   ###############
data_RRMS_Cat4 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==3)
data_RRMS_Cat4$Geschlecht<-as.factor(data_RRMS_Cat4$Geschlecht)

LMM_GCIPL_RRMS_Cat4 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_Cat4)
summary(LMM_GCIPL_RRMS_Cat4)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_GCIPL_RRMS_Cat4)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            13.4761  3.6710  
#ID        (Intercept)            39.0931  6.2524  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.2078  0.4559 

#0 variance of RE of DD --> Try without RE of DD
LMM_GCIPL_RRMS_Cat4_2 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                              data=data_RRMS_Cat4)
summary(LMM_GCIPL_RRMS_Cat4_2)
confint(LMM_GCIPL_RRMS_Cat4_2)

#add covariates to best fitting model
#a:age
LMM_GCIPL_RRMS_Cat4age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye), 
                               data=data_RRMS_Cat4)
summary(LMM_GCIPL_RRMS_Cat4age)
anova(LMM_GCIPL_RRMS_Cat4age, LMM_GCIPL_RRMS_Cat4_2)
#age doesn't improve model fit

#b: sex
LMM_GCIPL_RRMS_Cat4sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye), 
                               data=data_RRMS_Cat4)
summary(LMM_GCIPL_RRMS_Cat4sex)
anova(LMM_GCIPL_RRMS_Cat4sex, LMM_GCIPL_RRMS_Cat4_2)
#sex doesn't improve model fit

#b: both
LMM_GCIPL_RRMS_Cat4sexage <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                                  data=data_RRMS_Cat4)
summary(LMM_GCIPL_RRMS_Cat4sexage)
confint(LMM_GCIPL_RRMS_Cat4sexage)
anova(LMM_GCIPL_RRMS_Cat4sexage, LMM_GCIPL_RRMS_Cat4_2)
#better than model without covariates
#best model fit - presented in paper

#################      FIFTH INTERVAL 10,6-13,5 YEARS DISEASE DURATION   ################ 
data_RRMS_Cat5 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==4)
data_RRMS_Cat5$Geschlecht<-as.factor(data_RRMS_Cat5$Geschlecht)

LMM_GCIPL_RR_Cat5 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_RRMS_Cat5)
summary(LMM_GCIPL_RR_Cat5)
confint(LMM_GCIPL_RR_Cat5)
#best model fit - presented in paper

#a: age
LMM_GCIPL_RR_Cat5age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat5)
summary(LMM_GCIPL_RR_Cat5age)
anova(LMM_GCIPL_RR_Cat5, LMM_GCIPL_RR_Cat5age)
#age doesn't improve model fit

#b: sex
LMM_GCIPL_RR_Cat5sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat5)
summary(LMM_GCIPL_RR_Cat5sex)
anova(LMM_GCIPL_RR_Cat5, LMM_GCIPL_RR_Cat5sex)
#sex doesn't improve model fit

#c: both
LMM_GCIPL_RR_Cat5both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_RRMS_Cat5)
summary(LMM_GCIPL_RR_Cat5both)
#model failed to converge

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            1.399e+01 3.7400105
#ID        (Intercept)            4.677e-07 0.0006839
#ID.1      Disease_Duration_FINAL 6.003e-01 0.7747974
#Residual                         1.302e+00 1.1409059

#check if eliminating random intercept for ID (almost 0 variance) but including sex and age 
#sign. improves model compared to model with random intercept for ID but not sex and age

LMM_GCIPL_RR_Cat5both_2 <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_RRMS_Cat5)
summary(LMM_GCIPL_RR_Cat5both_2)
model_performance(LMM_GCIPL_RR_Cat5both_2)
#Fixed effects:
#Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)             65.4493     8.9882 28.0739   7.282 6.19e-08 ***
#Disease_Duration_FINAL  -1.1990     0.4194 21.7782  -2.859  0.00919 ** 
#Geschlecht1              2.2525     4.2797 23.0613   0.526  0.60368    
#Age                      0.2463     0.2279 24.7467   1.081  0.29032

#AIC     |    AICc |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
#302.464 | 305.196 | 315.707 |      0.988 |      0.056 | 0.987 | 0.560 | 1.141

summary(LMM_GCIPL_RR_Cat5)
model_performance(LMM_GCIPL_RR_Cat5)
#Fixed effects:
#Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)             72.5472     4.3028 18.4475   16.86 1.17e-12 ***
#Disease_Duration_FINAL  -0.8787     0.3755 18.4000   -2.34   0.0307 * 

#AIC     |    AICc |     BIC |  RMSE | Sigma
#304.238 | 306.238 | 315.589 | 0.502 | 1.048

#model results are quite similar, slightly in favor for model including covariates 
#but from a theoretical perspective, including a random intercept per subject is important
#compare model fit
plot(LMM_GCIPL_RR_Cat5both_2)
qqnorm(resid(LMM_GCIPL_RR_Cat5both_2))
qqline(resid(LMM_GCIPL_RR_Cat5both_2))

plot(LMM_GCIPL_RR_Cat5)
qqnorm(resid(LMM_GCIPL_RR_Cat5))
qqline(resid(LMM_GCIPL_RR_Cat5))
#qqplot in favor of LMM_GCIPL_RR_Cat5 --> stick with this model

##################     SIXTH INTERVAL 13,6-16,5 YEARS DISEASE DURATION   ################
data_RRMS_SPCat1 <- subset(Data_LMM_RRMS, DiseaseDurationSPMSCategories_LMM==1)
data_RRMS_SPCat1$Geschlecht<-as.factor(data_RRMS_SPCat1$Geschlecht)

LMM_GCIPL_RR_SPCat1 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_RRMS_SPCat1)
summary(LMM_GCIPL_RR_SPCat1)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_GCIPL_RR_SPCat1)
#TRUE

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             5.9882  2.4471  
#ID        (Intercept)            90.3101  9.5032  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.6467  0.8042

#RE of disease duration doesn't explain any variance --> Try without ID-specific Slope
LMM_GCIPL_RR_SPCat1_new <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                                data=data_RRMS_SPCat1)
summary(LMM_GCIPL_RR_SPCat1_new)
confint(LMM_GCIPL_RR_SPCat1_new)
#best model fit - presented in paper

#a: age as covariate
LMM_GCIPL_RR_SPCat1age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                               data=data_RRMS_SPCat1)
summary(LMM_GCIPL_RR_SPCat1age)
anova(LMM_GCIPL_RR_SPCat1_new, LMM_GCIPL_RR_SPCat1age)
#age doesn't improve model fit

#b: sex
LMM_GCIPL_RR_SPCat1_sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                                data=data_RRMS_SPCat1)
summary(LMM_GCIPL_RR_SPCat1_sex)
anova(LMM_GCIPL_RR_SPCat1_new, LMM_GCIPL_RR_SPCat1_sex)
#adding sex doesn't improve model fit

#c: both
LMM_GCIPL_RR_SPCat1both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                                data=data_RRMS_SPCat1)
summary(LMM_GCIPL_RR_SPCat1both)
anova(LMM_GCIPL_RR_SPCat1_new, LMM_GCIPL_RR_SPCat1both)
#adding both sex and age doesn't improve model fit

######################################################################################### 
###########################               PPMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 0-3,5 YEARS DISEASE DURATION      #################
data_PPMS_Cat1 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==0)
data_PPMS_Cat1$Geschlecht<-as.factor(data_PPMS_Cat1$Geschlecht)

LMM_GCIPL_PP_Cat1 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_Cat1)
summary(LMM_GCIPL_PP_Cat1)
confint(LMM_GCIPL_PP_Cat1)
#best model fit - presented in paper

#a: age
LMM_GCIPL_PP_Cat1age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat1)
summary(LMM_GCIPL_PP_Cat1age)
anova(LMM_GCIPL_PP_Cat1age, LMM_GCIPL_PP_Cat1)
#adding age doesn't improve model fit

#b: sex
LMM_GCIPL_PP_Cat1sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_Cat1)
summary(LMM_GCIPL_PP_Cat1sex)
anova(LMM_GCIPL_PP_Cat1sex, LMM_GCIPL_PP_Cat1)
#adding sex doesn't improve model fit

#c: both
LMM_GCIPL_PP_Cat1both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_Cat1)
summary(LMM_GCIPL_PP_Cat1both)
anova(LMM_GCIPL_PP_Cat1both, LMM_GCIPL_PP_Cat1)
#adding both sex and age doesn't improve model fit

##################      SECOND INTERVAL 3,6-5,5 YEARS DISEASE DURATION   ################
data_PPMS_Cat2 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==1)
data_PPMS_Cat2$Geschlecht<-as.factor(data_PPMS_Cat2$Geschlecht)

LMM_GCIPL_PP_Cat2 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_Cat2)
summary(LMM_GCIPL_PP_Cat2)
confint(LMM_GCIPL_PP_Cat2)
#best model fit - presented in paper

#a: age
LMM_GCIPL_PP_Cat2age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat2)
summary(LMM_GCIPL_PP_Cat2age)
anova(LMM_GCIPL_PP_Cat2age, LMM_GCIPL_PP_Cat2)
#adding age doesn't improve model fit

#b: sex
LMM_GCIPL_PP_Cat2sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat2)
summary(LMM_GCIPL_PP_Cat2sex)
anova(LMM_GCIPL_PP_Cat2sex, LMM_GCIPL_PP_Cat2)
#adding sex doesn't improve model fit

#c: both
LMM_GCIPL_PP_Cat2both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_PPMS_Cat2)
summary(LMM_GCIPL_PP_Cat2both)
anova(LMM_GCIPL_PP_Cat2both, LMM_GCIPL_PP_Cat2)
#adding both sex and age doesn't improve model fit

##################     THIRD INTERVAL 5,6-7,5 YEARS DISEASE DURATION    ################# 
data_PPMS_Cat3 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==2)
data_PPMS_Cat3$Geschlecht<-as.factor(data_PPMS_Cat3$Geschlecht)

LMM_GCIPL_PP_Cat3 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_Cat3)
summary(LMM_GCIPL_PP_Cat3)
confint(LMM_GCIPL_PP_Cat3)
#best model fit - presented in paper

#a: age
LMM_GCIPL_PP_Cat3age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat3)
summary(LMM_GCIPL_PP_Cat3age)
anova(LMM_GCIPL_PP_Cat3age, LMM_GCIPL_PP_Cat3)
#adding age doesn't improve model fit

#b: sex
LMM_GCIPL_PP_Cat3sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat3)
summary(LMM_GCIPL_PP_Cat3sex)
anova(LMM_GCIPL_PP_Cat3sex, LMM_GCIPL_PP_Cat3)
#adding sex doesn't improve model fit

#c: both
LMM_GCIPL_PP_Cat3both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_PPMS_Cat3)
summary(LMM_GCIPL_PP_Cat3both)
anova(LMM_GCIPL_PP_Cat3both, LMM_GCIPL_PP_Cat3)
#adding age and sex doesn't include model fit

##################      FOURTH INTERVAL 7,6-10,5 YEARS DISEASE DURATION   ###############
data_PPMS_Cat4 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==3)
data_PPMS_Cat4$Geschlecht<-as.factor(data_PPMS_Cat4$Geschlecht)

LMM_GCIPL_PP_Cat4 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_Cat4)
summary(LMM_GCIPL_PP_Cat4)
confint(LMM_GCIPL_PP_Cat4)
#best model fit - presented in paper

#a: age
LMM_GCIPL_PP_Cat4age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat4)
summary(LMM_GCIPL_PP_Cat4age)
anova(LMM_GCIPL_PP_Cat4age, LMM_GCIPL_PP_Cat4)
#adding age doesn't improve model fit

#b: sex
LMM_GCIPL_PP_Cat4sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat4)
summary(LMM_GCIPL_PP_Cat4sex)
anova(LMM_GCIPL_PP_Cat4sex, LMM_GCIPL_PP_Cat4)
#adding sex doesn't improve model fit

#c: both
LMM_GCIPL_PP_Cat4both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_PPMS_Cat4)
summary(LMM_GCIPL_PP_Cat4both)
anova(LMM_GCIPL_PP_Cat4both, LMM_GCIPL_PP_Cat4)
#adding both sex and age doesn't improve model fit

#################      FIFTH INTERVAL 10,6-13,5 YEARS DISEASE DURATION   ################ 
data_PPMS_Cat5 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==4)
data_PPMS_Cat5$Geschlecht<-as.factor(data_PPMS_Cat5$Geschlecht)

LMM_GCIPL_PP_Cat5 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_Cat5)
summary(LMM_GCIPL_PP_Cat5)
confint(LMM_GCIPL_PP_Cat5)
#best model fit - presented in paper

#a: age
LMM_GCIPL_PP_Cat5age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat5)
summary(LMM_GCIPL_PP_Cat5age)
anova(LMM_GCIPL_PP_Cat5age, LMM_GCIPL_PP_Cat5)
#age doesn't improve model fit

#b: sex
LMM_GCIPL_PP_Cat5sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_Cat5)
summary(LMM_GCIPL_PP_Cat5sex)
anova(LMM_GCIPL_PP_Cat5sex, LMM_GCIPL_PP_Cat5)
#sex doesn't improve model fit

#c: both
LMM_GCIPL_PP_Cat5both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_PPMS_Cat5)
summary(LMM_GCIPL_PP_Cat5both)
anova(LMM_GCIPL_PP_Cat5both, LMM_GCIPL_PP_Cat5)
#adding both sex and age doesn't improve model fit

##################     SIXTH INTERVAL 13,6-16,5 YEARS DISEASE DURATION   ################
data_PPMS_SPCat1 <- subset(Data_LMM_PPMS, DiseaseDurationSPMSCategories_LMM==1)
data_PPMS_SPCat1$Geschlecht<-as.factor(data_PPMS_SPCat1$Geschlecht)

LMM_GCIPL_PP_SPCat1 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_SPCat1)
summary(LMM_GCIPL_PP_SPCat1)
confint(LMM_GCIPL_PP_SPCat1)
#best model fit - presented in paper

#a: age
LMM_GCIPL_PP_SPCat1age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_PPMS_SPCat1)
summary(LMM_GCIPL_PP_SPCat1age)
anova(LMM_GCIPL_PP_SPCat1age, LMM_GCIPL_PP_SPCat1)
#age doesn't improve model fit

#b: sex
LMM_GCIPL_RR1sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=data_PPMS_SPCat1)
summary(LMM_GCIPL_RR1sex)
anova(LMM_GCIPL_RR1sex,LMM_GCIPL_PP_SPCat1)
#sex doesn't improve model fit

#c: both
LMM_GCIPL_RR1both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_SPCat1)
summary(LMM_GCIPL_RR1both)
anova(LMM_GCIPL_RR1both, LMM_GCIPL_PP_SPCat1)
#adding both sex and age doesn't improve model fit

#################      SEVENTH INTERVAL 16,6-20,5 YEARS DISEASE DURATION   ############## 
data_PPMS_SPCat2 <- subset(Data_LMM_PPMS, DiseaseDurationSPMSCategories_LMM==2)
data_PPMS_SPCat2$Geschlecht<-as.factor(data_PPMS_SPCat2$Geschlecht)

LMM_GCIPL_PP_SPCat2 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_SPCat2)
summary(LMM_GCIPL_PP_SPCat2)
confint(LMM_GCIPL_PP_SPCat2)
model_performance(LMM_GCIPL_PP_SPCat2)

#AIC     |    AICc |     BIC |  RMSE | Sigma
-------------------------------------------
#208.271 | 210.542 | 218.976 | 0.072 | 0.165

#best model fit - presented in paper

#a: age
LMM_GCIPL_PP_SPCat2age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_PPMS_SPCat2)
summary(LMM_GCIPL_PP_SPCat2age)
#optimizer (nloptwrap) convergence code: 0 (OK)
#unable to evaluate scaled gradient
#Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

#Random effects:
#Groups    Name                   Variance  Std.Dev.
#ID_Eye.ID (Intercept)            6.355e+00 2.520984
#ID        (Intercept)            1.145e-05 0.003384
#ID.1      Disease_Duration_FINAL 1.350e-01 0.367431
#Residual                         5.673e-02 0.238190

#adding age results in non-convergence; random intercept per subject explains ~0 variance
#check if model without this RE but with age as a covariate is better than model with RE but without age

LMM_GCIPL_PP_SPCat2age_2 <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_PPMS_SPCat2)
summary(LMM_GCIPL_PP_SPCat2age_2)
model_performance(LMM_GCIPL_PP_SPCat2age_2)
#Fixed effects:
#Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)             81.8775     6.3468 33.2488  12.901 1.71e-14 ***
#Disease_Duration_FINAL  -0.4393     0.1713 32.1386  -2.565   0.0152 *  
#Age                     -0.2022     0.1473 34.4174  -1.373   0.1786  

#AIC     |    AICc |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
#212.862 | 215.132 | 223.567 |      0.999 |      0.058 | 0.999 | 0.117 | 0.238

#based on AIC model LMM_GCIPL_PP_SPCat2 is better --> stick with this

#b: sex
LMM_GCIPL_PP_SPCat2sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_PPMS_SPCat2)
summary(LMM_GCIPL_PP_SPCat2sex)
anova(LMM_GCIPL_PP_SPCat2, LMM_GCIPL_PP_SPCat2sex)
#adding sex doesn't improve model performance

#c: both
LMM_GCIPL_PP_SPCat2both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_PPMS_SPCat2)
summary(LMM_GCIPL_PP_SPCat2both)
anova(LMM_GCIPL_PP_SPCat2, LMM_GCIPL_PP_SPCat2both)
#adding both sex and age doesn't improve model fit 

######################################################################################### 
###########################               SPMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 3,5-12,5 YEARS DISEASE DURATION      ##############
data_SPMS_SPCat1 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==0)
data_SPMS_SPCat1$Geschlecht<-as.factor(data_SPMS_SPCat1$Geschlecht)

LMM_GCIPL_SP_SPCat1 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_SPMS_SPCat1)
summary(LMM_GCIPL_SP_SPCat1)
confint(LMM_GCIPL_SP_SPCat1)
#best model fit - presented in paper

#a: age
LMM_GCIPL_SP_SPCat1age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_SPMS_SPCat1)
summary(LMM_GCIPL_SP_SPCat1age)
anova(LMM_GCIPL_SP_SPCat1age, LMM_GCIPL_SP_SPCat1)
#age doesn't improve model fit

#b: sex
LMM_GCIPL_SP_SPCat1sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_SPMS_SPCat1)
summary(LMM_GCIPL_SP_SPCat1sex)
anova(LMM_GCIPL_SP_SPCat1sex, LMM_GCIPL_SP_SPCat1)
#sex doesn't improve model fit

#c: both
LMM_GCIPL_SP_SPCat1both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=data_SPMS_SPCat1)
summary(LMM_GCIPL_SP_SPCat1both)
anova(LMM_GCIPL_SP_SPCat1both, LMM_GCIPL_SP_SPCat1)
#adding both sex and age doesn't improve model fit

##################     SECOND INTERVAL 12,6-16,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat2 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==1)
data_SPMS_SPCat2$Geschlecht<-as.factor(data_SPMS_SPCat2$Geschlecht)

LMM_GCIPL_SP_SPCat2 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_SPMS_SPCat2)
summary(LMM_GCIPL_SP_SPCat2)
#boundary (singular) fit: see help('isSingular')

isSingular(LMM_GCIPL_SP_SPCat2)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            11.13    3.336   
#ID        (Intercept)            50.86    7.131   
#ID.1      Disease_Duration_FINAL  0.00    0.000   
#Residual                          0.98    0.990

#random effect of disease duration doesn't explain variance
#try without this RE

LMM_GCIPL_SP_SPCat2_1 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat2)
summary(LMM_GCIPL_SP_SPCat2_1)
confint(LMM_GCIPL_SP_SPCat2_1)
#best model fit - presented in paper

#add covariates to best fitting model
#a: age
LMM_GCIPL_SP_SPCat2age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat2)
summary(LMM_GCIPL_SP_SPCat2age)
anova(LMM_GCIPL_SP_SPCat2age, LMM_GCIPL_SP_SPCat2_1)
#adding age doesn't improve model fit

#b: sex
LMM_GCIPL_SP_SPCat2sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                         data=data_SPMS_SPCat2)
summary(LMM_GCIPL_SP_SPCat2sex)
anova(LMM_GCIPL_SP_SPCat2sex, LMM_GCIPL_SP_SPCat2_1)
#adding sex doesn't improve model fit

#c: both
LMM_GCIPL_SP_SPCat2both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                          data=data_SPMS_SPCat2)
summary(LMM_GCIPL_SP_SPCat2both)
anova(LMM_GCIPL_SP_SPCat2both, LMM_GCIPL_SP_SPCat2_1)
#adding both sex and age doesn't improve model fit

##################     THIRD INTERVAL 16,6-20,5 YEARS DISEASE DURATION       ############
data_SPMS_SPCat3 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==2)
data_SPMS_SPCat3$Geschlecht<-as.factor(data_SPMS_SPCat3$Geschlecht)

LMM_GCIPL_SP_SPCat3 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_SPMS_SPCat3)
summary(LMM_GCIPL_SP_SPCat3)
confint(LMM_GCIPL_SP_SPCat3)
model_performance(LMM_GCIPL_SP_SPCat3)
plot(LMM_GCIPL_SP_SPCat3)
qqnorm(resid(LMM_GCIPL_SP_SPCat3))
qqline(resid(LMM_GCIPL_SP_SPCat3))
#AIC     |    AICc |     BIC |  RMSE | Sigma
#310.639 | 312.113 | 323.593 | 0.212 | 0.381


#a: age
LMM_GCIPL_SP_SPCat3age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_SPMS_SPCat3)
summary(LMM_GCIPL_SP_SPCat3age)
anova(LMM_GCIPL_SP_SPCat3age, LMM_GCIPL_SP_SPCat3)
#adding age doesn't improve model fit

#b: sex
LMM_GCIPL_SP_SPCat3sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_SPMS_SPCat3)
summary(LMM_GCIPL_SP_SPCat3sex)
#convergence error

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            14.76728 3.8428  
#ID        (Intercept)             8.57901 2.9290  
#ID.1      Disease_Duration_FINAL  0.04644 0.2155  
#Residual                          0.13585 0.3686 

#none of the RE with variance 0--> stick to model with all RE but not sex
#RE of disease duration explains least variance
#check if model without this RE but with sex as a covariate is better than model
#with RE of disease duration but without age

LMM_GCIPL_SP_SPCat3sex_1 <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat3)
summary(LMM_GCIPL_SP_SPCat3sex_1)
confint(LMM_GCIPL_SP_SPCat3sex_1)
model_performance(LMM_GCIPL_SP_SPCat3sex_1)
plot(LMM_GCIPL_SP_SPCat3sex_1)
qqnorm(resid(LMM_GCIPL_SP_SPCat3sex_1))
qqline(resid(LMM_GCIPL_SP_SPCat3sex_1))

#AIC     |    AICc |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
#305.544 | 307.017 | 318.497 |      0.995 |      0.116 | 0.995 | 0.260 | 0.442

#AIC and qqplot in favor of LMM_GCIPL_SP_SPCat3sex_1
#reduce complexity of the model and include sex instead of RE of disease duration 
#best model fit - presented in paper

#c: both
LMM_GCIPL_SP_SPCat3_both <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye), 
                                 data=data_SPMS_SPCat3)
summary(LMM_GCIPL_SP_SPCat3_both)
anova(LMM_GCIPL_SP_SPCat3_both, LMM_GCIPL_SP_SPCat3sex_1)
#adding age in addition to sex  doesn't improve model fit

##################     FOURTH INTERVAL 20,6-25,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat4 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==3)
data_SPMS_SPCat4$Geschlecht<-as.factor(data_SPMS_SPCat4$Geschlecht)

LMM_GCIPL_SP_SPCat4 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_SPMS_SPCat4)
summary(LMM_GCIPL_SP_SPCat4)
isSingular(LMM_GCIPL_SP_SPCat4)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             5.4089  2.3257  
#ID        (Intercept)            60.2157  7.7599  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.6584  0.8114 

#random effect of disease duration 0,0 variance --> try model without this RE
LMM_GCIPL_SP_SPCat4_new <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                                data=data_SPMS_SPCat4)
summary(LMM_GCIPL_SP_SPCat4_new)
confint(LMM_GCIPL_SP_SPCat4_new)
#best model fit - presented in paper

#add covariates to the best fitting model
#a: age
LMM_GCIPL_SP_SPCat4age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat4)
summary(LMM_GCIPL_SP_SPCat4age)
anova(LMM_GCIPL_SP_SPCat4age, LMM_GCIPL_SP_SPCat4_new)
#adding age doesn't improve model fit

#b: sex
LMM_GCIPL_SP_SPCat4sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat4)
summary(LMM_GCIPL_SP_SPCat4sex)
anova(LMM_GCIPL_SP_SPCat4sex, LMM_GCIPL_SP_SPCat4_new)
#adding sex doesn't improve model fit

#c: both
LMM_GCIPL_SP_SPCat4both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                                data=data_SPMS_SPCat4)
summary(LMM_GCIPL_SP_SPCat4both)
anova(LMM_GCIPL_SP_SPCat4both, LMM_GCIPL_SP_SPCat4_new)
#adding both sex and age doesn't improve model fit

##################     FIFTH INTERVAL 25,6-30,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat5 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==4)
data_SPMS_SPCat5$Geschlecht<-as.factor(data_SPMS_SPCat5$Geschlecht)

LMM_GCIPL_SP_SPCat5 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_SPMS_SPCat5)
summary(LMM_GCIPL_SP_SPCat5)
isSingular(LMM_GCIPL_SP_SPCat5)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)             5.9523  2.4397  
#ID        (Intercept)            87.1462  9.3352  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.3344  0.5783 

#RE of disease duration: variance 0,0 --> try model without this RED 

LMM_GCIPL_SP_SPCat5_new <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                                data=data_SPMS_SPCat5)
summary(LMM_GCIPL_SP_SPCat5_new)
confint(LMM_GCIPL_SP_SPCat5_new)

#add covariates to best fitting model
#a: age
LMM_GCIPL_SP_SPCat5age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat5)
summary(LMM_GCIPL_SP_SPCat5age)
anova(LMM_GCIPL_SP_SPCat5age, LMM_GCIPL_SP_SPCat5_new)
#adding age doesn't improve model fit

#b: sex
LMM_GCIPL_SP_SPCat5sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat5)
summary(LMM_GCIPL_SP_SPCat5sex)
anova(LMM_GCIPL_SP_SPCat5sex, LMM_GCIPL_SP_SPCat5_new)
#adding sex doesn't improve model fit

#c: both
LMM_GCIPL_SP_SPCat5both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                                data=data_SPMS_SPCat5)
summary(LMM_GCIPL_SP_SPCat5both)
anova(LMM_GCIPL_SP_SPCat5both, LMM_GCIPL_SP_SPCat5_new)
#adding both sex and age sign. improve model fit 
confint(LMM_GCIPL_SP_SPCat5both)
#best model fit - presented in paper

##################        SIXTH INTERVAL 30,6+ YEARS DISEASE DURATION        ############
data_SPMS_SPCat6 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==5)
data_SPMS_SPCat6$Geschlecht<-as.factor(data_SPMS_SPCat6$Geschlecht)

LMM_GCIPL_SP_SPCat6 <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_SPMS_SPCat6)
summary(LMM_GCIPL_SP_SPCat6)
confint(LMM_GCIPL_SP_SPCat6)
model_performance(LMM_GCIPL_SP_SPCat6)
plot(LMM_GCIPL_SP_SPCat6)
#AIC     |    AICc |     BIC |  RMSE | Sigma
#389.274 | 390.409 | 403.641 | 0.630 | 0.839

qqnorm(resid(LMM_GCIPL_SP_SPCat6))
qqline(resid(LMM_GCIPL_SP_SPCat6))
#best model fit - presented in paper

#a: age
LMM_GCIPL_SP_SPCat6age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_SPMS_SPCat6)
summary(LMM_GCIPL_SP_SPCat6age)
isSingular(LMM_GCIPL_SP_SPCat6age)
#true 

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            40.7273  6.3818  
#ID        (Intercept)            15.1148  3.8878  
#ID.1      Disease_Duration_FINAL  0.0000  0.0000  
#Residual                          0.7196  0.8483

#if age is added, variance of RE disease duration is reduced to 0 
#--> try model without this RE but with age as covariate

LMM_GCIPL_SP_SPCat6age_1 <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                               data=data_SPMS_SPCat6)
summary(LMM_GCIPL_SP_SPCat6age_1)
model_performance(LMM_GCIPL_SP_SPCat6age_1)
#AIC     |    AICc |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
#388.790 | 389.925 | 403.157 |      0.988 |      0.067 | 0.987 | 0.640 | 0.848

qqnorm(resid(LMM_GCIPL_SP_SPCat6age_1))
qqline(resid(LMM_GCIPL_SP_SPCat6age_1))
#aic slightly in favor of model with age
#qqplot in favor of model without age --> stick with model without age 

#b: sex
LMM_GCIPL_SP_SPCat6sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=data_SPMS_SPCat6)
summary(LMM_GCIPL_SP_SPCat6sex)
anova(LMM_GCIPL_SP_SPCat6sex, LMM_GCIPL_SP_SPCat6)
#adding sex doesn't improve model fit

#########################################################################################              
#########################################################################################              
#                                         INL
#########################################################################################              
######################################################################################### 

######################################################################################### 
###########################                HC           #################################
######################################################################################### 

###    FIRST (ONLY) INTERVAL 0-3,5 YEARS "DISEASE DURATION" (=Time since baseline)    ### 
data_HC_Cat1 <- subset(Data_LMM_HC, DiseaseDurationPPMSCategories_LMM==0)
data_HC_Cat1$Verlaufsform<-as.factor(data_HC_Cat1$Verlaufsform)

LMM_INL_HC_Cat1 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_HC_Cat1)
summary(LMM_INL_HC_Cat1)
confint(LMM_INL_HC_Cat1)
#best model fit - presented in paper

#Add Covariates
#a: Age
LMM_INL_HC_Cat1age <- lmer(INL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_HC_Cat1)
summary(LMM_INL_HC_Cat1age)
anova(LMM_INL_HC_Cat1age, LMM_INL_HC_Cat1)
#age doesn't improve model fit

#b: Sex
data_HC_Cat1$Geschlecht<-as.factor(data_HC_Cat1$Geschlecht)
LMM_INL_HC_Cat1sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_HC_Cat1)
summary(LMM_INL_HC_Cat1sex)
anova(LMM_INL_HC_Cat1sex, LMM_INL_HC_Cat1)
#sex doesn't improve model fit

#c: both
LMM_INL_HC_Cat1both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_HC_Cat1)
summary(LMM_INL_HC_Cat1both)
anova(LMM_INL_HC_Cat1both, LMM_INL_HC_Cat1)
#adding both sex and age doesn't improve model fit

######################################################################################### 
###########################               RRMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 0-3,5 YEARS DISEASE DURATION      #################
data_RRMS_Cat1 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==0)
data_RRMS_Cat1$Geschlecht<-as.factor(data_RRMS_Cat1$Geschlecht)

LMM_INL_RR_Cat1 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_RRMS_Cat1)
summary(LMM_INL_RR_Cat1)
confint(LMM_INL_RR_Cat1)

#add age as covariate
LMM_INL_RR_Cat1age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_RRMS_Cat1)
summary(LMM_INL_RR_Cat1age)
anova(LMM_INL_RR_Cat1, LMM_INL_RR_Cat1age)
#age doesn't improve the model

#b: sex
LMM_INL_RR1sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                       data=data_RRMS_Cat1)
summary(LMM_INL_RR1sex)
anova(LMM_INL_RR_Cat1, LMM_INL_RR1sex)
#better fit than without sex
confint(LMM_INL_RR1sex)
#best model fit - presented in paper

#c: both
LMM_INL_RR1both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_RRMS_Cat1)
summary(LMM_INL_RR1both)
anova(LMM_INL_RR1both, LMM_INL_RR1sex)
#adding age in addition to sex doesn't further improve model fit

##################      SECOND INTERVAL 3,6-5,5 YEARS DISEASE DURATION   ################
data_RRMS_Cat2 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==1)
data_RRMS_Cat2$Geschlecht<-as.factor(data_RRMS_Cat2$Geschlecht)

LMM_INL_RR_Cat2 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_RRMS_Cat2)
summary(LMM_INL_RR_Cat2)
isSingular(LMM_INL_RR_Cat2)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.3113   0.5579  
#ID        (Intercept)            4.7890   2.1884  
#ID.1      Disease_Duration_FINAL 0.0000   0.0000  
#Residual                         0.1604   0.4005

#RE of DD: variance=0 --> try without this RE

LMM_INL_RR_Cat2A <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye),
                         data=data_RRMS_Cat2)
summary(LMM_INL_RR_Cat2A)
confint(LMM_INL_RR_Cat2A)

#age
LMM_INL_RR_Cat2age <- lmer(INL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye), 
                           data=data_RRMS_Cat2)
summary(LMM_INL_RR_Cat2age)
anova(LMM_INL_RR_Cat2age,LMM_INL_RR_Cat2A)
#age doesn't improve model fit

#sex
LMM_INL_RR_Cat2sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye), 
                           data=data_RRMS_Cat2)
summary(LMM_INL_RR_Cat2sex)
anova(LMM_INL_RR_Cat2sex,LMM_INL_RR_Cat2A)
#better fit than without sex
confint(LMM_INL_RR_Cat2sex)
#best model fit - presented in paper

#both
LMM_INL_RR_Cat2both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                            data=data_RRMS_Cat2)
summary(LMM_INL_RR_Cat2both)
anova(LMM_INL_RR_Cat2both,LMM_INL_RR_Cat2sex)
#adding age in addition to sex doesn't further improve model fit

##################     THIRD INTERVAL 5,6-7,5 YEARS DISEASE DURATION    ################# 
data_RRMS_Cat3 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==2)
data_RRMS_Cat3$Geschlecht<-as.factor(data_RRMS_Cat3$Geschlecht)

LMM_INL_RR_Cat3 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_RRMS_Cat3)
summary(LMM_INL_RR_Cat3)
isSingular(LMM_INL_RR_Cat3)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.35539  0.5961  
#ID        (Intercept)            5.55517  2.3569  
#ID.1      Disease_Duration_FINAL 0.00000  0.0000  
#Residual                         0.09644  0.3106

#random effect of dd explains no variance --> try without it 

#random effect not necessary --> try without it
LMM_INL_RR_Cat3_new <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                            data=data_RRMS_Cat3)
summary(LMM_INL_RR_Cat3_new)
confint(LMM_INL_RR_Cat3_new)

#age
LMM_INL_RR_Cat3age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                           data=data_RRMS_Cat3)
summary(LMM_INL_RR_Cat3age)
anova(LMM_INL_RR_Cat3_new, LMM_INL_RR_Cat3age)
#age doesn't improve model fit

#sex
LMM_INL_RR_Cat3sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                           data=data_RRMS_Cat3)
summary(LMM_INL_RR_Cat3sex)
anova(LMM_INL_RR_Cat3sex, LMM_INL_RR_Cat3_new)
#adding sex sign. improves model fit
confint(LMM_INL_RR_Cat3sex)
#best model fit - presented in paper

#both
LMM_INL_RR_Cat3both <- lmer(INL_new~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye), 
                            data=data_RRMS_Cat3)
summary(LMM_INL_RR_Cat3both)
anova(LMM_INL_RR_Cat3both, LMM_INL_RR_Cat3sex)
#adding age in addition to sex doesn't further improve model fit

##################      FOURTH INTERVAL 7,6-10,5 YEARS DISEASE DURATION   ###############
data_RRMS_Cat4 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==3)
data_RRMS_Cat4$Geschlecht<-as.factor(data_RRMS_Cat4$Geschlecht)

LMM_INL_RRMS_Cat4 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_RRMS_Cat4)
summary(LMM_INL_RRMS_Cat4)
confint(LMM_INL_RRMS_Cat4)
#best model fit - presented in paper

#a:age
LMM_INL_RRMS_Cat4age <- lmer(INL_new~Disease_Duration_FINAL + Age +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat4)
summary(LMM_INL_RRMS_Cat4age)
anova(LMM_INL_RRMS_Cat4age, LMM_INL_RRMS_Cat4)
#age doesn't improve model fit

#b: sex
LMM_INL_RRMS_Cat4sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht +(1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_Cat4)
summary(LMM_INL_RRMS_Cat4sex)
anova(LMM_INL_RRMS_Cat4sex, LMM_INL_RRMS_Cat4)
#sex doesn't improve model fit

#b: both
LMM_INL_RRMS_Cat4sexage <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                data=data_RRMS_Cat4)
summary(LMM_INL_RRMS_Cat4sexage)
anova(LMM_INL_RRMS_Cat4, LMM_INL_RRMS_Cat4sexage)
#adding both sex and age doesn't improve model fit

#################      FIFTH INTERVAL 10,6-13,5 YEARS DISEASE DURATION   ################ 
data_RRMS_Cat5 <- subset(Data_LMM_RRMS, DiseaseDurationPPMSCategories_LMM==4)
data_RRMS_Cat5$Geschlecht<-as.factor(data_RRMS_Cat5$Geschlecht)

LMM_INL_RR_Cat5 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_RRMS_Cat5)
summary(LMM_INL_RR_Cat5)
isSingular(LMM_INL_RR_Cat5)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.2667   0.5164  
#ID        (Intercept)            7.5194   2.7421  
#ID.1      Disease_Duration_FINAL 0.0000   0.0000  
#Residual                         0.4389   0.6625

#RE of disease duration doesn't explain any variance --> eliminate this unnecessary RE

LMM_INL_RR_Cat5_2 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                          data=data_RRMS_Cat5)
summary(LMM_INL_RR_Cat5_2)
confint(LMM_INL_RR_Cat5_2)

#a: age
LMM_INL_RR_Cat5age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                           data=data_RRMS_Cat5)
summary(LMM_INL_RR_Cat5age)
anova(LMM_INL_RR_Cat5_2, LMM_INL_RR_Cat5age)
#doesn't improve model

#b: sex
LMM_INL_RR_Cat5sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                           data=data_RRMS_Cat5)
summary(LMM_INL_RR_Cat5sex)
anova(LMM_INL_RR_Cat5_2, LMM_INL_RR_Cat5sex)
#improves model
confint(LMM_INL_RR_Cat5sex)
#best model fit - presented in paper

#b: both
LMM_INL_RR_Cat5sexage <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                              data=data_RRMS_Cat5)
summary(LMM_INL_RR_Cat5sexage)
anova(LMM_INL_RR_Cat5sexage, LMM_INL_RR_Cat5sex)
#adding age in addition to sex doesn't further improve model fit

##################     SIXTH INTERVAL 13,6-16,5 YEARS DISEASE DURATION   ################
data_RRMS_SPCat1 <- subset(Data_LMM_RRMS, DiseaseDurationSPMSCategories_LMM==1)
data_RRMS_SPCat1$Geschlecht<-as.factor(data_RRMS_SPCat1$Geschlecht)

LMM_INL_RR_SPCat1 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_RRMS_SPCat1)
summary(LMM_INL_RR_SPCat1)
confint(LMM_INL_RR_SPCat1)
#best model fit - presented in paper

#add age as covariate
LMM_INL_RR_SPCat1age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_RRMS_SPCat1)
summary(LMM_INL_RR_SPCat1age)
anova(LMM_INL_RR_SPCat1, LMM_INL_RR_SPCat1age)
#age doesn't improve model fit

#b: sex
LMM_INL_RR_SPCat1_sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_RRMS_SPCat1)
summary(LMM_INL_RR_SPCat1_sex)
anova(LMM_INL_RR_SPCat1, LMM_INL_RR_SPCat1_sex)
#sex doesn't improve model fit

#c: both
LMM_INL_RR_SPCat1_both <- lmer(INL_new~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_RRMS_SPCat1)
summary(LMM_INL_RR_SPCat1_both)
anova(LMM_INL_RR_SPCat1, LMM_INL_RR_SPCat1_both)
#adding both age and sex doesn't improve model fit

######################################################################################### 
###########################               PPMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 0-3,5 YEARS DISEASE DURATION      #################
data_PPMS_Cat1 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==0)
data_PPMS_Cat1$Geschlecht<-as.factor(data_PPMS_Cat1$Geschlecht)

LMM_INL_PP_Cat1 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat1)
summary(LMM_INL_PP_Cat1)
confint(LMM_INL_PP_Cat1)
#best model fit - presented in paper

#a: age
LMM_INL_PP_Cat1age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_PPMS_Cat1)
summary(LMM_INL_PP_Cat1age)
anova(LMM_INL_PP_Cat1age, LMM_INL_PP_Cat1)
#age doesn't improve model fit

#b: sex
LMM_INL_PP1sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                       data=data_PPMS_Cat1)
summary(LMM_INL_PP1sex)
anova(LMM_INL_PP1sex, LMM_INL_PP_Cat1)
#sex doesn't improve model fit

#c: both
LMM_INL_PP1both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat1)
summary(LMM_INL_PP1both)
anova(LMM_INL_PP1both, LMM_INL_PP_Cat1)
#adding both age and sex doesn't improve model fit

##################      SECOND INTERVAL 3,6-5,5 YEARS DISEASE DURATION   ################
data_PPMS_Cat2 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==1)
data_PPMS_Cat2$Geschlecht<-as.factor(data_PPMS_Cat2$Geschlecht)

LMM_INL_PP_Cat2 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat2)
summary(LMM_INL_PP_Cat2)
isSingular(LMM_INL_PP_Cat2)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.0000   0.0000  
#ID        (Intercept)            3.0049   1.7335  
#ID.1      Disease_Duration_FINAL 0.0596   0.2441  
#Residual                         0.2947   0.5429 

#RE of Eye within subject explains 0 variance --> try without this RE

LMM_INL_PP_Cat2_1 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat2)
summary(LMM_INL_PP_Cat2_1)
confint(LMM_INL_PP_Cat2_1)
#best model fit - presented in paper

#a: age
LMM_INL_PP_Cat2age <- lmer(INL_new~Disease_Duration_FINAL +  Age+ (1|ID) + (0+Disease_Duration_FINAL|ID), 
                           data=data_PPMS_Cat2)
summary(LMM_INL_PP_Cat2age)
anova(LMM_INL_PP_Cat2age, LMM_INL_PP_Cat2_1)
#age doesn't improve model fit

#b: sex
LMM_INL_PP_Cat2sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID) + (0+Disease_Duration_FINAL|ID), 
                       data=data_PPMS_Cat2)
summary(LMM_INL_PP_Cat2sex)
anova(LMM_INL_PP_Cat2sex, LMM_INL_PP_Cat2_1)
#sex doesn't improve model fit

#c: both
LMM_INL_PP_Cat2both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat2)
summary(LMM_INL_PP_Cat2both)
anova(LMM_INL_PP_Cat2both, LMM_INL_PP_Cat2_1)
#adding both sex and age doesn't improve model fit

##################     THIRD INTERVAL 5,6-7,5 YEARS DISEASE DURATION    ################# 
data_PPMS_Cat3 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==2)
data_PPMS_Cat3$Geschlecht<-as.factor(data_PPMS_Cat3$Geschlecht)

LMM_INL_PP_Cat3 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat3)
summary(LMM_INL_PP_Cat3)
confint(LMM_INL_PP_Cat3)
#best model fit - presented in paper

#a: age
LMM_INL_PP_Cat3age <- lmer(INL_new~Disease_Duration_FINAL + Age +  (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_PPMS_Cat3)
summary(LMM_INL_PP_Cat3age)
anova(LMM_INL_PP_Cat3age, LMM_INL_PP_Cat3)
#adding age doesn't improve model fit

#b: sex
LMM_INL_PP_Cat3sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht +  (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=data_PPMS_Cat3)
summary(LMM_INL_PP_Cat3sex)
anova(LMM_INL_PP_Cat3sex, LMM_INL_PP_Cat3)
#adding sex doesn't improve model fit

#c: both
LMM_INL_PP_Cat3both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_Cat3)
summary(LMM_INL_PP_Cat3both)
anova(LMM_INL_PP_Cat3both, LMM_INL_PP_Cat3)
#adding sex and age doesn't improve model fit

##################      FOURTH INTERVAL 7,6-10,5 YEARS DISEASE DURATION   ###############
data_PPMS_Cat4 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==3)
data_PPMS_Cat4$Geschlecht<-as.factor(data_PPMS_Cat4$Geschlecht)

LMM_INL_PP_Cat4 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat4)
summary(LMM_INL_PP_Cat4)
isSingular(LMM_INL_PP_Cat4)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.2943   0.5425  
#ID        (Intercept)            4.3774   2.0922  
#ID.1      Disease_Duration_FINAL 0.0000   0.0000  
#Residual                         0.1890   0.4348  

#RE of Disease duration doesn't explain any variance --> try model without this RE

LMM_INL_PP_Cat4_new <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat4)
summary(LMM_INL_PP_Cat4_new)
confint(LMM_INL_PP_Cat4_new)
#best model fit - presented in paper

#a: age
LMM_INL_PP_Cat4age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                           data=data_PPMS_Cat4)
summary(LMM_INL_PP_Cat4age)
anova(LMM_INL_PP_Cat4age, LMM_INL_PP_Cat4_new)
#age doesn't improve model fit

#b: sex
LMM_INL_PP_Cat4sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                           data=data_PPMS_Cat4)
summary(LMM_INL_PP_Cat4sex)
anova(LMM_INL_PP_Cat4sex, LMM_INL_PP_Cat4_new)
#age doesn't improve model fit

#c: both
LMM_INL_PP_Cat4both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                            data=data_PPMS_Cat4)
summary(LMM_INL_PP_Cat4both)
anova(LMM_INL_PP_Cat4both, LMM_INL_PP_Cat4_new)
#adding sex and age doesn't improve model fit

#################      FIFTH INTERVAL 10,6-13,5 YEARS DISEASE DURATION   ################ 
data_PPMS_Cat5 <- subset(Data_LMM_PPMS, DiseaseDurationPPMSCategories_LMM==4)
data_PPMS_Cat5$Geschlecht<-as.factor(data_PPMS_Cat5$Geschlecht)

LMM_INL_PP_Cat5 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat5)
summary(LMM_INL_PP_Cat5)
isSingular(LMM_INL_PP_Cat5)
#singular fit

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            8.709e-09 9.332e-05
#ID        (Intercept)            0.000e+00 0.000e+00
#ID.1      Disease_Duration_FINAL 3.350e-02 1.830e-01
#Residual                         6.095e-01 7.807e-01

#random intercept per subject explains 0 variance --> try model without this RE

LMM_INL_PP_Cat5_new <- lmer(INL_new~Disease_Duration_FINAL + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_Cat5)
summary(LMM_INL_PP_Cat5_new)
#still singular
#Random effects:
#Groups   Name                   Variance  Std.Dev. 
#ID_Eye   (Intercept)            1.052e-10 1.026e-05
#ID       Disease_Duration_FINAL 3.350e-02 1.830e-01
#Residual                        6.095e-01 7.807e-01

#random intercept per eye doesn't explain any variance 

LMM_INL_PP_Cat5_new2 <- lmer(INL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_Cat5)
summary(LMM_INL_PP_Cat5_new2)
confint(LMM_INL_PP_Cat5_new2)
#best model fit - presented in paper

#a: age
LMM_INL_PP_Cat5age <- lmer(INL_new~Disease_Duration_FINAL + Age +(0+Disease_Duration_FINAL|ID), 
                           data=data_PPMS_Cat5)
summary(LMM_INL_PP_Cat5age)
anova(LMM_INL_PP_Cat5age, LMM_INL_PP_Cat5_new2)
#adding age doesn't improve model fit

#b: sex
LMM_INL_PP_Cat5sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                           data=data_PPMS_Cat5)
summary(LMM_INL_PP_Cat5sex)
anova(LMM_INL_PP_Cat5sex, LMM_INL_PP_Cat5_new2)
#adding sex doesn't improve model fit

#c: both
LMM_INL_PP_Cat5both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age+ (0+Disease_Duration_FINAL|ID), 
                            data=data_PPMS_Cat5)
summary(LMM_INL_PP_Cat5both)
anova(LMM_INL_PP_Cat5both, LMM_INL_PP_Cat5_new2)
#adding both sex and age doesn't improve model fit

##################     SIXTH INTERVAL 13,6-16,5 YEARS DISEASE DURATION   ################
data_PPMS_SPCat1 <- subset(Data_LMM_PPMS, DiseaseDurationSPMSCategories_LMM==1)
data_PPMS_SPCat1$Geschlecht<-as.factor(data_PPMS_SPCat1$Geschlecht)

LMM_INL_PP_SPCat1 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_SPCat1)
summary(LMM_INL_PP_SPCat1)
confint(LMM_INL_PP_SPCat1)
#best model fit - presented in paper

#a: age
LMM_INL_PP_SPCat1age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_PPMS_SPCat1)
summary(LMM_INL_PP_SPCat1age)
anova(LMM_INL_PP_SPCat1age, LMM_INL_PP_SPCat1)
#age doesn't improve model fit

#b: sex
LMM_INL_PP_SPCat1sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                       data=data_PPMS_SPCat1)
summary(LMM_INL_PP_SPCat1sex)
anova(LMM_INL_PP_SPCat1sex, LMM_INL_PP_SPCat1)
#sex doesn't improve model fit

#c: both
LMM_INL_PP_SPCat1both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_PPMS_SPCat1)
summary(LMM_INL_PP_SPCat1both)
anova(LMM_INL_PP_SPCat1both, LMM_INL_PP_SPCat1)
#adding both sex and age doesn't improve model fit

#################      SEVENTH INTERVAL 16,6-20,5 YEARS DISEASE DURATION   ##############
data_PPMS_SPCat2 <- subset(Data_LMM_PPMS, DiseaseDurationSPMSCategories_LMM==2)
data_PPMS_SPCat2$Geschlecht<-as.factor(data_PPMS_SPCat2$Geschlecht)

LMM_INL_PP_SPCat2 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_PPMS_SPCat2)
summary(LMM_INL_PP_SPCat2)
#Model failed to converge with max|grad| = 0.0025786 (tol = 0.002, component 1)

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.217106 0.46595 
#ID        (Intercept)            1.463014 1.20955 
#ID.1      Disease_Duration_FINAL 0.003561 0.05968 
#Residual                         0.192062 0.43825 

#RE of disease duration explains ~0 variance --> eliminate this RE

LMM_INL_PP_SPCat2_new <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                          data=data_PPMS_SPCat2)
summary(LMM_INL_PP_SPCat2_new)
confint(LMM_INL_PP_SPCat2_new)
#best model fit - presented in paper

#a: age
LMM_INL_PP_SPCat2age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                             data=data_PPMS_SPCat2)
summary(LMM_INL_PP_SPCat2age)
anova(LMM_INL_PP_SPCat2age, LMM_INL_PP_SPCat2_new)
#adding age doesn't improve model fit

#b: sex
LMM_INL_PP_SPCat2sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                             data=data_PPMS_SPCat2)
summary(LMM_INL_PP_SPCat2sex)
anova(LMM_INL_PP_SPCat2_new, LMM_INL_PP_SPCat2sex)
#adding sex doesn't improve model fit

#c: both
LMM_INL_PP_SPCat2both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                              data=data_PPMS_SPCat2)
summary(LMM_INL_PP_SPCat2both)
anova(LMM_INL_PP_SPCat2_new, LMM_INL_PP_SPCat2both)
#adding both age and sex doesn't improve model fit

######################################################################################### 
###########################               SPMS          #################################
######################################################################################### 

##################     FIRST INTERVAL 3,5-12,5 YEARS DISEASE DURATION      ##############
data_SPMS_SPCat1 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==0)
data_SPMS_SPCat1$Geschlecht<-as.factor(data_SPMS_SPCat1$Geschlecht)

LMM_INL_SP_SPCat1 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_SPMS_SPCat1)
summary(LMM_INL_SP_SPCat1)
confint(LMM_INL_SP_SPCat1)
#best model fit - presented in paper

#a: age
LMM_INL_SP_SPCat1age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat1)
summary(LMM_INL_SP_SPCat1age)
anova(LMM_INL_SP_SPCat1age, LMM_INL_SP_SPCat1)
#adding age doesn't improve model 

#b: sex
LMM_INL_SP_SPCat1sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat1)
summary(LMM_INL_SP_SPCat1sex)
anova(LMM_INL_SP_SPCat1sex, LMM_INL_SP_SPCat1)
#adding sex doesn't improve model fit

#c: both
LMM_INL_SP_SPCat1both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_SPMS_SPCat1)
summary(LMM_INL_SP_SPCat1both)
anova(LMM_INL_SP_SPCat1both, LMM_INL_SP_SPCat1)
#adding both sex and age doesn't improve model fit

##################     SECOND INTERVAL 12,6-16,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat2 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==1)
data_SPMS_SPCat2$Geschlecht<-as.factor(data_SPMS_SPCat2$Geschlecht)

LMM_INL_SP_SPCat2 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_SPMS_SPCat2)
summary(LMM_INL_SP_SPCat2)
confint(LMM_INL_SP_SPCat2)
model_performance(LMM_INL_SP_SPCat2)

#AIC     |    AICc |     BIC |  RMSE | Sigma
#203.788 | 205.212 | 216.926 | 0.290 | 0.440

#a: age
LMM_INL_SP_SPCat2age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat2)
summary(LMM_INL_SP_SPCat2age)
isSingular(LMM_INL_SP_SPCat2age)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.32125  0.5668  
#ID        (Intercept)            0.00000  0.0000  
#ID.1      Disease_Duration_FINAL 0.01453  0.1206  
#Residual                         0.18394  0.4289  

#random intercept per subject explains no variance if age is included
#--> try model wihtout random intercept per subject but with age as a covariate
#and check which model is better

LMM_INL_SP_SPCat2age_2 <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat2)
summary(LMM_INL_SP_SPCat2age_2)
model_performance(LMM_INL_SP_SPCat2age_2)
#AIC     |    AICc |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
#207.117 | 208.541 | 220.255 |      0.951 |      0.040 | 0.949 | 0.281 | 0.429

#AIC is in favor of the model with a random intercept and without age

#b: sex
LMM_INL_SP_SPCat2sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                       data=data_SPMS_SPCat2)
summary(LMM_INL_SP_SPCat2sex)
#singular fit

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            3.203e-01 0.5659222
#ID        (Intercept)            1.122e-09 0.0000335
#ID.1      Disease_Duration_FINAL 9.723e-03 0.0986037
#Residual                         2.195e-01 0.468508

##random intercept per subject explains no variance if sex is included
#--> try model wihtout random intercept per subject but with sex as a covariate
#and check which model is better

LMM_INL_SP_SPCat2sex_2 <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat2)
summary(LMM_INL_SP_SPCat2sex_2)
model_performance(LMM_INL_SP_SPCat2sex_2)
#AIC     |    AICc |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
#198.086 | 199.510 | 211.224 |      0.931 |      0.181 | 0.916 | 0.314 | 0.469

#aic in favor of model with sex but without random intercept per subject

confint(LMM_INL_SP_SPCat2sex_2)
#best model fit - presented in paper

#c: both
LMM_INL_SP2both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                        data=data_SPMS_SPCat2)
summary(LMM_INL_SP2both)
anova(LMM_INL_SP2both, LMM_INL_SP_SPCat2sex_2)
#adding age in addition to sex doesn't further improve model fit

##################     THIRD INTERVAL 16,6-20,5 YEARS DISEASE DURATION       ############
data_SPMS_SPCat3 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==2)
data_SPMS_SPCat3$Geschlecht<-as.factor(data_SPMS_SPCat3$Geschlecht)

LMM_INL_SP_SPCat3 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_SPMS_SPCat3)
summary(LMM_INL_SP_SPCat3)
isSingular(LMM_INL_SP_SPCat3)
#singular fit

#Random effects:
#Groups    Name                   Variance Std.Dev.
#ID_Eye.ID (Intercept)            0.055291 0.23514 
#ID        (Intercept)            0.000000 0.00000 
#ID.1      Disease_Duration_FINAL 0.009975 0.09987 
#Residual                         0.226933 0.47637 

#random intercept per subject doesn't explain any variance --> try without this RE
LMM_INL_SP_SPCat3_new <- lmer(INL_new~Disease_Duration_FINAL + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_SPMS_SPCat3)
summary(LMM_INL_SP_SPCat3_new)
confint(LMM_INL_SP_SPCat3_new)
#best model fit - presented in paper

#a: age
LMM_INL_SP_SPCat3age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat3)
summary(LMM_INL_SP_SPCat3age)
anova(LMM_INL_SP_SPCat3age, LMM_INL_SP_SPCat3_new)
#adding age doesn't improve model fit

#b: sex
LMM_INL_SP_SPCat3sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat3)
summary(LMM_INL_SP_SPCat3sex)
anova(LMM_INL_SP_SPCat3sex, LMM_INL_SP_SPCat3_new)
#adding sex doesn't improve model fit

#c: Both
LMM_INL_SP_SPCat3both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht +  Age + (1|ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_SPMS_SPCat3)
summary(LMM_INL_SP_SPCat3both)
anova(LMM_INL_SP_SPCat3both, LMM_INL_SP_SPCat3_new)
#adding both age and sex doesn't improve model fit

##################     FOURTH INTERVAL 20,6-25,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat4 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==3)
data_SPMS_SPCat4$Geschlecht<-as.factor(data_SPMS_SPCat4$Geschlecht)

LMM_INL_SP_SPCat4 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_SPMS_SPCat4)
summary(LMM_INL_SP_SPCat4)
#Model failed to converge with max|grad| = 0.0166431 (tol = 0.002, component 1)

#Random effects:
#Groups    Name                   Variance  Std.Dev. 
#ID_Eye.ID (Intercept)            9.190e-02 0.3031579
#ID        (Intercept)            5.275e+00 2.2966585
#ID.1      Disease_Duration_FINAL 1.063e-08 0.0001031
#Residual                         2.316e-01 0.4812384

#RE of diseae duration explains ~0 variance--> try model without this RE
LMM_INL_SP_SPCat4_new <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat4)
summary(LMM_INL_SP_SPCat4_new)
confint(LMM_INL_SP_SPCat4_new)
#best model fit - presented in paper

#a: age
LMM_INL_SP_SPCat4age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye), 
                             data=data_SPMS_SPCat4)
summary(LMM_INL_SP_SPCat4age)
anova(LMM_INL_SP_SPCat4age, LMM_INL_SP_SPCat4_new)
#age doesn't improve model fit

#b: sex
LMM_INL_SP_SPCat4sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye), 
                             data=data_SPMS_SPCat4)
summary(LMM_INL_SP_SPCat4sex)
anova(LMM_INL_SP_SPCat4sex, LMM_INL_SP_SPCat4_new)
#sex doesn't improve model fit

#c: both
LMM_INL_SP_SPCat4both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye), 
                              data=data_SPMS_SPCat4)
summary(LMM_INL_SP_SPCat4both)
anova(LMM_INL_SP_SPCat4both, LMM_INL_SP_SPCat4_new)
#adding both age and sex doesn't improve model fit

##################     FIFTH INTERVAL 25,6-30,5 YEARS DISEASE DURATION      ############
data_SPMS_SPCat5 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==4)
data_SPMS_SPCat5$Geschlecht<-as.factor(data_SPMS_SPCat5$Geschlecht)

LMM_INL_SP_SPCat5 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_SPMS_SPCat5)
summary(LMM_INL_SP_SPCat5)
confint(LMM_INL_SP_SPCat5)
#best model fit - presented in paper

#a: age
LMM_INL_SP_SPCat5age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat5)
summary(LMM_INL_SP_SPCat5age)
anova(LMM_INL_SP_SPCat5age, LMM_INL_SP_SPCat5)
#doesn't improve model fit

#b: sex
LMM_INL_SP_SPCat5sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat5)
summary(LMM_INL_SP_SPCat5sex)
anova(LMM_INL_SP_SPCat5sex, LMM_INL_SP_SPCat5)
#doesn't improve model fit

#c: both
LMM_INL_SP_SPCat5both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=data_SPMS_SPCat5)
summary(LMM_INL_SP_SPCat5both)
anova(LMM_INL_SP_SPCat5both, LMM_INL_SP_SPCat5)
#doesn't improve model fit

##################        SIXTH INTERVAL 30,6+ YEARS DISEASE DURATION        ############
data_SPMS_SPCat6 <- subset(Data_LMM_SPMS, DiseaseDurationSPMSCategories_LMM==5)
data_SPMS_SPCat6$Geschlecht<-as.factor(data_SPMS_SPCat6$Geschlecht)

LMM_INL_SP_SPCat6 <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=data_SPMS_SPCat6)
summary(LMM_INL_SP_SPCat6)
confint(LMM_INL_SP_SPCat6)
#best model fit - presented in paper

#a: age
LMM_INL_SP_SPCat6age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat6)
summary(LMM_INL_SP_SPCat6age)
anova(LMM_INL_SP_SPCat6age, LMM_INL_SP_SPCat6)
#doesn't improve model fit

#b: sex
LMM_INL_SP_SPCat6sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat6)
summary(LMM_INL_SP_SPCat6sex)
anova(LMM_INL_SP_SPCat6sex, LMM_INL_SP_SPCat6)
#doesn't improve model fit

#c: both
LMM_INL_SP_SPCat6sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=data_SPMS_SPCat6)
summary(LMM_INL_SP_SPCat6sex)
anova(LMM_INL_SP_SPCat6sex, LMM_INL_SP_SPCat6)
#doesn't improve model fit
