#----------------Longitudinal linear mixed-effects models (LMM)---------------#

# This script computes LMM analyzing the overall annualized rate of thickness 
# change (globalpRNFL, macular RNFL, GCIPL, INL) over time. To account for 
# confounding influences, age and/or sex were subsequently added to the model 
# if the model fit could be improved. 
# The results are presented in Supplementary Table 3. 
# Usual runtime ~= 5 minutes 

#---------------------------------R packages---------------------------------#  

library(lme4)
library(lmerTest) #needed for p-values
library(haven)
library(performance)
library(car)

#-----------------------Data loading & data management-----------------------# 

OCT_data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)
Data_LMM <- subset(OCT_data, OCT_included==1)
Data_LMM_HC <- subset(OCT_data, OCT_included==1 & Verlaufsform==0)
Data_LMM_RRMS <- subset(OCT_data, OCT_included==1 & Verlaufsform==1)
Data_LMM_PPMS <- subset(OCT_data, OCT_included==1 & Verlaufsform==2)
Data_LMM_SPMS <- subset(OCT_data, OCT_included==1 & Verlaufsform==3)
#Verlaufsform=0 --> HC
#Verlaufsform=1 --> RRMS
#Verlaufsform=2 --> PPMS
#Verlaufsform=3 --> SPMS


#-------------------------------------GCIPL-----------------------------------#    
#---RRMS---#
Data_LMM_RRMS$Geschlecht<-as.factor(Data_LMM_RRMS$Geschlecht)

LMM_GCIPL_RR_total <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_RRMS)
summary(LMM_GCIPL_RR_total)
confint(LMM_GCIPL_RR_total)
#best model fit - presented in paper

#add age as covariate
LMM_GCIPL_RR_total_age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_RRMS)
summary(LMM_GCIPL_RR_total_age)
anova(LMM_GCIPL_RR_total_age, LMM_GCIPL_RR_total)
#age doesn't improve model fit

#b: sex
LMM_GCIPL_RR_total_sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_RRMS)
summary(LMM_GCIPL_RR_total_sex)
anova(LMM_GCIPL_RR_total, LMM_GCIPL_RR_total_sex)
#sex doesn't sign. improve model 

#c: both
LMM_GCIPL_RR_total_both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye)+ (0+Disease_Duration_FINAL|ID), 
                                data=Data_LMM_RRMS)
summary(LMM_GCIPL_RR_total_both)
anova(LMM_GCIPL_RR_total, LMM_GCIPL_RR_total_both)
#adding both sex and age doesn't improve model fit

#---PPMS---#
Data_LMM_PPMS$Geschlecht<-as.factor(Data_LMM_PPMS$Geschlecht)
LMM_GCIPL_PP_total <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_PPMS)
summary(LMM_GCIPL_PP_total)
confint(LMM_GCIPL_PP_total)
#best model fit - presented in paper

#a: age
LMM_GCIPL_PP_totalage <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_PPMS)
summary(LMM_GCIPL_PP_totalage)
anova(LMM_GCIPL_PP_totalage, LMM_GCIPL_PP_total)
#age doesn't sign. improve model fit

#b: sex
LMM_GCIPL_PP_total_sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_PPMS)
summary(LMM_GCIPL_PP_total_sex)
anova(LMM_GCIPL_PP_total_sex, LMM_GCIPL_PP_total)
#sex doesn't sign. improve model fit

#c: both
LMM_GCIPL_PP_total_both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=Data_LMM_PPMS)
summary(LMM_GCIPL_PP_total_both)
anova(LMM_GCIPL_PP_total_both, LMM_GCIPL_PP_total)
#adding both sex and age doesn't improve model fit

#---SPMS---#
Data_LMM_SPMS$Geschlecht<-as.factor(Data_LMM_SPMS$Geschlecht)
LMM_GCIPL_SP_total <- lmer(GCIPL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_SPMS)
summary(LMM_GCIPL_SP_total)
confint(LMM_GCIPL_SP_total)
#best model fit - presented in paper

#a: age
LMM_GCIPL_SP_total_age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_SPMS)
summary(LMM_GCIPL_SP_total_age)
anova(LMM_GCIPL_SP_total_age, LMM_GCIPL_SP_total)
#age doesn't sign. improve model fit

#b: sex
LMM_GCIPL_SP_total_sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_SPMS)
summary(LMM_GCIPL_SP_total_sex)
anova(LMM_GCIPL_SP_total_sex, LMM_GCIPL_SP_total)
#sex doesn't sign. improve model fit

#c: both
LMM_GCIPL_SP_total_both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                data=Data_LMM_SPMS)
summary(LMM_GCIPL_SP_total_both)
anova(LMM_GCIPL_SP_total_both, LMM_GCIPL_SP_total)
#adding both sex and age doesn't sign. improve model fit


#-------------------------------------pRNFL-----------------------------------#    
#---RRMS---#
Data_LMM_RRMS$Geschlecht<-as.factor(Data_LMM_RRMS$Geschlecht)
LMM_globalpRNFL_RR_total <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                 data=Data_LMM_RRMS)
summary(LMM_globalpRNFL_RR_total)
confint(LMM_globalpRNFL_RR_total)
#best model fit - presented in paper

#add age as covariate
LMM_globalpRNFL_RR_totalage <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                    data=Data_LMM_RRMS)
summary(LMM_globalpRNFL_RR_totalage)
anova(LMM_globalpRNFL_RR_total, LMM_globalpRNFL_RR_totalage)
#age doesn't sign. improve model fit

#b: sex
LMM_globalpRNFL_RR_total_sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=Data_LMM_RRMS)
summary(LMM_globalpRNFL_RR_total_sex)
anova(LMM_globalpRNFL_RR_total, LMM_globalpRNFL_RR_total_sex)
#sex doesn't sign. improve model fit

#c: both
LMM_globalpRNFL_RR_total_both <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=Data_LMM_RRMS)
summary(LMM_globalpRNFL_RR_total_both)
anova(LMM_globalpRNFL_RR_total, LMM_globalpRNFL_RR_total_both)
#adding both sex and age doesn't improve model fit


#----PPMS---#
Data_LMM_PPMS$Geschlecht<-as.factor(Data_LMM_PPMS$Geschlecht)
LMM_globalpRNFL_PP_total <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                 data=Data_LMM_PPMS)
summary(LMM_globalpRNFL_PP_total)

#a: age
LMM_globalpRNFL_PP_total_age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=Data_LMM_PPMS)
summary(LMM_globalpRNFL_PP_total_age)
confint(LMM_globalpRNFL_PP_total_age)
anova(LMM_globalpRNFL_PP_total_age, LMM_globalpRNFL_PP_total)
#age sign. improves model fit
#best model fit - presented in paper

#b: sex
LMM_globalpRNFL_PP_total_sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=Data_LMM_PPMS)
summary(LMM_globalpRNFL_PP_total_sex)
anova(LMM_globalpRNFL_PP_total_sex, LMM_globalpRNFL_PP_total)
#sex doesn't sign. improve model fit

#c: both
LMM_globalpRNFL_SP_total_both <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=Data_LMM_PPMS)
summary(LMM_globalpRNFL_SP_total_both)
anova(LMM_globalpRNFL_SP_total_both, LMM_globalpRNFL_PP_total_age)
#adding sex in addition to age doesn't further improve model fit

#---SPMS---#
Data_LMM_SPMS$Geschlecht<-as.factor(Data_LMM_SPMS$Geschlecht)
LMM_globalpRNFL_SP_total <- lmer(globalpRNFL~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                 data=Data_LMM_SPMS)
summary(LMM_globalpRNFL_SP_total)
confint(LMM_globalpRNFL_SP_total)
#best model fit - presented in paper

#a: age
LMM_globalpRNFL_SP_total_age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=Data_LMM_SPMS)
summary(LMM_globalpRNFL_SP_total_age)
anova(LMM_globalpRNFL_SP_total_age, LMM_globalpRNFL_SP_total)
#age doesn't sign. improve model fit

#b: sex
LMM_globalpRNFL_SP_total_sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                     data=Data_LMM_SPMS)
summary(LMM_globalpRNFL_SP_total_sex)
anova(LMM_globalpRNFL_SP_total_sex, LMM_globalpRNFL_SP_total)
#sex doesn't sign. improve model fit

#c: both
LMM_globalpRNFL_SP_total_both <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                                      data=Data_LMM_SPMS)
summary(LMM_globalpRNFL_SP_total_both)
anova(LMM_globalpRNFL_SP_total_both, LMM_globalpRNFL_SP_total)
#adding both sex and age doesn't improve model fit


#-------------------------------------mRNFL-----------------------------------#    
#---RRMS---#
Data_LMM_RRMS$Geschlecht<-as.factor(Data_LMM_RRMS$Geschlecht)
LMM_RNFL_RR_total <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_RRMS)
summary(LMM_RNFL_RR_total)
confint(LMM_RNFL_RR_total)
#best model fit - presented in paper

#add age as covariate
LMM_RNFL_RR_total_age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_RRMS)
summary(LMM_RNFL_RR_total_age)
anova(LMM_RNFL_RR_total, LMM_RNFL_RR_total_age)
#age doesn't sign. improve model fit

#b: sex
LMM_RNFL_RR_total_sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_RRMS)
summary(LMM_RNFL_RR_total_sex)
anova(LMM_RNFL_RR_total, LMM_RNFL_RR_total_sex)
#sex doesn't sign. improve model fit

#c: both
LMM_RNFL_RR_total_both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_RRMS)
summary(LMM_RNFL_RR_total_both)
anova(LMM_RNFL_RR_total_both, LMM_RNFL_RR_total)
#adding both sex and age doesn't improve model fit

#---PPMS---#
Data_LMM_PPMS$Geschlecht<-as.factor(Data_LMM_PPMS$Geschlecht)
LMM_RNFL_PP_total <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_PPMS)
summary(LMM_RNFL_PP_total)
confint(LMM_RNFL_PP_total)
#best model fit - presented in paper

#a: age
LMM_RNFL_PP_total_age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye)  + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_PPMS)
summary(LMM_RNFL_PP_total_age)
anova(LMM_RNFL_PP_total, LMM_RNFL_PP_total_age)
#age doesn't sign. improve model fit

#b: sex
LMM_RNFL_PP_total_sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_PPMS)
summary(LMM_RNFL_PP_total_sex)
anova(LMM_RNFL_PP_total, LMM_RNFL_PP_total_sex)
#sex sign. improves model fit
confint(LMM_RNFL_PP_total_sex)

#c: both
LMM_RNFL_PP_total_both <- lmer(RNFL_new~Disease_Duration_FINAL + Age+ Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_PPMS)
summary(LMM_RNFL_PP_total_both)
anova(LMM_RNFL_PP_total_sex, LMM_RNFL_PP_total_both)
#adding age in addition to sex doesn't further improve model fit


#---SPMS---#
Data_LMM_SPMS$Geschlecht<-as.factor(Data_LMM_SPMS$Geschlecht)
LMM_RNFL_SP_total <- lmer(RNFL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_SPMS)
summary(LMM_RNFL_SP_total)

#a: age
LMM_RNFL_SP_total_age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_SPMS)
summary(LMM_RNFL_SP_total_age)
anova(LMM_RNFL_SP_total, LMM_RNFL_SP_total_age)
#age sign. improves model fit
confint(LMM_RNFL_SP_total_age)
#best model fit - presented in paper

#b: sex
LMM_RNFL_SP_total_sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_SPMS)
summary(LMM_RNFL_SP_total_sex)
anova(LMM_RNFL_SP_total, LMM_RNFL_SP_total_sex)
#sex doesn't improve model fit

#c: both
LMM_RNFL_SP_total_both <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_SPMS)
summary(LMM_RNFL_SP_total_both)
anova(LMM_RNFL_SP_total_both, LMM_RNFL_SP_total_age)
#adding sex in addition to age doesn't further improve model 


#--------------------------------------INL------------------------------------#    
#---RRMS---#
Data_LMM_RRMS$Geschlecht<-as.factor(Data_LMM_RRMS$Geschlecht)
LMM_INL_RR_total <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_RRMS)
summary(LMM_INL_RR_total)
confint(LMM_INL_RR_total)

#a: age
LMM_INL_RR_total_age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_RRMS)
summary(LMM_INL_RR_total_age)
anova(LMM_INL_RR_total_age, LMM_INL_RR_total)
#age doesn't sign. improve model fit

#b: sex
LMM_INL_RR_total_sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_RRMS)
summary(LMM_INL_RR_total_sex)
anova(LMM_INL_RR_total_sex, LMM_INL_RR_total)
#adding sex sign. improves model fit
confint(LMM_INL_RR_total_sex)
#best model fit - presented in paper

#c: both
LMM_INL_RR_total_both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_RRMS)
summary(LMM_INL_RR_total_both)
anova(LMM_INL_RR_total_both, LMM_INL_RR_total_sex)
#adding age in addition to sex doesn't further improve model fit

#---PPMS---#
Data_LMM_PPMS$Geschlecht<-as.factor(Data_LMM_PPMS$Geschlecht)
LMM_INL_PP_total <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_PPMS)
summary(LMM_INL_PP_total)
confint(LMM_INL_PP_total)
#best model fit - presented in paper

#a: age
LMM_INL_PP_total_age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_PPMS)
summary(LMM_INL_PP_total_age)
anova(LMM_INL_PP_total_age, LMM_INL_PP_total)
#age doesn't sign. improve model fit

#b: sex
LMM_INL_PP_total_sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_PPMS)
summary(LMM_INL_PP_total_sex)
anova(LMM_INL_PP_total_sex, LMM_INL_PP_total)
#sex doesn't sign. improve model fit

#c: both
LMM_INL_PP_total_both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_PPMS)
summary(LMM_INL_PP_total_both)
anova(LMM_INL_PP_total_both, LMM_INL_PP_total)
#adding both sex and age doesn't sign. improve model fit

#---SPMS---#
Data_LMM_SPMS$Geschlecht<-as.factor(Data_LMM_SPMS$Geschlecht)
LMM_INL_SP_total <- lmer(INL_new~Disease_Duration_FINAL + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_SPMS)
summary(LMM_INL_SP_total)
confint(LMM_INL_SP_total)
#best model fit - presented in paper

#a: age
LMM_INL_SP_total_age <- lmer(INL_new~Disease_Duration_FINAL + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_SPMS)
summary(LMM_INL_SP_total_age)
anova(LMM_INL_SP_total_age, LMM_INL_SP_total)
#age doesn't sign. improve model fit

#b: sex
LMM_INL_SP_total_sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_SPMS)
summary(LMM_INL_SP_total_sex)
anova(LMM_INL_SP_total_sex, LMM_INL_SP_total)
#sex doesn't sign. improve model fit

#c: both
LMM_INL_SP_total_both <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + Age + (1|ID/ID_Eye) + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_SPMS)
summary(LMM_INL_SP_total_both)
anova(LMM_INL_SP_total_both, LMM_INL_SP_total)
#adding both sex and age doesn't sign. improve model fit 
