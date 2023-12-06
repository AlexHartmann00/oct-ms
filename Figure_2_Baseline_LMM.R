
#----------------Linear mixed-effects models (LMM) at baseline---------------#

# This script computes LMMs analyzing the influence of disease duration on retinal 
# layer thickness (globalpRNFL, macular RNFL, GCIPL, INL), visual acuity, and VEP 
# latency at baseline. To account for confounding influences, age and/or sex were
# subsequently added to the model if the model fit could be improved.
# The results are presented in Figure 2. 
# Usual runtime ~= 2 minutes


#---------------------------------R packages---------------------------------#  

library(lme4)
library(lmerTest) #needed for p-values
library(haven)
library(performance)
library(car)
library(openxlsx)
library(model_performance)


#-----------------------Data loading & data management-----------------------# 

OCT_data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)

#Correct disease duration, if VEP dates do not match OCT dates
library(lubridate)
OCT_data$Disease_Duration_Latency <- as.numeric((as_date(OCT_data$VEP) - as_date(OCT_data$OCTDatum))/365 + OCT_data$Disease_Duration_FINAL)
OCT_data$Age_Latency <- as.numeric((as_date(OCT_data$VEP) - as_date(OCT_data$OCTDatum))/365 + OCT_data$Age)


Data_Baseline <- subset(OCT_data, OCT_included==1 & Sequenz==1)
Data_Baseline$Verlaufsform <- as.factor(Data_Baseline$Verlaufsform)
Data_Baseline$Geschlecht <- as.factor(Data_Baseline$Geschlecht)


Data_LMM_RRMS_BL <- subset(OCT_data, OCT_included==1 & Sequenz==1 &  Verlaufsform==1)
Data_LMM_RRMS_BL$Geschlecht <- as.factor(Data_LMM_RRMS_BL$Geschlecht)
Data_LMM_PPMS_BL <- subset(OCT_data, OCT_included==1 & Sequenz==1 & Verlaufsform==2)
Data_LMM_PPMS_BL$Geschlecht <- as.factor(Data_LMM_PPMS_BL$Geschlecht)
Data_LMM_SPMS_BL <- subset(OCT_data, OCT_included==1 & Sequenz==1 & Verlaufsform==3)
Data_LMM_SPMS_BL$Geschlecht <- as.factor(Data_LMM_SPMS_BL$Geschlecht)

#Sequenz=1 --> Baseline (BL)
#Verlaufsform=1 --> RRMS
#Verlaufsform=2 --> PPMS
#Verlaufsform=3 --> SPMS
#OCT_included=1 --> only OCT measurements included in the study

#One random effect for disease duration per person, not per eye, because disease 
#duration doesn't differ between eyes but only between subjects.
#Random intercept not necessary, because intercept would represent 
#the baseline values, which are the only values included in the model. 


#---------------------------------globalpRNFL---------------------------------#    
#---RRMS---#
LMM_globalpRNFL_DD_RR <- lmer(globalpRNFL~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_RRMS_BL)
summary(LMM_globalpRNFL_DD_RR)
model_performance(LMM_globalpRNFL_DD_RR)
#best model fit - presented in paper

#add age
LMM_globalpRNFL_DD_RR_age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_RRMS_BL)
summary(LMM_globalpRNFL_DD_RR_age)
anova(LMM_globalpRNFL_DD_RR, LMM_globalpRNFL_DD_RR_age)
#age does not improve model fit

#add sex
LMM_globalpRNFL_DD_RR_sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_RRMS_BL)
summary(LMM_globalpRNFL_DD_RR_sex)
anova(LMM_globalpRNFL_DD_RR, LMM_globalpRNFL_DD_RR_sex)
#sex does not improve model fit

#add age+sex
LMM_globalpRNFL_DD_RR_both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                                data=Data_LMM_RRMS_BL)
summary(LMM_globalpRNFL_DD_RR_both)
anova(LMM_globalpRNFL_DD_RR_both, LMM_globalpRNFL_DD_RR)
#adding both age and sex doesn't improve model fit

#---PPMS---#
LMM_globalpRNFL_DD_PP <- lmer(globalpRNFL~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_PPMS_BL)
summary(LMM_globalpRNFL_DD_PP)
model_performance(LMM_globalpRNFL_DD_PP)

#add age
LMM_globalpRNFL_DD_PP_age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                                  data=Data_LMM_PPMS_BL)
summary(LMM_globalpRNFL_DD_PP_age)
anova(LMM_globalpRNFL_DD_PP, LMM_globalpRNFL_DD_PP_age)
#age doesn't improve model fit

#add sex
LMM_globalpRNFL_DD_PP_sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                                  data=Data_LMM_PPMS_BL)
summary(LMM_globalpRNFL_DD_PP_sex)
anova(LMM_globalpRNFL_DD_PP, LMM_globalpRNFL_DD_PP_sex)
#sex sign. improves the model fit
#best model fit - presented in paper

#add age+sex
LMM_globalpRNFL_DD_PP_both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                                   data=Data_LMM_PPMS_BL)
summary(LMM_globalpRNFL_DD_PP_both)
anova(LMM_globalpRNFL_DD_PP_sex, LMM_globalpRNFL_DD_PP_both)
#adding age doesn't further improve the model compared to the model with sex only

#---SPMS---#
LMM_globalpRNFL_DD_SP <- lmer(globalpRNFL~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_SPMS_BL)
summary(LMM_globalpRNFL_DD_SP)
model_performance(LMM_globalpRNFL_DD_SP)
#best model fit - presented in paper

#add age
LMM_globalpRNFL_DD_SP_age <- lmer(globalpRNFL~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                                  data=Data_LMM_SPMS_BL)
summary(LMM_globalpRNFL_DD_SP_age)
anova(LMM_globalpRNFL_DD_SP, LMM_globalpRNFL_DD_SP_age)
#age doesn't improve model fit

#add sex
LMM_globalpRNFL_DD_SP_sex <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                                  data=Data_LMM_SPMS_BL)
summary(LMM_globalpRNFL_DD_SP_sex)
anova(LMM_globalpRNFL_DD_SP_sex, LMM_globalpRNFL_DD_SP)
#sex doesn't improve model fit

#add sex+age
LMM_globalpRNFL_DD_SP_both <- lmer(globalpRNFL~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                                   data=Data_LMM_SPMS_BL)
summary(LMM_globalpRNFL_DD_SP_both)
anova(LMM_globalpRNFL_DD_SP_both, LMM_globalpRNFL_DD_SP)
#adding both age and sex doesn't improve model fit


#-----------------------------------GCIPL------------------------------------# 
#---RRMS---#
LMM_GCIPL_DD_RR <- lmer(GCIPL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                        data=Data_LMM_RRMS_BL)
summary(LMM_GCIPL_DD_RR)
model_performance(LMM_GCIPL_DD_RR)
#best model fit - presented in paper

#add age
LMM_GCIPL_DD_RR_age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_RRMS_BL)
summary(LMM_GCIPL_DD_RR_age)
anova(LMM_GCIPL_DD_RR_age, LMM_GCIPL_DD_RR)
#age doesn't improve the model fit

#add sex
LMM_GCIPL_DD_RR_sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_RRMS_BL)
summary(LMM_GCIPL_DD_RR_sex)
anova(LMM_GCIPL_DD_RR_sex, LMM_GCIPL_DD_RR)
#sex doesn't improve model fit

#add sex+age
LMM_GCIPL_DD_RR_both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_RRMS_BL)
summary(LMM_GCIPL_DD_RR_both)
anova(LMM_GCIPL_DD_RR_both, LMM_GCIPL_DD_RR)
#adding both sex+age doesn't improve model fit

#---PPMS---#
LMM_GCIPL_DD_PP <- lmer(GCIPL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                        data=Data_LMM_PPMS_BL)
summary(LMM_GCIPL_DD_PP)
model_performance(LMM_GCIPL_DD_PP)

#add age
LMM_GCIPL_DD_PP_age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                            data=Data_LMM_PPMS_BL)
summary(LMM_GCIPL_DD_PP_age)
anova(LMM_GCIPL_DD_PP_age, LMM_GCIPL_DD_PP)
#age significantly improves model
#best model fit - presented in paper

#add sex
LMM_GCIPL_DD_PP_sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                            data=Data_LMM_PPMS_BL)
summary(LMM_GCIPL_DD_PP_sex)
anova(LMM_GCIPL_DD_PP_sex, LMM_GCIPL_DD_PP)
#sex doesn't improve model

#add sex+age
LMM_GCIPL_DD_PP_both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_PPMS_BL)
summary(LMM_GCIPL_DD_PP_both)
anova(LMM_GCIPL_DD_PP_both, LMM_GCIPL_DD_PP_age)
#adding sex to the model which already includes age doesn't further improve the model fit
#model including age is the best fitting model

#---SPMS---#
LMM_GCIPL_DD_SP <- lmer(GCIPL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                        data=Data_LMM_SPMS_BL)
summary(LMM_GCIPL_DD_SP)
model_performance(LMM_GCIPL_DD_SP)
#best model fit - presented in paper

#add age
LMM_GCIPL_DD_SP_age <- lmer(GCIPL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                            data=Data_LMM_SPMS_BL)
summary(LMM_GCIPL_DD_SP_age)
anova(LMM_GCIPL_DD_SP_age, LMM_GCIPL_DD_SP)
#age doesn't improve model fit

#add sex
LMM_GCIPL_DD_SP_sex <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                            data=Data_LMM_SPMS_BL)
summary(LMM_GCIPL_DD_SP_sex)
anova(LMM_GCIPL_DD_SP_sex, LMM_GCIPL_DD_SP)
#sex doesn't improve model fit

#add both
LMM_GCIPL_DD_SP_both <- lmer(GCIPL_new~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_SPMS_BL)
summary(LMM_GCIPL_DD_SP_both)
anova(LMM_GCIPL_DD_SP_both, LMM_GCIPL_DD_SP)
#adding sex + age doesn't improve the model fit


#------------------------------------INL-------------------------------------# 
#---RRMS---#
LMM_INL_DD_RR <- lmer(INL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                      data=Data_LMM_RRMS_BL)
summary(LMM_INL_DD_RR)
model_performance(LMM_INL_DD_RR)

#add age
LMM_INL_DD_RR_age <- lmer(INL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                       data=Data_LMM_RRMS_BL)
summary(LMM_INL_DD_RR_age)
anova(LMM_INL_DD_RR_age,LMM_INL_DD_RR)
#age doesn't improve model fit

#add sex
LMM_INL_DD_RR_sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                       data=Data_LMM_RRMS_BL)
summary(LMM_INL_DD_RR_sex)
anova(LMM_INL_DD_RR_sex,LMM_INL_DD_RR)
#sex sign. improves model fit.
#best model fit - presented in paper

#add age+sex
LMM_INL_DD_RR_both <- lmer(INL_new~Disease_Duration_FINAL + Age + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                        data=Data_LMM_RRMS_BL)
summary(LMM_INL_DD_RR_both)
anova(LMM_INL_DD_RR_sex, LMM_INL_DD_RR_both)
#adding age in addition to sex doesn't further improve model fit

#---PPMS---#
LMM_INL_DD_PP <- lmer(INL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                      data=Data_LMM_PPMS_BL)
summary(LMM_INL_DD_PP)
model_performance(LMM_INL_DD_PP)

#add age
LMM_INL_DD_PP_age <- lmer(INL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_PPMS_BL)
summary(LMM_INL_DD_PP_age)
anova(LMM_INL_DD_PP_age, LMM_INL_DD_PP)
#age sign. improves model fit
#best model fit - presented in paper

#add sex
LMM_INL_DD_PP_sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_PPMS_BL)
summary(LMM_INL_DD_PP_sex)
#sex sign. improves model fit 

#add age+sex
LMM_INL_DD_PP_both <- lmer(INL_new~Disease_Duration_FINAL + Age + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_PPMS_BL)
summary(LMM_INL_DD_PP_both)
anova(LMM_INL_DD_PP_both, LMM_INL_DD_PP)
#adding both sex and age sign. improves model fit 
#check if age or sex alone is enough
anova(LMM_INL_DD_PP_both, LMM_INL_DD_PP_sex)
anova(LMM_INL_DD_PP_both, LMM_INL_DD_PP_age)
#age is enough, sex doesn't further improve model 

#---SPMS---#
LMM_INL_DD_SP <- lmer(INL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                      data=Data_LMM_SPMS_BL)
summary(LMM_INL_DD_SP)
model_performance(LMM_INL_DD_SP)

#add age
LMM_INL_DD_SP_age <- lmer(INL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_SPMS_BL)
summary(LMM_INL_DD_SP_age)
anova(LMM_INL_DD_SP_age, LMM_INL_DD_SP)
#age doesn't improve model fit

#add sex
LMM_INL_DD_SP_sex <- lmer(INL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_SPMS_BL)
summary(LMM_INL_DD_SP_sex)
anova(LMM_INL_DD_SP_sex, LMM_INL_DD_SP)
#sex improves model fit
#best model fit - presented in paper

#add both
LMM_INL_DD_SP_both <- lmer(INL_new~Disease_Duration_FINAL + Age + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_SPMS_BL)
summary(LMM_INL_DD_SP_both)
anova(LMM_INL_DD_SP_sex, LMM_INL_DD_SP_both)
#adding age to the model already including sex doesn't improve model fit


#------------------------------------mRNFL-------------------------------------# 
#---RRMS---#
LMM_RNFL_DD_RR <- lmer(RNFL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                       data=Data_LMM_RRMS_BL)
summary(LMM_RNFL_DD_RR)
model_performance(LMM_RNFL_DD_RR)
#best model fit - presented in paper

#add age
LMM_RNFL_DD_RR_age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                        data=Data_LMM_RRMS_BL)
summary(LMM_RNFL_DD_RR_age)
anova(LMM_RNFL_DD_RR_age, LMM_RNFL_DD_RR)
#age doesn't improve model fi

#add sex
LMM_RNFL_DD_RR_sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                        data=Data_LMM_RRMS_BL)
summary(LMM_RNFL_DD_RR_sex)
anova(LMM_RNFL_DD_RR_sex, LMM_RNFL_DD_RR)
#sex doesn't improve model fit

#add both
LMM_RNFL_DD_RR_both <- lmer(RNFL_new~Disease_Duration_FINAL + Age + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_RRMS_BL)
summary(LMM_RNFL_DD_RR_both)
anova(LMM_RNFL_DD_RR_both, LMM_RNFL_DD_RR)
#adding both sex and age doesn't improve model fit

#---PPMS---#
LMM_RNFL_DD_PP <- lmer(RNFL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                       data=Data_LMM_PPMS_BL)
summary(LMM_RNFL_DD_PP)
model_performance(LMM_RNFL_DD_PP)

#add age
LMM_RNFL_DD_PP_age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_PPMS_BL)
summary(LMM_RNFL_DD_PP_age)
anova(LMM_RNFL_DD_PP_age, LMM_RNFL_DD_PP)
#adding age doesn't improve the model

#add sex
LMM_RNFL_DD_PP_sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_PPMS_BL)
summary(LMM_RNFL_DD_PP_sex)
anova(LMM_RNFL_DD_PP, LMM_RNFL_DD_PP_sex)
#sex improves the model fit
#best model fit - presented in paper

#add both
LMM_RNFL_DD_PP_both <- lmer(RNFL_new~Disease_Duration_FINAL + Age + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                            data=Data_LMM_PPMS_BL)
summary(LMM_RNFL_DD_PP_both)
anova(LMM_RNFL_DD_PP_both, LMM_RNFL_DD_PP_sex)
#adding age in addition to sex doesn't further improve the model fit

#---SPMS---#
LMM_RNFL_DD_SP <- lmer(RNFL_new~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                       data=Data_LMM_SPMS_BL)
summary(LMM_RNFL_DD_SP)
model_performance(LMM_RNFL_DD_SP)
#best model fit - presented in paper

#add age
LMM_RNFL_DD_SP_age <- lmer(RNFL_new~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_SPMS_BL)
summary(LMM_RNFL_DD_SP_age)
anova(LMM_RNFL_DD_SP_age, LMM_RNFL_DD_SP)
#age doesn't improve model fit

#add sex
LMM_RNFL_DD_SP_sex <- lmer(RNFL_new~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                           data=Data_LMM_SPMS_BL)
summary(LMM_RNFL_DD_SP_sex)
anova(LMM_RNFL_DD_SP_sex, LMM_RNFL_DD_SP)
#sex doesn't improve model fit

#add both
LMM_RNFL_DD_SP_both <- lmer(RNFL_new~Disease_Duration_FINAL + Age + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                            data=Data_LMM_SPMS_BL)
summary(LMM_RNFL_DD_SP_both)
anova(LMM_RNFL_DD_SP_both, LMM_RNFL_DD_SP)
#adding both sex and age doesn't improve model fit.


#--------------------------------VEP Latency---------------------------------#
#---RRMS---#
LMM_Latenz_DD_RR <- lmer(Latenz~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_RRMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_RR)
#best model fit - presented in paper

#add age
LMM_Latenz_DD_RR_age <- lmer(Latenz~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_RRMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_RR_age)
anova(LMM_Latenz_DD_RR_age, LMM_Latenz_DD_RR)
#adding age does not improve the model fit

#add sex
LMM_Latenz_DD_RR_sex <- lmer(Latenz~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_RRMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_RR_sex)
anova(LMM_Latenz_DD_RR_sex, LMM_Latenz_DD_RR)
#adding sex does not improve model fit

#add sex+age
LMM_Latenz_DD_RR_both <- lmer(Latenz~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_RRMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_RR_both)
anova(LMM_Latenz_DD_RR_both, LMM_Latenz_DD_RR)
#adding both sex and age doesn't improve model fit

#---PPMS---# 
LMM_Latenz_DD_PP <- lmer(Latenz~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_PPMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_PP)
#best model fit - presented in paper

#add age
LMM_Latenz_DD_PP_age <- lmer(Latenz~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_PPMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_PP_age)
anova(LMM_Latenz_DD_PP_age, LMM_Latenz_DD_PP)
#adding age doesn't improve model fit

#add sex
LMM_Latenz_DD_PP_sex <- lmer(Latenz~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_PPMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_PP_sex)
anova(LMM_Latenz_DD_PP_sex, LMM_Latenz_DD_PP)
#adding sex doesn't improve model fit

#add age+sex
LMM_Latenz_DD_PP_both <- lmer(Latenz~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_PPMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_PP_both)
anova(LMM_Latenz_DD_PP_both, LMM_Latenz_DD_PP)
#adding both sex and age doesn't improve model fit

#---SPMS---# 
LMM_Latenz_DD_SP <- lmer(Latenz~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                         data=Data_LMM_SPMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_SP)
#best model fit - presented in paper

#add age
LMM_Latenz_DD_SP_age <- lmer(Latenz~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_SPMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_SP_age)
anova(LMM_Latenz_DD_SP_age, LMM_Latenz_DD_SP)
#adding age doesn't improve model fit

#add sex
LMM_Latenz_DD_SP_sex <- lmer(Latenz~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                             data=Data_LMM_SPMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_SP_sex)
anova(LMM_Latenz_DD_SP_sex, LMM_Latenz_DD_SP)
#adding sex doesn't improve model fit

#add age+sex
LMM_Latenz_DD_SP_both <- lmer(Latenz~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_SPMS_BL,subset = VEPMatch == 1)
summary(LMM_Latenz_DD_SP_both)
anova(LMM_Latenz_DD_SP_both, LMM_Latenz_DD_SP)
#adding both sex and age doesn't improve model fit


#-----------------------------------Visus ------------------------------------#
#---RRMS---#
LMM_VisusHC_DD_RR <- lmer(VisusHC~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_RRMS_BL)
summary(LMM_VisusHC_DD_RR)

#add age
LMM_VisusHC_DD_RR_age <- lmer(VisusHC~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_RRMS_BL)
summary(LMM_VisusHC_DD_RR_age)
anova(LMM_VisusHC_DD_RR,LMM_VisusHC_DD_RR_age)
#adding age doesn't improve model fit

#add sex
LMM_VisusHC_DD_RR_sex <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_RRMS_BL)
summary(LMM_VisusHC_DD_RR_sex)
anova(LMM_VisusHC_DD_RR_sex, LMM_VisusHC_DD_RR)
#adding sex improves model fit
#best model fit - presented in paper

#add age+sex
LMM_VisusHC_DD_RR_both <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_RRMS_BL)
summary(LMM_VisusHC_DD_RR_both)
anova(LMM_VisusHC_DD_RR_both, LMM_VisusHC_DD_RR_sex)
#adding age in addition to sex doesn't further improve model fit

#---PPMS---#
LMM_VisusHC_DD_PP <- lmer(VisusHC~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_PPMS_BL)
summary(LMM_VisusHC_DD_PP)
#best model fit - presented in paper

#add age
LMM_VisusHC_DD_PP_age <- lmer(VisusHC~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_PPMS_BL)
summary(LMM_VisusHC_DD_PP_age)
anova(LMM_VisusHC_DD_PP_age, LMM_VisusHC_DD_PP)
#adding age doesn't improve model fit

#add sex
LMM_VisusHC_DD_PP_sex <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_PPMS_BL)
summary(LMM_VisusHC_DD_PP_sex)
anova(LMM_VisusHC_DD_PP_sex, LMM_VisusHC_DD_PP)
#adding sex doesn't improve model fit

#add age+sex
LMM_VisusHC_DD_PP_both <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_PPMS_BL)
summary(LMM_VisusHC_DD_PP_both)
anova(LMM_VisusHC_DD_PP_both, LMM_VisusHC_DD_PP)
#adding both sex and age doesn't improve model fit

#---SPMS---#
LMM_VisusHC_DD_SP <- lmer(VisusHC~Disease_Duration_FINAL + (0+Disease_Duration_FINAL|ID), 
                          data=Data_LMM_SPMS_BL)
summary(LMM_VisusHC_DD_SP)
#best model fit - presented in paper

#add age
LMM_VisusHC_DD_SP_age <- lmer(VisusHC~Disease_Duration_FINAL + Age + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_SPMS_BL)
summary(LMM_VisusHC_DD_SP_age)
anova(LMM_VisusHC_DD_SP_age, LMM_VisusHC_DD_SP)
#adding age doesn't improve model fit

#add sex
LMM_VisusHC_DD_SP_sex <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + (0+Disease_Duration_FINAL|ID), 
                              data=Data_LMM_SPMS_BL)
summary(LMM_VisusHC_DD_SP_sex)
anova(LMM_VisusHC_DD_SP_sex, LMM_VisusHC_DD_SP)
#adding sex doesn't improve model fit

#add age+sex
LMM_VisusHC_DD_SP_both <- lmer(VisusHC~Disease_Duration_FINAL + Geschlecht + Age + (0+Disease_Duration_FINAL|ID), 
                               data=Data_LMM_SPMS_BL)
summary(LMM_VisusHC_DD_SP_both)
anova(LMM_VisusHC_DD_SP_both, LMM_VisusHC_DD_SP)
#adding both sex and age doesn't improve model fit