
#----------------Longitudinal linear mixed-effects models (LMM)---------------#

# This script computes LMMs analyzing the influence of retinal layer thickness 
# (globalpRNFL, macular RNFL, GCIPL, INL) on visual acuity and VEP 
# latency. To account for confounding influences, age and/or sex 
# were subsequently added to the model if the model fit could be improved. 
# The results are presented in Figure 4. 
# Usual runtime ~= 2 minutes 


#---------------------------------R packages---------------------------------# 

library(lme4)
library(lmerTest) #needed for p-values
library(haven)
library(performance)
library(car)


#-----------------------Data loading & data management-----------------------# 

OCT_data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)
Data_Visus <- subset(OCT_data, OCT_included==1)
Data_Visus_HC <- subset(Data_Visus, Verlaufsform==0)
Data_Visus_RRMS <- subset(Data_Visus, Verlaufsform==1)
Data_Visus_PPMS <- subset(Data_Visus, Verlaufsform==2)
Data_Visus_SPMS <- subset(Data_Visus, Verlaufsform==3)
Data_Visus_HC$Geschlecht <- as.factor(Data_Visus_HC$Geschlecht)
Data_Visus_RRMS$Geschlecht <- as.factor(Data_Visus_RRMS$Geschlecht)
Data_Visus_PPMS$Geschlecht <- as.factor(Data_Visus_PPMS$Geschlecht)
Data_Visus_SPMS$Geschlecht <- as.factor(Data_Visus_SPMS$Geschlecht)
#Verlaufsform=0 --> HC
#Verlaufsform=1 --> RRMS
#Verlaufsform=2 --> PPMS
#Verlaufsform=3 --> SPMS


#--------------------Visus HC (high contrast visual acuity)-------------------#
#------------------------------------pRNFL------------------------------------# 
#---RRMS---#
LMM_pRNFL_RR <- lmer(VisusHC~globalpRNFL + (globalpRNFL|ID/ID_Eye), 
                     data=Data_Visus_RRMS)
summary(LMM_pRNFL_RR)
#boundary (singular) fit: see help('isSingular')
isSingular(LMM_pRNFL_RR)
#singular fit --> reduce to less complex model 

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 9.135e-02 0.302244      
#globalpRNFL 1.129e-05 0.003359 -1.00
#ID        (Intercept) 3.065e-01 0.553620      
#globalpRNFL 3.144e-05 0.005607 -0.99
#Residual              3.598e-03 0.059981

#variation in baseline pRNFL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of pRNFL = RE with lowest variance

LMM_pRNFL_RR1 <- lmer(VisusHC~globalpRNFL + (1|ID/ID_Eye), 
                      data=Data_Visus_RRMS)
summary(LMM_pRNFL_RR1)
#best model fit - presented in paper

#add age
LMM_pRNFL_RR1_Age <- lmer(VisusHC~globalpRNFL + Age + (1|ID/ID_Eye), 
                          data=Data_Visus_RRMS)
summary(LMM_pRNFL_RR1_Age)
anova(LMM_pRNFL_RR1_Age, LMM_pRNFL_RR1)
#age doesn't sign. improve model fit

#add sex
LMM_pRNFL_RR1_sex <- lmer(VisusHC~globalpRNFL + Geschlecht + (1|ID/ID_Eye), 
                          data=Data_Visus_RRMS)
summary(LMM_pRNFL_RR1_sex)
anova(LMM_pRNFL_RR1_sex, LMM_pRNFL_RR1)
#sex doesn't sign. improve model fit

#add both
LMM_pRNFL_RR_both <- lmer(VisusHC~globalpRNFL + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Visus_RRMS)
summary(LMM_pRNFL_RR_both)
anova(LMM_pRNFL_RR_both, LMM_pRNFL_RR1)
#adding both sex and age doesn't improve model fit 


#---PPMS---#
LMM_pRNFL_PP <- lmer(VisusHC~globalpRNFL + (globalpRNFL|ID/ID_Eye), 
                     data=Data_Visus_PPMS)
summary(LMM_pRNFL_PP)
#Error: number of observations (=315) <= number of random effects (=318) for term (globalpRNFL | ID_Eye:ID);
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#reduce to less complex model without random slop vor pRNFL

LMM_pRNFL_PP <- lmer(VisusHC~globalpRNFL + (1|ID/ID_Eye), 
                     data=Data_Visus_PPMS)
summary(LMM_pRNFL_PP)
#best model fit - presented in paper

#add age
LMM_pRNFL_PP_age <- lmer(VisusHC~globalpRNFL + Age + (1|ID/ID_Eye), 
                         data=Data_Visus_PPMS)
summary(LMM_pRNFL_PP_age)
anova(LMM_pRNFL_PP_age, LMM_pRNFL_PP)
#adding age doesn't improve model fit

#add sex
LMM_pRNFL_PP_sex <- lmer(VisusHC~globalpRNFL + Geschlecht + (1|ID/ID_Eye), 
                         data=Data_Visus_PPMS)
summary(LMM_pRNFL_PP_sex)
anova(LMM_pRNFL_PP_sex, LMM_pRNFL_PP)
#adding sex doesn't improve model fit

#add both
LMM_pRNFL_PP_both <- lmer(VisusHC~globalpRNFL + Geschlecht + Age +(1|ID/ID_Eye), 
                          data=Data_Visus_PPMS)
summary(LMM_pRNFL_PP_both)
anova(LMM_pRNFL_PP_both, LMM_pRNFL_PP)
#adding both sex and age doesn't improve model fit


#---SPMS---#
LMM_pRNFL_SP <- lmer(VisusHC~globalpRNFL + (globalpRNFL|ID/ID_Eye), 
                     data=Data_Visus_SPMS)
summary(LMM_pRNFL_SP)
#Error: number of observations (=236) <= number of random effects (=252) for term (globalpRNFL | ID_Eye:ID); 
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#reduce to less complex model without random slop vor pRNFL

LMM_pRNFL_SP <- lmer(VisusHC~globalpRNFL + (1|ID/ID_Eye), 
                     data=Data_Visus_SPMS)
summary(LMM_pRNFL_SP)
#best model fit - presented in paper

#add age
LMM_pRNFL_SP_age <- lmer(VisusHC~globalpRNFL + Age + (1|ID/ID_Eye), 
                         data=Data_Visus_SPMS)
summary(LMM_pRNFL_SP_age)
anova(LMM_pRNFL_SP_age, LMM_pRNFL_SP)
#adding age doesn't sign. improve model fit

#add sex
LMM_pRNFL_SP_sex <- lmer(VisusHC~globalpRNFL + Geschlecht +  (1|ID/ID_Eye), 
                         data=Data_Visus_SPMS)
summary(LMM_pRNFL_SP_sex)
anova(LMM_pRNFL_SP_sex, LMM_pRNFL_SP)
#adding sex doesn't sign. improve model fit

#add both
LMM_pRNFL_SP_both <- lmer(VisusHC~globalpRNFL + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Visus_SPMS)
summary(LMM_pRNFL_SP_both)
anova(LMM_pRNFL_SP_both, LMM_pRNFL_SP)
#adding both sex and age doesn't improve model fit


#------------------------------------GCIPL------------------------------------#
#---HC---#
LMM_GCIPL_HC <- lmer(VisusHC~GCIPL_new + (GCIPL_new|ID/ID_Eye), 
                     data=Data_Visus_HC)
summary(LMM_GCIPL_HC)
#Error: number of observations (=222) <= number of random effects (=226) for term (GCIPL_new | ID_Eye:ID);
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#reduce to less complex model without random slop vor GCIPL

LMM_GCIPL_HC1 <- lmer(VisusHC~GCIPL_new + (1|ID/ID_Eye), 
                      data=Data_Visus_HC)
summary(LMM_GCIPL_HC1)
#Model failed to converge with 1 negative eigenvalue: -7.3e-01 

#Random effects:
#Groups    Name        Variance  Std.Dev. 
#ID_Eye:ID (Intercept) 1.365e-10 1.168e-05
#ID        (Intercept) 9.650e-03 9.823e-02
#Residual              4.995e-03 7.068e-02

#random intercept per eye nested in subject explains almost no variance --> reduce to less complex model

LMM_GCIPL_HC2 <- lmer(VisusHC~GCIPL_new + (1|ID_Eye), 
                      data=Data_Visus_HC)
summary(LMM_GCIPL_HC2)
#best model fit - presented in paper

#add age
LMM_GCIPL_HC2_Age <- lmer(VisusHC~GCIPL_new + Age + (1|ID), 
                          data=Data_Visus_HC)
summary(LMM_GCIPL_HC2_Age)
anova(LMM_GCIPL_HC2_Age, LMM_GCIPL_HC2)
#adding age doesn't improve model fit

#add sex
LMM_GCIPL_HC2_sex <- lmer(VisusHC~GCIPL_new + Geschlecht + (1|ID), 
                          data=Data_Visus_HC)
summary(LMM_GCIPL_HC2_sex)
anova(LMM_GCIPL_HC2_sex, LMM_GCIPL_HC2)
#adding sex doesn't improve model fit

#add both
LMM_GCIPL_HC2_both <- lmer(VisusHC~GCIPL_new + Geschlecht + Age +  (1|ID), 
                          data=Data_Visus_HC)
summary(LMM_GCIPL_HC2_both)
anova(LMM_GCIPL_HC2_both, LMM_GCIPL_HC2)
#adding both sex and age doesn't improve model fit

#---RRMS---#
LMM_GCIPL_RR <- lmer(VisusHC~GCIPL_new + (GCIPL_new|ID/ID_Eye), 
                     data=Data_Visus_RRMS)
summary(LMM_GCIPL_RR)
#Warning messages:
#1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#unable to evaluate scaled gradient
#2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
#3: Model failed to converge with 1 negative eigenvalue: -8.2e-01

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 4.291e-01 0.655054      
#GCIPL_new   7.845e-05 0.008857 -1.00
#ID        (Intercept) 3.186e-03 0.056446      
#GCIPL_new   7.413e-07 0.000861 0.20 
#Residual              3.600e-03 0.060004     

#variation in baseline GCIPL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of GCIPL = RE with lowest variance

LMM_GCIPL_RR1 <- lmer(VisusHC~GCIPL_new + (1|ID/ID_Eye), 
                      data=Data_Visus_RRMS)
summary(LMM_GCIPL_RR1)
#best model fit - presented in paper

#add age
LMM_GCIPL_RR1_Age <- lmer(VisusHC~GCIPL_new + Age + (1|ID/ID_Eye), 
                          data=Data_Visus_RRMS)
summary(LMM_GCIPL_RR1_Age)
anova(LMM_GCIPL_RR1_Age, LMM_GCIPL_RR1)
#adding age doesn't sign. improve model fit

#add sex
LMM_GCIPL_RR1_sex <- lmer(VisusHC~GCIPL_new + Geschlecht + (1|ID/ID_Eye), 
                          data=Data_Visus_RRMS)
summary(LMM_GCIPL_RR1_sex)
anova(LMM_GCIPL_RR1_sex, LMM_GCIPL_RR1)
#adding sex doesn't sign. improve model fit

#add both
LMM_GCIPL_RR_both <- lmer(VisusHC~GCIPL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Visus_RRMS)
summary(LMM_GCIPL_RR_both)
anova(LMM_GCIPL_RR_both, LMM_GCIPL_RR1_Age)
#adding both sex and age doesn't sign. improve model fit

#---PPMS---#
LMM_GCIPL_PP <- lmer(VisusHC~GCIPL_new + (GCIPL_new|ID/ID_Eye), 
                     data=Data_Visus_PPMS)
summary(LMM_GCIPL_PP)
isSingular(LMM_GCIPL_PP)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 5.651e-02 0.237714      
#GCIPL_new   5.894e-06 0.002428 -1.00
#ID        (Intercept) 8.390e-02 0.289648      
#GCIPL_new   6.152e-06 0.002480 -1.00
#Residual              8.387e-03 0.091578 

#variation in baseline GCIPL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of GCIPL = RE with lowest variance

LMM_GCIPL_PP <- lmer(VisusHC~GCIPL_new + (1|ID/ID_Eye), 
                     data=Data_Visus_PPMS)
summary(LMM_GCIPL_PP)
#best model fit - presented in paper

#add age
LMM_GCIPL_PP_age <- lmer(VisusHC~GCIPL_new + Age + (1|ID/ID_Eye), 
                         data=Data_Visus_PPMS)
summary(LMM_GCIPL_PP_age)
anova(LMM_GCIPL_PP_age, LMM_GCIPL_PP)
#adding age doesn't improve model fit

#add sex
LMM_GCIPL_PP_sex <- lmer(VisusHC~GCIPL_new + Geschlecht + (1|ID/ID_Eye), 
                         data=Data_Visus_PPMS)
summary(LMM_GCIPL_PP_sex)
anova(LMM_GCIPL_PP_sex, LMM_GCIPL_PP)
#adding sex doesn't improve model fit

#add both
LMM_GCIPL_PP_both <- lmer(VisusHC~GCIPL_new + Geschlecht + Age +(1|ID/ID_Eye), 
                          data=Data_Visus_PPMS)
summary(LMM_GCIPL_PP_both)
anova(LMM_GCIPL_PP_both, LMM_GCIPL_PP)
#adding both sex and age doesn't improve model fit

#---SPMS---#
LMM_GCIPL_SP <- lmer(VisusHC~GCIPL_new + (GCIPL_new|ID/ID_Eye), 
                     data=Data_Visus_SPMS)
#Error: number of observations (=262) <= number of random effects (=270) for term (GCIPL_new | ID_Eye:ID);
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope 

LMM_GCIPL_SP <- lmer(VisusHC~GCIPL_new + (1|ID/ID_Eye), 
                     data=Data_Visus_SPMS)
summary(LMM_GCIPL_SP)
#best model fit - presented in paper

#add age
LMM_GCIPL_SP_age <- lmer(VisusHC~GCIPL_new + Age + (1|ID/ID_Eye), 
                         data=Data_Visus_SPMS)
summary(LMM_GCIPL_SP_age)
anova(LMM_GCIPL_SP_age, LMM_GCIPL_SP)
#adding age doesn't improve model fit

#add sex
LMM_GCIPL_SP_sex <- lmer(VisusHC~GCIPL_new + Geschlecht +  (1|ID/ID_Eye), 
                         data=Data_Visus_SPMS)
summary(LMM_GCIPL_SP_sex)
anova(LMM_GCIPL_SP_sex, LMM_GCIPL_SP)
#adding sex doesn't improve model fit

#add both
LMM_GCIPL_SP_both <- lmer(VisusHC~GCIPL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Visus_SPMS)
summary(LMM_GCIPL_SP_both)
anova(LMM_GCIPL_SP_both, LMM_GCIPL_SP)
#adding both sex and age doesn't improve model fit


#-------------------------------------INL-------------------------------------#
#---HC---#

LMM_INL_HC <- lmer(VisusHC~INL_new + (INL_new|ID/ID_Eye), 
                   data=Data_Visus_HC)
#Error: number of observations (=222) <= number of random effects (=226) for term (INL_new | ID_Eye:ID); 
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope (in line with approach for the other retinal layers)

LMM_INL_HC1 <- lmer(VisusHC~INL_new + (1|ID/ID_Eye), 
                    data=Data_Visus_HC)
summary(LMM_INL_HC1)
#best model fit - presented in paper

#add age
LMM_INL_HC1_Age <- lmer(VisusHC~INL_new + Age + (1|ID/ID_Eye), 
                        data=Data_Visus_HC)
summary(LMM_INL_HC1_Age)
anova(LMM_INL_HC1_Age, LMM_INL_HC1)
#adding age doesn't improve model fit

#add sex
LMM_INL_HC1_sex <- lmer(VisusHC~INL_new + Geschlecht + (1|ID/ID_Eye), 
                        data=Data_Visus_HC)
summary(LMM_INL_HC1_sex)
anova(LMM_INL_HC1_sex, LMM_INL_HC1)
#adding sex doesn't improve model fit

#add both
LMM_INL_HC_both <- lmer(VisusHC~INL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                        data=Data_Visus_HC)
summary(LMM_INL_HC_both)
anova(LMM_INL_HC_both, LMM_INL_HC1)
#adding both sex and age doesn't improve model fit

#---RRMS---#
LMM_INL_RR <- lmer(VisusHC~INL_new + (INL_new|ID/ID_Eye), 
                   data=Data_Visus_RRMS)
summary(LMM_INL_RR)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 4.599e-01 0.678166      
#INL_new               5.611e-04 0.023687 -1.00
#ID        (Intercept) 1.894e-02 0.137608      
#INL_new               3.634e-05 0.006028 -0.92
#Residual              3.545e-03 0.059539 

#variation in baseline INL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of INL = RE with lowest variance

LMM_INL_RR1 <- lmer(VisusHC~INL_new + (1|ID/ID_Eye), 
                    data=Data_Visus_RRMS)
summary(LMM_INL_RR1)
#best model fit - presented in paper

#add age
LMM_INL_RR1_Age <- lmer(VisusHC~INL_new + Age + (1|ID/ID_Eye), 
                        data=Data_Visus_RRMS)
summary(LMM_INL_RR1_Age)
anova(LMM_INL_RR1_Age, LMM_INL_RR1)
#adding age doesn't improve model fit

#add sex
LMM_INL_RR1_sex <- lmer(VisusHC~INL_new + Geschlecht + (1|ID/ID_Eye), 
                        data=Data_Visus_RRMS)
summary(LMM_INL_RR1_sex)
anova(LMM_INL_RR1_sex, LMM_INL_RR1)
#adding sex doesn't improve model fit

#add both
LMM_INL_RR_both <- lmer(VisusHC~INL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                        data=Data_Visus_RRMS)
summary(LMM_INL_RR_both)
anova(LMM_INL_RR_both, LMM_INL_RR1)
#adding both sex and age doesn't improve model fit

#---PPMS---#
LMM_INL_PP <- lmer(VisusHC~INL_new + (INL_new|ID/ID_Eye), 
                   data=Data_Visus_PPMS)
summary(LMM_INL_PP)
#Model failed to converge with 1 negative eigenvalue: -5.3e-02 

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 0.4987952 0.70625       
#INL_new               0.0003462 0.01861  -1.00
#ID        (Intercept) 0.0609043 0.24679       
#INL_new               0.0001340 0.01158  -1.00
#Residual              0.0079085 0.08893  

#variation in baseline INL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of INL = RE with lowest variance

LMM_INL_PP <- lmer(VisusHC~INL_new + (1|ID/ID_Eye), 
                   data=Data_Visus_PPMS)
summary(LMM_INL_PP)
#best model fit - presented in paper

#add age
LMM_INL_PP_age <- lmer(VisusHC~INL_new + Age + (1|ID/ID_Eye), 
                       data=Data_Visus_PPMS)
summary(LMM_INL_PP_age)
anova(LMM_INL_PP_age, LMM_INL_PP)
#adding age doesn't improve model fit

#add sex
LMM_INL_PP_sex <- lmer(VisusHC~INL_new + Geschlecht + (1|ID/ID_Eye), 
                       data=Data_Visus_PPMS)
summary(LMM_INL_PP_sex)
anova(LMM_INL_PP_sex, LMM_INL_PP)
#adding sex doesn't improve model fit

#add both
LMM_INL_PP_both <- lmer(VisusHC~INL_new + Geschlecht + Age +(1|ID/ID_Eye), 
                        data=Data_Visus_PPMS)
summary(LMM_INL_PP_both)
anova(LMM_INL_PP_both, LMM_INL_PP)
#adding both sex and age doesn't improve model fit

#---SPMS---#
LMM_INL_SP <- lmer(VisusHC~INL_new + (INL_new|ID/ID_Eye), 
                   data=Data_Visus_SPMS)
#Error: number of observations (=262) <= number of random effects (=270) for term (INL_new | ID_Eye:ID);
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope (in line with approach for the other retinal layers)

LMM_INL_SP <- lmer(VisusHC~INL_new + (1|ID/ID_Eye), 
                   data=Data_Visus_SPMS)
summary(LMM_INL_SP)
#best model fit - presented in paper

#add age
LMM_INL_SP_age <- lmer(VisusHC~INL_new + Age + (1|ID/ID_Eye), 
                       data=Data_Visus_SPMS)
summary(LMM_INL_SP_age)
anova(LMM_INL_SP_age, LMM_INL_SP)
#adding age doesn't improve model fit

#add sex
LMM_INL_SP_sex <- lmer(VisusHC~INL_new + Geschlecht +  (1|ID/ID_Eye), 
                       data=Data_Visus_SPMS)
summary(LMM_INL_SP_sex)
anova(LMM_INL_SP_sex, LMM_INL_SP)
#adding sex doesn't improve model fit

#add both
LMM_INL_SP_both <- lmer(VisusHC~INL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                        data=Data_Visus_SPMS)
summary(LMM_INL_SP_both)
anova(LMM_INL_SP_both, LMM_INL_SP)
#adding both sex and age doesn't improve model fit


#------------------------------------ mRNFL-----------------------------------#
#---HC---#
LMM_RNFL_HC <- lmer(VisusHC~RNFL_new + (RNFL_new|ID/ID_Eye), 
                     data=Data_Visus_HC)
summary(LMM_RNFL_HC)
#Error: number of observations (=222) <= number of random effects (=226) for term (RNFL_new | ID_Eye:ID); 
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
##--> reduce to less complex model without random slope (in line with approach for the other retinal layers)

LMM_RNFL_HC1 <- lmer(VisusHC~RNFL_new + (1|ID/ID_Eye), 
                     data=Data_Visus_HC)
summary(LMM_RNFL_HC1)
#  Model failed to converge with max|grad| = 0.00367405 (tol = 0.002, component 1)

#Random effects:
#Groups    Name        Variance  Std.Dev.
#ID_Eye:ID (Intercept) 3.260e-06 0.001806
#ID        (Intercept) 9.543e-03 0.097690
#Residual              4.905e-03 0.070036

#variation in baseline RNFL thickness per subject seems to be more important than variation per eye and subject
#--> eliminate random intercept per eye = RE with lowest variance

LMM_RNFL_HC2 <- lmer(VisusHC~RNFL_new + (1|ID), 
                     data=Data_Visus_HC)
summary(LMM_RNFL_HC2)
#best model fit - presented in paper

#add age
LMM_RNFL_HC1_Age <- lmer(VisusHC~RNFL_new + Age + (1|ID), 
                         data=Data_Visus_HC)
summary(LMM_RNFL_HC1_Age)
anova(LMM_RNFL_HC1_Age, LMM_RNFL_HC2)
#adding age doesn't improve model fit

#add sex
LMM_RNFL_HC1_sex <- lmer(VisusHC~RNFL_new + Geschlecht + (1|ID), 
                         data=Data_Visus_HC)
summary(LMM_RNFL_HC1_sex)
anova(LMM_RNFL_HC1_sex, LMM_RNFL_HC2)
#adding sex doesn't improve model fit

#add both
LMM_RNFL_HC_both <- lmer(VisusHC~RNFL_new + Geschlecht + Age +  (1|ID), 
                         data=Data_Visus_HC)
summary(LMM_RNFL_HC_both)
anova(LMM_RNFL_HC_both, LMM_RNFL_HC2)
#adding both sex and age doesn't improve model fit

#---RRMS---#
LMM_RNFL_RR <- lmer(VisusHC~RNFL_new + (RNFL_new|ID/ID_Eye), 
                    data=Data_Visus_RRMS)
summary(LMM_RNFL_RR)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID  (Intercept) 9.456e-02 0.307504      
#RNFL_new               7.064e-05 0.008405 -1.00
#ID        (Intercept)  2.293e-01 0.478858      
#RNFL_new               1.647e-04 0.012834 -0.99
#Residual               3.408e-03 0.058380

#variation in baseline RNFL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of RNFL = RE with lowest variance

LMM_RNFL_RR1 <- lmer(VisusHC~RNFL_new + (1|ID/ID_Eye), 
                     data=Data_Visus_RRMS)
summary(LMM_RNFL_RR1)
#best model fit - presented in paper

#add age
LMM_RNFL_RR1_Age <- lmer(VisusHC~RNFL_new + Age + (1|ID/ID_Eye), 
                         data=Data_Visus_RRMS)
summary(LMM_RNFL_RR1_Age)
anova(LMM_RNFL_RR1_Age, LMM_RNFL_RR1)
#adding age doesn't improve model fit

#add sex
LMM_RNFL_RR1_sex <- lmer(VisusHC~RNFL_new + Geschlecht + (1|ID/ID_Eye), 
                         data=Data_Visus_RRMS)
summary(LMM_RNFL_RR1_sex)
anova(LMM_RNFL_RR1_sex, LMM_RNFL_RR1_Age)
#adding sex doesn't improve model fit

#add both
LMM_RNFL_RR_both <- lmer(VisusHC~RNFL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                         data=Data_Visus_RRMS)
summary(LMM_RNFL_RR_both)
anova(LMM_RNFL_RR_both, LMM_RNFL_RR1)
#adding both sex and age doesn't improve model fit

#---PPMS---#
LMM_RNFL_PP <- lmer(VisusHC~RNFL_new + (RNFL_new|ID/ID_Eye), 
                    data=Data_Visus_PPMS)
summary(LMM_RNFL_PP)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 2.162e-02 0.147040      
#RNFL_new              5.554e-06 0.002357 -1.00
#ID        (Intercept) 1.255e-01 0.354312      
#RNFL_new              6.984e-05 0.008357 -0.97
#Residual              8.332e-03 0.091282  

#variation in baseline RNFL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of RNFL = RE with lowest variance

LMM_RNFL_PP <- lmer(VisusHC~RNFL_new + (1|ID/ID_Eye), 
                    data=Data_Visus_PPMS)
summary(LMM_RNFL_PP)
#best model fit - presented in paper

#add age
LMM_RNFL_PP_age <- lmer(VisusHC~RNFL_new + Age + (1|ID/ID_Eye), 
                        data=Data_Visus_PPMS)
summary(LMM_RNFL_PP_age)
anova(LMM_RNFL_PP_age, LMM_RNFL_PP)
#adding age doesn't improve model fit

#add sex
LMM_RNFL_PP_sex <- lmer(VisusHC~RNFL_new + Geschlecht + (1|ID/ID_Eye), 
                        data=Data_Visus_PPMS)
summary(LMM_RNFL_PP_sex)
anova(LMM_RNFL_PP_sex, LMM_RNFL_PP)
#adding sex doesn't improve model fit

#add both
LMM_RNFL_PP_both <- lmer(VisusHC~RNFL_new + Geschlecht + Age +(1|ID/ID_Eye), 
                         data=Data_Visus_PPMS)
summary(LMM_RNFL_PP_both)
anova(LMM_RNFL_PP_both, LMM_RNFL_PP)
#adding both sex and age doesn't improve model fit

#---SPMS---#
LMM_RNFL_SP <- lmer(VisusHC~RNFL_new + (RNFL_new|ID/ID_Eye), 
                    data=Data_Visus_SPMS)
summary(LMM_RNFL_SP)

#Error: number of observations (=262) <= number of random effects (=270) for term (RNFL_new | ID_Eye:ID);
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope (in line with approach for the other retinal layers)

LMM_RNFL_SP <- lmer(VisusHC~RNFL_new + (1|ID/ID_Eye), 
                    data=Data_Visus_SPMS)
summary(LMM_RNFL_SP)
#best model fit - presented in paper

#add age
LMM_RNFL_SP_age <- lmer(VisusHC~RNFL_new + Age + (1|ID/ID_Eye), 
                        data=Data_Visus_SPMS)
summary(LMM_RNFL_SP_age)
anova(LMM_RNFL_SP_age, LMM_RNFL_SP)
#adding age doesn't improve model fit

#add sex
LMM_RNFL_SP_sex <- lmer(VisusHC~RNFL_new + Geschlecht +  (1|ID/ID_Eye), 
                        data=Data_Visus_SPMS)
summary(LMM_RNFL_SP_sex)
anova(LMM_RNFL_SP_sex, LMM_RNFL_SP)
#adding both sex and age doesn't improve model fit

#add both
LMM_RNFL_SP_both <- lmer(VisusHC~RNFL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                         data=Data_Visus_SPMS)
summary(LMM_RNFL_SP_both)
anova(LMM_RNFL_SP_both, LMM_RNFL_SP)
#adding both sex and age doesn't improve model fit


#-----------------------Data loading & data management-----------------------# 

rm(list = ls())
OCT_data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)
OCT_data <- read_sav("20210310_OCTstudy_ACTUAL.sav")
Data_Latency <- subset(OCT_data, OCT_included==1 & VEPMatch==1)
Data_Latency_HC <- subset(Data_Latency, Verlaufsform==0)
Data_Latency_RRMS <- subset(Data_Latency, Verlaufsform==1)
Data_Latency_PPMS <- subset(Data_Latency, Verlaufsform==2)
Data_Latency_SPMS <- subset(Data_Latency, Verlaufsform==3)
Data_Latency_HC$Geschlecht <- as.factor(Data_Latency_HC$Geschlecht)
Data_Latency_RRMS$Geschlecht <- as.factor(Data_Latency_RRMS$Geschlecht)
Data_Latency_PPMS$Geschlecht <- as.factor(Data_Latency_PPMS$Geschlecht)
Data_Latency_SPMS$Geschlecht <- as.factor(Data_Latency_SPMS$Geschlecht)
#Verlaufsform=0 --> HC
#Verlaufsform=1 --> RRMS
#Verlaufsform=2 --> PPMS
#Verlaufsform=3 --> SPMS


#-----------------------------------LATENCY-----------------------------------#
#------------------------------------mRNFL------------------------------------# 
#---HC---#
LMM_RNFL_HC <- lmer(Latenz~RNFL_new + (RNFL_new|ID/ID_Eye), 
                    data=Data_Latency_HC)
summary(LMM_RNFL_HC)
#Model failed to converge with 1 negative eigenvalue: -6.7e+01 

#Groups    Name        Variance Std.Dev. Corr 
#ID_Eye:ID (Intercept) 16.30611 4.0381        
#           RNFL_new    0.01461 0.1209   -1.00
#ID        (Intercept)  0.00000 0.0000        
#           RNFL_new    0.01387 0.1178    NaN 

#random intercept of ID doesn't explain any variance --> reduce model complexity

LMM_RNFL_HC1 <- lmer(Latenz~RNFL_new + (0+RNFL_new|ID/ID_Eye), 
                     data=Data_Latency_HC)
summary(LMM_RNFL_HC1)
#still singular fit 

#Random effects:
#Groups    Name     Variance Std.Dev.
#ID_Eye:ID RNFL_new  0.00000 0.0000  
#ID        RNFL_new  0.01413 0.1189  
#Residual           13.69635 3.7009

#random intercept of eye per subject doesn't explain any variance --> reduce model complexity

LMM_RNFL_HC2 <- lmer(Latenz~RNFL_new + (0+RNFL_new|ID), 
                     data=Data_Latency_HC)
summary(LMM_RNFL_HC2)

#add age
LMM_RNFL_HC2_Age <- lmer(Latenz~RNFL_new + Age + (0+RNFL_new|ID), 
                         data=Data_Latency_HC)
summary(LMM_RNFL_HC2_Age)
anova(LMM_RNFL_HC2_Age, LMM_RNFL_HC2)
#adding age doesn't sign. improve model fit

#add sex
LMM_RNFL_HC2_sex <- lmer(Latenz~RNFL_new + Geschlecht + (0+RNFL_new|ID), 
                         data=Data_Latency_HC)
summary(LMM_RNFL_HC2_sex)
anova(LMM_RNFL_HC2_sex, LMM_RNFL_HC2)
#adding sex sign. improve model fit
#best model fit - presented in paper

#add both
LMM_RNFL_HC2_both <- lmer(Latenz~RNFL_new + Geschlecht + Age + (0+RNFL_new|ID), 
                         data=Data_Latency_HC)
summary(LMM_RNFL_HC2_both)
anova(LMM_RNFL_HC2_both, LMM_RNFL_HC2_sex)
#adding age in addition to sex doesn't further improve model fit

#---RRMS---#
LMM_RNFL_RR <- lmer(Latenz~RNFL_new + (RNFL_new|ID/ID_Eye), 
                    data=Data_Latency_RRMS)
summary(LMM_RNFL_RR)
#boundary (singular) fit: see help('isSingular')
#Warning message:
#Model failed to converge with 1 negative eigenvalue: -6.6e+02 

#Random effects:
#Groups    Name        Variance Std.Dev. Corr 
#ID_Eye:ID (Intercept) 252.2573 15.8826       
#RNFL_new              0.1727  0.4156  -1.00
#ID        (Intercept) 622.0436 24.9408       
#RNFL_new              0.1852  0.4304  -1.00
#Residual              37.8117  6.1491

#variation in baseline RNFL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of RNFL = RE with lowest variance

LMM_RNFL_RR1 <- lmer(Latenz~RNFL_new + (1|ID/ID_Eye), 
                     data=Data_Latency_RRMS)
summary(LMM_RNFL_RR1)

#add age
LMM_RNFL_RR1_Age <- lmer(Latenz~RNFL_new + Age + (1|ID/ID_Eye), 
                         data=Data_Latency_RRMS)
summary(LMM_RNFL_RR1_Age)
anova(LMM_RNFL_RR1_Age, LMM_RNFL_RR1)
#adding age sign. improves model fit
#best model fit - presented in paper

#add sex
LMM_RNFL_RR1_sex <- lmer(Latenz~RNFL_new + Geschlecht + (1|ID/ID_Eye), 
                         data=Data_Latency_RRMS)
summary(LMM_RNFL_RR1_sex)
anova(LMM_RNFL_RR1_sex, LMM_RNFL_RR1)
#adding sex doesn't improve model fit

#add both
LMM_RNFL_RR_both <- lmer(Latenz~RNFL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                         data=Data_Latency_RRMS)
summary(LMM_RNFL_RR_both)
anova(LMM_RNFL_RR_both, LMM_RNFL_RR1_Age)
#adding sex in addition to age doesn't further improve model 

#---PPMS---#
LMM_RNFL_PP <- lmer(Latenz~RNFL_new + (RNFL_new|ID/ID_Eye), 
                    data=Data_Latency_PPMS)
summary(LMM_RNFL_PP)

#Warning message:
#Model failed to converge with 1 negative eigenvalue: -3.1e-01 

#Random effects:
#Groups    Name        Variance Std.Dev. Corr 
#ID_Eye:ID (Intercept) 477.8764 21.8604       
#RNFL_new              0.2592  0.5092  -1.00
#ID        (Intercept) 175.8353 13.2603       
#RNFL_new              0.1721  0.4148  -0.79
#Residual              70.2115  8.3792

#variation in baseline RNFL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of RNFL = RE with lowest variance

LMM_RNFL_PP <- lmer(Latenz~RNFL_new + (1|ID/ID_Eye), 
                    data=Data_Latency_PPMS)
summary(LMM_RNFL_PP)
#best model fit - presented in paper

#add age
LMM_RNFL_PP_age <- lmer(Latenz~RNFL_new + Age + (1|ID/ID_Eye), 
                        data=Data_Latency_PPMS)
summary(LMM_RNFL_PP_age)
anova(LMM_RNFL_PP_age, LMM_RNFL_PP)
#adding age doesn't improve model fit

#add sex
LMM_RNFL_PP_sex <- lmer(Latenz~RNFL_new + Geschlecht + (1|ID/ID_Eye), 
                        data=Data_Latency_PPMS)
summary(LMM_RNFL_PP_sex)
anova(LMM_RNFL_PP_sex, LMM_RNFL_PP)
#adding sex doesn't improve model fit

#add both
LMM_RNFL_PP_both <- lmer(Latenz~RNFL_new + Geschlecht + Age +(1|ID/ID_Eye), 
                         data=Data_Latency_PPMS)
summary(LMM_RNFL_PP_both)
anova(LMM_RNFL_PP_both, LMM_RNFL_PP)
#adding both sex and age doesn't improve model fit

#---SPMS---#
LMM_RNFL_SP <- lmer(Latenz~RNFL_new + (RNFL_new|ID/ID_Eye), 
                    data=Data_Latency_SPMS)
summary(LMM_RNFL_SP)

#Error: number of observations (=132) <= number of random effects (=140) for term (RNFL_new | ID_Eye:ID); 
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope (in line with approach for VA)

LMM_RNFL_SP <- lmer(Latenz~RNFL_new + (1|ID/ID_Eye), 
                    data=Data_Latency_SPMS)
summary(LMM_RNFL_SP)
#best model fit - presented in paper

#add age
LMM_RNFL_SP_age <- lmer(Latenz~RNFL_new + Age + (1|ID/ID_Eye), 
                        data=Data_Latency_SPMS)
summary(LMM_RNFL_SP_age)
anova(LMM_RNFL_SP_age, LMM_RNFL_SP)
#adding age doesn't improve model fit

#add sex
LMM_RNFL_SP_sex <- lmer(Latenz~RNFL_new + Geschlecht +  (1|ID/ID_Eye), 
                        data=Data_Latency_SPMS)
summary(LMM_RNFL_SP_sex)
anova(LMM_RNFL_SP_sex, LMM_RNFL_SP)
#adding sex doesn't improve model fit

#add both
LMM_RNFL_SP_both <- lmer(Latenz~RNFL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                         data=Data_Latency_SPMS)
summary(LMM_RNFL_SP_both)
anova(LMM_RNFL_SP_both, LMM_RNFL_SP)
#adding both age and sex doesn't improve model fit


#------------------------------------GCIPL------------------------------------# 
#---HC---#

LMM_GCIPL_HC <- lmer(Latenz~GCIPL_new + (GCIPL_new|ID/ID_Eye), 
                     data=Data_Latency_HC)
summary(LMM_GCIPL_HC)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev.  Corr 
#ID_Eye:ID (Intercept) 1.361e+01 3.6893620      
#GCIPL_new             2.568e-03 0.0506732 -1.00
#ID        (Intercept) 1.435e+01 3.7876771      
#GCIPL_new             9.189e-07 0.0009586 -1.00
#Residual              1.416e+01 3.7631417

##variation in baseline GCIPL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of GCIPL = RE with lowest variance

LMM_GCIPL_HC1 <- lmer(Latenz~GCIPL_new + (1|ID/ID_Eye), 
                      data=Data_Latency_HC)
summary(LMM_GCIPL_HC1)
#still singular fit

#Random effects:
#Groups    Name        Variance Std.Dev.
#ID_Eye:ID (Intercept)  0.00    0.000   
#ID        (Intercept) 14.64    3.826   
#Residual              14.10    3.755

#random intercept of eye per subject doesn't explain any variance --> reduce model complexity

LMM_GCIPL_HC2 <- lmer(Latenz~GCIPL_new + (1|ID), 
                      data=Data_Latency_HC)
summary(LMM_GCIPL_HC2)

#add age
LMM_GCIPL_HC1_Age <- lmer(Latenz~GCIPL_new + Age + (1|ID), 
                          data=Data_Latency_HC)
summary(LMM_GCIPL_HC1_Age)
anova(LMM_GCIPL_HC1_Age, LMM_GCIPL_HC2)
#adding age doesn't improve model fit

#add sex
LMM_GCIPL_HC1_sex <- lmer(Latenz~GCIPL_new + Geschlecht + (1|ID), 
                          data=Data_Latency_HC)
summary(LMM_GCIPL_HC1_sex)
anova(LMM_GCIPL_HC1_sex, LMM_GCIPL_HC2)
#adding sex sign. improves model fit
#best model fit - presented in paper

#add both
LMM_GCIPL_HC_both <- lmer(Latenz~GCIPL_new + Geschlecht + Age +  (1|ID), 
                          data=Data_Latency_HC)
summary(LMM_GCIPL_HC_both)
anova(LMM_GCIPL_HC_both, LMM_GCIPL_HC1_sex)
#adding age in addition to sex doesn't further improve model fit

#---RRMS---#
LMM_GCIPL_RR <- lmer(Latenz~GCIPL_new + (GCIPL_new|ID/ID_Eye), 
                     data=Data_Latency_RRMS)
summary(LMM_GCIPL_RR)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev.  Corr 
#ID_Eye:ID (Intercept) 1.080e+02 10.393026      
#GCIPL_new             1.220e-02  0.110432 -1.00
#ID        (Intercept) 1.053e+02 10.259758      
#GCIPL_new             6.279e-05  0.007924 1.00 
#Residual              3.910e+01  6.253148

##variation in baseline GCIPL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of GCIPL = RE with lowest variance

LMM_GCIPL_RR1 <- lmer(Latenz~GCIPL_new + (1|ID/ID_Eye), 
                      data=Data_Latency_RRMS)
summary(LMM_GCIPL_RR1)

#add age
LMM_GCIPL_RR1_Age <- lmer(Latenz~GCIPL_new + Age + (1|ID/ID_Eye), 
                          data=Data_Latency_RRMS)
summary(LMM_GCIPL_RR1_Age)
anova(LMM_GCIPL_RR1_Age, LMM_GCIPL_RR1)
#age sign. improves model fit
#best model fit - presented in paper

#add sex
LMM_GCIPL_RR1_sex <- lmer(Latenz~GCIPL_new + Geschlecht + (1|ID/ID_Eye), 
                          data=Data_Latency_RRMS)
summary(LMM_GCIPL_RR1_sex)
anova(LMM_GCIPL_RR1_sex, LMM_GCIPL_RR1)
#age doesn't sign. improve model fit

#add both
LMM_GCIPL_RR_both <- lmer(Latenz~GCIPL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Latency_RRMS)
summary(LMM_GCIPL_RR_both)
anova(LMM_GCIPL_RR_both, LMM_GCIPL_RR1_Age)
#adding sex in addition to age doesn't further improve model fit

#---PPMS---#
LMM_GCIPL_PP <- lmer(Latenz~GCIPL_new + (GCIPL_new|ID/ID_Eye), 
                     data=Data_Latency_PPMS)
summary(LMM_GCIPL_PP)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 9.436e+01  9.71398      
#GCIPL_new             3.891e-03  0.06238 -1.00
#ID        (Intercept) 1.305e+02 11.42301      
#GCIPL_new             1.529e-03  0.03910 -1.00
#Residual              7.260e+01  8.52048

#variation in baseline GCIPL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of GCIPL = RE with lowest variance

LMM_GCIPL_PP <- lmer(Latenz~GCIPL_new + (1|ID/ID_Eye), 
                     data=Data_Latency_PPMS)
summary(LMM_GCIPL_PP)
#best model fit - presented in paper

#add age
LMM_GCIPL_PP_age <- lmer(Latenz~GCIPL_new + Age + (1|ID/ID_Eye), 
                         data=Data_Latency_PPMS)
summary(LMM_GCIPL_PP_age)
anova(LMM_GCIPL_PP_age, LMM_GCIPL_PP)
#age doesn't improve model fit

#add sex
LMM_GCIPL_PP_sex <- lmer(Latenz~GCIPL_new + Geschlecht + (1|ID/ID_Eye), 
                         data=Data_Latency_PPMS)
summary(LMM_GCIPL_PP_sex)
anova(LMM_GCIPL_PP_sex, LMM_GCIPL_PP)
#sex doesn't improve model fit

#add both
LMM_GCIPL_PP_both <- lmer(Latenz~GCIPL_new + Geschlecht + Age +(1|ID/ID_Eye), 
                          data=Data_Latency_PPMS)
summary(LMM_GCIPL_PP_both)
anova(LMM_GCIPL_PP_both, LMM_GCIPL_PP)
#adding both sex and age doesn't improve model fit

#---SPMS---#
LMM_GCIPL_SP <- lmer(Latenz~GCIPL_new + (GCIPL_new|ID/ID_Eye), 
                     data=Data_Latency_SPMS)
summary(LMM_GCIPL_SP)

#Error: number of observations (=132) <= number of random effects (=140) for term (GCIPL_new | ID_Eye:ID);
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope (in line with approach for VA and other retinal layers)

LMM_GCIPL_SP <- lmer(Latenz~GCIPL_new + (1|ID/ID_Eye), 
                     data=Data_Latency_SPMS)
summary(LMM_GCIPL_SP)
#best model fit - presented in paper

#add age
LMM_GCIPL_SP_age <- lmer(Latenz~GCIPL_new + Age + (1|ID/ID_Eye), 
                         data=Data_Latency_SPMS)
summary(LMM_GCIPL_SP_age)
anova(LMM_GCIPL_SP_age, LMM_GCIPL_SP)
#adding age doesn't improve model fit

#add sex
LMM_GCIPL_SP_sex <- lmer(Latenz~GCIPL_new + Geschlecht +  (1|ID/ID_Eye), 
                         data=Data_Latency_SPMS)
summary(LMM_GCIPL_SP_sex)
anova(LMM_GCIPL_SP_sex, LMM_GCIPL_SP)
#adding sex doesn't improve model fit

#add both
LMM_GCIPL_SP_both <- lmer(Latenz~GCIPL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Latency_SPMS)
summary(LMM_GCIPL_SP_both)
anova(LMM_GCIPL_SP_both, LMM_GCIPL_SP)
#adding both sex and age doesn't improve model fit


#-------------------------------------INL-------------------------------------# 
#---HC---#
LMM_INL_HC <- lmer(Latenz~INL_new + (INL_new|ID/ID_Eye), 
        data=Data_Latency_HC)
summary(LMM_INL_HC)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 2.265e-01  0.47594      
#INL_new               1.892e-04  0.01376 -1.00
#ID        (Intercept) 2.211e+02 14.86869      
#INL_new               1.106e-01  0.33258 -1.00
#Residual              1.437e+01  3.79103

##variation in baseline INL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of INL = RE with lowest variance

LMM_INL_HC1 <- lmer(Latenz~INL_new + (1|ID/ID_Eye), 
                    data=Data_Latency_HC)
summary(LMM_INL_HC1)
#still singular fit

#Random effects:
#Groups    Name        Variance Std.Dev.
#ID_Eye:ID (Intercept)  0.00    0.000   
#ID        (Intercept) 13.80    3.715   
#Residual              14.37    3.791 

#random intercept of eye per subject doesn't explain any variance --> reduce model complexity

LMM_INL_HC2 <- lmer(Latenz~INL_new + (1|ID), 
                    data=Data_Latency_HC)
summary(LMM_INL_HC2)

#add age
LMM_INL_HC1_Age <- lmer(Latenz~INL_new + Age + (1|ID), 
                        data=Data_Latency_HC)
summary(LMM_INL_HC1_Age)
anova(LMM_INL_HC1_Age, LMM_INL_HC2)
#age sign. improves model fit

#add sex
LMM_INL_HC1_sex <- lmer(Latenz~INL_new + Geschlecht + (1|ID), 
                        data=Data_Latency_HC)
summary(LMM_INL_HC1_sex)
anova(LMM_INL_HC1_sex, LMM_INL_HC2)
#sex sign. improves model fit
#best model fit - presented in paper

#add both
LMM_INL_HC_both <- lmer(Latenz~INL_new + Geschlecht + Age +  (1|ID), 
                        data=Data_Latency_HC)
summary(LMM_INL_HC_both)
anova(LMM_INL_HC_both, LMM_INL_HC1_sex)
#adding age in addition to sex doesn't further improve model fit

#---RRMS---#
LMM_INL_RR <- lmer(Latenz~INL_new + (INL_new|ID/ID_Eye), 
                   data=Data_Latency_RRMS)
summary(LMM_INL_RR)

#Warning messages:
#1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#unable to evaluate scaled gradient
#2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
#3: Model failed to converge with 1 negative eigenvalue: -7.6e+02

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 3.060e+01 5.531687      
#INL_new               1.984e-05 0.004455 -0.18
#ID        (Intercept) 3.252e+01 5.703066      
#INL_new               3.057e-01 0.552889 -1.00
#Residual              3.784e+01 6.151473 

#variation in baseline INL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of INL = RE with lowest variance

LMM_INL_RR1 <- lmer(Latenz~INL_new + (1|ID/ID_Eye), 
                    data=Data_Latency_RRMS)
summary(LMM_INL_RR1)

#add age
LMM_INL_RR1_Age <- lmer(Latenz~INL_new + Age + (1|ID/ID_Eye), 
                        data=Data_Latency_RRMS)
summary(LMM_INL_RR1_Age)
anova(LMM_INL_RR1_Age, LMM_INL_RR1)
#age sign. improves model fit
#best model fit - presented in paper

#add sex
LMM_INL_RR1_sex <- lmer(Latenz~INL_new + Geschlecht + (1|ID/ID_Eye), 
                        data=Data_Latency_RRMS)
summary(LMM_INL_RR1_sex)
anova(LMM_INL_RR1_sex, LMM_INL_RR1)
#sex doesn't improve model fit

#add both
LMM_INL_RR_both <- lmer(Latenz~INL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                        data=Data_Latency_RRMS)
summary(LMM_INL_RR_both)
anova(LMM_INL_RR_both, LMM_INL_RR1_Age)
#adding sex in addition to age doesn't further improve model fit

#---PPMS---#
LMM_INL_PP <- lmer(Latenz~INL_new + (INL_new|ID/ID_Eye), 
                   data=Data_Latency_PPMS)
summary(LMM_INL_PP)
#Singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 8.701e+02 29.49729      
#INL_new               1.116e+00  1.05639 -1.00
#ID        (Intercept) 1.635e+02 12.78689      
#INL_new               6.302e-03  0.07939 -1.00
#Residual              7.179e+01  8.47287

#variation in baseline INL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of INL = RE with lowest variance

LMM_INL_PP <- lmer(Latenz~INL_new + (1|ID/ID_Eye), 
                   data=Data_Latency_PPMS)
summary(LMM_INL_PP)
#best model fit - presented in paper

#add age
LMM_INL_PP_age <- lmer(Latenz~INL_new + Age + (1|ID/ID_Eye), 
                       data=Data_Latency_PPMS)
summary(LMM_INL_PP_age)
anova(LMM_INL_PP_age, LMM_INL_PP)
#adding age doesn't improve model fit

#add sex
LMM_INL_PP_sex <- lmer(Latenz~INL_new + Geschlecht + (1|ID/ID_Eye), 
                       data=Data_Latency_PPMS)
summary(LMM_INL_PP_sex)
anova(LMM_INL_PP_sex, LMM_INL_PP)
#adding sex doesn't improve model fit

#add both
LMM_INL_PP_both <- lmer(Latenz~INL_new + Geschlecht + Age +(1|ID/ID_Eye), 
                        data=Data_Latency_PPMS)
summary(LMM_INL_PP_both)
anova(LMM_INL_PP_both, LMM_INL_PP)
#adding both sex and age doesn't improve model fit


#---SPMS---#
LMM_INL_SP <- lmer(Latenz~INL_new + (INL_new|ID/ID_Eye), 
                   data=Data_Latency_SPMS)
summary(LMM_INL_SP)

#Error: number of observations (=132) <= number of random effects (=140) for term (INL_new | ID_Eye:ID);
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope (in line with approach for VA and other retinal layers)

LMM_INL_SP <- lmer(Latenz~INL_new + (1|ID/ID_Eye), 
                   data=Data_Latency_SPMS)
summary(LMM_INL_SP)
#best model fit - presented in paper

#add age
LMM_INL_SP_age <- lmer(Latenz~INL_new + Age + (1|ID/ID_Eye), 
                       data=Data_Latency_SPMS)
summary(LMM_INL_SP_age)
anova(LMM_INL_SP_age, LMM_INL_SP)
#age doesn't improve model fit

#add sex
LMM_INL_SP_sex <- lmer(Latenz~INL_new + Geschlecht +  (1|ID/ID_Eye), 
                       data=Data_Latency_SPMS)
summary(LMM_INL_SP_sex)
anova(LMM_INL_SP_sex, LMM_INL_SP)
#adding sex doesn't improve model fit

#add both
LMM_INL_SP_both <- lmer(Latenz~INL_new + Geschlecht + Age +  (1|ID/ID_Eye), 
                        data=Data_Latency_SPMS)
summary(LMM_INL_SP_both)
anova(LMM_INL_SP_both, LMM_INL_SP)
#adding both sex and age doesn't improve model fit


#------------------------------------pRNFL------------------------------------# 
#---RRMS---#

LMM_pRNFL_RR <- lmer(Latenz~globalpRNFL + (globalpRNFL|ID/ID_Eye), 
                     data=Data_Latency_RRMS)
summary(LMM_pRNFL_RR)
#singular fit

#Random effects:
#Groups    Name        Variance  Std.Dev. Corr 
#ID_Eye:ID (Intercept) 58.392330 7.64149       
#globalpRNFL           0.001761 0.04196  -1.00
#ID        (Intercept) 36.692117 6.05740       
#globalpRNFL           0.002333 0.04830  1.00 
#Residual              41.588855 6.44894

#variation in baseline pRNFL thickness seems to be more important than variation in degeneration 
#or regeneration over time per eye and subject --> eliminate random slope of pRNFL = RE with lowest variance

LMM_pRNFL_RR1 <- lmer(Latenz~globalpRNFL + (1|ID/ID_Eye), 
                      data=Data_Latency_RRMS)
summary(LMM_pRNFL_RR1)
#best model fit - presented in paper

#add age
LMM_pRNFL_RR1_Age <- lmer(Latenz~globalpRNFL + Age + (1|ID/ID_Eye), 
                          data=Data_Latency_RRMS)
summary(LMM_pRNFL_RR1_Age)
anova(LMM_pRNFL_RR1_Age, LMM_pRNFL_RR1)
#adding age doesn't sign. improve model fit

#add sex
LMM_pRNFL_RR1_sex <- lmer(Latenz~globalpRNFL + Geschlecht + (1|ID/ID_Eye), 
                          data=Data_Latency_RRMS)
summary(LMM_pRNFL_RR1_sex)
anova(LMM_pRNFL_RR1_sex, LMM_pRNFL_RR1)
#adding sex doesn't sign. improve model fit

#add both
LMM_pRNFL_RR_both <- lmer(Latenz~globalpRNFL + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Latency_RRMS)
summary(LMM_pRNFL_RR_both)
anova(LMM_pRNFL_RR_both, LMM_pRNFL_RR1)
#adding both sex and age doesn't sign. improve model fit

#---PPMS---#
LMM_pRNFL_PP <- lmer(Latenz~globalpRNFL + (globalpRNFL|ID/ID_Eye), 
                     data=Data_Latency_PPMS)
summary(LMM_pRNFL_PP)

#Error: number of observations (=172) <= number of random effects (=172) for term (globalpRNFL | ID_Eye:ID);
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope (in line with approach for VA and other retinal layers)

LMM_pRNFL_PP <- lmer(Latenz~globalpRNFL + (1|ID/ID_Eye), 
                     data=Data_Latency_PPMS)
summary(LMM_pRNFL_PP)
#best model fit - presented in paper

#add age
LMM_pRNFL_PP_age <- lmer(Latenz~globalpRNFL + Age + (1|ID/ID_Eye), 
                         data=Data_Latency_PPMS)
summary(LMM_pRNFL_PP_age)
anova(LMM_pRNFL_PP_age, LMM_pRNFL_PP)
#adding age doesn't sign. improve model fit

#add sex
LMM_pRNFL_PP_sex <- lmer(Latenz~globalpRNFL + Geschlecht + (1|ID/ID_Eye), 
                         data=Data_Latency_PPMS)
summary(LMM_pRNFL_PP_sex)
anova(LMM_pRNFL_PP_sex, LMM_pRNFL_PP)
#adding sex doesn't sign. improve model fit

#add both
LMM_pRNFL_PP_both <- lmer(Latenz~globalpRNFL + Geschlecht + Age +(1|ID/ID_Eye), 
                          data=Data_Latency_PPMS)
summary(LMM_pRNFL_PP_both)
anova(LMM_pRNFL_PP_both, LMM_pRNFL_PP)
#adding sex and age doesn't sign. improve model fit

#---SPMS---#
LMM_pRNFL_SP <- lmer(Latenz~globalpRNFL + (globalpRNFL|ID/ID_Eye), 
                     data=Data_Latency_SPMS)
summary(LMM_pRNFL_SP)

#Error: number of observations (=119) <= number of random effects (=126) for term (globalpRNFL | ID_Eye:ID); 
#the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
#--> reduce to less complex model without random slope (in line with approach for VA and other retinal layers)

LMM_pRNFL_SP <- lmer(Latenz~globalpRNFL + (1|ID/ID_Eye), 
                     data=Data_Latency_SPMS)
summary(LMM_pRNFL_SP)
#best model fit - presented in paper

#add age
LMM_pRNFL_SP_age <- lmer(Latenz~globalpRNFL + Age + (1|ID/ID_Eye), 
                         data=Data_Latency_SPMS)
summary(LMM_pRNFL_SP_age)
anova(LMM_pRNFL_SP_age, LMM_pRNFL_SP)
#adding age doesn't improve model fit

#add sex
LMM_pRNFL_SP_sex <- lmer(Latenz~globalpRNFL + Geschlecht +  (1|ID/ID_Eye), 
                         data=Data_Latency_SPMS)
summary(LMM_pRNFL_SP_sex)
anova(LMM_pRNFL_SP_sex, LMM_pRNFL_SP)
#adding sex and age doesn't improve model fit

#add both
LMM_pRNFL_SP_both <- lmer(Latenz~globalpRNFL + Geschlecht + Age +  (1|ID/ID_Eye), 
                          data=Data_Latency_SPMS)
summary(LMM_pRNFL_SP_both)
anova(LMM_pRNFL_SP_both, LMM_pRNFL_SP)
#adding both sex and age doesn't improve model fit
