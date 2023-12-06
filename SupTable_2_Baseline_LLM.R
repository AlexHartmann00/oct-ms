

#----------------Linear mixed-effects models (LMM) at baseline---------------#

# This script computes LMMs analyzing the influence of disease duration on retinal 
# layer thickness (globalpRNFL, macular RNFL, GCIPL, INL), visual acuity, and VEP 
# latency at baseline. To account for confounding influences, age, sex, and/or ON 
# were subsequently added to the model if the model fit could be improved. 
# The results are presented in Supplementary Table 2. 


#---------------------------------R packages---------------------------------#  

library(lme4)
library(lmerTest) #needed for p-values
library(haven)
library(performance)
library(car)
library(openxlsx)

#-----------------------Data loading & data management-----------------------# 

OCT_data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)

OCT_data$ONHistory <- ifelse(OCT_data$RBNinVergangenheit == 1,
                                 "RBN",
                                 ifelse(
                                   OCT_data$HistoryofON_GCIPL == 1,
                                   "Any",
                                   "None"))

#Correct disease duration, if VEP dates do not match OCT dates
library(lubridate)
OCT_data$Disease_Duration_Latency <- as.numeric((as_date(OCT_data$VEP) - as_date(OCT_data$OCTDatum))/365 + OCT_data$Disease_Duration_FINAL)
OCT_data$Age_Latency <- as.numeric((as_date(OCT_data$VEP) - as_date(OCT_data$OCTDatum))/365 + OCT_data$Age)


OCT_data <- subset(OCT_data, OCT_included==1)
OCT_data$Verlaufsform <- as.factor(OCT_data$Verlaufsform)
Data_Baseline <- subset(OCT_data, Sequenz==1)
Data_Baseline$Geschlecht <- as.factor(Data_Baseline$Geschlecht)

Data_LMM_RRMS_BL <- subset(OCT_data, OCT_included==1 & Sequenz==1 &  Verlaufsform==1)
Data_LMM_RRMS_BL$Geschlecht <- as.factor(Data_LMM_RRMS_BL$Geschlecht)
Data_LMM_PPMS_BL <- subset(OCT_data, OCT_included==1 & Sequenz==1 & Verlaufsform==2)
Data_LMM_PPMS_BL$Geschlecht <- as.factor(Data_LMM_PPMS_BL$Geschlecht)
Data_LMM_SPMS_BL <- subset(OCT_data, OCT_included==1 & Sequenz==1 & Verlaufsform==3)
Data_LMM_SPMS_BL$Geschlecht <- as.factor(Data_LMM_SPMS_BL$Geschlecht)

covariates_to_consider <- c("Age","Geschlecht","ONHistory")

#Sequenz=1 --> Baseline (BL)
#Verlaufsform=1 --> RRMS
#Verlaufsform=2 --> PPMS
#Verlaufsform=3 --> SPMS
#OCT_included=1 --> only OCT measurements included in the study

#Influence of disease duration on baseline scores.
#One random effect for disease duration per person, not per eye, because disease 
#duration doesn't differ between eyes but only between subjects.
#Random intercept not necessary, because intercept would represent 
#the baseline values, which are the only values included in the model. 

#Define function for model fitting
build_covariates_list <- function(covariates){
  # Build a list of all combinations of a set of covariates
  # To be used for feature selection
  lst <- list()
  factor_list <- list()
  i <- 1
  for(covariate in covariates){
    factor_list[[i]] <- factor(0:1)
    i <- i + 1
  }
  design_mat <- expand.grid(factor_list)
  for(i in 1:nrow(design_mat)){
    tmp <- c()
    for(j in 1:ncol(design_mat)){
      if(design_mat[i,j] == 1){
        tmp <- c(tmp,covariates[j])
      }
    }
    lst[[i]] <- tmp
  }
  return(lst[2:length(lst)])
}

step.Lmer <- function(X,depvar,idvar){
  #Function to perform full grid search for best model configuration
  #Naively brute forces the search across all possible covariate combinations
  X <- as.data.frame(X)[c(depvar, idvar, covariates_to_consider, "ID", "Age_Latency")]

  if(depvar == "Latenz"){
    necessary <- c(depvar, idvar, covariates_to_consider, "ID", "Age_Latency")
  }else{
    necessary <- c(depvar, idvar, covariates_to_consider, "ID")
  }

  X <- X[complete.cases(X[,necessary]),]
  covariate_combinations <- build_covariates_list(covariates_to_consider)

  
  #Keep counter for list appending
  i <- 1
  
  #Get participant ID column
  ids <- X$ID
  #Get focal predictor
  IV <- X[,idvar]
  #Select input data
  data_ <- X[,covariates_to_consider]
  #Build final input DF
  data_$disease_duration <- IV
  data_$depvar <- X[,depvar]
  data_$ids <- ids
  if(depvar == "Latenz"){
    data_$Age <- X$Age_Latency
  }
  
  #Initialize with a null model
  bestmod <- lmer("depvar ~ disease_duration + (0+disease_duration|ids)",data=data_)
  bestaic <- AIC(bestmod)
  bestcovariates <- c()
  
  # If there should only be positive / negative outcomes, skip computation
  if(length(table(X[,depvar])) == 1){
    return(list(Model=bestmod,Covariates = bestcovariates,AIC=bestaic,Data=X,DV=dv))
  }
  
  
  for(covariates in covariate_combinations){
    #Iterate over all combinations of covariates
    #Build model formula by appending covariates
    formula_string <- "depvar ~ disease_duration + "
    for(covariate in covariates){
      formula_string <- paste(formula_string, covariate, "+")
    }

    formula_string <- paste(formula_string, " (0+disease_duration|ids)")
    
    #Fit model
    mod <- lmer(formula(formula_string),data=data_)
    convergence_error = F
    tryCatch(
      {mod <- lmer(formula(formula_string),data=data_)},
      error=function(cond){
        improvement <- F
        convergence_error <- T
      }
    )
    if(i > 1 && !convergence_error){
      #If i == 1, this is always an improvement
      #Check chi2 for significance
      p_value <- anova(mod, bestmod,refit=FALSE)$`Pr(>Chisq)`[2] 

      threshold <- 0.05
      
      improvement <- p_value < threshold
      
      if(is.na(improvement)){
        improvement <- FALSE
      }
    }else{
      improvement <- TRUE
    }
    
    #If the covariates list provided improvement, save this as best model
    if(improvement){#} && !lme4::isSingular(mod)){
      bestcovariates <- covariates
      bestmod <- mod
      bestaic <- AIC(mod)
      dv <- insight::get_response(mod)
    }
    
    i <- i+ 1
  }
  return(list(Model=bestmod,Covariates = bestcovariates,AIC=bestaic,Data=X,DV=dv))
}

custom.summary <- function(model){
  df <- as.data.frame(summary(model)$coefficients)
  new_df <- data.frame(row.names = row.names(df))
  new_df$b <- df$Estimate
  new_df$LLCI <- new_df$b - 1.96 * df$`Std. Error`
  new_df$ULCI <- new_df$b + 1.96 * df$`Std. Error`
  new_df$p <- df$`Pr(>|t|)`
  return(new_df)
}

#---------------------------------globalpRNFL---------------------------------#    
#---RRMS---#
best_fit <- step.Lmer(
  Data_LMM_RRMS_BL,
  "globalpRNFL",
  "Disease_Duration_FINAL"
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---PPMS---#
best_fit <- step.Lmer(
  Data_LMM_PPMS_BL,
  "globalpRNFL",
  "Disease_Duration_FINAL"
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---SPMS---#
best_fit <- step.Lmer(
  Data_LMM_SPMS_BL,
  "globalpRNFL",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#-----------------------------------GCIPL------------------------------------# 
#---RRMS---#
best_fit <- step.Lmer(
  Data_LMM_RRMS_BL,
  "GCIPL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---PPMS---#
best_fit <- step.Lmer(
  Data_LMM_PPMS_BL,
  "GCIPL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---SPMS---#
best_fit <- step.Lmer(
  Data_LMM_SPMS_BL,
  "GCIPL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#------------------------------------INL-------------------------------------# 
#---RRMS---#
best_fit <- step.Lmer(
  Data_LMM_RRMS_BL,
  "INL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---PPMS---#
best_fit <- step.Lmer(
  Data_LMM_PPMS_BL,
  "INL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---SPMS---#
best_fit <- step.Lmer(
  Data_LMM_SPMS_BL,
  "INL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)


#------------------------------------mRNFL-----------------------------------# 
#---RRMS---#
best_fit <- step.Lmer(
  Data_LMM_RRMS_BL,
  "RNFL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---PPMS---#
best_fit <- step.Lmer(
  Data_LMM_PPMS_BL,
  "RNFL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---SPMS---#
best_fit <- step.Lmer(
  Data_LMM_SPMS_BL,
  "RNFL_new",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#--------------------------------VEP Latency----------------------------------#
Data_LMM_RRMS_BL$Disease_Duration_FINAL <- Data_LMM_RRMS_BL$Disease_Duration_Latency

#---RRMS---#
best_fit <- step.Lmer(
  Data_LMM_RRMS_BL,
  "Latenz",
  "Disease_Duration_Latency"
)



custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---PPMS---#
best_fit <- step.Lmer(
  Data_LMM_PPMS_BL,
  "Latenz",
  "Disease_Duration_Latency"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---SPMS---#
best_fit <- step.Lmer(
  Data_LMM_SPMS_BL,
  "Latenz",
  "Disease_Duration_Latency"
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#--------------------Visus HC (high contrast visual acuity)-------------------#
#---RRMS---#
best_fit <- step.Lmer(
  Data_LMM_RRMS_BL,
  "VisusHC",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---PPMS---#
best_fit <- step.Lmer(
  Data_LMM_PPMS_BL,
  "VisusHC",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)

#---SPMS---#
best_fit <- step.Lmer(
  Data_LMM_SPMS_BL,
  "VisusHC",
  "Disease_Duration_FINAL"
  
)


custom.summary(best_fit$Model)
#Chosen covariates
print(best_fit$Covariates)


#------------------------------Group Comparisons------------------------------#

step.Lmer.Interaction <- function(X,depvar,idvar){
  #Function to perform full grid search for best model configuration
  #Naively brute forces the search across all possible covariate combinations
  X <- as.data.frame(X)[c(depvar, idvar, covariates_to_consider, "ID", "Verlaufsform","Age_Latency")]
  #print(X[1:5,])
  if(depvar == "Latenz"){
    necessary <- c(depvar, idvar, covariates_to_consider, "ID", "Verlaufsform", "Age_Latency")
  }else{
    necessary <- c(depvar, idvar, covariates_to_consider, "ID", "Verlaufsform")
  }
  
  X <- X[complete.cases(X[,necessary]),]
  covariate_combinations <- build_covariates_list(covariates_to_consider)
  
  
  #Keep counter for list appending
  i <- 1
  
  #Get participant ID column
  ids <- X$ID
  #Get focal predictor
  IV <- X[,idvar]
  #Select input data
  data_ <- X[,covariates_to_consider]
  #Build final input DF
  data_$disease_duration <- IV
  data_$depvar <- X[,depvar]
  data_$ids <- ids
  data_$disease_course <- X$Verlaufsform
  if(depvar == "Latenz"){
    data_$Age <- X$Age_Latency
  }
  
  #Initialize with a null model
  bestmod <- lmer(formula("depvar ~ disease_duration*disease_course + (0+disease_duration|ids)"),data=data_)
  bestaic <- AIC(bestmod)
  bestcovariates <- c()
  
  # If there should only be positive / negative outcomes, skip computation
  if(length(table(X[,depvar])) == 1){
    return(list(Model=bestmod,Covariates = bestcovariates,AIC=bestaic,Data=X,DV=dv))
  }
  
  
  for(covariates in covariate_combinations){
    #Iterate over all combinations of covariates
    #Build model formula by appending covariates
    formula_string <- "depvar ~ disease_duration*disease_course + "
    for(covariate in covariates){
      formula_string <- paste(formula_string, covariate, "+")
    }
    
    formula_string <- paste(formula_string, " (0+disease_duration|ids)")
    
    #Fit model
    convergence_error = F
    tryCatch(
      {mod <- lmer(formula(formula_string),data=data_)},
      error=function(cond){
        improvement <- F
        convergence_error <- T
      }
    )
    
    if(i > 1 && !convergence_error){
      #If i == 1, this is always an improvement
      #Check chi2 for significance
      p_value <- anova(mod, bestmod,refit=FALSE)$`Pr(>Chisq)`[2] 
      threshold <- 0.05 
      
      improvement <- p_value < threshold
      
      if(is.na(improvement)){
        improvement <- FALSE
      }
    }else{
      improvement <- TRUE
    }
    
    #If the covariates list provided improvement, save this as best model
    if(improvement && !lme4::isSingular(mod)){
      bestcovariates <- covariates
      bestmod <- mod
      bestaic <- AIC(mod)
      dv <- insight::get_response(mod)
    }
    
    i <- i+ 1
  }
  return(list(Model=bestmod,Covariates = bestcovariates,AIC=bestaic,Data=X,DV=dv))
}

#---------------------------------globalpRNFL---------------------------------#    

best_fit <- step.Lmer.Interaction(
  OCT_data,
  "globalpRNFL",
  "Disease_Duration_FINAL"
  
)

custom.summary(best_fit$Model)
anova(best_fit$Model)

#------------------------------------GCIPL------------------------------------#

best_fit <- step.Lmer.Interaction(
  OCT_data,
  "GCIPL_new",
  "Disease_Duration_FINAL"
  
)

custom.summary(best_fit$Model)
anova(best_fit$Model)

#--------------------------------------INL------------------------------------#

best_fit <- step.Lmer.Interaction(
  OCT_data,
  "INL_new",
  "Disease_Duration_FINAL"
  
)

custom.summary(best_fit$Model)
anova(best_fit$Model)

#------------------------------------mRNFL------------------------------------#

best_fit <- step.Lmer.Interaction(
  OCT_data,
  "RNFL_new",
  "Disease_Duration_FINAL"
  
)

custom.summary(best_fit$Model)
anova(best_fit$Model)

#--------------------------------VEP Latency----------------------------------#

best_fit <- step.Lmer.Interaction(
  OCT_data,
  "Latenz",
  "Disease_Duration_Latency"
  
)

custom.summary(best_fit$Model)
anova(best_fit$Model)

#-----------------------------------Visus ------------------------------------#

best_fit <- step.Lmer.Interaction(
  OCT_data,
  "VisusHC",
  "Disease_Duration_FINAL"
  
)

custom.summary(best_fit$Model)
anova(best_fit$Model)



############################################
##############Plotting######################
############################################

library(ggplot2)

pdf("test123.pdf")
ggplot(data=Data_Baseline,aes(x=Disease_Duration_Latency, y=Latenz, group = Verlaufsform, color=Verlaufsform, fill=Verlaufsform)) +
  geom_smooth(method = "lm", formula = y ~ x, alpha = 0.15) +
  geom_point(aes(shape=Verlaufsform)) +
  theme_classic()
dev.off()

custom.summary(lme4::lmer(globalpRNFL ~ Disease_Duration_FINAL + ONHistory +  (1 | ID), data= Data_Baseline, subset =  Verlaufsform == 3))
