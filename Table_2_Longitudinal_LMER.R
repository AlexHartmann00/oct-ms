#------------------Mixed models analyses for prediction----------------------#

# This script computes logistic mixed-effects regression (LMER) to analyze the 
# potential of retinal layer thicknesses to predict disability progression 
# (EDSS worsening), MRI progression/activity (new or enlarging T2-weighted/
# gadolinium-enhancing lesions), and relapses.
# The results are presented in Table 2. 
# Usual runtime ~= 12 minutes on 16 CPU cores

#---------------------------------R packages---------------------------------# 
#Set up parallelization
library(doParallel)


#-----------------------Data loading & data management-----------------------# 

data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)

#Filter OCTs
data <- data[data$OCT_included == 1,]
data <- data[data$Verlaufsform > 0,]

indvars <- c("RNFL_new","globalpRNFL","GCIPL_new","INL_new")





cl <- makeCluster(16)
registerDoParallel(cl)

scale.ifNumeric <- function(X,covs){
  # Function to scale (standardize) all numeric columns of a data frame "X"
  X <- as.data.frame(X)
  X <- as.matrix(X,ncol=length(covs))
  xnew <- data.frame(X)
  #Iterate over columns. If numeric, standardize it
  factors <- c()
  for(i in 1:ncol(X)){
    if(!all(is.na(as.numeric(X[,i])))){
      xnew[,i] <- scale(as.numeric(X[,i]))
    }
  }
  colnames(xnew) <- covs
  return(xnew)
}

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

step.Lmer <- function(X,depvar,idvar,idcolumn,covs){
  #Function to perform full grid search for best model configuration
  #Naively brute forces the search across all possible covariate combinations
  X <- as.data.frame(X)[c(depvar, idvar, covs, "ID")]
  #print(X[1:5,])
  X <- X[complete.cases(X),]
  covariate_combinations <- build_covariates_list(covs)
  
  
  #Keep counter for list appending
  i <- 1
  
  #Get participant ID column
  ids <- X$ID#[,idcolumn]
  #Get focal predictor
  IV <- scale(X[,idvar])
  #Scale input data
  data_ <- suppressWarnings(scale.ifNumeric(X[,covs],covs))
  #Build final input DF
  data_$IV <- IV
  data_$depvar <- X[,depvar]
  print(table(data_$depvar))
  data_$ids <- ids
  
  #Initialize with a null model
  bestmod <- lme4::glmer("depvar ~ IV + (1|ids)",data=data_,family="binomial",control = lme4::glmerControl(optimizer="bobyqa"))
  bestaic <- AIC(bestmod)
  bestcovariates <- c()
  
  # If there should only be positive / negative outcomes, skip computation
  if(length(table(X[,depvar])) == 1){
    return(list(Model=bestmod,Covariates = bestcovariates,AIC=bestaic,Data=X,DV=dv))
  }
  
  
  for(covariates in covariate_combinations){
    #Iterate over all combinations of covariates
    
    
    #Build model formula by appending covariates
    formula_string <- "depvar ~ IV + "
    for(covariate in covariates){
      formula_string <- paste(formula_string, covariate, "+")
    }
    
    formula_string <- paste(formula_string, " (1|ids)")
    
    #Fit model
    convergence_error = F
    tryCatch(
      {mod <- lme4::glmer(formula(formula_string),data=data_,family="binomial",control = lme4::glmerControl(optimizer="bobyqa"))},
      error=function(cond){
        improvement <- F
        convergence_error <- T
      }
    )
    
    if(i > 1 && !convergence_error){
      #If i == 1, this is always an improvement
      #Check chi2 for significance
      p_value <- anova(mod, bestmod)$`Pr(>Chisq)`[2] 
      #Bonferroni correction
      threshold <- 0.05 #/ length(covariate_combinations)
      
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

format_result <- function(stepmod,depvar="",indvar=""){
  #Small function to format results for table generation
  mod <- stepmod$Model
  X <- stepmod$Data
  dv <- stepmod$DV
  if("glmerMod" %in% class(mod)){
    hr <- lme4::fixef(mod)[2]
  }
  else{
    hr <- coef(mod)[2]
  }
  se <- summary(mod)$coefficients[2,2]
  hrlow <- hr - 1.96*se
  hrhigh <- hr + 1.96*se
  p <- summary(mod)$coefficients[2,4]
  covs <- stepmod$Covariates
  cell1 <- paste(round(exp(hr),3)," (",round(exp(hrlow),3),"-",round(exp(hrhigh),3),")",sep="")
  cell2 <- round(p,3)
  cell3 <- paste(round(100*mean(dv,na.rm=T),2),"% (", sum(dv,na.rm=T),"/",length(dv),")",sep="")
  covs[covs %in% c("mri_time_to_assessment","relapse_time_to_assessment","edss_time_to_assessment")] <- "TTA"
  covs[covs == "ONHistory"] <- "ON"
  covs[covs == "Disease_Duration_FINAL"] <- "DD"
  cell4 <- paste(covs,collapse=", ")
  
  unit_size <- sd(X[,indvar], na.rm=T)
  print(unit_size)
  print(hr)
  cell5 <- exp(hr)^(1/unit_size)
  cell6 <- (1/exp(hr))^(1/unit_size)
  c(cell1,cell2,cell3,cell4,cell5,cell6)
}

#EDSS worsening definition
library(tidyverse)

last_edss_df <- data %>% group_by(ID) %>% 
  summarize(LastEDSS = last(EDSS, na_rm=T), lastOCT = max(Sequenz))

max_following_edss <- function(id, sequenz){
  subdata <- data[data$ID == id[1] & data$Sequenz > sequenz[1],]
  if(nrow(subdata[!is.na(subdata$EDSS),]) == 0){
    return(NA)
  }
  return(max(subdata$EDSS,na.rm=T))
}

is_worsening <- function(edss_old, edss_new){
  if(edss_old < 6){
    return(edss_new >= edss_old + 1)
  }
    
  return(edss_new >= edss_old + 0.5)
    
  
}

date_of_first_worsening <- function(id, sequenz){
  subdata <- data[data$ID == id[1] & data$Sequenz >= sequenz[1],]
  if(nrow(subdata[!is.na(subdata$EDSS),]) == 0){
    return(NA)
  }
  worsenings <- is_worsening(first(subdata$EDSS, na_rm=T), subdata$EDSS)

  if(!any(worsenings,na.rm=T)){
    return(NA)
  }
  date <- min(subdata$OCTDatum[worsenings],na.rm=T)
  if(is.infinite(date)){
    return(NA)
  }
  return(date)
}

max_oct_dates <- data %>% group_by(ID) %>% 
  summarize(Followup_date = max(OCTDatum,na.rm=T))

max_following_edss_df <- data %>% 
  inner_join(
    max_oct_dates,
    by="ID") %>%
  group_by(ID, Sequenz) %>% 
  summarize(
    MaxEDSS = max_following_edss(ID, Sequenz),
    EDSS_Worsening_date = date_of_first_worsening(ID, Sequenz),
    EDSS_Worsening_tta = as.numeric(EDSS_Worsening_date - first(OCTDatum)),
    Followup_tta = as.numeric(Followup_date - first(OCTDatum)))

max_following_edss_df <- max_following_edss_df[!duplicated(max_following_edss_df$ID),]

data <- data %>%
  left_join(last_edss_df,by="ID") %>% 
    left_join(max_following_edss_df, by=c("ID","Sequenz"))

data$EDSS_Worsening_Until_Followup <- ifelse(data$EDSS < 6,
                                             data$LastEDSS - data$EDSS >= 1,
                                             data$LastEDSS - data$EDSS >= 0.5)

data$EDSS_Worsening_Until_Followup[data$lastOCT == data$Sequenz | is.na(data$Followup_tta)] <- NA

data$EDSS_Worsening_Any_Future <- ifelse(data$EDSS < 6,
                                             data$MaxEDSS - data$EDSS >= 1,
                                             data$MaxEDSS - data$EDSS >= 0.5)

#Use follow up tta for cases without worsening
no_worsening <- !is.na(data$EDSS_Worsening_Any_Future) & !data$EDSS_Worsening_Any_Future
data$EDSS_Worsening_tta[no_worsening] <- data$Followup_tta[no_worsening]

#############################################
################## MRI ######################
#############################################

mri_data_long <- tidyr::pivot_longer(data,c("MRTAktivität1", "MRTAktivität2", "MRTAktivität3", "MRTAktivität4"),names_to = "MRI_Type",values_to = "MRIActivity")

mri_data_long$MRIActivity <- ifelse(is.na(mri_data_long$MRIActivity),
                                    NA,
                                    ifelse(mri_data_long$MRIActivity %in% c(1,2,3),
                                           1,
                                           0))

#Compute times to assessment

mri_time_to_assessment <- matrix(nrow=nrow(data),ncol=4)
mridatecols <- c("DatumMRT1",  "DatumMRT2",  "DatumMRT3",  "DatumMRT4")
for(i in 1:nrow(mri_time_to_assessment)){
  for(j in 1:length(mridatecols)){
    mri_time_to_assessment[i,j] <- as.Date(data[i,mridatecols[j]]) - data$OCTDatum[i]
  }
}
mri_time_to_assessment <- tidyr::pivot_longer(as.data.frame(mri_time_to_assessment),1:4)
mri_time_to_assessment <- mri_time_to_assessment$value/365

mri_data_long$mri_time_to_assessment <- mri_time_to_assessment

#Initialize some extra variables

mri_data_long$ONHistory <- ifelse(mri_data_long$RBNinVergangenheit == 1,
                                    "RBN",
                                      ifelse(mri_data_long$HistoryofON_GCIPL == 1, "Any","None"))

#Category therapies with respect to their efficacy
mri_data_long$DMT <- ifelse(mri_data_long$Therapie %in% c(0,3,13,14,15),0,
                            ifelse(mri_data_long$Therapie %in% c(10,16,6,5),1,
                                   ifelse(mri_data_long$Therapie %in% c(26,23) | is.na(mri_data_long$Therapie),NA,2)))

#Classification of DMT: 
#DMT=0 --> no therapy: 0, 3, 13, 14, 15
#DMT=1 --> low efficacy: 10, 16, 6, 5
#DMT=2 --> high efficacy: 1, 2, 4, 7, 8, 9, 11, 12, 17, 19, 20, 21, 22, 24, 25


#Create baseline EDSS variable
if(!("BaseEDSS" %in% colnames(mri_data_long))){
  BaseEDSS <- numeric(nrow(mri_data_long))
  for(i in 1:nrow(mri_data_long)){
    id <- mri_data_long$ID[i]
    subdata <- mri_data_long[mri_data_long$ID == id,]
    BaseEDSS[i] <- subdata$EDSS[subdata$OCTDatum == min(subdata$OCTDatum,na.rm=T)][1]
  }
  mri_data_long$BaseEDSS <- BaseEDSS
}

#Create subsets of participants

Data_NoRbn_long <- subset(mri_data_long, OCT_included == 1 & ONHistory == "None")
Data_LMM_RRMS_long <- subset(mri_data_long, OCT_included==1 & Verlaufsform==1)
Data_LMM_PPMS_long <- subset(mri_data_long, OCT_included==1 & Verlaufsform==2)
Data_LMM_SPMS_long <- subset(mri_data_long, OCT_included==1 & Verlaufsform==3)
Data_LMM_RRMSnoRBN_long <- subset(mri_data_long, OCT_included==1 & Verlaufsform == 1 & ONHistory == "None")
Data_LMM_SPMSnoRBN_long <- subset(mri_data_long, OCT_included==1 & Verlaufsform == 3 & ONHistory == "None")

datasets <- list(subset(mri_data_long, OCT_included == 1)
                 ,Data_NoRbn_long,
                 Data_LMM_RRMS_long,
                 Data_LMM_RRMSnoRBN_long,
                 Data_LMM_SPMS_long,
                 Data_LMM_SPMSnoRBN_long,
                 Data_LMM_PPMS_long)

#Get best model for all "participant group" ~ "focal predicor"-combinations
i <- 1
combination_list <- list()
for(indvar in indvars){
  for(dataset in datasets){
    combination_list[[i]] <- list(var=indvar, data=dataset)
    i <- i + 1
  }
}

mriresults <- foreach(element=combination_list) %dopar% {
  format_result(step.Lmer(
    element$data,
    "MRIActivity",
    element$var,
    5,
    c("Age","Geschlecht","mri_time_to_assessment", "EDSS") #"ONHistory",
  ),"MRIActivity",element$var)
}

mriresults_mixmods <- do.call(rbind.data.frame,mriresults)

#The process as descibed for MRI is now repeated for relapse and EDSS.


###########################################
################# Relapse ###################
###########################################

relapse_data <- data

relapsedatecols <- c("DatumSchub1", "DatumSchub2", "DatumSchub3", "DatumSchub4")

relapse <- numeric(nrow(relapse_data))
relapse_time_to_assessment <- numeric(nrow(relapse_data))

for(i in 1:nrow(relapse_data)){
  #For each OCT...
  id <- relapse_data$ID[i]
  octdate <- relapse_data$OCTDatum[i]
  #get all OCTs related to the current patient
  subdata <- relapse_data[relapse_data$ID == id,]
  
  if(all(is.na(relapse_data[i,relapsedatecols]))){
    # If all relapse dates for this OCT are NA, there was no relapse related to this OCT
    relapse[i] <- 0
    #Time to assessment is then the next OCT after the current one minus the current OCT
    upcoming_octs <- subdata$OCTDatum[subdata$OCTDatum > octdate]
    if(length(upcoming_octs) == 0){
      nextoct <- NA
    }else{
      nextoct <- min(upcoming_octs,na.rm=T)
    }
    
    relapse_time_to_assessment[i] <- ifelse(is.na(nextoct),NA, nextoct - octdate)
  }
  else{
    # If there are non-missing relapse dates, a relapse must have occurred related to this OCT
    max_relapse_date <- octdate
    for(datecol in relapsedatecols){
      # Iterate over all four possible relapse dates
      if(!is.na(relapse_data[i,datecol]) & relapse_data[i,datecol] > max_relapse_date){
        # If relapse date is not NA, a relapse happened at this date. If the date is later than the current max relapse date, overwrite it.
        max_relapse_date <- relapse_data[i,datecol]
      }
    }
    #A relapse happened, and the time to assessment is date - oct date
    relapse[i] <- 1
    relapse_time_to_assessment[i] <- max_relapse_date - octdate
  }
}
#Replace infinities and assign calculated variables
relapse_time_to_assessment[is.infinite(relapse_time_to_assessment)] <- NA
relapse_data$Relapse <- relapse
relapse_data$relapse_time_to_assessment <- relapse_time_to_assessment/365

#Initialize some extra variables
relapse_data$ONHistory <- ifelse(relapse_data$RBNinVergangenheit == 1,
                                 "RBN",
                                 ifelse(
                                   relapse_data$HistoryofON_GCIPL == 1,
                                   "Any",
                                   "None"))

# Extract information whether patient was on high-efficacy DMT
relapse_data$DMT <- ifelse(relapse_data$Therapie %in% c(0,3,13,14,15),0,
                           ifelse(relapse_data$Therapie %in% c(10,16,6,5),1,
                                  ifelse(relapse_data$Therapie %in% c(26,23) | is.na(relapse_data$Therapie),NA,2)))

#Classification of DMT: 
#DMT=0 --> no therapy: 0, 3, 13, 14, 15
#DMT=1 --> low efficacy: 10, 16, 6, 5
#DMT=2 --> high efficacy: 1, 2, 4, 7, 8, 9, 11, 12, 17, 19, 20, 21, 22, 24, 25

# Generate baseline EDSS column
if(!("BaseEDSS" %in% colnames(relapse_data))){
  baseline_edss <- numeric(nrow(relapse_data))
  for(i in 1:nrow(relapse_data)){
    #For each row in the dataset, get the first EDSS measurement for the corresponding patients
    id <- relapse_data$ID[i]
    subdata <- relapse_data[relapse_data$ID == id,]
    baseline_edss[i] <- subdata$EDSS[subdata$OCTDatum == min(subdata$OCTDatum,na.rm=T)][1]
  }
  relapse_data$BaseEDSS <- baseline_edss
}

#Create subsets of participants

Data_NoRbn <- subset(relapse_data, OCT_included == 1 & ONHistory == "None" & Verlaufsform != 2)
Data_LMM_RRMS <- subset(relapse_data, OCT_included==1 & Verlaufsform==1)
Data_LMM_SPMS <- subset(relapse_data, OCT_included==1 & Verlaufsform==3)
Data_LMM_RRMSnoRBN <- subset(relapse_data, OCT_included==1 & Verlaufsform==1 & ONHistory == "None")
Data_LMM_SPMSnoRBN <- subset(relapse_data, OCT_included==1 & Verlaufsform==3 & ONHistory == "None")

datasets <- list(subset(relapse_data, OCT_included==1),
                 Data_NoRbn,
                 Data_LMM_RRMS,
                 Data_LMM_RRMSnoRBN,
                 Data_LMM_SPMS,
                 Data_LMM_SPMSnoRBN)

#Get best model for all "participant group" ~ "focal predictor"-combinations
i <- 1
combination_list <- list()
for(indvar in indvars){
  for(dataset in datasets){
    combination_list[[i]] <- list(var=indvar, data=dataset)
    i <- i + 1
  }
}

relapseresults <- foreach(element=combination_list) %dopar% {
  format_result(step.Lmer(
    element$data,
    "Relapse",
    element$var,
    5,
    c("Age","Geschlecht","relapse_time_to_assessment","EDSS","DMT") #"ONHistory", 
  ),
  "Relapse",element$var)
}

relapseresults_mixmods <- do.call(rbind.data.frame,relapseresults)

###########################################
################## EDSS ###################
###########################################

edss_data <- relapse_data

#Rename TTA variable:
edss_data$edss_time_to_assessment <- edss_data$EDDSSchubOCTTimeDiff

Data_NoRbn <- subset(edss_data, OCT_included == 1 & ONHistory == "None")
Data_LMM_RRMS <- subset(edss_data, OCT_included==1 & Verlaufsform==1)
Data_LMM_PPMS <- subset(edss_data, OCT_included==1 & Verlaufsform==2)
Data_LMM_SPMS <- subset(edss_data, OCT_included==1 & Verlaufsform==3)
Data_LMM_RRMSnoRBN <- subset(edss_data, OCT_included==1 & Verlaufsform==1 & ONHistory == "None")
Data_LMM_SPMSnoRBN <- subset(edss_data, OCT_included==1 & Verlaufsform==3 & ONHistory == "None")

datasets <- list(subset(edss_data,OCT_included==1)
                 ,Data_NoRbn,
                 Data_LMM_RRMS,
                 Data_LMM_RRMSnoRBN,
                 Data_LMM_SPMS,
                 Data_LMM_SPMSnoRBN,
                 Data_LMM_PPMS)

# Extract information whether patient was on high-efficacy DMT
edss_data$DMT <- ifelse(edss_data$Therapie %in% c(0,3,13,14,15),0,
                           ifelse(edss_data$Therapie %in% c(10,16,6,5),1,
                                  ifelse(edss_data$Therapie %in% c(26,23) | is.na(edss_data$Therapie),NA,2)))

#Classification of DMT: 
#DMT=0 --> no therapy: 0, 3, 13, 14, 15
#DMT=1 --> low efficacy: 10, 16, 6, 5
#DMT=2 --> high efficacy: 1, 2, 4, 7, 8, 9, 11, 12, 17, 19, 20, 21, 22, 24, 25


#Get best model for all "participant group" ~ "focal predictor"-combinations
i <- 1
combination_list <- list()
for(indvar in indvars){
  for(dataset in datasets){
    combination_list[[i]] <- list(var=indvar, data=dataset)
    i <- i + 1
  }
}

edssresults <- foreach(element=combination_list) %dopar% {
  format_result(step.Lmer(
    element$data,
    "EDSS_Worsening_Until_Followup",
    element$var,
    5,
    c("Age","Geschlecht","Followup_tta","DMT") #"ONHistory", 
  ),"EDSS_Worsening_Until_Followup",element$var)
}

edssresults_mixmods <- do.call(rbind.data.frame,edssresults)

###############################################################
##################EDSS with any worsening in future###########
###############################################################


#Get best model for all "participant group" ~ "focal predictor"-combinations
i <- 1
combination_list <- list()
for(indvar in indvars){
  for(dataset in datasets){
    combination_list[[i]] <- list(var=indvar, data=dataset)
    i <- i + 1
  }
}

edssresults_any <- foreach(element=combination_list) %dopar% {
  format_result(step.Lmer(
    element$data,
    "EDSS_Worsening_Any_Future",
    element$var,
    5,
    c("Age","Geschlecht","EDSS_Worsening_tta", "DMT") #"ONHistory", 
  ),"EDSS_Worsening_Any_Future",element$var)
}

edssresults_any_mixmods <- do.call(rbind.data.frame, edssresults_any)

#############################################################
#########################Format tables ######################
#############################################################

colnames(relapseresults_mixmods) <- c("Odds Ratio","p-value","Prob. Event","Covariates","unit_increase","unit_decrease")
colnames(edssresults_mixmods) <- c("Odds Ratio","p-value","Prob. Event","Covariates","unit_increase","unit_decrease")
colnames(edssresults_any_mixmods) <- c("Odds Ratio","p-value","Prob. Event","Covariates","unit_increase","unit_decrease")
colnames(mriresults_mixmods) <- c("Odds Ratio","p-value","Prob. Event","Covariates","unit_increase","unit_decrease")

View(relapseresults_mixmods)
View(edssresults_mixmods)
View(edssresults_any_mixmods)
View(mriresults_mixmods)


# Descriptives for table 2 input data
library(tidyverse)
mri_data_long %>% group_by(ID) %>% 
  summarize(
    disease_course = first(Verlaufsform),
    MRI = ifelse(max(MRIActivity, na.rm=T) > 0,1,0)) %>% 
  group_by(disease_course) %>% 
  summarize(
    activity = sum(MRI),
    total = n(),
    ratio = activity / total
  )

relapse_data %>% group_by(ID) %>% 
  summarize(
    disease_course = first(Verlaufsform),
    Relapse = ifelse(max(Relapse, na.rm=T) > 0,1,0)) %>% 
  group_by(disease_course) %>% 
  summarize(
    activity = sum(Relapse),
    total = n(),
    ratio = activity / total
  )

edss_data %>% group_by(ID) %>% 
  summarize(
    disease_course = first(Verlaufsform),
    EDSS = ifelse(max(EDSS_Worsening_Until_Followup, na.rm=T) > 0,1,0)) %>% 
  group_by(disease_course) %>% 
  summarize(
    activity = sum(EDSS),
    total = n(),
    ratio = activity / total
  )

median(edss_data$EDDSSchubOCTTimeDiff/365, na.rm=T)

edss_data %>% group_by(ID) %>% 
  summarize(
    disease_course = first(Verlaufsform),
    EDSS = ifelse(max(EDSS_Worsening_Any_Future, na.rm=T) > 0,1,0)) %>% 
  group_by(disease_course) %>% 
  summarize(
    activity = sum(EDSS),
    total = n(),
    ratio = activity / total
  )








