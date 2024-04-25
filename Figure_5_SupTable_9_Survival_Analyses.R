#------------------Cox regressions and Kaplan Meier analyses-----------------#

# This script computes Cox regressions and Kaplan Meier analyses to analyze the 
# potential of retinal layer thicknesses to predict disability progression 
# (EDSS worsening), MRI progression/activity (new or enlarging T2-weighted/
# gadolinium-enhancing lesions), and relapses.
# The results are presented in Supplementary Table 8 (Cox regressions) and 
# Figure 5 (Kaplan Meier analyses)
# Usual runtime ~= 3 minutes


#setwd("C:/octms5")

#---------------------------------R packages---------------------------------# 

library(survival)
library(tidyverse)
library(patchwork)


#-----------------------Data loading & data management-----------------------# 

data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)
data <- data[data$OCT_included == 1,]
data <- data[data$Verlaufsform != 0,]

data %>% group_by(Verlaufsform) %>% 
  summarize(missings = sum(is.na(Therapie)))

# ON was categorized in the categories 
# a) ON based on medical history ("RBN")
# b) ON based on inter-eye GCIPL thickness difference ("Any")
# c) no ON ("None")

data$ONHistory <- ifelse(
  data$RBNinVergangenheit == 1,
  "RBN",
  ifelse(
    data$HistoryofON_GCIPL == 1,
    "Any",
    "None"))

data_b <- data[data$Sequenz == 1,]

#Thresholds for retinal layer thickness tertiles at baseline
gciplthresh <- quantile(data_b$GCIPL_new,1/3,na.rm=T)
prnflthresh <- quantile(data_b$globalpRNFL,1/3,na.rm=T)
mrnflthresh <- quantile(data_b$RNFL_new,1/3,na.rm=T)
inlthresh <- quantile(data_b$INL_new,1/3,na.rm=T)

#Print
print(gciplthresh)
print(prnflthresh)
print(mrnflthresh)
print(inlthresh)

#-- Define general functions
kaplan.meier <- function(survdf,groupvar){
  #Function to compute Kaplan-Meier statistic
  survdf <- survdf[order(survdf$Time),]
  survdf <- survdf[!is.na(survdf[,groupvar]),]
  grouping_variable <- survdf[,groupvar]
  S <- numeric(nrow(survdf))
  for(group in unique(grouping_variable)){
    survdf_g <- survdf[grouping_variable == group,]
    atrisk <- nrow(survdf_g)
    props <- c()
    for(i in 1:nrow(survdf_g)){
      if(survdf_g$Status[i] == 1){
        p <- (atrisk - 1)/atrisk
      }
      else{
        p <- 1
      }
      props <- c(props,p)
      atrisk <- atrisk - 1
      S[survdf$ID == survdf_g$ID[i]] <- cumprod(props)[length(props)]
    }
    propsprod <- cumprod(props)
  }
  
  kmdf <- cbind(S=S,Time=survdf$Time,KMGroup = survdf[,groupvar])#survdf,
  kmdf <- as.data.frame(kmdf)
  
  kmdf <- rbind.data.frame(
    kmdf,
    data.frame(
      S=c(1,1),
      Time = c(0,0),
      KMGroup = c(0,1)
    )
  )
  
  colnames(kmdf)[ncol(kmdf)] <- stringr::str_replace(groupvar,"_Group","")
  return(kmdf)
}


KMCurve <- function(dataset, groupvar){
  #Function to display Kaplan-Meier-Curve
  kmd <- kaplan.meier(dataset,groupvar)
  
  kaplan.meier.curve <- function(kmdf){
    groupname <- colnames(kmdf)[ncol(kmdf)]
    kmdf[,groupname] <- ifelse(kmdf[,groupname] == 1, "low","high")
    ggplot(data=kmdf,aes(Time/365,1-S,color=factor(kmdf[,groupname]),shape=factor(kmdf[,groupname])))+
      geom_line()+
      geom_point()+
      scale_x_continuous(breaks=0:9)+
      scale_color_manual(values=c("low"="#bf1755","high"="#2cced4"))+
      scale_shape_manual(values=c("low"=16,"high"=17))+
      labs(x="Follow-up time in years",shape=paste(groupname,"thickness"),color=paste(groupname,"thickness"),y="Cumulative Worsening %")+
      scale_y_continuous(labels=scales::percent)+
      theme_classic()
  }
  
  kaplan.meier.curve(kmd)
}

build_covariates_list <- function(covariates){
  #Function to collect all combinations of a list of covariates
  #Basically builds the grid to search through
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

step.CoxReg <- function(surv,idvar,cov_list){
  # Function to perform grid search for best features to use for model fit
  # Selects model by likelihood ratio
  covariate_combinations <- build_covariates_list(cov_list)#c("Sex","Age","DMT")) # ,"ON", "DD",
  current_covariates <- c()
  surv <- surv[complete.cases(surv[,c("Time","Status",cov_list)]),] #,"ON", "DD",
  response <- Surv(surv$Time, surv$Status)

  #Initialize baseline model
  basemod <- coxph(response ~ surv[,idvar])
  aic <- AIC(basemod)
  bestmod <- basemod
  bestaic <- aic
  event_num <- sum(surv$Status)
  
  #Perform grid search
  for(covariate in covariate_combinations){
    #print(covariate)
    if(event_num/(1+length(covariate)) >= 10){
      formula_string = paste("response ~ ",idvar)
      for(cov in covariate){
        formula_string <- paste(formula_string,"+",cov)
      }
      
      mod <- coxph(formula(formula_string),data=surv)
      
      lr_test <- anova(bestmod,mod)
      p <- lr_test$`Pr(>|Chi|)`[2]
      threshold <- 0.05 

      if(p < threshold){
        bestaic <- AIC(mod)
        bestmod <- mod
        current_covariates <- covariate
      }
      
    }
  }
  return(list(Model=bestmod,
              AIC = bestaic,
              Covariates = current_covariates,
              Data = surv))
}

format_result <- function(stepmod){
  coxmod <- stepmod$Model
  X <- stepmod$Data
  X <- X[!is.na(X$Status),]
  hr <- coxmod$coefficients[1]
  se <- summary(coxmod)$coefficients[1,3]
  hrlow <- hr - 1.96*se
  hrhigh <- hr + 1.96*se
  p <- summary(coxmod)$coefficients[1,5]
  covs <- paste(stepmod$Covariates,collapse=", ")
  cell1 <- paste(round(exp(hr),3)," (",round(exp(hrlow),3),"-",round(exp(hrhigh),3),")",sep="")
  cell2 <- round(p,3)
  cell3 <- paste(round(100*mean(X$Status),2),"% (", sum(X$Status),"/",nrow(X),")",sep="")
  cell4 <- covs
  return_vec <- c(cell1,cell2,cell3,cell4)
  names(return_vec) <- c("b(CI)","p","#event","covariates")
  return(return_vec)
}



# Function to generate survival data frame
get_surv_df <- function(input_data, tta, outcome){
   #Takes in data frame, time to assessment and outcome variable names
  N <- length(unique(input_data$ID))
  
  if(!is.null(tta)){
    input_data$TTA <- input_data[,tta]    
  }else{
    input_data$TTA <- 0
  }

  input_data$Outcome <- input_data[,outcome]
  
  ids <- character(N)
  time <- numeric(N)
  status <- numeric(N)
  group <- numeric(N)
  uniqueids <- unique(input_data$ID)
  gcipls <- numeric(N)
  prnfls <- numeric(N)
  mrnfls <- numeric(N)
  inls <- numeric(N)
  ons <- character(N)
  ages <- numeric(N)
  dds <- numeric(N)
  sex <- numeric(N)
  baseedss <- numeric(N)
  dmt <- numeric(N)
  i <- 1

  
  for(id in uniqueids){

    ids[i] <- id
    subdata <- input_data[input_data$ID == id,]

    group[i] <- first(subdata$Verlaufsform, na_rm=T)
    ages[i] <- min(subdata$Age,na.rm=T)
    dds[i] <- min(subdata$Disease_Duration_FINAL,na.rm=T)
    sex[i] <- first(subdata$Geschlecht, na_rm=T)
    baseedss[i] <- first(subdata$BaseEDSS, na_rm=T)
    dmt[i] <- first(subdata$DMT,na_rm=T)
    gcipls[i] <- mean(subdata$GCIPL_new[1:2],na.rm=T)
    prnfls[i] <- mean(subdata$globalpRNFL[1:2],na.rm=T)
    mrnfls[i] <- mean(subdata$RNFL_new[1:2],na.rm=T)
    inls[i] <- mean(subdata$INL_new[1:2],na.rm=T)
    ons[i] <- first(subdata$ONHistory, na_rm=T)
    
    
    
    if(all(is.na(subdata$Outcome))){
      time[i] <- NA
      status[i] <- 0
    }
    else{
      subdata <- subdata[!is.na(subdata$Outcome) & !is.na(subdata$TTA),]
      if(nrow(subdata) == 0){
        time[i] <- NA
        status[i] <- 0
      }
      else{
        if(all(subdata$Outcome == 0)){
          time[i] <- max(365*subdata$TTA + subdata$OCTDatum) - min(subdata$OCTDatum)
          status[i] <- 0
        }
        else{
          subsubdata <- subdata[subdata$Outcome > 0,]
          time[i] <- min(365*subsubdata$TTA + subsubdata$OCTDatum) - min(subdata$OCTDatum)
          status[i] <- 1
        }
      }
      
    }
    
    i <- i+1
  }
  
  survdata <- data.frame(ID = ids, Time = time, Status = status, Group = group,
                         GCIPL_raw = gcipls, GCIPL_Group = gcipls < gciplthresh,
                         pRNFL_raw = prnfls, pRNFL_Group = prnfls < prnflthresh,
                         mRNFL_raw = mrnfls, mRNFL_Group = mrnfls < mrnflthresh,
                         INL_raw = inls, INL_Group = inls < inlthresh,
                         Sex = sex,
                         ON = ons, Age=ages, DD=dds,
                         BaseEDSS = baseedss, DMT = dmt) 
  
  survdata <- survdata[!is.na(survdata$Group),]
  return(survdata)
}


#Compute baseline EDSS scores
BaseEDSS <- numeric(nrow(data))
for(i in 1:nrow(data)){
  id <- data$ID[i]
  subdata <- data[data$ID == id,]
  BaseEDSS[i] <- subdata$EDSS[subdata$OCTDatum == min(subdata$OCTDatum,na.rm=T)][1]
}
data$BaseEDSS <- BaseEDSS

#Category therapies with respect to their efficacy
data$DMT <- ifelse(data$Therapie %in% c(0,3,13,14,15),0,
                            ifelse(data$Therapie %in% c(10,16,17,6,5),1,
                                   ifelse(data$Therapie %in% c(26,23) | is.na(data$Therapie),NA,2)))



# Preprocessing for MRI
# Pivot in order to consider each MRI measurement after one OCT and before the next OCT
Data_MRI_long <- tidyr::pivot_longer(data,
                                     c("MRTAktivität1", "MRTAktivität2", "MRTAktivität3", "MRTAktivität4"),
                                     names_to = "MRI_Type",
                                     values_to = "MRIActivity")

# Compute MRI times to assessment 
mri_time_to_assessment <- matrix(nrow=nrow(data),ncol=4)
mridatecols <- c("DatumMRT1",  "DatumMRT2",  "DatumMRT3",  "DatumMRT4")
for(i in 1:nrow(mri_time_to_assessment)){
  for(j in 1:length(mridatecols)){
    mri_time_to_assessment[i,j] <- as.Date(data[i,mridatecols[j]]) - as.Date(data$OCTDatum[i])
  }
}
mri_time_to_assessment <- tidyr::pivot_longer(as.data.frame(mri_time_to_assessment),1:4)
mri_time_to_assessment <- mri_time_to_assessment$value/365

Data_MRI_long <- cbind(Data_MRI_long,mri_time_to_assessment)

#Get binary MRI activity 
Data_MRI_long$MRIActivity <- ifelse(Data_MRI_long$MRIActivity > 0,1,0)
Data_MRI_long <- Data_MRI_long[!is.na(Data_MRI_long$MRIActivity),]

#------------------------------------MRI--------------------------------------# 
survdata <- get_surv_df(Data_MRI_long, "mri_time_to_assessment", "MRIActivity")

#----------------------------MRI Cox-Regressions------------------------------#
focal_predictors <- c("pRNFL_Group","mRNFL_Group","GCIPL_Group","INL_Group")
Data_NoRbn_long <- subset(Data_MRI_long, OCT_included == 1 & ONHistory == "None")
Data_LMM_RRMS_long <- subset(Data_MRI_long, OCT_included==1 & Verlaufsform==1)
Data_LMM_PPMS_long <- subset(Data_MRI_long, OCT_included==1 & Verlaufsform==2)
Data_LMM_SPMS_long <- subset(Data_MRI_long, OCT_included==1 & Verlaufsform==3)
Data_LMM_RRMSnoRBN_long <- subset(Data_MRI_long, OCT_included==1 & Verlaufsform == 1 & ONHistory == "None")
Data_LMM_SPMSnoRBN_long <- subset(Data_MRI_long, OCT_included==1 & Verlaufsform == 3 & ONHistory == "None")

datasets <- list(Data_MRI_long
                 ,Data_NoRbn_long,
                 Data_LMM_RRMS_long,
                 Data_LMM_RRMSnoRBN_long,
                 Data_LMM_SPMS_long,
                 Data_LMM_SPMSnoRBN_long,
                 Data_LMM_PPMS_long)

mriresults <- list()
i <- 1
for(focal_predictor in focal_predictors){
  for(dataset in datasets){
    survdat <- survdata[survdata$ID %in% dataset$ID,]
    print(nrow(survdat))
    mriresults[[i]] <- format_result(step.CoxReg(survdat,
                                                 focal_predictor,
                                                 c("Sex","Age"))) #"DMT",
    i <- i + 1
  }
}
mriresults_coxregs <- do.call(rbind.data.frame,mriresults)

#----------------------------Kaplan Meier analyses---------------------------#

#------------------------------MRI all patients------------------------------#

prnflmri <- KMCurve(survdata,"pRNFL_Group")
mrnflmri <- KMCurve(survdata,"mRNFL_Group")
gciplmri <- KMCurve(survdata,"GCIPL_Group")
inlmri <- KMCurve(survdata,"INL_Group")

current_plot <- (prnflmri + mrnflmri) / (gciplmri + inlmri)

#-------------------MRI all patients without history of ON-------------------#
# without history of ON means without ON based on medical history and without ON 
# based on inter-eye GCIPL thickness difference 


survdata2 <- survdata[survdata$ON == "None",]
prnflmri <- KMCurve(survdata2,"pRNFL_Group")
mrnflmri <- KMCurve(survdata2,"mRNFL_Group")
gciplmri <- KMCurve(survdata2,"GCIPL_Group")
inlmri <- KMCurve(survdata2,"INL_Group")

current_plot <- (prnflmri + mrnflmri) / (gciplmri + inlmri)
current_plot


#----------------------------MRI only RRMS-----------------------------------#

survdatarrms <- survdata[survdata$Group == 1,]
prnflmri <- KMCurve(survdatarrms,"pRNFL_Group")
mrnflmri <- KMCurve(survdatarrms,"mRNFL_Group")
gciplmri <- KMCurve(survdatarrms,"GCIPL_Group")
inlmri <- KMCurve(survdatarrms,"INL_Group")

current_plot <- (prnflmri + mrnflmri) / (gciplmri + inlmri)
current_plot


#---------------------MRI only RRMS without history of ON-------------------#
# without history of ON means without ON based on medical history and without ON 
# based on inter-eye GCIPL thickness difference 

survdatarrms <- survdata[survdata$Group == 1 & survdata$ON == "None",]
prnflmri <- KMCurve(survdatarrms,"pRNFL_Group")
mrnflmri <- KMCurve(survdatarrms,"mRNFL_Group")
gciplmri <- KMCurve(survdatarrms,"GCIPL_Group")
inlmri <- KMCurve(survdatarrms,"INL_Group")

current_plot <- (prnflmri + mrnflmri) / (gciplmri + inlmri)
current_plot

#----------------------------MRI only PPMS------------------------------------#

survdatappms <- survdata[survdata$Group == 2,]
prnflmri <- KMCurve(survdatappms,"pRNFL_Group")
mrnflmri <- KMCurve(survdatappms,"mRNFL_Group")
gciplmri <- KMCurve(survdatappms,"GCIPL_Group")
inlmri <- KMCurve(survdatappms,"INL_Group")

current_plot <- (prnflmri + mrnflmri) / (gciplmri + inlmri)
current_plot


#----------------------------MRI only SPMS------------------------------------#

survdataspms <- survdata[survdata$Group == 3,]
prnflmri <- KMCurve(survdataspms,"pRNFL_Group")
mrnflmri <- KMCurve(survdataspms,"mRNFL_Group")
gciplmri <- KMCurve(survdataspms,"GCIPL_Group")
inlmri <- KMCurve(survdataspms,"INL_Group")

current_plot <- (prnflmri + mrnflmri) / (gciplmri + inlmri)
current_plot


#---------------------MRI only SPMS without history of ON-------------------#
# without history of ON means without ON based on medical history and without ON 
# based on inter-eye GCIPL thickness difference 

survdataspms <- survdata[survdata$Group == 3 & survdata$ON == "None",]
prnflmri <- KMCurve(survdataspms,"pRNFL_Group")
mrnflmri <- KMCurve(survdataspms,"mRNFL_Group")
gciplmri <- KMCurve(survdataspms,"GCIPL_Group")
inlmri <- KMCurve(survdataspms,"INL_Group")

current_plot <- (prnflmri + mrnflmri) / (gciplmri + inlmri)
current_plot


#----------------Number of at-risk participants per year--------------------#
survdatappms$years <- survdatappms$Time %/% 365

yrs <- sort(unique(survdatappms$years))
vars <- c("pRNFL_Group","mRNFL_Group","GCIPL_Group","INL_Group")

#Initialize matrix. We need 2 (high / low groups) * #variables * #years rows.
N <- length(vars)*2*length(yrs)

retmat <- data.frame(vars=character(N),
                     year=numeric(N),
                     value=character(N),
                     atrisk=numeric(N))

i <- 1
for(yearindex in 1:length(yrs)){
  for(varindex in 1:length(vars)){
    subfirst <- survdatappms[!is.na(survdatappms$years),]
    subdat <- subfirst[subfirst$years >= yrs[yearindex],]
    #high:
    atrisk <- sum(complete.cases(subdat[subdat[,vars[varindex]] == F,]))
    retmat$vars[i] <- vars[varindex]
    retmat$year[i] <- yrs[yearindex]
    retmat$value[i] <- "high"
    retmat$atrisk[i] <- atrisk
    i <- i + 1
    #low:
    atrisk <- nrow(subdat[subdat[,vars[varindex]] == T,])
    retmat$vars[i] <- vars[varindex]
    retmat$year[i] <- yrs[yearindex]
    retmat$value[i] <- "low"
    retmat$atrisk[i] <- atrisk
    i <- i + 1
  }
}

print(retmat)

#--------------------------------EDSS--------------------------------------#
data$edss_diff <- data$EDSS - data$BaseEDSS
data$edss_worsening <- as.numeric(ifelse(data$BaseEDSS < 6, data$edss_diff >= 1, data$edss_diff >= 0.5))

survdata <- get_surv_df(data, NULL, "edss_worsening")

#--------------------------EDSS Cox-Regressions----------------------------#
focal_predictors <- c("pRNFL_Group","mRNFL_Group","GCIPL_Group","INL_Group") ##Group
Data_NoRbn <- subset(data, OCT_included == 1 & ONHistory == "None")
Data_LMM_RRMS <- subset(data, OCT_included==1 & Verlaufsform==1)
Data_LMM_PPMS <- subset(data, OCT_included==1 & Verlaufsform==2)
Data_LMM_SPMS <- subset(data, OCT_included==1 & Verlaufsform==3)
Data_LMM_RRMSnoRBN <- subset(data, OCT_included==1 & Verlaufsform == 1 & ONHistory == "None")
Data_LMM_SPMSnoRBN <- subset(data, OCT_included==1 & Verlaufsform == 3 & ONHistory == "None")

datasets <- list(data
                 ,Data_NoRbn,
                 Data_LMM_RRMS,
                 Data_LMM_RRMSnoRBN,
                 Data_LMM_SPMS,
                 Data_LMM_SPMSnoRBN,
                 Data_LMM_PPMS)

edssresults <- list()
i <- 1
for(focal_predictor in focal_predictors){
  for(dataset in datasets){
    survdat <- survdata[survdata$ID %in% dataset$ID,]
    edssresults[[i]] <- format_result(step.CoxReg(survdat,
                                                  focal_predictor,
                                                  c("Sex","Age"))) #"BaseEDSS", "DMT"
    i <- i + 1
  }
}
edssresults_coxregs <- do.call(rbind.data.frame,edssresults)

#----------------------------Kaplan Meier analyses---------------------------#

#------------------------------EDSS all patients-----------------------------#

prnfledss <- KMCurve(survdata,"pRNFL_Group")
mrnfledss <- KMCurve(survdata,"mRNFL_Group")
gcipledss <- KMCurve(survdata,"GCIPL_Group")
inledss <- KMCurve(survdata,"INL_Group")

current_plot <- (prnfledss + mrnfledss) / (gcipledss + inledss)
current_plot


#-------------------EDSS all patients without history of ON------------------#
# without history of ON means without ON based on medical history and without ON 
# based on inter-eye GCIPL thickness difference 

survdata2 <- survdata[survdata$ON == "None",]
prnfledss <- KMCurve(survdata2,"pRNFL_Group")
mrnfledss <- KMCurve(survdata2,"mRNFL_Group")
gcipledss <- KMCurve(survdata2,"GCIPL_Group")
inledss <- KMCurve(survdata2,"INL_Group")

current_plot <- (prnfledss + mrnfledss) / (gcipledss + inledss)
current_plot


#-----------------------------EDSS only RRMS--------------------------------#

survdatarrms <- survdata[survdata$Group == 1,]
prnfledss <- KMCurve(survdatarrms,"pRNFL_Group")
mrnfledss <- KMCurve(survdatarrms,"mRNFL_Group")
gcipledss <- KMCurve(survdatarrms,"GCIPL_Group")
inledss <- KMCurve(survdatarrms,"INL_Group")

current_plot <- (prnfledss + mrnfledss) / (gcipledss + inledss)
current_plot


#------------------- EDSS only RRMS without history of ON-------------------#
# without history of ON means without ON based on medical history and without ON 
# based on inter-eye GCIPL thickness difference 

survdatarrms <- survdata[survdata$Group == 1 & survdata$ON == "None",]
prnfledss <- KMCurve(survdatarrms,"pRNFL_Group")
mrnfledss <- KMCurve(survdatarrms,"mRNFL_Group")
gcipledss <- KMCurve(survdatarrms,"GCIPL_Group")
inledss <- KMCurve(survdatarrms,"INL_Group")

current_plot <- (prnfledss + mrnfledss) / (gcipledss + inledss)
current_plot


#----------------------------EDSS only PPMS---------------------------------#

survdatappms <- survdata[survdata$Group == 2,]
prnfledss <- KMCurve(survdatappms,"pRNFL_Group")
mrnfledss <- KMCurve(survdatappms,"mRNFL_Group")
gcipledss <- KMCurve(survdatappms,"GCIPL_Group")
inledss <- KMCurve(survdatappms,"INL_Group")

current_plot <- (prnfledss + mrnfledss) / (gcipledss + inledss)
current_plot


#-----------------------------EDSS only SPMS--------------------------------#

survdataspms <- survdata[survdata$Group == 3,]
prnfledss <- KMCurve(survdataspms,"pRNFL_Group")
mrnfledss <- KMCurve(survdataspms,"mRNFL_Group")
gcipledss <- KMCurve(survdataspms,"GCIPL_Group")
inledss <- KMCurve(survdataspms,"INL_Group")

current_plot <- (prnfledss + mrnfledss) / (gcipledss + inledss)
current_plot


#------------------- EDSS only SPMS without history of ON-------------------#
# without history of ON means without ON based on medical history and without ON 
# based on inter-eye GCIPL thickness difference 


survdataspms <- survdata[survdata$Group == 3 & survdata$ON == "None",]
prnfledss <- KMCurve(survdataspms,"pRNFL_Group")
mrnfledss <- KMCurve(survdataspms,"mRNFL_Group")
gcipledss <- KMCurve(survdataspms,"GCIPL_Group")
inledss <- KMCurve(survdataspms,"INL_Group")

current_plot <- (prnfledss + mrnfledss) / (gcipledss + inledss)
current_plot


#---------------------------------Relapse----------------------------------#

#Add features to data

relapsedatecols <- c("DatumSchub1","DatumSchub2","DatumSchub3","DatumSchub4")

relapse <- numeric(nrow(data))
relapse_time_to_assessment <- numeric(nrow(data))

for(i in 1:nrow(data)){
  id <- data$ID[i]
  octdate <- data$OCTDatum[i]
  subdata <- data[data$ID == id,]
  if(all(is.na(data[i,relapsedatecols]))){
    relapse[i] <- 0
    minoct <- min(subdata$OCTDatum[subdata$OCTDatum > octdate],na.rm=T)
    relapse_time_to_assessment[i] <- ifelse(is.na(minoct),NA, minoct - octdate)
  }
  else{
    maxd <- octdate
    for(datecol in relapsedatecols){
      if(!is.na(data[i,datecol]) & data[i,datecol] > maxd){
        maxd <- data[i,datecol]
      }
    }
    relapse[i] <- 1
    relapse_time_to_assessment[i] <- maxd - octdate
  }
}
relapse_time_to_assessment[is.infinite(relapse_time_to_assessment)] <- NA
data$Relapse <- relapse
data$relapse_time_to_assessment <- relapse_time_to_assessment/365

# # Generate survival dataframe for relapse
survdata <- get_surv_df(data, "relapse_time_to_assessment", "Relapse")

#------------------------Relapse Cox-Regressions----------------------------#
Data_NoRbn <- subset(data, OCT_included == 1 & ONHistory == "None" & Verlaufsform != 2)
Data_LMM_RRMS <- subset(data, OCT_included==1 & Verlaufsform==1)
Data_LMM_SPMS <- subset(data, OCT_included==1 & Verlaufsform==3)
Data_LMM_RRMSnoRBN <- subset(data, OCT_included==1 & Verlaufsform == 1 & ONHistory == "None")
Data_LMM_SPMSnoRBN <- subset(data, OCT_included==1 & Verlaufsform == 3 & ONHistory == "None")

datasets <- list(as.data.frame(data)
                 ,as.data.frame(Data_NoRbn_long),
                 as.data.frame(Data_LMM_RRMS_long),
                 as.data.frame(Data_LMM_RRMSnoRBN_long),
                 as.data.frame(Data_LMM_SPMS_long),
                 as.data.frame(Data_LMM_SPMSnoRBN_long))

relapseresults <- list()
i <- 1
for(focal_predictor in focal_predictors){
  for(dataset in datasets){
    survdat <- survdata[survdata$ID %in% dataset$ID,]
    relapseresults[[i]] <- format_result(step.CoxReg(survdat,
                                                     focal_predictor,
                                                     c("Sex","Age"))) #,"BaseEDSS"
    i <- i + 1
  }
}
relapseresults_coxregs <- do.call(rbind.data.frame,relapseresults)


#----------------------------Kaplan Meier analyses---------------------------#

#---------------------------Relapse all patients-----------------------------#

prnflrelapse <- KMCurve(survdata, "pRNFL_Group")
mrnflrelapse <- KMCurve(survdata,"mRNFL_Group")
gciplrelapse <- KMCurve(survdata,"GCIPL_Group")
inlrelapse <- KMCurve(survdata,"INL_Group")

current_plot <- (prnflrelapse + mrnflrelapse) / (gciplrelapse + inlrelapse)
current_plot



#---------------Relapse all patients without history of ON-------------------#
# without history of ON means without ON based on medical history and without ON 
# based on inter-eye GCIPL thickness difference 

survdata2 <- survdata[survdata$ON == "None",]
prnflrelapse <- KMCurve(survdata2,"pRNFL_Group")
mrnflrelapse <- KMCurve(survdata2,"mRNFL_Group")
gciplrelapse <- KMCurve(survdata2,"GCIPL_Group")
inlrelapse <- KMCurve(survdata2,"INL_Group")

current_plot <- (prnflrelapse + mrnflrelapse) / (gciplrelapse + inlrelapse)
current_plot


#---------------------------Relapse only RRMS-------------------------------#

survdatarrms <- survdata[survdata$Group == 1,]
prnflrelapse <- KMCurve(survdatarrms,"pRNFL_Group")
mrnflrelapse <- KMCurve(survdatarrms,"mRNFL_Group")
gciplrelapse <- KMCurve(survdatarrms,"GCIPL_Group")
inlrelapse <- KMCurve(survdatarrms,"INL_Group")

current_plot <- (prnflrelapse + mrnflrelapse) / (gciplrelapse + inlrelapse)
current_plot


#-----------------Relapse only RRMS without history of ON-------------------#
# without history of ON means without ON based on medical history and without ON 
# based on inter-eye GCIPL thickness difference 

survdata2 <- survdata[survdata$RBN != 1 & survdata$Group == 1,]
prnflrelapse <- KMCurve(survdata2,"pRNFL_Group")
mrnflrelapse <- KMCurve(survdata2,"mRNFL_Group")
gciplrelapse <- KMCurve(survdata2,"GCIPL_Group")
inlrelapse <- KMCurve(survdata2,"INL_Group")

current_plot <- (prnflrelapse + mrnflrelapse) / (gciplrelapse + inlrelapse)
current_plot


#----------------------------  Relapse only SPMS----------------------------#

survdata2 <- survdata[survdata$Group == 3,]
prnflrelapse <- KMCurve(survdata2,"pRNFL_Group")
mrnflrelapse <- KMCurve(survdata2,"mRNFL_Group")
gciplrelapse <- KMCurve(survdata2,"GCIPL_Group")
inlrelapse <- KMCurve(survdata2,"INL_Group")

current_plot <- (prnflrelapse + mrnflrelapse) / (gciplrelapse + inlrelapse)
current_plot


#Format

#############################################################
#########################Format tables ######################
#############################################################

colnames(relapseresults_coxregs) <- c("Hazard Ratio","p-value","Prob. Event","Covariates")
colnames(edssresults_coxregs) <- c("Hazard Ratio","p-value","Prob. Event","Covariates")
colnames(mriresults_coxregs) <- c("Hazard Ratio","p-value","Prob. Event","Covariates")

View(relapseresults_coxregs)
View(edssresults_coxregs)
View(mriresults_coxregs)
