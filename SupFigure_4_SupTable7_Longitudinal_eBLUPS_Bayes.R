#------------------empirical Bayes estimates of BLUPs------------------------#
#Script to compute empirical Bayes estimates of BLUPs of annualized thickness
# change rates and group difference statistics.
#To decrease the high variability associated with the small number of observations 
#in each stratum of disease duration, coefficients were estimated over the entire 
#disease duration.
# The results are presented in Supplementary Table 6 and Supplementary Figure 4. 
# Usual runtime ~= 1 minute

#---------------------------------R packages---------------------------------# 

library(Matrix)
library(tidyverse)
library(lme4)
library(ggplot2)


#-----------------------Data loading & data management-----------------------# 

#Load dataset
data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)

#Filter OCTs
data <- data[data$OCT_included==1,]

#Get baseline disease duration
baseline_dds <- data %>% group_by(ID) %>% 
  summarize(bldd = min(Disease_Duration_FINAL,na.rm=T))

row.names(baseline_dds) <- baseline_dds$ID

#Get person base line corrected disease duration
data$pbl_DD <- data$Disease_Duration_FINAL - baseline_dds[data$ID,]$bldd

data$bl_DD <- baseline_dds[data$ID,]$bldd


#---------------------------------Estimation model---------------------------# 

vars <- c("globalpRNFL", "RNFL_new", "GCIPL_new", "INL_new")


get_blups <- function(mod){
  ranef_vals <- ranef(mod)
  
  return(ranef_vals$ID)
}

for(var in vars){
  print("------------------------------")
  print(var)
  print("------------------------------")

  fm <- formula(paste(var,"~ pbl_DD + Age + Geschlecht +
                      (pbl_DD + 1|ID) + (1|ID_Eye)"))
  model <- lme4::lmer(fm,
                      data=data)
  eb_blups <- get_blups(model)
  
  dd_estimates <- eb_blups$pbl_DD + lme4::fixef(model)[2]
  
  ids <- row.names(get_blups(model))
  
  dcs <- data$Verlaufsform[match(row.names(get_blups(model)), data$ID)]
  
  disease_course <- as.factor(dcs)
  
  df_gc <- data.frame(dd_estimates = dd_estimates, disease_course = disease_course)[disease_course != 0,]

  #Group comparisons
  print(summary(aov(df_gc$dd_estimates ~ df_gc$disease_course)))
  print(TukeyHSD(aov(df_gc$dd_estimates ~ df_gc$disease_course)))
  
  
  df <- data.frame(b=dd_estimates, disease_course=disease_course)
  #Significance tests for coefficients, separated by disease course
  
  labels <- c("HC","RRMS","PPMS","SPMS")
  
  for(i in 1:4){
    v <- wilcox.test(df$b[df$disease_course == i-1])$p.value
    b <- mean(df$b[df$disease_course == i-1])
    print(paste(labels[i],": B =", b, ", p =",v))  
    

  }
  
  #Plotting
  value_range <- c(-max(abs(df$b)), max(abs(df$b)))
  binwidth <- max(abs(df$b)) / 10
  
  hc <- ggplot(df[df$disease_course == 0,])+
    geom_histogram(aes(b),binwidth = binwidth, color="grey", fill="#F8766D")+
    lims(x=value_range)+
    labs(x="BLUPs for disease duration", y="Count")+
    theme_classic()
  
  rrms <- ggplot(df[df$disease_course == 1,])+
    geom_histogram(aes(b),binwidth = binwidth, color="grey", fill="#7CAE00")+
    lims(x=value_range)+
    labs(x="BLUPs for disease duration", y="Count")+
    theme_classic()
  
  ppms <- ggplot(df[df$disease_course == 2,])+
    geom_histogram(aes(b),binwidth = binwidth, color="grey", fill="#00BFC4")+
    lims(x=value_range)+
    labs(x="BLUPs for disease duration", y="Count")+
    theme_classic()
  
  spms <- ggplot(df[df$disease_course == 3,])+
    geom_histogram(aes(b),binwidth = binwidth, color="grey", fill="#C77CFF")+
    lims(x=value_range)+
    labs(x="BLUPs for disease duration", y="Count")+
    theme_classic()
  
  library(patchwork)
  plt <- (hc + rrms) / (ppms + spms)
  pdf(paste("blup_figs/",var,"_blup.pdf",sep=""),width = 12,height = 6)
  print(plt)
  dev.off()
}


#------------------------------------Histograms------------------------------# 

ggplot(df[df$disease_course == 0,],aes(x=b)) +
  geom_histogram(bins=15)



