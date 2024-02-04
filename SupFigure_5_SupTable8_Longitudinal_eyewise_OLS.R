#---------------------------------Eye-wise OLS-------------------------------#
# Used to compute more robust estimates of change rates
# Usual runtime ~= 1 minute
# The results are presented in Supplementary Table 7 and Supplementary Figure 5
#---------------------------------R packages---------------------------------# 

library(lsmeans)
library(emmeans)
library(tidyverse)


#-----------------------Data loading & data management-----------------------# 

data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)

data <- data[data$OCT_included == 1,]


#-------------------------------Estimation model and test---------------------# 

#Select dependent variable to use 
# # (dependent_variable %in% )




for(dependent_variable in c("globalpRNFL", "RNFL_new", "GCIPL_new", "INL_new")){
  
  print(paste("------------",dependent_variable,"------------"))
  
  dat <- data[c("Disease_Duration_FINAL",dependent_variable,"Verlaufsform","ID_Eye","ID")]
  get_estimates <- function(id){
    #Function to get regression estimate of disease duration for a single eye
    sub <- dat[dat$ID_Eye == id,]
    sub$disease_duration <- sub$Disease_Duration_FINAL - min(sub$Disease_Duration_FINAL,na.rm=T)
    pid <- sub$ID[1]
    n <- nrow(sub)
    v <- sub$Verlaufsform[1]
    sub <- sub[complete.cases(sub),]
    if (nrow(sub) < 3){
      # To further reduce the variability of OLS estimates, we discarded eyes with less than three OCTs
      return(data.frame(person = pid,eye=id, b=NA, se=NA, disease_course=v, num_oct=n))
    }
    #Create formula
    formula_string <- paste(dependent_variable,"~ disease_duration")
    model <- lm(formula(formula_string),data=sub)
    c <- suppressWarnings(summary(model)$coefficients)
    #Extract coefficient and standard error for disease duration
    b <- c[2,1]
    se <- c[2,2]
    return(data.frame(person = pid,eye=id, b=b, se=se, disease_course=v, num_oct=n))
  }
  
  #Get all eye IDs to iterate over
  eyeids <- unique(dat$ID_Eye)
  
  resdfs <- list()
  
  for(eyeid in eyeids){
    r <- get_estimates(eyeid)
    resdfs[[eyeid]] <- r
  }
  
  df <- do.call(rbind.data.frame, resdfs)
  
  df <- df[complete.cases(df),]
  
  df.test <- df %>% group_by(disease_course) %>% 
    summarize(mean = mean(b,na.rm=T), t = t.test(b)$statistic, p = t.test(b)$p.value)
  
  print(df.test)
  #To test whether there are significant differences between disease courses, 
  #we performed an ANOVA with random intercepts per subject.
  #In cases where the ANOVA revealed significant group differences, 
  #we conducted a Tukey post-hoc test to identify which groups differed significantly.
  
  value_range <- c(-max(abs(df$b)), max(abs(df$b)))
  binwidth <- max(abs(df$b)) / 10
  
  #Plotting
  hc <- ggplot(df[df$disease_course == 0,])+
    geom_histogram(aes(b),binwidth = binwidth, color="grey", fill="#F8766D")+
    lims(x=value_range)+
    labs(x="OLS coefficient for disease duration", y="Count")+
    theme_classic()
  
  rrms <- ggplot(df[df$disease_course == 1,])+
    geom_histogram(aes(b),binwidth = binwidth, color="grey", fill="#7CAE00")+
    lims(x=value_range)+
    labs(x="OLS coefficient for disease duration", y="Count")+
    theme_classic()
  
  ppms <- ggplot(df[df$disease_course == 2,])+
    geom_histogram(aes(b),binwidth = binwidth, color="grey", fill="#00BFC4")+
    lims(x=value_range)+
    labs(x="OLS coefficient for disease duration", y="Count")+
    theme_classic()
  
  spms <- ggplot(df[df$disease_course == 3,])+
    geom_histogram(aes(b),binwidth = binwidth, color="grey", fill="#C77CFF")+
    lims(x=value_range)+
    labs(x="OLS coefficient for disease duration", y="Count")+
    theme_classic()
  
  library(patchwork)
  plt <- (hc + rrms) / (ppms + spms)
  
  # Without HC
  model <- lmerTest::lmer(b ~ as.factor(disease_course) + (1|person), data = df[df$disease_course != 0,])
  print(anova(model,type=3,ddf="Kenward-Roger"))
  print(lsmeans::lsmeans(model, pairwise~as.factor(disease_course), adjust="tukey"))
  

}

