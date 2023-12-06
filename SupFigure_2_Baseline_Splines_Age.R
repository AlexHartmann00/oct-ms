#-----------------Baseline spline fits depending on age---------------------#
# Script to compute spline fits for baseline retinal layer thickness
# depending on age.
# The results are presented in Supplementary Figure 2.
# Usual runtime ~= 20 seconds


#---------------------------------R packages---------------------------------# 
library(ggplot2)
library(ggforce)


#-----------------------Data loading & data management-----------------------#
data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)

#Correct diesease durations, if VEP dates do not exactly match OCT dates
library(lubridate)
data$Age_Latency <- (as_date(data$VEP) - as_date(data$OCTDatum))/365 + data$Age

# Select only baseline observations (from first OCT per patient)
data_first <- data[data$Sequenz == 1 & !is.na(data$Age) & data$OCT_included == 1,]

#Function to generate jackknife bounds, given a spline fit
jackknife.bounds <- function(splinefit){
  #Get residuals and divide by leverage
  res <- (splinefit$yin - splinefit$y)/(1-splinefit$lev)  
  # Compute lower and upper bounds from residual standard deviation
  sigma <- sd(res)      
  upper <- splinefit$y + 1.96*sigma*sqrt(splinefit$lev)
  lower <- splinefit$y - 1.96*sigma*sqrt(splinefit$lev)
  
  return(cbind(x=splinefit$x,lower,upper))
}

spline.plot.baseline <- function(depvar,df){
  # Function to generate spline fit plot for
  # # a given dependent variable
  # # a given number of degrees of freedom
  if(depvar == "Latenz"){
    # Remove some latency outliers
    data <- data[data$Latenz < 200,]
    data$Age <- data$Age_Latency
  }
  # Initialize vectors
  ind <- c()
  ind_conf <- c()
  # Get only data with available measurements of the target layer thickness
  data_base <- data[!is.na(data[,depvar]) & !is.na(data$Age),]
  #Only keep first OCT per person
  data_base <- data_base[!duplicated(data_base$ID),]
  
  
  # RRMS
  dat_rrms <- data_base[data_base$Verlaufsform == 1,]
  spline_RRMS <- smooth.spline(dat_rrms$Age,dat_rrms[,depvar],df = df)
  rrms <- broom::augment(spline_RRMS)[,c("x",".fitted")]
  
  rrms_conf <- jackknife.bounds(spline_RRMS)
  
  ind_conf <- c(ind_conf,rep(1,nrow(rrms_conf)))
  ind <- c(ind, rep(1,nrow(rrms)))
  
  #PPMS
  dat_ppms <- data_base[data_base$Verlaufsform == 2,]
  spline_PPMS <- smooth.spline(dat_ppms$Age,dat_ppms[,depvar],df = df)
  ppms <- broom::augment(spline_PPMS)[,c("x",".fitted")]
  
  ppms_conf <- jackknife.bounds(spline_PPMS)
  
  ind_conf <- c(ind_conf,rep(2,nrow(ppms_conf)))
  ind <- c(ind, rep(2,nrow(ppms)))
  
  # SPMS
  dat_spms <- data_base[data_base$Verlaufsform == 3,]
  spline_SPMS <- smooth.spline(dat_spms$Age,dat_spms[,depvar],df = df)
  spms <- broom::augment(spline_SPMS)[,c("x",".fitted")]
  
  spms_conf <- jackknife.bounds(spline_SPMS)
  
  ind_conf <- c(ind_conf,rep(3,nrow(spms_conf)))
  ind <- c(ind, rep(3,nrow(spms)))
  plotdat <- cbind(as.data.frame(rbind(rrms,ppms,spms)),ind)
  
  #Combine results
  splinedat1 <- as.data.frame(rbind(rrms_conf,ppms_conf,spms_conf))
  
  splinedat <- cbind(splinedat1,ind_conf)
  
  #return Plot object
  ggplot()+
    geom_point(aes(data_base$Age,data_base[,depvar],fill=factor(data_base$Verlaufsform),color=factor(data_base$Verlaufsform),shape=factor(data_base$Verlaufsform)),size=0.5)+
    geom_line(aes(x,.fitted,color=factor(ind)),plotdat,size=0.45)+
    geom_ribbon(aes(x=x,y=rowMeans(splinedat[,c("lower","upper")]),ymin=lower,ymax=upper,fill=factor(ind_conf)),data.frame(splinedat),alpha=0.15)+
    scale_color_discrete(labels=c("HC","RRMS","PPMS","SPMS"))+
    scale_fill_discrete(labels=c("HC","RRMS","PPMS","SPMS"))+
    scale_shape_discrete(labels=c("HC","RRMS","PPMS","SPMS"))+
    scale_x_continuous(breaks=seq(15,85,5))+
    labs(fill="",color="",shape="",x="Age in years",y=depvar)+
    theme_classic()
}

library(patchwork)
prnfl <- spline.plot.baseline("globalpRNFL",3)
rnfl <- spline.plot.baseline("RNFL_new",3)
gcipl <- spline.plot.baseline("GCIPL_new",3)
inl <- spline.plot.baseline("INL_new",3)
va <- spline.plot.baseline("VisusHC",3)
latency <- spline.plot.baseline("Latenz",3)

pdf("plots.pdf",width=18,height=6)
(prnfl + rnfl + gcipl) / (inl + va + latency)

dev.off()
