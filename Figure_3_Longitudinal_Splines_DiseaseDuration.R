#------------Longitudinal spline fits depending on disease duration----------#
# Script to compute spline fits for retinal layer thickness
# depending on disease duration.
#Script to fit splines to data using linear mixed models with cross validation
#Run 10-fold CV in order to find the best degree / knot configurations.
# The results are presented in Figure 3.
# Usual runtime ~= 20 seconds

#---------------------------------R packages---------------------------------# 
library(lme4)
library(lmerTest)
library(haven)
library(ggplot2)

#-----------------------Data loading & data management-----------------------# 
data <- openxlsx::read.xlsx("data_04112023_JK.xlsx",detectDates = T)
data <- data[data$OCT_included == 1,]

#-- Define general functions

spline_fits_lmer <- function(x,y,sex,id,groups,knots=3,degree = 3){
  # Function to compute a single spline fit
  # Is used in next function, to iteratively optimize fit
  library(splines)
  #Build input data frame
  df <- data.frame(x,sex,id,groups,y)
  df <- df[complete.cases(df),]
  models <- list()
  xs <- c()
  ys <- c()
  grps <- c()
  #For all disease courses (1 - 4)
  for(i in 1:4){
    if((i-1) %in% unique(df$groups)){
      #Set knots at evenly spaced quantiles
      knots_val <- quantile(df$x[df$groups == i-1],(1:knots)/(knots+1))
      #Fit mixed model with splines
      models[[i]] <- lme4::lmer(y~bs(x,knots = knots_val,degree = degree)+sex+(1|id),data=df[df$groups == i-1,])
      preddf <- df[df$groups == i-1,]
      modmat <- model.matrix(models[[i]])
      #Compute predictions
      pred <- modmat[,1:(ncol(modmat)-1)] %*% lme4::fixef(models[[i]])[1:(ncol(modmat)-1)]
      xs <- c(xs,df$x[df$groups == i-1])
      ys <- c(ys,pred)
      grps <- c(grps,rep(i-1,length(pred)))
    }
  }
  return(list(frame=data.frame(x=xs,Group=grps,.fitted=ys),models=models))
}

spline_cv_fits_lmer <- function(x,y,sex,id,groups,cv = 10,knotrange = c(1:4),degreerange = c(2:3)){
  #Function to compute best spline fit, given covariate values and ranges to search for best configuration in

  df <- data.frame(x=x,y=y,sex=sex,id=id,groups=groups)
  df <- df[complete.cases(df),]
  xs <- c()
  preds <- c()
  grps <- c()
  for(group in unique(df$groups)){
    x_cv <- c()
    df_ <- df[df$groups == group,]
    newindices <- sample(1:nrow(df_),nrow(df_),replace=F)
    split_points <- ceiling(seq(1,nrow(df_),length.out=cv+1))
    #Iterate over all combinations
    bestknots = 1
    bestdegree = 1
    besterror = Inf
    bestpred <- c()
    for(knots in knotrange){
      for(degree in degreerange){
        #Cross-Validate and keep track of best fit
        mse <- 0
        prediction_cv <- c()
        for(i in 1:cv){
          # Generate cross-validation splits
          split <- split_points[i]
          split_upper <- split_points[i+1]
          df_test <- df_[newindices[split:split_upper],]
          df_train <- df_[-newindices[split:split_upper],]
          
          #Compute fit
          ret <- spline_fits_lmer(df_train$x,df_train$y,df_train$sex,df_train$id,df_train$groups,knots,degree)
          x <- df_test$x
          
          knots_val <- quantile(x,(1:knots)/(knots+1),na.rm=T)
          
          dftest <- as.matrix(data.frame(
            one <- rep(1,length(df_test$sex)),
            bs(x,knots=knots_val,degree=degree),
            sex = df_test$sex
          ))
          
          #Compute MSEs
          if(ncol(dftest) == length(lme4::fixef(ret$models[[group+1]]))){
            pred <- dftest %*% lme4::fixef(ret$models[[group+1]])
            prediction_cv <- c(prediction_cv,dftest[,1:(ncol(dftest-1))] %*% lme4::fixef(ret$models[[group+1]])[1:(ncol(dftest-1))])
            x_cv <- c(x_cv,x)
            mse <- mse + mean((df_test$y-pred)^2)
          }
          else{
            mse <- mse + mean((df_test$y - mean(df_test$y))^2)
          }
          
        }
        #Update current best spline fit
        if(mse < besterror){
          besterror <- mse
          bestknots <- knots
          bestdegree <- degree
          knots_val <- quantile(df_$x,(1:bestknots)/(bestknots+1),na.rm=T)
          # train model on
          model <- suppressWarnings(lme4::lmer(y ~ bs(df_$x,knots=knots_val,degree=degree) + sex + (1|id),data=df_))
          #Build input data for prediction
          dffin <- as.matrix(data.frame(
            one <- rep(1,length(df_$sex)),
            bs(df_$x,knots=knots_val,degree=degree)
          ))
          bestpred <- dffin %*% lme4::fixef(model)[1:ncol(dffin)]
        }
      }
    }
    
    preds <- c(preds,bestpred)
    
    xs <- c(xs,df_$x)
    
    grps <- c(grps,df_$groups)
  }
  # Return best fit
  return(data.frame(
    x = xs,
    Group = grps,
    .fitted = preds
  ))
}

show_plot <- function(outcome){
  if(outcome == "Latenz"){
    best_fit <- spline_cv_fits_lmer(data$Disease_Duration_Latency,
                                    data[,outcome],
                                    data$Geschlecht,
                                    data$ID,
                                    data$Verlaufsform)
  }else{
    best_fit <- spline_cv_fits_lmer(data$Disease_Duration_FINAL,
                                    data[,outcome],
                                    data$Geschlecht,
                                    data$ID,
                                    data$Verlaufsform)
  }

  
  
  #-- Plot the result
  plt <- ggplot() +
    geom_point(aes(data$Disease_Duration_FINAL,data[,outcome],shape=as.factor(data$Verlaufsform),color=as.factor(data$Verlaufsform)),size=0.5)+
    geom_line(aes(data$Disease_Duration_FINAL,data[,outcome],group=as.factor(data$ID_Eye),color=as.factor(data$Verlaufsform)),size=0.45)+
    geom_line(aes(best_fit$x,best_fit$.fitted, group=best_fit$Group, color=as.factor(best_fit$Group)), size=1) +
    labs(x = "Disease Duration in years", y=outcome, color="Disease Course", group="Disease Course", shape="Disease Course")+
    theme_classic()
 
  return(plt)
     
}

#Correct diesease durations, if VEP dates do not match OCT dates
library(lubridate)
data$Disease_Duration_Latency <- (as_date(data$VEP) - as_date(data$OCTDatum))/365 + data$Disease_Duration_FINAL


library(patchwork)
prnfl <- show_plot("globalpRNFL")
mrnfl <- show_plot("RNFL_new")
gcipl <- show_plot("GCIPL_new")
inl <- show_plot("INL_new")
va <- show_plot("VisusHC")
latency <- show_plot("Latenz")

pdf("plots.pdf",width=24,height=14.5)
(prnfl + mrnfl) / (gcipl + inl) / (va + latency)
dev.off()


