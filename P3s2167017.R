###### Practical 3: Maeve LI(Minqing LI)(s2167017) ######
## Overview: This R code creates 4 functions, among which linmod fits a linear model, print.linmod prints out 
## a simple summary of the fitted model containing the model formula and the parameter estimates and their standard deviations,
## plot.linmod plots the model residuals against the model fitted values, and predict.linmod gives prediction values of 
## the response variable based on the fitted model when given a new set of values of predictor variables. 


linmod <- function(formula,dat){
  ## A function that fits a linear model using the given formula (y~X) and data using the QR decomposition method
  ## and outputs an object of class "linmod", a list that contains: 
  ## "beta", the least squares parameter estimates of the model; 
  ## "V", the estimated covariance matrix of the least squares estimators; "mu", the fitted values;
  ## "y", the vector containing the response variable; "yname", the name of the response variable;
  ## "formula", the model formula; "flev", a named list containing the factor variables 
  ## used in the model and their respective levels; "sigma", the estimated standard deviation of the response.
  
  #### Data Preparation
  ## y
  yname = all.vars(formula)[1] ## get y name
  y <- model.frame(formula,dat)[[1]] ## take out the y values
  names(y) = 1:length(y) ## labeling y
  ## X
  X = model.matrix(formula,dat)
  ## getting dat's factor variable names
  flev = list()  ## initialize an empty list to store the factor names and levels
  j=0
  ifac <- c() ## initialize an empty vector to store factors' indices
  for (i in 1:ncol(dat)){  ## a loop is then created to get each factor's levels
    if (is.factor(dat[,i])){
      j = j + 1
      flev[[j]] <- unique(dat[,i]) ## store the levels in the jth element of the list
      ifac[j] <- i ## store the index of this particular factor
    }
  }
  if (j) 
    ## if dat contains at least one factor variable (flev is not empty)
    names(flev)=names(dat[ifac]) ## name the list's elements with the factor variables' names
  if (all.vars(formula)[2] != "."){ 
    ## if the formula is not a regression on all variables in dat, 
    ## we need to just keep the variables used in the formula
    xnames <- all.vars(formula)[-1] ## select x's names in the formula
    truefalse <- names(flev) %in% xnames ## check if factor variables in dat is in the formula
    usedfac <- names(flev)[truefalse]
    flev <- flev[usedfac] ## Keep the factors used in the formula in list flev
  }
  
  #### Fit linear model using QR decomposition method
  qrx <- qr(X)
  p <- ncol(X)
  beta <- backsolve(qr.R(qrx),qr.qty(qrx,y)[1:p])
  names(beta)=colnames(X) ## name the elements of beta the same as the corresponding model components
  
  #### Several computations
  ## fitted values
  mu = X %*% beta
  mu = drop(mu)  ## drop mu's dimension and turn it into a vector
  names(mu) = 1:length(mu) ## labeling mu
  ## standard deviation of the response (residuals)
  sigma = sqrt(sum((y-mu)^2)/(nrow(X)-p))
  ## covariance matrix of the least square estimators
  V <- backsolve(qr.R(qrx),forwardsolve(t(qr.R(qrx)),diag(nrow(qr.R(qrx)))*(sigma^2)))
  
  #### Return a list
  l = list(beta=beta,V=V,mu=mu,y=y,yname=yname,formula=formula,flev=flev,sigma=sigma)
  class(l)<-"linmod" ## define a class of the object to be returned as "linmod"
  return(l)
}

print.linmod <- function(x){
  ## A print method function that gives the model formula defining the model,
  ## and report the parameter estimates and their standard deviations in a specified format
  ## for x, which is an object of class linmod
  print(x$formula)
  cat('\n')
  b_st_error <- sqrt(diag(x$V)) ## calculate each estimate's standard error
  M = cbind(x$beta,b_st_error) ## combine two vectors to form the output matrix
  colnames(M) = c("Estimate","s.e.") ## name the output's columns
  print(M)
}

plot.linmod <- function(x){
  ## A plot method function that plots a Residuals vs Fitted values scatter plot 
  ## for x, which is an object of class linmod
  resids <- x$y - x$mu  ## compute residuals
  plot(x$mu,resids,main="Residuals vs Fitted",xlab="fitted values",ylab="residuals")
  abline(h=0, col = "dimgray", lwd=1, lty=2) ## form a horizontal dashed line passing through 0
}

predict.linmod <- function(x,newdata){
  ## A predict method funtion for an object of class linmod, which uses the built model fit to predict a new
  ## vector of values when given a new dataframe "newdata". It returns a vector of prediction values,
  ## each of which is a corresponding prediction for each row of newdata. It deals with the situations where 
  ## factor variables are provided as character variables, and where factors with the wrong levels are provided.
  
  ## When fewer levels are provided, the function will still work. But when more levels (i.e. new levels) are provided,
  ## this function will output error, for it's essentially impossible to predict with new levels 
  ## that didn't exist in the original model fit.
  
  ## Any variable that is a factor variable in the original model fit but isn't one in newdata
  ## is transformed into a factor variable, and any factor variable with the wrong number of levels 
  ## is ensured that it has the right number of levels using this if function and two loops inside of it
  if(length(x$flev) != 0){  ## for fitted models that contain at least 1 factor variable
    for (i in 1:length(x$flev)){
      a = names(x$flev)[i]
      if(class(newdata[,a]) != "factor")  ## if the factor variable in the model fit but isn't one in newdata
        newdata[,a] <- factor(newdata[,a]) ## transform it into a factor variable
    }
    for (i in 1:length(x$flev)){
      a = names(x$flev)[i]
      if (length(levels((x$flev)[[i]])) > length(levels(newdata[,a]))){ ## if fewer levels are provided in newdata
        newdata[,a] <- factor(as.character(newdata[,a]),levels=levels((x$flev)[[i]])) ## ensure it has enough levels as the one in the model fit
      }
    }
  }
  newdata$dummy <- rep(0,nrow(newdata))  ## creates a dummy variable in new data
  names(newdata)[length(newdata)] <- x$yname ## name the dummy variable by "yname"
  newdat <- model.matrix(x$formula,newdata)
  prediction = newdat %*% (x$beta)
  return(prediction)
}