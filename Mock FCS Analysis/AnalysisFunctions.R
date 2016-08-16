#===========================================================
# The version of TMLE used in the analysis
#===========================================================
hc.tmle <- function(Y, # outcome
                    X, # data frame of observed data including column for trt
                    trt="trt", # name of treatment variable
                    fm, # regression model for prediction (SuperLearner or other 
                    # method as long as predict works)
                    X0, # data frame with trt = 0
                    X1, # data frame with trt = 1
                    onlySL=FALSE,
                    ... # additional arguments passed to predict
){
  # upper and lower bounds on Y
  a <- min(Y); b <- max(Y)
  
  # predicted values at observed treatment
  # predicted values at set trt
  if(class(fm)!="SuperLearner"){
    Qobs <- predict(fm, ...)
    Q1 <- predict(fm, newdata=X1,...)
    Q0 <- predict(fm, newdata=X0,...)
  }else{
    Qobs <- predict(fm,onlySL=onlySL)
    Q1 <- predict(fm, newdata=X1,onlySL=onlySL)
    Q0 <- predict(fm, newdata=X0,onlySL=onlySL)
  }
  X$H0 <- as.numeric(X[,trt]==0)/(1-mean(X[,trt]))
  X$H1 <- as.numeric(X[,trt]==1)/mean(X[,trt])
  X0$H0 <- 1/(1-mean(X[,trt])); X0$H1 <- 0
  X1$H0 <- 0; X1$H1 <- 1/mean(X$trt)
  
  getTMLE <- function(Qobs,Q1,Q0){
    Ystar <- (Y-a)/(b-a)
    X$logitQobs <- SuperLearner:::trimLogit((Qobs - a)/(b-a))
    X1$logitQobs <- SuperLearner:::trimLogit((Q1 - a)/(b-a))
    X0$logitQobs <- SuperLearner:::trimLogit((Q0 - a)/(b-a))
    
    suppressWarnings(
      fmFluc <- glm(Ystar ~ -1 + offset(logitQobs) + H0 + H1, data=X,
                    family="binomial")
    )
    QstarObs <- predict(fmFluc, type="response")*(b-a) + a
    Qstar1 <- predict(fmFluc, newdata=X1, type="response")*(b-a) + a
    Qstar0 <- predict(fmFluc, newdata=X0, type="response")*(b-a) + a
    
    psi1Star <- mean(Qstar1)
    psi0Star <- mean(Qstar0)
    D1Star <- X$H1 * (Y - Qstar1) + Qstar1 - psi1Star
    D0Star <- X$H0 * (Y - Qstar0) + Qstar0 - psi0Star
    n <- length(Y)
    v1 <- mean(D1Star^2); se1 <- sqrt(v1/n)
    v0 <- mean(D0Star^2); se0 <- sqrt(v0/n)
    vDiff <- mean((D1Star - D0Star)^2)
    seDiff <- sqrt(vDiff/n)
    
    out <- list(round(as.numeric(psi1Star + c(0,-1.96, 1.96)%o%se1),1),
                round(as.numeric(psi0Star + c(0,-1.96, 1.96)%o%se0),1),
                round(as.numeric(psi1Star - psi0Star + c(0,-1.96, 1.96)%o%seDiff),1),                
                round(2*pnorm(-abs((psi1Star - psi0Star)/seDiff)),5))
    names(out) <- c("trt","cntrl","diff","pval")
    out
  }
  
  if(class(fm)=="SuperLearner"){
    slTMLE <- getTMLE(Qobs = Qobs[[1]], Q1=Q1[[1]], Q0=Q0[[1]])
    # discrete Super Learner
    dslCol <- which(fm$cvRisk == min(fm$cvRisk))
    dslTMLE <- getTMLE(Qobs = Qobs[[2]][,dslCol], Q1 = Q1[[2]][,dslCol], Q0 = Q0[[2]][,dslCol])
    # for all candidates
    if(!onlySL){
      J <- ncol(Qobs[[2]])
      allCand <- lapply(split(1:J, 1:J), FUN=function(j){
        getTMLE(Qobs=Qobs[[2]][,j],Q1=Q1[[2]][,j],Q0=Q0[[2]][,j])
      })
      names(allCand) <- fm$libraryNames
    }else{
      allCand <- NULL
    }
    out <- list(slTMLE=slTMLE, dslTMLE=dslTMLE, allCandTMLE = allCand)
    class(out) <- "hc.tmle"
  }else{
    slTMLE <- dslTMLE <- NULL
    allCand <- getTMLE(Qobs=Qobs,Q1=Q1,Q0=Q0)
    out <- list(slTMLE=NULL, dslTMLE=NULL, allCandTMLE = allCand)
  }
  out
}

print.hc.tmle <- function(x,...){
  cat("Mean in trt = 1 arm (95% CI):  ", x$slTMLE$trt[1], " (", x$slTMLE$trt[2], ",", x$slTMLE$trt[3],") \n")
  cat("Mean in trt = 0 arm (95% CI):  ", x$slTMLE$cntrl[1], " (", x$slTMLE$cntrl[2], ",", x$slTMLE$cntrl[3],") \n")
  cat("Difference (95% CI):  ", x$slTMLE$diff[1], " (", x$slTMLE$diff[2], ",", x$slTMLE$diff[3],") \n")
  cat("p-value Difference = ", x$slTMLE$pval)
}


#=========================================================
# Screening functions for analysis
#=========================================================
demoScreen <- function (Y, X, family, obsWeights, id){ 
  return(colnames(X) %in% c("trt", "age", "female", "race"))
}

medScreen <- function(Y, X, family, obsWeights, id){
  return(colnames(X) %in% c("trt","sofa","score"))
}

demoScreen_ageInt <- function(Y, X, family, obsWeights, id){ 
  return(colnames(X) %in% c("trt", "age", "female", "race", "ageInt"))
}

demoScreen_femaleInt <- function(Y, X, family, obsWeights, id){ 
  return(colnames(X) %in% c("trt", "age", "female", "race", "femaleInt"))
}

demoScreen_raceInt <- function(Y, X, family, obsWeights, id){ 
  return(colnames(X) %in% c("trt", "age", "female", "race", "raceInt"))
}

medScreen_sofaInt <- function(Y, X, family, obsWeights, id){
  return(colnames(X) %in% c("trt","sofa","sofaInt","score"))
}
medScreen_charlInt <- function(Y, X, family, obsWeights, id){
  return(colnames(X) %in% c("trt","sofa","score","scoreInt"))
}

allScreen_noInt <- function(Y, X, family, obsWeights, id){
  intVars <- grep("Int",colnames(X))
  allVars <- rep.int(TRUE, ncol(X))
  allVars[intVars] <- FALSE
  return(allVars)
}

