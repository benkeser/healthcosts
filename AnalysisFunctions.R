#===========================================================
# The version of TMLE used in the analysis
#===========================================================
hc.tmle <- function(Y, # outcome used to fit fm
                    X, # data frame of trt + covariates used to fit fm
                    X0, # data frame with all relevant trt variables set to 1
                    X1, # data frame with all relevant trt variables set to 0
                    fm, # an object that predict can be called on, possibly a Super Learner object
                    trt, # name of the treatment variable
                    onlySL=FALSE, # return predictions for only Super Learner or all candidate methods too?
                    divider=1, # factor to divide results by
                    digits=1 # digits to round to after divider is used
){
  if(length(unique(X[,"trt"]))>2 | any(!(unique(X[,"trt"]) %in% c(0,1)))){
    stop("trt variable should only have two unique values - (0,1)")
  }
  
  if(all(X1==X0)){
    warning("X1 and X0 are exactly the same. Are you sure there is a trt variable?")
  }
  
  if(class(fm)!="SuperLearner"){
    stop("fm is not of class SuperLearner")
  }
  
  allPred <- predict(fm)
  allPred1 <- predict(fm, newdata=X1,onlySL=FALSE)
  allPred0 <- predict(fm, newdata=X0,onlySL=FALSE)
  
  out.SL <- get.tmle(Y=Y,X0=X0,X1=X1,X=X,divider=divider,digits=digits,trt=trt,
                     allPred=allPred$pred, allPred1=allPred1$pred,
                     allPred0=allPred0$pred)
  
  # get inference for other methods as well
  out.All <- NULL
  if(!onlySL){
    listMethods <- alply(1:ncol(allPred$library.predict), 1, function(x){
      cbind(allPred$library.predict[,x],allPred0$library.predict[,x],
            allPred1$library.predict[,x])
    })
    out.All <- lapply(listMethods, FUN=function(x){
      get.tmle(Y=Y,X0=X0,X1=X1,X=X,divider=divider,digits=digits,trt=trt,
               allPred=x[,1], allPred1=x[,3],
               allPred0=x[,2])
    })
    names(out.All) <- fm$libraryNames
  }
  
  # format output
  out <- list(SL=out.SL, All=invisible(out.All))
  class(out) <- "hc.tmle"
  out
}


get.tmle <- function(allPred,allPred1,allPred0,X,X0,X1,Y,trt,
                     divider, digits){
  X$H0 <- as.numeric(X[,trt]==0)/(1-mean(X[,trt]))
  X$H1 <- as.numeric(X[,trt]==1)/mean(X[,trt])
  X0$H0 <- 1/(1-mean(X[,trt])); X0$H1 <- 0
  X1$H0 <- 0; X1$H1 <- 1/mean(X[,trt])
  
  # predicted values
  X$QobsT <- allPred
  X1$QobsT <- allPred1
  X0$QobsT <- allPred0
  
  fmFluc <- glm(Y ~ -1 + offset(QobsT) + H0 + H1, data=X)
  Q1StarT <- predict(fmFluc, newdata=X1)
  Q0StarT <- predict(fmFluc, newdata=X0)
  psi1T <- mean(Q1StarT)
  psi0T <- mean(Q0StarT)
  D1StarT <- X$H1 * (Y - Q1StarT) + Q1StarT - mean(Q1StarT)
  D0StarT <- X$H0 * (Y - Q0StarT) + Q0StarT - mean(Q0StarT)
  DT <- matrix(rbind(D1StarT, D0StarT),nrow=2)
  covMatT <- DT%*%t(DT)/length(Y)
  a <- matrix(c(1, -1),nrow=2)
  psiTDiff <- psi1T - psi0T
  seTDiff <- sqrt(t(a)%*%covMatT%*%a / length(Y))
  
  #confidence intervals
  out <- list(round(as.numeric(psi1T + c(0,-1.96, 1.96)%o%sqrt(covMatT[1,1]/length(Y)))/divider,digits),
              round(as.numeric(psi0T + c(0,-1.96, 1.96)%o%sqrt(covMatT[2,2]/length(Y)))/divider,digits),
              round(as.numeric(psi1T - psi0T + c(0,-1.96, 1.96)%o%seTDiff)/divider,digits),
              round(2*pnorm(-abs(psiTDiff/seTDiff)),5))
  names(out) <- c("trt1","trt0","diff","pval")
  return(out)
}

print.hc.tmle <- function(x,...){
  print(x$SL)
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

