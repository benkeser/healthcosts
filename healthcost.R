# semiparametric transformation method of Welsh and Zhou (2006) and Zhou et al (2009)**
# discrete conditional density estimation **

# finite mixture models
# SL.fmm <- function(Y, X, newX, family, obsWeights, 
#                    k=2, concomitantFormula= as.formula("~ 1"),
#                    iter=500,
#                    ...){
#   if(family$family=="gaussian"){
#     require(flexmix)
#     
#     fit.flexmix <- flexmix(as.formula(paste0("Y~",paste(colnames(X),collapse="+"))), data=data.frame(Y, X), k=k, 
#                            model=FLXMRglm(family="Gamma"),
#                            control=list(iter=iter),
#                            concomitant=FLXPmultinom(concomitantFormula))
#     
#     
#     fit <- list(object = list(fit.flexmix, k))
#     class(fit) <- "SL.fmm"
#     
#     cHat <- clusters(fit.flexmix, newdata=newX)
#     muList <- predict(fit.flexmix, newdata=newX)
#     
#     pred <- rep(NA, times=nrow(newX))
#     for(i in 1:k){ pred[cHat==i] <- muList[[i]] }
#   
#   }else{
#     stop("SL.fmm not written for binomial family")
#   }
#   
#   out <- list(fit=fit, pred=pred)
#   out
# }
# 
# predict.SL.fmm <- function(object, newdata, ...){
#   cHat <- clusters(object$object[[1]])
#   muList <- predict(object$object[[1]], newdata=newdata)
#   pred <- rep(NA, times=nrow(newdata))
#   for(i in 1:object$object[[2]]){ pred[cHat==i] <- muList[[i]][cHat==i] }
#   pred
# }
#   
#   

# gilleskie and mroz 
SL.gilleskie <- function(Y, X, newX, family, obsWeights,
                         kValues=c(5,15,25), # number of intervals
                         yBoundaries, # boundaries on the y-variable
                         maxPoly=2,
                         ...){
  library(sandwich)
  
  #==========================
  # loop this over K values
  #==========================
  if(any(kValues < maxPoly)){ 
    warning("kValue specified that is less than maxPoly. These kValues will be ignored")
    kValues <- kValues[kValues>maxPoly]
  }
  outList <- lapply(split(kValues, 1:length(kValues)), FUN=function(K){
    Ytilde <- cut(Y, breaks=quantile(Y, p=seq(0,1,length=K+1)), labels=FALSE,
                  include.lowest=TRUE)
    longDat <- data.frame(Ytilde,X,id=1:length(Y))[rep(1:length(Y),Ytilde),]
    row.names(longDat)[row.names(longDat) %in% paste(row.names(data.frame(Ytilde,X)))] <- paste(row.names(data.frame(Ytilde,X)),".0",sep="")  
    longDat$k <- as.numeric(paste(unlist(strsplit(row.names(longDat),".",fixed=T))[seq(2,nrow(longDat)*2,2)]))+1
    longDat$indY <- as.numeric(longDat$k==longDat$Ytilde)
  
    pVal <- Inf
    d <- maxPoly
    while(pVal > 0.05 & d>=1){
    rhs <- NULL
      for(i in 1:(ncol(X)-1)){
        rhs <- c(rhs, paste0("poly(",colnames(X)[i],",",ifelse(length(unique(X[,i]))>d, d, length(unique(X[,i]))-1),")*poly(k,",d,")*",colnames(X)[(i+1):(ncol(X))],collapse="+"))
      }
      rhs <- c(rhs, paste0("poly(",colnames(X)[ncol(X)],",",ifelse(length(unique(X[,i]))>d, d, length(unique(X[,i]))-1),")*poly(k,",d,")"))
      
      suppressWarnings(
      fm <- glm(as.formula(paste0("indY ~ ",paste0(rhs,collapse="+"))),
                data=longDat, family="binomial")
      )
      dropNum <- NULL
      for(cn in colnames(X)){
        dropNum <- c(dropNum, grep(paste0(cn,", ",d,")",d), names(fm$coef[!is.na(fm$coef)])))
      }
      dropCoef <- fm$coef[!is.na(fm$coef)][dropNum]
      fullCov <- vcovHC(fm,type="HC0")
      dropCov <- fullCov[dropNum,dropNum]
      chi2Stat <- tryCatch({
        t(dropCoef)%*%solve(dropCov)%*%dropCoef
      },error=function(e){ return(0) }
      )
      pVal <- pchisq(chi2Stat, lower.tail=FALSE, df=length(dropCoef))
      d <- d-1
    }
    suppressWarnings(
      longDat$haz <- predict(fm, newdata=longDat,type="response")
    )
    
    # calculate likelihood
    tmp <- by(longDat, factor(longDat$id), FUN=function(x){
      prod((c(1,cumprod(1-x$haz[x$k < x$Ytilde])) * x$haz)^x$indY)
    })
    LKR <- sum(log(as.numeric(tmp))) + length(Y)*log(K)
    
    return(list(LKR=LKR, fm=fm))
  })
  
  LKRs <- unlist(lapply(outList, function(x){x[[1]]}))
  maxLKR <- which(LKRs==max(LKRs))
  maxK <- kValues[maxLKR]
  thisOut <- outList[[maxLKR]]
  
  # calculate mean
  # get mean in each bin
  Ytilde <- cut(Y, breaks=quantile(Y, p=seq(0,1,length=maxK+1)), labels=FALSE,
                include.lowest=TRUE)
  hK <- apply(matrix(1:maxK), 1, function(x){
    mean(Y[Ytilde==x])
  })
  
  # calculate mean
  pred <- apply(matrix(1:length(newX[,1])), 1, FUN=function(x){
    suppressWarnings(
      haz <- predict(thisOut$fm, 
                     newdata=data.frame(newX[x,],k=1:maxK),
                     type="response")
    )
    dens <- c(1,cumprod(1-haz[1:(maxK-1)])) * haz
    sum(hK*dens)
  })
  
  fit <- list(object=thisOut$fm, maxK=maxK, hK=hK)
  class(fit) <- "SL.gilleskie"
  out <- list(pred=pred, fit=fit)
  out
}

predict.SL.gilleskie <- function(object, newdata,...){
  pred <- apply(matrix(1:length(newdata[,1])), 1, FUN=function(x){
    suppressWarnings(
      haz <- predict(object$object,newdata=data.frame(newdata[rep(x,object$maxK),],k=1:object$maxK),
                     type="response")
    )
    dens <- c(1,cumprod(1-haz[1:(object$maxK-1)])) * haz
    sum(object$hK*dens)
  })
  pred
}

# wang and zhou quantile method
SL.wangZhou <- function(Y, X, newX, family, obsWeights, 
                        g="log",m=length(Y), c=0.2, b=0.05,
                        ...){
  require(quantreg)
  if(family$family=="gaussian"){
    n <- length(Y)
    # calculate alpha_n
    alpha <- c*n^(-1/(1+4*b))
    tau <- seq(alpha, 1-alpha, length=m)
    
    # transform Y
    if(g=="log"){
      thisY <- log(Y)
      ginv <- function(x){ exp(x) }
    }else{
      stop("SL.wangZhou only implemented for log transforms")
    }
    
    # get quantile regressions
    suppressWarnings(
    fm <- rq(formula=as.formula("thisY~."), tau=tau, weights=obsWeights,
             data=data.frame(thisY,X))
    )
    QhatList <- predict(fm, newdata=newX, stepfun=TRUE, type="Qhat")
    QhatRearrange <- lapply(QhatList, rearrange)
    pred <- unlist(lapply(QhatRearrange, FUN=function(Q){
      Qw <- ginv(environment(Q)$y[-which(duplicated(environment(Q)$x))])
      1/(1-2*alpha) * sum(Qw * diff(c(0,tau)))
    }))
  }else{
    stop("SL.wangZhou not written for binomial family")
  }
  fit <- list(object=fm, alpha=alpha, ginv=ginv)
  class(fit) <- "SL.wangZhou"
  out <- list(pred=pred, fit=fit)
  return(out)
}

predict.SL.wangZhou <- function(object, newdata, ...){
  require(quantreg)
  QhatList <- predict(object$object, newdata=newdata, stepfun=TRUE, type="Qhat")
  QhatRearrange <- lapply(QhatList, rearrange)
  pred <- mapply(Q=QhatRearrange, dt=diff(c(0,object$object$tau)), function(Q,dt){
    Qw <- do.call(object$ginv,args=list(x=environment(Q)$y[-which(duplicated(environment(Q)$x))]))
    1/(1-2*object$alpha) * sum(Qw * dt)
  })    
  pred
}

SL.gammaLogGLM <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    fit.glm <- glm(Y ~ ., data=X, family=Gamma(link='log'), weights=obsWeights,
                   control=list(maxit=100))
    pred <- predict(fit.glm, newdata=newX, type="response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.logGLM not written for binomial family")
  }
}

SL.gammaIdentityGLM <- function(Y, X, newX, family, obsWeights,...){
  if(family$family=="gaussian"){
    fit.glm <- glm(Y ~ ., data=X, family=Gamma(link='identity'), 
                   weights=obsWeights,
                   control=list(maxit=100), start=c(mean(Y),rep(0,ncol(X))))
    pred <- predict(fit.glm, newdata=newX, type="response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.gammaIdentityGLM not written for binomial family")
  }
}


# log transformed glm
SL.gaussianLogGLM <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    myFamily <- gaussian(link="log")
    startValues <- c(log(mean(Y)),rep(0,ncol(X)))
    fit.glm <- glm(Y ~ ., data=X, family=myFamily, weights=obsWeights,start=startValues,
                   control=list(maxit=100))
    pred <- predict(fit.glm, newdata=newX, type="response")
    
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.logGLM not written for binomial family")
  }
}

# log transformed glm + Duan (1983) correction
SL.logOLS.smear <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
  logY <- log(Y)
  fit.logGLM <- glm(logY ~ ., data=X, family=family, weights=obsWeights)
  
  mu <- predict(fit.logGLM, type="response", newdata=X)
  resid <- logY - mu
  
  pred <- exp(predict(fit.logGLM, type="response",newdata=newX))*mean(exp(resid))
  fit <- list(object=fit.logGLM, mean(exp(resid)))
  class(fit) <- "SL.logOLS.smear"
  }else{
    stop("SL.logGLM.smear not written for binomial family")
  }
  
  out <- list(fit=fit, pred=pred)
  return(out)
}

predict.SL.logOLS.smear <- function(object, newdata, ...){
  mu <- predict(object$object, newdata=newdata, type="response")
  correction <- object[[2]]
  
  return(exp(mu)*correction) 
}

# manning 2001 algorithm
SL.manningGLM <- function(Y, X, newX, family, obsWeights, 
                          kCut = 3, lambdaCut = c(0.5,1.5,2.5), 
                          startNLS=0,
                          ...){
  if(family$family=="gaussian"){
    require(moments)
    # ols on log scale
    logY <- log(Y)
    fit.logGLM <- glm(logY ~ ., data=X, family=family, weights=obsWeights)
    
    mu <- predict(fit.logGLM, type="response", newdata=X)
    resid <- logY - mu
    
    # check kurtosis of residuals
    k <- kurtosis(resid)
    
    pred <- exp(predict(fit.logGLM, type="response",newdata=newX))*mean(exp(resid))
    fit <- list(object=fit.logGLM, mean(exp(resid)))
    class(fit) <- "SL.logOLS.smear"
    try({
    if(k < kCut){
      # park test
      fit.initGLM <- glm(Y ~ ., data=X, weights=obsWeights, family="gaussian")
      muPark <- predict(fit.initGLM, type="response", newdata=X)
      resid2Park <- (Y - muPark)^2
      fit.parkGLM <- glm(resid2Park ~ muPark, family="gaussian")
      lambda1 <- fit.parkGLM$coef[2]
      if(lambda1 < lambdaCut[1]){
        xNames <- colnames(X)
        d <- length(xNames)
        bNames <- paste0("b",1:d)
        form <- apply(matrix(1:d), 1, function(i){
          paste(c(bNames[i],xNames[i]),collapse="*")
        })
        formula <- paste(form,collapse=" + ")
        try({
        fit.nls <- nls(as.formula(paste0("Y ~ exp(b0 +",formula,")")), data= data.frame(Y, X),
                       start=eval(parse(text=paste0(
                         "list(b0=0.5,",paste(paste0(bNames, "=", startNLS),collapse=","),")"))))
        })
        pred <- predict(fit.nls, newdata=newX)
        fit <- list(object=fit.nls)
        class(fit) <- "SL.manningGLM"
      }else if(lambda1 < lambdaCut[2] & lambda1 >= lambdaCut[1]){
        fit.poisGLM <- suppressWarnings(
          glm(Y ~ ., data=X, weights=obsWeights, family="poisson",control=list(maxit=100))
        )
        pred <- predict(fit.poisGLM, newdata=newX, type="response")
        fit <- list(object=fit.poisGLM)
        class(fit) <- "SL.manningGLM"
      }else if(lambda1 < lambdaCut[3] & lambda1 >= lambdaCut[2]){
        fit.gammaGLM <- glm(Y ~ ., data=X, weights=obsWeights, family=Gamma(link='log'),control=list(maxit=100))
        pred <- predict(fit.gammaGLM, newdata=newX, type="response")
        fit <- list(object=fit.gammaGLM)
        class(fit) <- "SL.manningGLM"
      }else if(lambda1 > lambdaCut[3]){
        fit.invGaussianGLM <- glm(Y ~ ., data=X,weights=obsWeights, family=inverse.gaussian(link="log"),control=list(maxit=100))
        pred <- predict(fit.invGaussianGLM, newdata=newX, type="response")
        fit <- list(object=fit.invGaussianGLM)
        class(fit) <- "SL.manningGLM"
      }
    }
  }, silent=TRUE)
  }else{
    stop("SL.manningGLM doesn't work with binomial family.")
  }
  out <- list(pred = pred, fit=fit)
  return(out)
}

predict.SL.manningGLM <- function(object, newdata,...){
  if(!is.list(object$object)){
    cat("I'm trying")
    pred <- predict(object=object$object, newdata=newdata, type="response")
  }else{
    pred <- predict(object=object$object, newdata=newdata, type="response")
  }
  pred
}

# generalized gamma regression (or other 'survival' regression models) *
SL.flexsurvreg <- function(Y, X, newX, family, obsWeights,
                           dist="gengamma",...){
  require(flexsurv)
  
  if(family$family=="gaussian"){
    fit.flexSurv <- flexsurvreg(
      as.formula(paste0("Surv(Y, rep(1, length(Y))) ~", paste(colnames(X),collapse="+"))) ,
      data=X, dist=dist
    )
    
    pred <- predict.SL.flexsurvreg(object=list(object=fit.flexSurv), newdata=newX, type="mean")
    
    fit <- list(object=fit.flexSurv)
    class(fit) <- "SL.flexsurvreg"
    out <- list(fit=fit, pred=pred)
  }else{
    stop("SL.genGamma not implemented for binominal family")
  }
  
  out
}

SL.gengamma <- function(..., dist="gengamma"){
  SL.flexsurvreg(...,dist=dist)
}

SL.weibull <- function(..., dist="weibull"){
  SL.flexsurvreg(...,dist=dist)
}

SL.lognormalsurv <- function(..., dist="lognormal"){
  SL.flexsurvreg(...,dist=dist)
}

predict.SL.flexsurvreg <- function(object, newdata, type="mean", ...){
  .getSurv <- function(x, fit, thisnewdata){
    summary(fit, t=x, B=0, newdata=thisnewdata)[[1]][,2]
  }
  pred <- as.numeric(apply(matrix(1:nrow(newdata)), 1, function(i){
    upper <- Inf
    out <- NA; class(out) <- "try-error"
    while(class(out)=="try-error" & upper > 1e6){
      out <- try(integrate(.getSurv, fit=object$object, thisnewdata=newdata[i,],lower=0,upper=upper)$value, silent=TRUE)
      if(upper==Inf){
        upper <- 1e8
      }else{
        upper <- upper/2
      }
    }
    if(class(out)=="try-error"){
      warning("Unable to integrate survival function. Returning random number between 0 and 100k")
      out <- runif(1,0,100000)
    }
    out
  }))
  pred
}
SL.coxph  <- function(Y, X, newX, family, obsWeights,
                      dist="gengamma",...){
  if(family$family=="gaussian"){
    library(survival)
    fit.coxph <- coxph(Surv(Y,rep(1,length(Y)))~., data=X)
    fit <- list(object=fit.coxph)
    class(fit) <- "SL.coxph"
    pred <- predict.SL.coxph(object=list(object=fit.coxph), newdata=newX)
  }else{
    stop("SL.coxph not implemented for binominal family")
  }
  return(list(fit=fit,pred=pred))
}
predict.SL.coxph <- function(object,newdata,type="mean",...){
  surv.fit <- survfit(object$object, newdata=newdata)
  pred <- colSums(diff(c(0,surv.fit$time))*rbind(rep(1,dim(surv.fit$surv)[2]),surv.fit$surv[1:(dim(surv.fit$surv)[1]-1),]))
  pred
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

# other SL functions
SL.caret1 <- function (Y, X, newX, family, obsWeights, method = "rf", tuneLength = 3, 
                       trControl = trainControl(method = "cv", number = 20, verboseIter = FALSE), 
                       metric,...) 
{
  if (length(unique(Y))>2){
    if(is.matrix(Y)) Y <- as.numeric(Y)
    metric <- "RMSE"
    if(method=="gbm"){
      suppressWarnings(
    fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                              metric = metric, method = method, 
                              tuneLength = tuneLength, 
                              trControl = trControl,verbose=FALSE)
      )
    }else{
      suppressWarnings(
      fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                metric = metric, method = method, 
                                tuneLength = tuneLength, 
                                trControl = trControl)
      )
    }
    pred <- predict(fit.train, newdata = newX, type = "raw")
  }
  if (length(unique(Y))<=2) {
    metric <- "Accuracy"
    Y.f <- as.factor(Y)
    levels(Y.f) <- c("A0", "A1")
    fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                              metric = metric, method = method, 
                              tuneLength = tuneLength, 
                              trControl = trControl)
    pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
  }
  fit <- list(object = fit.train)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.caret")
  return(out)
}

SL.rpart.caret1 <- function(...,method="rpart",tuneLength = 8){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

SL.rf.caret1 <- function(...,method="rf",tuneLength=8){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

SL.gbm.caret1 <- function(...,method="gbm",tuneLength=8){
  SL.caret1(...,method=method,tuneLength=tuneLength)
}

my.tmle <- function(Y,X,fm,X0,X1){
  X$H0 <- as.numeric(X$trt==0)/(1-mean(X$trt))
  X$H1 <- as.numeric(X$trt==1)/mean(X$trt)
  X0$H0 <- 1/(1-mean(X$trt)); X0$H1 <- 0
  X1$H0 <- 0; X1$H1 <- 1/mean(X$trt)
  
  # predicted values
  X$QobsT <- predict(fm)[[1]]
  Q1T <- predict(fm, newdata=X1,onlySL=TRUE)$pred
  Q0T <- predict(fm, newdata=X0,onlySL=TRUE)$pred
  X1$QobsT <- Q1T
  X0$QobsT <- Q0T
  
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
  out <- list(round(as.numeric(psi1T + c(0,-1.96, 1.96)%o%sqrt(covMatT[1,1]/length(dat$totalcost)))/1000,1),
              round(as.numeric(psi0T + c(0,-1.96, 1.96)%o%sqrt(covMatT[2,2]/length(dat$totalcost)))/1000,1),
              round(as.numeric(psi1T - psi0T + c(0,-1.96, 1.96)%o%seTDiff)/1000,1),
              round(2*pnorm(-abs(psiTDiff/seTDiff)),5))
  
  names(out) <- c("trt","cntrl","diff","pval")
  out
  # direct costs
  # QobsD <- predict(fmD)[[1]]
  # Q1D <- predict(fmD, newdata=X1,onlySL=TRUE)$pred
  # Q0D <- predict(fmD, newdata=X0,onlySL=TRUE)$pred
  # X1$QobsD <- Q1D
  # X0$QobsD <- Q0D
  # 
  # fmFlucD <- glm(dat$directcost ~ -1 + offset(QobsD) + H0 + H1, data=X)
  # Q1StarD <- predict(fmFlucD, newdata=X1)
  # Q0StarD <- predict(fmFlucD, newdata=X0)
  # psi1D <- mean(Q1StarD)
  # psi0D <- mean(Q0StarD)
  # D1StarD <- X$H1 * (dat$directcost - Q1StarD) + Q1StarD - mean(Q1StarD)
  # D0StarD <- X$H0 * (dat$directcost - Q0StarD) + Q0StarD - mean(Q0StarD)
  # DD <- matrix(rbind(D1StarD, D0StarD),nrow=2)
  # covMatD <- DD%*%t(DD)/length(dat$directcost)
  # a <- matrix(c(1, -1),nrow=2)
  # psiDDiff <- psi1D - psi0D
  # seDDiff <- sqrt(t(a)%*%covMatD%*%a / length(dat$directcost))
  # 
  # #confidence intervals
  # round(as.numeric(psi1D + c(0,-1.96, 1.96)%o%sqrt(covMatD[1,1]/length(dat$totalcost)))/1000,1)
  # round(as.numeric(psi0D + c(0,-1.96, 1.96)%o%sqrt(covMatD[2,2]/length(dat$totalcost)))/1000,1)
  # round(as.numeric(psi1D - psi0D + c(0,-1.96, 1.96)%o%seDDiff)/1000,1)
  # round(pnorm(psiDDiff/seDDiff),2)
  
}


genLabels <- function(sl){
  out <- lapply(sl, function(x){
    ind <- FALSE
    if(length(grep("gbm",x[1]))>0 &
       length(grep("caret",x[1]))==0){
      meth <- "GBM"
      ind <- TRUE
    }else if(length(grep("SL.gbm.caret1",x[1]))>0){
      meth <- "CV-GBM"
      ind <- TRUE
    }else if(length(grep("rf",x[1]))>0){
      meth <- "CV-Random Forest"
      ind <- TRUE
    }else if(length(grep("randomForest",x[1]))>0){
      meth <- "Random Forest"
      ind <- TRUE
    }else if(length(grep("SL.rpart.",x[1]))>0){
      meth <- "CV-Regression Tree"
      ind <- TRUE
    }else if(length(grep("SL.rpart",x[1]))>0 &
             length(grep("caret",x[1]))==0){
      meth <- "Regression Tree"
      ind <- TRUE
    }else if(length(grep("SL.rpart",x[1]))>0 &
             length(grep("caret",x[1]))>0){
      meth <- "CV-Regression Tree"
      ind <- TRUE
    }else if(length(grep("coxph",x[1]))>0){
      meth <- "Cox PH"
    }else if(length(grep("gilleskie",x[1]))>0){
      meth <- "Gilleskie (2004)"
      ind <- TRUE
    }else if(length(grep("wangZhou",x[1]))>0){
      meth <- "Wang (2009)"
    }else if(length(grep("lognormalsurv",x[1]))>0){
      meth <- "AFT-LogNormal"
    }else if(length(grep("weibull",x[1]))>0){
      meth <- "AFT-Weibull"
    }else if(length(grep("gengamma",x[1]))>0){
      meth <- "AFT-Gen.Gamma"
    }else if(length(grep("manning",x[1]))>0){
      meth <- "Manning (2001)"
    }else if(length(grep("logOLS.smear",x[1]))>0){
      meth <- "GLM-Smearing"
    }else if(length(grep("gammaId",x[1]))>0){
      meth <- "GLM-Gamma/Identity"
    }else if(length(grep("gammaLog",x[1]))>0){
      meth <- "GLM-Gamma/Log"
    }else if(length(grep("SL.glm",x[1]))>0){
      meth <- "GLM-Normal/Identity"
    }
    if(length(grep("All",x[2]))>0){
      cov <-"Demographics + Medical - all"
    }else if(length(grep("demoScreen",x[2]))>0 &
             length(grep("_",x[2]))==0){
      if(!ind)  cov <- "Demographics - none"
      if(ind) cov <-"Demographics - all"
    }else if(length(grep("allScreen_noInt",x[2]))>0){
      if(!ind) cov <-"Demographics + Medical - none"
      if(ind) cov <-"Demographics + Medical - all"
    }else if(length(grep("demoScreen_ageInt",x[2]))>0){
      cov <-"Demographics - age"
    }else if(length(grep("demoScreen_raceInt",x[2]))>0){
      cov <-"Demographics - race"
    }else if(length(grep("demoScreen_femaleInt",x[2]))>0){
      cov <-"Demographics - gender"
    }else if(length(grep("medScreen",x[2]))>0 &
             length(grep("_",x[2]))==0){
      if(!ind) cov <-"Medical - none"
      if(ind) cov <-"Medical - all"
    }else if(length(grep("medScreen_sofaInt",x[2]))>0){
      cov <-"Medical - sofa"
    }else if(length(grep("medScreen_charlInt",x[2]))>0){
      cov <-"Medical - charlton"
    }
    return(c(meth, cov))
  })
}

plotInferenceResults <- function(allOut, total=TRUE, xl, SLrslt, SLlibrary, t1,
                                 lmar=9.1){
  lab <- genLabels(sl=SLlibrary)
  
  par(mar=c(4.1,lmar,0.1,0.3))
  plot(0,0,pch="", ylim=c(0,length(SLlibrary)-3),xaxt="n",xlab="ATE (95% CI)",ylab="",
       xlim=xl,yaxt="n",bty="n")
  axis(side=1)
  abline(v=0,lty=3)
  ct <- length(SLlibrary)
  ct2 <- 0
  for(i in 1:length(allOut)){
    ct2 <- ct2+1
    ct <- ct-1
    #   if(x$diff[1]>0) cat(ct2, " " ,paste0(mySLlibrary[[ct2]],collapse="_"),"\n")
    points(y=ct,x=allOut[[i]]$diff[1],col=1)
    segments(y0=ct,x0=allOut[[i]]$diff[2],x1=allOut[[i]]$diff[3], col=1)
    text(x=par()$usr[1], y=ct, cex=0.4, pos=2, lab[[i]][2], xpd=TRUE)
    text(x=par()$usr[1]-t1, y=ct, cex=0.4, pos=4, lab[[i]][1],xpd=TRUE)
  }
  points(y=-2,x=SLrslt[1],pch=18,cex=1.2, col=1)
  segments(y0=-2,x0=SLrslt[2], x1=SLrslt[3])
  text(x=par()$usr[1]-t1,pos=4,y=-2,xpd=TRUE, "Super Learner",cex=0.4)
  
  text(x=par()$usr[1]-t1,pos=4,y=length(SLlibrary),xpd=TRUE, expression(underline(paste("Method"))),cex=0.4)
  text(x=par()$usr[1],pos=2,y=length(SLlibrary),xpd=TRUE, expression(underline(paste("Covariates-interactions"))),cex=0.4)
}

