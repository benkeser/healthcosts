#==================================================#
# Adaptive Hazard method of Gilleskie and Mroz 
#==================================================#
SL.gilleskie <- function(Y, X, newX, family, obsWeights,
                         kValues=c(5,15,25), # number of intervals
                         yBoundaries, # boundaries on the y-variable
                         maxPoly=2, # maximum polynomial in hazard regressions
                         ...){
  
  # need the sandwich package for covariate selection algorithm
  library(sandwich)
  
  # maxPoly describes the polynomial used for the partition variable
  # in the hazard regression. Choosing the number of partitions to be
  # less than this value leads to rank definiciency in glm()
  if(any(kValues < maxPoly)){ 
    warning("kValue specified that is less than maxPoly. These kValues will be ignored")
    kValues <- kValues[kValues>maxPoly]
  }
  
  #====================================================
  # get hazard fit over different partitions of data
  #====================================================
  outList <- lapply(split(kValues, 1:length(kValues)), FUN=function(K){
    # break up Y into K+1 partitions
    Ytilde <- cut(Y, breaks=quantile(Y, p=seq(0,1,length=K+1)), labels=FALSE,
                  include.lowest=TRUE)
    # make a long versions data set
    longDat <- data.frame(Ytilde,X,id=1:length(Y))[rep(1:length(Y),Ytilde),]
    # assign parition number variable
    row.names(longDat)[row.names(longDat) %in% paste(row.names(data.frame(Ytilde,X)))] <- paste(row.names(data.frame(Ytilde,X)),".0",sep="")  
    longDat$k <- as.numeric(paste(unlist(strsplit(row.names(longDat),".",fixed=T))[seq(2,nrow(longDat)*2,2)]))+1
    # indicator of falling in a particular partition
    longDat$indY <- as.numeric(longDat$k==longDat$Ytilde)
  
    # loop to do covariate selection
    pVal <- Inf
    d <- maxPoly
    while(pVal > 0.05 & d>=1){
      # generate the regression equation
      rhs <- NULL
      for(i in 1:(ncol(X)-1)){
        rhs <- c(rhs, paste0("poly(",colnames(X)[i],",",ifelse(length(unique(X[,i]))>d, d, length(unique(X[,i]))-1),")*poly(k,",d,")*",colnames(X)[(i+1):(ncol(X))],collapse="+"))
      }
      rhs <- c(rhs, paste0("poly(",colnames(X)[ncol(X)],",",ifelse(length(unique(X[,i]))>d, d, length(unique(X[,i]))-1),")*poly(k,",d,")"))
      
      # fit the hazard regression
      suppressWarnings(
      fm <- glm(as.formula(paste0("indY ~ ",paste0(rhs,collapse="+"))),
                data=longDat, family="binomial")
      )
      # get coefficients of degree d
      dropNum <- NULL
      for(cn in colnames(X)){
        dropNum <- c(dropNum, grep(paste0(cn,", ",d,")",d), names(fm$coef[!is.na(fm$coef)])))
      }
      dropCoef <- fm$coef[!is.na(fm$coef)][dropNum]
      # get covariance matrix for all thos ecoefficients
      fullCov <- vcovHC(fm,type="HC0")
      dropCov <- fullCov[dropNum,dropNum]
      # test significance of those coefficients
      chi2Stat <- tryCatch({
        t(dropCoef)%*%solve(dropCov)%*%dropCoef
      },error=function(e){ return(0) }
      )
      pVal <- pchisq(chi2Stat, lower.tail=FALSE, df=length(dropCoef))
      d <- d-1
    }
    # after done dropping polynomial terms, get hazard predictions
    suppressWarnings(
      longDat$haz <- predict(fm, newdata=longDat,type="response")
    )
    
    # calculate likelihood
    tmp <- by(longDat, factor(longDat$id), FUN=function(x){
      prod((c(1,cumprod(1-x$haz[x$k < x$Ytilde])) * x$haz)^x$indY)
    })
    LKR <- sum(log(as.numeric(tmp))) + length(Y)*log(K)
    # return the likelihood ratio and estimate hazard regression
    return(list(LKR=LKR, fm=fm))
  })
  
  # figure out which one had highest likelihood ratio
  LKRs <- unlist(lapply(outList, function(x){x[[1]]}))
  maxLKR <- which(LKRs==max(LKRs))
  maxK <- kValues[maxLKR]
  thisOut <- outList[[maxLKR]]
  
  # get mean in each partition for transforming back to mean-scale
  Ytilde <- cut(Y, breaks=quantile(Y, p=seq(0,1,length=maxK+1)), labels=FALSE,
                include.lowest=TRUE)
  hK <- apply(matrix(1:maxK), 1, function(x){
    mean(Y[Ytilde==x])
  })
  
  # calculate mean by calculating density of each partition 
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

# predict function for SL.gilleskie
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

#=======================================================
# Quantile regression method of Wang and Zhou (2009)
#=======================================================
SL.wangZhou <- function(Y, X, newX, family, obsWeights, 
                        g="log", # transformation of Y
                        m=length(Y), # number of quantiles
                        c=0.2, # for calculating truncated mean
                        b=0.05,# for calculating truncated mean
                        ...){
  require(quantreg)
  if(family$family=="gaussian"){
    n <- length(Y)
    # calculate alpha_n for calculating truncated mean
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
    # transform to means
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

# predict function for SL.wangZhou
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

#======================================================
# GLM with Gamma family and log link
#======================================================
SL.gammaLogGLM <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    fit.glm <- glm(Y ~ ., data=X, family=Gamma(link='log'), weights=obsWeights,
                   control=list(maxit=100))
    pred <- predict(fit.glm, newdata=newX, type="response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm" # can use predict.SL.glm
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.logGLM not written for binomial family")
  }
}
#======================================================
# GLM with Gamma family and identity link
#======================================================
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

#===========================================================
# GLM with Gaussian family and log link -- not very stable
#===========================================================
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

#======================================================
# GLM with Gaussian family and id link on log(Y) 
# + Duan (1983) correction
#======================================================
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
# predict function for SL.logOLS.smear
predict.SL.logOLS.smear <- function(object, newdata, ...){
  mu <- predict(object$object, newdata=newdata, type="response")
  correction <- object[[2]]
  
  return(exp(mu)*correction) 
}

#======================================================
# Adaptive GLM algorithm of Manning (2001)
#======================================================
SL.manningGLM <- function(Y, X, newX, family, obsWeights, 
                          kCut = 3, # kurtosis cutpoint
                          lambdaCut = c(0.5,1.5,2.5), # skew cutpoint
                          startNLS=0, # starting values for NLS?
                          ...){
  if(family$family=="gaussian"){
    require(moments)
    # first do ols on log scale
    logY <- log(Y)
    fit.logGLM <- glm(logY ~ ., data=X, family=family, weights=obsWeights)
    
    mu <- predict(fit.logGLM, type="response", newdata=X)
    resid <- logY - mu
    
    # check kurtosis of residuals
    k <- kurtosis(resid)
    
    # by default use these methods
    # some of the other GLMs are unstable and if they fail, this 
    # algorithm returns log OLS + smearning estimate
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
      # use nls
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
        # use poisson glm
        fit.poisGLM <- suppressWarnings(
          glm(Y ~ ., data=X, weights=obsWeights, family="poisson",control=list(maxit=100))
        )
        pred <- predict(fit.poisGLM, newdata=newX, type="response")
        fit <- list(object=fit.poisGLM)
        class(fit) <- "SL.manningGLM"
      }else if(lambda1 < lambdaCut[3] & lambda1 >= lambdaCut[2]){
        # use gamma glm
        fit.gammaGLM <- glm(Y ~ ., data=X, weights=obsWeights, family=Gamma(link='log'),control=list(maxit=100))
        pred <- predict(fit.gammaGLM, newdata=newX, type="response")
        fit <- list(object=fit.gammaGLM)
        class(fit) <- "SL.manningGLM"
      }else if(lambda1 > lambdaCut[3]){
        # use inverse gaussian glm -- not very stable
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
# predict function
predict.SL.manningGLM <- function(object, newdata,...){
  if(!is.list(object$object)){
    pred <- predict(object=object$object, newdata=newdata, type="response")
  }else{
    pred <- predict(object=object$object, newdata=newdata, type="response")
  }
  pred
}

#======================================================
# Accelerated Failure Time Models
#======================================================
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
  # function to return survival probability based on flexsurv object
  .getSurv <- function(x, fit, thisnewdata){
    summary(fit, t=x, B=0, newdata=thisnewdata)[[1]][,2]
  }
  pred <- as.numeric(apply(matrix(1:nrow(newdata)), 1, function(i){
    upper <- Inf
    out <- NA; class(out) <- "try-error"
    # integrate can be finnicky, so for stability, we first try to integrate with
    # upper limit = Inf, but if that fails move to 1e8, which sometimes is able to 
    # provide a sane answer when upper limit=Inf fails. Keep trying smaller and smaller
    # values, but don't go smaller than 1e6. If you try, then it just returns a random 
    # number between 0 and 1e6, which prevents Super Learner from crashing. 
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


#===========================================================
# Cox Proportional Hazard for estimating conditional mean
#===========================================================
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
  # use surv.fit to get survival estimate and because by default it uses
  # nelson-aalen hazard, easy to convert back to an estimate of the mean
  surv.fit <- survfit(object$object, newdata=newdata)
  pred <- colSums(
    diff(c(0,surv.fit$time))*rbind(
      rep(1,dim(surv.fit$surv)[2]),
      surv.fit$surv[1:(dim(surv.fit$surv)[1]-1),]
    )
  )
  pred
}

#===================================================================
# Modified version of SL.caret that prints less annoying GBM output
#===================================================================
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


#====================================================================
# The Highly Adaptive Lasso (HAL) estimator and associated functions
#====================================================================
SL.hal <- function(Y, X, newX, family=gaussian(), 
                   verbose=FALSE,
                   obsWeights=rep(1,length(Y)), 
                   sparseMat=TRUE,
                   nfolds = ifelse(length(Y)<=100, 20, 10), 
                   nlambda = 100, useMin = TRUE,...){
  require(Matrix)
  d <- ncol(X)
  if(!sparseMat){
    uniList <- alply(as.matrix(X),2,function(x){
      # myX <- matrix(x,ncol=length(unique(x)), nrow=length(x)) - 
      #   matrix(unique(x), ncol=length(unique(x)), nrow=length(x), byrow=TRUE)
      myX <- matrix(x,ncol=length(x), nrow=length(x)) -
        matrix(x, ncol=length(x), nrow=length(x), byrow=TRUE)
      myX <- ifelse(myX < 0, 0, 1)
      myX
    })
    
    if(d >=2){
      highDList <- alply(matrix(2:d),1,function(k){
        thisList <- alply(combn(d,k),2,function(x){
          Reduce("*",uniList[x])
        })
        Reduce("cbind",thisList)
      })
      initX <- cbind(Reduce("cbind",uniList), Reduce("cbind",highDList))
      dup <- duplicated(t(initX))
      designX <- initX[,!dup]
    }else{
      initX <- Reduce("cbind",uniList)
      dup <- duplicated(t(initX))
      designX <- initX[,!dup]
    }
    
    fitCV <- glmnet::cv.glmnet(x = designX, y = Y, weights = obsWeights, 
                               lambda.min.ratio=0.001,
                               lambda = NULL, type.measure = "deviance", nfolds = nfolds, 
                               family = family$family, alpha = 1, nlambda = nlambda)
    
    
    ## get predictions back
    mynewX <- matrix(newX[,1],ncol=length(X[,1]), nrow=length(newX[,1])) - 
      matrix(X[,1], ncol=length(X[,1]), nrow=length(newX[,1]), byrow=TRUE)
    mynewX <- ifelse(mynewX < 0, 0, 1)
    
    makeNewDesignX <- TRUE
    if(all(dim(X)==dim(newX)))
      makeNewDesignX <- !all(X==newX)
    
    if(makeNewDesignX){
      uniList <- alply(matrix(1:ncol(X)),1,function(x){
        myX <- matrix(newX[,x],ncol=length(X[,x]), nrow=length(newX[,x])) - 
          matrix(X[,x], ncol=length(X[,x]), nrow=length(newX[,x]), byrow=TRUE)
        myX <- ifelse(myX < 0, 0, 1)
        myX
      })
      
      if(d >=2){
        highDList <- alply(matrix(2:d),1,function(k){
          thisList <- alply(combn(d,k),2,function(x){
            Reduce("*",uniList[x])
          })
          Reduce("cbind",thisList)
        })
        
        initX <- cbind(Reduce("cbind",uniList), Reduce("cbind",highDList))
        designNewX <- initX[,!dup]
      }else{
        initX <- Reduce("cbind",uniList)
        designNewX <- initX[,!dup]
      }
    }else{
      designNewX <- designX
    }
    
    pred <- predict(fitCV$glmnet.fit, newx = designNewX, 
                    s = ifelse(useMin,fitCV$lambda.min, fitCV$lambda.1se), type = "response")
    fit <- list(object = fitCV, useMin = useMin, X=X, dup=dup, sparseMat=sparseMat)
  }else{
    SuperLearner:::.SL.require("plyr")
    SuperLearner:::.SL.require("data.table")
    SuperLearner:::.SL.require("stringr")
    SuperLearner:::.SL.require("bit")
    
    if(is.vector(X)) X <- matrix(X, ncol=1)
    if(is.vector(newX)) newX <- matrix(newX, ncol=1)
    n <- length(X[,1])
    d <- ncol(X)
    
    if(verbose) cat("Making sparse matrix \n")
    X.init <- makeSparseMat(X=X,newX=X,verbose=TRUE)
    
    ## find duplicated columns
    if(verbose) cat("Finding duplicate columns \n")
    # Number of columns will become the new number of observations in the data.table
    nIndCols <- ncol(X.init)
    # Pre-allocate a data.table with one column, each row will store a single column from X.init
    datDT <- data.table(ID = 1:nIndCols, bit_to_int_to_str = rep.int("0", nIndCols))
    # Each column in X.init will be represented by a unique vector of integers.
    # Each indicator column in X.init will be converted to a row of integers or a string of cat'ed integers in data.table
    # The number of integers needed to represent a single column is determined automatically by package "bit" and it depends on nrow(X.init)
    nbits <- nrow(X.init) # number of bits (0/1) used by each column in X.init
    bitvals <- bit(length = nbits) # initial allocation (all 0/FALSE)
    nints_used <- length(unclass(bitvals)) # number of integers needed to represent each column
    # For loop over columns of X.init
    ID_withNA <- NULL # track which results gave NA in one of the integers
    for (i in 1:nIndCols) {
      bitvals <- bit(length = nbits) # initial allocation (all 0/FALSE)
      Fidx_base0 <- (X.init@p[i]) : (X.init@p[i + 1]-1) # zero-base indices of indices of non-zero rows for column i=1
      nonzero_rows <- X.init@i[Fidx_base0 + 1] + 1 # actual row numbers of non-zero elements in column i=1
      # print(i); print(nonzero_rows)
      # X.init@i[i:X.init@p[i]]+1 # row numbers of non-zero elements in first col
      bitvals[nonzero_rows] <- TRUE
      # str(bitwhich(bitvals))
      intval <- unclass(bitvals) # integer representation of the bit vector
      # stringval <- str_c(intval, collapse = "")
      if (any(is.na(intval))) ID_withNA <- c(ID_withNA, i)
      set(datDT, i, 2L, value = str_c(str_replace_na(intval), collapse = ""))
    }
    # create a hash-key on the string representation of the column, 
    # sorts it by bit_to_int_to_str using radix sort:
    setkey(datDT, bit_to_int_to_str)
    # add logical column indicating duplicates, 
    # following the first non-duplicate element
    datDT[, duplicates := duplicated(datDT)]
    # just get the column IDs and duplicate indicators:
    datDT[, .(ID, duplicates)]
    
    dupInds <- datDT[,ID][which(datDT[,duplicates])]
    uniqDup <- unique(datDT[duplicates==TRUE,bit_to_int_to_str])
    
    colDups <- alply(uniqDup, 1, function(l){
      datDT[,ID][which(datDT[,bit_to_int_to_str] == l)]
    })
    
    if(verbose) cat("Fitting lasso \n")
    if(length(dupInds)>0){
      notDupInds <- (1:ncol(X.init))[-unlist(colDups, use.names = FALSE)]
      keepDupInds <- unlist(lapply(colDups, function(x){ x[[1]] }), use.names=FALSE)
      
      fitCV <- glmnet::cv.glmnet(x = X.init[,c(keepDupInds,notDupInds)], y = Y, weights = obsWeights, 
                                 lambda = NULL, lambda.min.ratio=0.001, type.measure = "deviance", nfolds = nfolds, 
                                 family = family$family, alpha = 1, nlambda = nlambda)
    }else{
      fitCV <- glmnet::cv.glmnet(x = X.init, y = Y, weights = obsWeights, 
                                 lambda = NULL, lambda.min.ratio=0.001, type.measure = "deviance", nfolds = nfolds, 
                                 family = family$family, alpha = 1, nlambda = nlambda)
    }
    
    fit <- list(object = fitCV, useMin = useMin, X=X, dupInds = dupInds, colDups=colDups, sparseMat=sparseMat)
    class(fit) <- "SL.hal"
    
    if(identical(X,newX)){
      if(length(dupInds) > 0){
        pred <- predict(fitCV, newx = X.init[,c(keepDupInds,notDupInds)], s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), 
                        type = "response")
      }else{
        pred <- predict(fitCV, newx = X.init, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), 
                        type = "response")
      }
    }else{
      pred <- predict(fit, newdata=newX, bigDesign=FALSE, chunks=10000)          
    }
  }
  
  out <- list(pred = pred, fit = fit)
  cat("Done with SL.hal")
  return(out)
}



predict.SL.hal <- function (object, newdata, bigDesign=FALSE, verbose=TRUE, 
                            chunks=10000,
                            s = ifelse(object$useMin, object$object$lambda.min, object$object$lambda.1se),...)
{
  if(!object$sparseMat){
    d <- ncol(object$X)
    # if you want to get predictions all at once (smaller newdata)
    if(bigDesign){
      uniList <- alply(matrix(1:ncol(object$X)),1,function(x){
        myX <- matrix(newdata[,x],ncol=length(object$X[,x]), nrow=length(newdata[,x])) - 
          matrix(object$X[,x], ncol=length(object$X[,x]), nrow=length(newdata[,x]), byrow=TRUE)
        myX <- ifelse(myX < 0, 0, 1)
        myX
      })
      
      if(d >=2){
        highDList <- alply(matrix(2:d),1,function(k){
          thisList <- alply(combn(d,k),2,function(x){
            Reduce("*",uniList[x])
          })
          Reduce("cbind",thisList)
        })
        
        initX <- cbind(Reduce("cbind",uniList), Reduce("cbind",highDList))
        designNewX <- initX[,!object$dup]
      }else{
        initX <- Reduce("cbind",uniList)
        designNewX <- initX[,!object$dup]
      }
      # get predictions
      pred <- predict(object$object$glmnet.fit, newx = designNewX, 
                      s = s, 
                      type = "response")
      
    }else{
      # get row by row predictions, so you never have to store a big design matrix
      # for newdata
      pred <- apply(as.matrix(newdata),1,function(i){
        uniList <- alply(matrix(1:ncol(object$X)),1,function(x){
          myX <- matrix(i[x],ncol=length(object$X[,x]), nrow=length(i[x])) - 
            matrix(object$X[,x], ncol=length(object$X[,x]), nrow=length(i[x]), byrow=TRUE)
          myX <- ifelse(myX < 0, 0, 1)
          myX
        })
        
        if(d >=2){
          highDList <- alply(matrix(2:d),1,function(k){
            thisList <- alply(combn(d,k),2,function(x){
              Reduce("*",uniList[x])
            })
            Reduce("cbind",thisList)
          })
          
          initX <- cbind(Reduce("cbind",uniList), Reduce("cbind",highDList))
          designNewX <- initX[!object$dup]
        }else{
          initX <- Reduce("cbind",uniList)
          designNewX <- initX[!object$dup]
        }
        # get predictions
        thispred <- predict(object$object$glmnet.fit, newx = matrix(designNewX,nrow=1), s=s,
                            type = "response")
        thispred
      })
    }
  }else{
    # all predictions at once
    if(bigDesign){
      pred <- doPred(object=object,newdata=newdata,verbose=verbose)
    }else{
      nNew <- length(newdata[,1])
      nChunks <- floor(nNew/chunks) + ifelse(nNew%%chunks==0, 0, 1)
      pred <- rep(NA, length(nNew))
      for(i in 1:nChunks){
        minC <- (i-1)*chunks + 1
        maxC <- ifelse(i==nChunks, nNew,i*chunks)
        pred[minC:maxC] <- doPred(object=object,s=s,newdata=newdata[minC:maxC,],verbose=verbose)
      }
    }
  }
  return(as.numeric(pred))
}

doPred <- function(object,newdata,verbose=TRUE,s){
  if(is.vector(newdata)) newdata <- matrix(newdata)
  
  if(verbose) cat("Making initial sparse matrix for predictions \n")
  tmp <- makeSparseMat(X=object$X, newX=newdata, verbose=verbose)
  
  if(length(object$dupInds) > 0){
    if(verbose) cat("Correcting duplicate columns in sparse matrix \n")
    # get vector of duplicate columns
    dupVec <- unlist(object$colDups,use.names=FALSE)
    # number of each duplicate
    nperDup <- unlist(lapply(object$colDups,length),use.names=FALSE)
    # number of duplicates to roll through
    K <- length(nperDup)
    # start and ending index
    startInd <- c(0, cumsum(nperDup)[1:(K-1)])
    endInd <- cumsum(nperDup)
    # not duplicate colums
    notdupVec <- (1:ncol(tmp))[-dupVec]
    # put all the duplicated guys first
    tmp <- tmp[,c(dupVec,notdupVec)]
    
    uniqRowsList <- list()
    myp <- c(0,rep(NA, K))
    # look at the i associatiated with 
    for(k in 1:K){
      # this condition ensures that not all the values of a given set of duplciates
      # are equal to zero.
      if(tmp@p[startInd[k]+1] != tmp@p[endInd[k]+1]){
        Fidx_base0 <- (tmp@p[startInd[k] + 1]) : (tmp@p[endInd[k] + 1] - 1)
        nonzero_rows <- tmp@i[Fidx_base0 + 1] + 1 # actual row numbers of non-zero elements in column i=1
        # unique nonzero_rows
        theseRows <- sort(unique(nonzero_rows))
        uniqRowsList[[k]] <- theseRows
        # a new p vector for duplicated columns
        myp[k+1] <- myp[k] + length(theseRows)
      }else{ # the whole block for this set of duplicates is 0
        uniqRowsList[[k]] <- NULL
        myp[k+1] <- myp[k]
      }
    }
    
    # look at this sparse matrix
    myi <- unlist(uniqRowsList)
    # check if it came out right
    # grbg1 <- sparseMatrix(i=myi, p=myp, x=1)
    
    # fix up p with nondup columns
    ## for this example every non-duplicated column in the new design
    ## matrix is 0, which is causing this to break. I think. 
    if(tmp@p[endInd[K] + 1] != tmp@p[length(tmp@p)]){
      fulli <- c(myi, tmp@i[(tmp@p[endInd[K] + 1] + 1) : tmp@p[length(tmp@p)]] + 1)
      fullp <- c(myp, 
                 tmp@p[((endInd[K]+1) + 1) : length(tmp@p)] - 
                   tmp@p[(endInd[K]+1)] + myp[K+1])
    }else{
      fulli <- myi
      fullp <- myp
    }
    # 
    tmp <- sparseMatrix(i=fulli, p=fullp, x=1, 
                        dims=c(length(newdata[,1]),
                               length(notdupVec) + length(object$colDup)))
  }
  pred <- predict(object$object$glmnet.fit, newx=tmp, 
                  s = s)
  pred
}

makeSparseMat <- function(X,newX=X, verbose=FALSE){
  if(is.vector(X)) X <- matrix(X, ncol=1)
  if(is.vector(newX)) newX <- matrix(newX, ncol=1)
  
  d <- ncol(X)
  stopifnot(ncol(X)==ncol(newX))
  
  nX <- length(X[,1])
  n <- length(newX[,1])
  
  # numbers used to correct column indices later
  colStart <- 1
  colEnd <- d 
  
  # start by creating a list of univariate indicators
  # length of the list is d and the entries are matrices
  # of row and column indices for a design matrix based 
  # only on that covariate, i.e. columns in each list entry 
  # run from 1:n, so we can use intersect later for higher
  # order terms. 
  if(verbose) cat("Making ", d," basis fns of dimension  1 \n")
  uni <- alply(matrix(1:d),1,function(x){
    j <- alply(matrix(newX[,x]),1,function(y){ 
      which(X[,x] <= y) 
    })
    i <- rep(1:n,unlist(lapply(j, length), use.names=FALSE))
    cbind(unlist(i, use.names=FALSE),unlist(j, use.names=FALSE))
  })
  # number of 1's for each variable -- for variables with
  # length(unique(x)) == length(x) will be n*(n+1)/2, but if there
  # are ties, the length will be different
  nperuni <- lapply(uni,nrow)
  
  # slap them all together
  uni.ij <- Reduce("rbind",uni)
  
  # fix up the column indices
  uni.ij[,2] <- uni.ij[,2] + 
    rep.int((colStart:colEnd)-1, times=unlist(nperuni, use.names=FALSE))*nX
  
  # i = row indices, j = column indices
  i <- uni.ij[,1]
  j <- uni.ij[,2]
  
  # functions used for higher order terms
  .myIntersect <- function(...){
    Reduce(intersect,list(...))
  }
  .getIntersect <- function(...){
    tmp <- lapply(..., function(b){
      split(b[,2],b[,1])
    })
    tmpNames <- lapply(tmp,function(l){
      as.numeric(names(l))
    })
    overlap <- Reduce(intersect,tmpNames)
    
    # indices of tmp that overlap
    newtmp <- lapply(tmp,function(b){
      b[paste(overlap)]
    })
    
    # get intersection
    out <- eval(parse(text=paste0(paste0(
      "mapply(.myIntersect,"),paste0("newtmp[[",1:length(tmp),"]]",collapse=","),",SIMPLIFY=FALSE)")
    ))
    out
  }
  
  # loop over higher order terms
  if(d > 1){
    for(k in 2:d){
      # matrix of all d choose k combinations
      combos <- combn(d,k)
      
      if(verbose) cat("Making ", ncol(combos), " basis fns of dimension ", k, "\n" )
      # adjust column indicators for column indices
      colStart <- colEnd + 1
      colEnd <- (colStart-1) + ncol(combos)
      
      # list of length d choose k, each entry 
      # containing n indices of columns corresponding to subjects
      j.list <- alply(combos, 2, function(a){
        .getIntersect(uni[a])
      })
      
      # list of length d choose k, each entry containing
      # n indices of rows corresponding to subjects
      i.list <- llply(j.list, function(x){
        rep(as.numeric(names(x)), unlist(lapply(x, length), use.names=FALSE))
      })
      
      # number of 1's for each combination
      nper <- lapply(i.list, length)
      
      # unlist column numbers
      j.list <- lapply(j.list, unlist, use.names=FALSE)
      
      # unlist rows and columns
      thisi <- unlist(i.list,use.names=FALSE)
      thisj <- unlist(j.list,use.names=FALSE)
      
      # fix up the column number 
      thisj <- thisj + 
        rep.int((colStart:colEnd)-1, times=unlist(nper, use.names=FALSE))*nX
      
      # put it together
      i <- c(i, thisi)
      j <- c(j, thisj)
    }
  }
  
  # make the sparseMatrix
  grbg <- sparseMatrix(i=i[order(i)],j=j[order(i)],x=1, dims=c(n, nX*(2^d - 1)))
  return(grbg)
}

#============================================================
# list the wrappers included in GitHub repository
#============================================================
hc.listWrappers <- function(what="both"){
  if (what == "both") {
    message("All prediction algorithm wrappers in healthcosts:\n")
    message("GLMs:\n")
    print(c("SL.gammaIdentityGLM","SL.gammaLogGLM", "SL.gaussianLogGLM",
            "SL.logOLS.smear", "SL.manningGLM"))
    message("Survival-like methods:\n")
    print(c("SL.coxph","SL.flexsurvreg", "SL.lognormalsurv", "SL.gengamma","SL.weibull",
            "SL.gilleskie"))
    message("Quantile-based methods:\n")
    print(c("SL.wangZhou"))
    message("  Nonparametric methods:\n")
    print(c("SL.caret1",  "SL.gbm.caret1","SL.rf.caret1", "SL.rpart.caret1"))

    message("\nAll screening algorithm wrappers in healthcosts:\n")
    print("No additional screening wrappers")
  }
  else if (what == "SL") {
    message("All prediction algorithm wrappers in healthcosts:\n")
    message("   GLMs:\n")
    print(c("SL.gammaIdentityGLM","SL.gammaLogGLM", "SL.gaussianLogGLM",
            "SL.logOLS.smear", "SL.manningGLM"))
    message("   Survival-like methods:\n")
    print(c("SL.coxph","SL.flexsurvreg", "SL.lognormalsurv", "SL.gengamma","SL.weibull",
            "SL.gilleskie"))
    message("   Quantile-based methods:\n")
    print(c("SL.wangZhou"))
    message("  Nonparametric methods:\n")
    print(c("SL.caret1",  "SL.gbm.caret1","SL.rf.caret1", "SL.rpart.caret1"))
  }
  else if (what == "screen") {
    message("\nAll screening algorithm wrappers in healthcosts:\n")
    print("No additional screening wrappers")
  }
  else if (what == "method") {
    message("All methods in SuperLearner package:\n")
    print("No additional methods")
  }
  else {
    stop("Please specify what = 'both', 'SL', or 'screen'")
  }
}

#============================================================
# make method object to use with Super Learner for 
# bounded log-likelihood loss
#============================================================
makeBoundedMethod <- function(upperBound,lowerBound,
                              name="method.CC_nloglik.Bounded"){
  eval(parse(text=paste0(name, "<<- list(
                         require='nloptr',
                         computeCoef= function (Z, Y, libraryNames, obsWeights, control, verbose, 
                         lowerBound=",lowerBound,", upperBound=",upperBound,", ...){
                         transZ <- (Z - lowerBound) / (upperBound - lowerBound)
                         transY <- (Y - lowerBound) / (upperBound - lowerBound)
                         logitZ <- trimLogit(transZ, control$trimLogit)
                         cvRisk <- apply(logitZ, 2, function(x){
                         -mean(obsWeights * (transY * plogis(x, log.p = TRUE) + 
                         (1-transY) * plogis(x, log.p = TRUE,lower.tail = FALSE)))
                         })
                         names(cvRisk) <- libraryNames
                         obj_and_grad <- function(y, x, w = NULL) {
                         y <- y
                         x <- x
                         function(beta) {
                         xB <- x %*% cbind(beta)
                         loglik <- y * plogis(xB, log.p = TRUE) + (1 - y) * 
                         plogis(xB, log.p = TRUE, lower.tail = FALSE)
                         if (!is.null(w)) 
                         loglik <- loglik * w
                         obj <- - sum(loglik)
                         p <- plogis(xB)
                         grad <- if (is.null(w)) 
                         crossprod(x, cbind(p - y))
                         else crossprod(x, w * cbind(p - y))
                         list(objective = obj, gradient = grad)
                         }
                         }
                         r <- nloptr::nloptr(x0 = rep(1/ncol(Z), ncol(Z)), 
                         eval_f = obj_and_grad(transY, logitZ), 
                         lb = rep(0, ncol(Z)), ub = rep(1, ncol(Z)), 
                         eval_g_eq = function(beta) (sum(beta) - 1), 
                         eval_jac_g_eq = function(beta) rep(1,length(beta)), 
                         opts = list(algorithm = 'NLOPT_LD_SLSQP',xtol_abs = 1e-08))
                         if (r$status < 1 || r$status > 4) {
                         warning(r$message)
                         }
                         coef <- r$solution
                         if (any(is.na(coef))) {
                         warning('Some algorithms have weights of NA, setting to 0.')
                         coef[is.na(coef)] <- 0
                         }
                         coef[coef < 1e-04] <- 0
                         coef <- coef/sum(coef)
                         out <- list(cvRisk = cvRisk, coef = coef)
                         return(out)
                         },
                         computePred=function (predY, coef, control, 
                         lowerBound=",lowerBound,", upperBound=",upperBound,", ...){
                         plogis(trimLogit((predY-lowerBound)/(upperBound - lowerBound), trim = control$trimLogit) %*% matrix(coef))*(upperBound - lowerBound) + lowerBound
                         }
  )
                         ")))
  cat(name," (logit ensemble) read into Global environment using upperBound =",upperBound," and lowerBound =",lowerBound)
}

#======================================================================
# Summary function for SuperLearner objects based on bounded
# log-likelihood loss function, adopted from summary.CV.SuperLearner
#======================================================================
my.summary.CV.SuperLearner <- function (object, obsWeights = NULL, ...) 
{
  method <- ifelse(is.null(as.list(object$call)[["method"]]), 
                   "method.NNLS", as.list(object$call)[["method"]])
  library.names <- colnames(coef(object))
  V <- object$V
  n <- length(object$SL.predict)
  if (is.null(obsWeights)) {
    obsWeights <- rep(1, length(object$Y))
  }
  folds <- object$folds
  SL.predict <- object$SL.predict
  discreteSL.predict <- object$discreteSL.predict
  library.predict <- object$library.predict
  Y <- object$Y
  Risk.SL <- rep(NA, length = V)
  Risk.dSL <- rep(NA, length = V)
  Risk.library <- matrix(NA, nrow = length(library.names), 
                         ncol = V)
  rownames(Risk.library) <- library.names
  if (method %in% c("method.NNLS", "method.NNLS2", "method.CC_LS")) {
    for (ii in seq_len(V)) {
      Risk.SL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] - 
                                                       SL.predict[folds[[ii]]])^2)
      Risk.dSL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] - 
                                                        discreteSL.predict[folds[[ii]]])^2)
      Risk.library[, ii] <- apply(library.predict[folds[[ii]], 
                                                  , drop = FALSE], 2, function(x) mean(obsWeights[folds[[ii]]] * 
                                                                                         (Y[folds[[ii]]] - x)^2))
    }
    se <- (1/sqrt(n)) * c(sd(obsWeights * (Y - SL.predict)^2), 
                          sd(obsWeights * (Y - discreteSL.predict)^2), apply(library.predict, 
                                                                             2, function(x) sd(obsWeights * (Y - x)^2)))
  }
  else if (method %in% c("method.NNloglik", "method.CC_nloglik")) {
    for (ii in seq_len(V)) {
      Risk.SL[ii] <- -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]], 
                                                            log(SL.predict[folds[[ii]]]), 
                                                            log(1 - SL.predict[folds[[ii]]])))
      Risk.dSL[ii] <- -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]], 
                                                             log(discreteSL.predict[folds[[ii]]]), 
                                                             log(1 - discreteSL.predict[folds[[ii]]])))
      Risk.library[, ii] <- apply(library.predict[folds[[ii]], , drop = FALSE], 2, function(x) {
        -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]], 
                                               log(x), log(1 - x)))
      })
    }
    se <- rep.int(NA, (length(library.names) + 2))
  }
  else if (method %in% c("method.AUC")) {
    requireNamespace("cvAUC")
    for (ii in seq_len(V)) {
      Risk.SL[ii] <- cvAUC::cvAUC(predictions = SL.predict[folds[[ii]]], 
                                  labels = Y[folds[[ii]]], folds = NULL)$cvAUC
      Risk.dSL[ii] <- cvAUC::cvAUC(predictions = discreteSL.predict[folds[[ii]]], 
                                   labels = Y[folds[[ii]]], folds = NULL)$cvAUC
      Risk.library[, ii] <- apply(library.predict[folds[[ii]], 
                                                  , drop = FALSE], 2, function(x) cvAUC::cvAUC(predictions = x, 
                                                                                               labels = Y[folds[[ii]]], folds = NULL)$cvAUC)
    }
    se <- rep.int(NA, (length(library.names) + 2))
  }
  else {
    # get the upper and lower bounds -- not sure this is fool proof...
    grbg <- strsplit(paste(object$method)[2],split="=")[[1]][2:3] 
    a <- as.numeric(strsplit(grbg[1],",")[[1]][1])
    b <- as.numeric(strsplit(grbg[2],",")[[1]][1])
    for (ii in seq_len(V)) {
      # replace extrapolated predictions
      SL.pred <- SL.predict[folds[[ii]]]
      dSL.pred <- discreteSL.predict[folds[[ii]]]
      SL.pred[SL.pred <= a] <- SL.predict[folds[[ii]]][SL.pred <= a] <- a + 1e-5
      SL.pred[SL.pred >= b] <- SL.predict[folds[[ii]]][SL.pred >= b] <- b - 1e-5
      dSL.pred[dSL.pred <= a] <- discreteSL.predict[folds[[ii]]][dSL.pred <= a] <- a + 1e-5
      dSL.pred[dSL.pred >= b] <- discreteSL.predict[folds[[ii]]][dSL.pred >= b] <- b - 1e-5
      librarypred <- library.predict[folds[[ii]], , drop = FALSE]
      librarypred[librarypred <= a] <- library.predict[folds[[ii]], ][librarypred <= a] <- a + 1e-5
      librarypred[librarypred >= b] <- library.predict[folds[[ii]], ][librarypred >= b] <- b - 1e-5
      
      Risk.SL[ii] <- -mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] - a)/(b-a) * log((SL.pred-a)/(b-a))
                           + (1 - (Y[folds[[ii]]] - a)/(b-a)) * log(1 - (SL.pred-a)/(b-a)))
      
      Risk.dSL[ii] <- -mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] - a)/(b-a) * log((dSL.pred-a)/(b-a))
                            + (1 - (Y[folds[[ii]]] - a)/(b-a)) * log(1 - (dSL.pred-a)/(b-a)))
      Risk.library[, ii] <- apply(librarypred, 2, function(x) {
        -mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] - a)/(b-a) * log((x-a)/(b-a))
              + (1 - (Y[folds[[ii]]] - a)/(b-a)) * log(1 - (x-a)/(b-a)))
      })
    }
    se <- (1/sqrt(n)) * c(sd(obsWeights * (Y - a)/(b-a) * log((SL.predict-a)/(b-a))
                             + (1 - (Y - a)/(b-a)) * log(1 - (SL.predict-a)/(b-a))), 
                          sd(obsWeights * (Y - a)/(b-a) * log((discreteSL.predict-a)/(b-a))
                             + (1 - (Y - a)/(b-a)) * log(1 - (discreteSL.predict-a)/(b-a))), 
                          apply(library.predict, 2, function(x) sd(obsWeights * (Y - a)/(b-a) * log((x-a)/(b-a))
                                                                   + (1 - (Y - a)/(b-a)) * log(1 - (x-a)/(b-a)))))
  }
  Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL", 
                                    library.names), Ave = c(mean(Risk.SL), mean(Risk.dSL), 
                                                            apply(Risk.library, 1, mean)), se = se, Min = c(min(Risk.SL), 
                                                                                                            min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(max(Risk.SL), 
                                                                                                                                                                 max(Risk.dSL), apply(Risk.library, 1, max)))
  out <- list(call = object$call, method = method, V = V, Risk.SL = Risk.SL, 
              Risk.dSL = Risk.dSL, Risk.library = Risk.library, Table = Table)
  class(out) <- "my.summary.CV.SuperLearner"
  return(out)
}

#======================================================================
# Plotting function for my.summary.CV.SuperLearner based on 
# plot.CV.SuperLearner
#======================================================================
plot.my.CV.SuperLearner <- function (x, package = "ggplot2", constant = qnorm(0.975), sort = TRUE, 
                                     ...) 
{
  sumx <- my.summary.CV.SuperLearner(x, ...)
  if (sort) 
    sumx$Table$Algorithm <- reorder(sumx$Table$Algorithm, 
                                    -sumx$Table$Ave)
  Mean <- sumx$Table$Ave
  se <- sumx$Table$se
  Lower <- Mean - constant * se
  Upper <- Mean + constant * se
  assign("d", data.frame(Y = Mean, X = sumx$Table$Algorithm, 
                         Lower = Lower, Upper = Upper))
  if (package == "lattice") {
    SuperLearner:::.SL.require("lattice")
    p <- lattice::dotplot(X ~ Y, data = d, xlim = c(min(d$Lower) - 
                                                      0.02, max(d$Upper) + 0.02), xlab = "V-fold CV Risk Estimate", 
                          ylab = "Method", panel = function(x, y) {
                            lattice::panel.xyplot(x, y, pch = 16, cex = 1)
                            lattice::panel.segments(d$Lower, y, d$Upper, 
                                                    y, lty = 1)
                          })
  }
  if (package == "ggplot2") {
    SuperLearner:::.SL.require("ggplot2")
    p <- ggplot2::ggplot(d, ggplot2::aes_string(x = "X", 
                                                y = "Y", ymin = "Lower", ymax = "Upper")) + ggplot2::geom_pointrange() + 
      ggplot2::coord_flip() + ggplot2::ylab("V-fold CV Risk Estimate") + 
      ggplot2::xlab("Method")
  }
  return(p)
}


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

