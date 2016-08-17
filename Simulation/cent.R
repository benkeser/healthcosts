#! /usr/bin/env Rscript


# This file was used to submit the simulation files to 
# a slurm-based Unix system. Using the sce.sh shell script
# one can submit each simulation in sequence. First, data files are
# created for each simulation. Those data files are then analyzed 
# in the 'run' execution. Then the results are collated in the 'merge'
# execution. 

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

ns <- c(100,250,500,750,1000,2000)
bigB <- 1000

# 
# # 
# # simulation parameters
parm <- expand.grid(seed=1:bigB,
                    n=ns)
# source in functions 
source("~/hc/sim/healthcosts.R")
source("~/hc/sim/makeData.R")

# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  for(i in 1:nrow(parm)){
    set.seed(parm$seed[i])
    dat <- makeFCSData(n=parm$n[i])
    save(dat, file=paste0("~/hc/sim/scratch/inFile_n=",parm$n[i],"_seed=",parm$seed[i],".RData"))
  }
  print(paste0('initial datasets saved to: ~/hc/sim/inFile ... .RData'))
}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id <- as.numeric(args[2])
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    print(paste(Sys.time(), "i:" , i))
    print(parm[i,])
    
    # load data
    load(paste0("~/hc/sim/scratch/inFile_n=",parm$n[i],"_seed=",parm$seed[i],".RData"))
    
    # load libraries
    library(SuperLearner)
    library(caret)
    
    algo <- c("SL.glm","SL.gammaLogGLM","SL.gammaIdentityGLM",        
                  "SL.logOLS.smear","SL.manningGLM","SL.gengamma",
                  "SL.weibull","SL.lognormalsurv",
                  "SL.wangZhou","SL.coxph","SL.gilleskie",
                  "SL.randomForest","SL.rpart","SL.gbm",
                  "SL.rpart.caret1","SL.rf.caret1","SL.gbm.caret1",
                  "SL.mean","SL.hal")
    screen <- c(rep("All",10),rep("allScreen_noInt",9))
    
    learners <- split(cbind(algo,screen),1:19)
    
    ## for sample splitting
    # spread outliers amongst data with snaking
    ordTC <- order(-dat$totalcost)
    v <- rep(c(1:10,10:1),1000)[1:nrow(dat)]
    folds <- split(ordTC, factor(v))
    
    Y <- dat$totalcost
    X <- dat[,which(!(colnames(dat) %in% c("pid","totalcost","directcost")))]
    
    ## first a super learner with random sample splitting and MSE
    set.seed(parm$seed[i])
    sl.mse.nosplit <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=learners,
                                   verbose=TRUE,method=method.NNLS)
    
    ## now a super learner with random sample splitting and MSE and splitting outliers up
    sl.mse.split <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=learners,
                                 verbose=TRUE,method=method.NNLS,
                                 cvControl=list(V=10L, stratify=FALSE, shuffle=FALSE, 
                                                validRows=folds))
    
    ## now a super learner with log-likelihood loss
    makeBoundedMethod(upperBound=max(dat$totalcost), lowerBound=min(dat$totalcost), name="methodTotal")
    set.seed(parm$seed[i])
    sl.llikB.nosplit <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=learners,
                                     verbose=FALSE,method=methodTotal)
    
    sl.llikB.split <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=learners,
                                   verbose=FALSE,method=methodTotal,
                                   cvControl=list(V=10L, stratify=FALSE, shuffle=FALSE, validRows=folds))
    
    ## get inference
    X1 <- X0 <- X
    X1$trt <- 1; X0$trt <- 0
    X1$ageInt <- X1$age
    X1$femaleInt <- X1$female
    X1$raceInt <- X1$race
    X1$sofaInt <- X1$sofa
    X1$score <- X1$score
    X0[,grep("Int",names(X0))] <- 0
    
    tmleOut1 <- hc.tmle(Y=Y,X=X,X0=X0,X1=X1,fm=sl.mse.nosplit,onlySL=FALSE)
    tmleOut2 <- hc.tmle(Y=Y,X=X,X0=X0,X1=X1,fm=sl.mse.split,onlySL=TRUE)
    tmleOut3 <- hc.tmle(Y=Y,X=X,X0=X0,X1=X1,fm=sl.llikB.nosplit,onlySL=TRUE)
    tmleOut4 <- hc.tmle(Y=Y,X=X,X0=X0,X1=X1,fm=sl.llikB.split,onlySL=TRUE)
    out <- list(tmleOut1,tmleOut2,tmleOut3,tmleOut4)
    save(out, file=paste0("~/hc/sim/out/outInf_n=",parm$n[i],"_seed=",parm$seed[i],".RData"))
    }
  
  
}

# merge job ###########################
if (args[1] == 'merge') {
  n <- sl.mse.ns <- sl.mse.s <- sl.ll.ns <- sl.ll.s <- cand <- NULL
  truth <- -13796.87
  allf <- list.files("~/hc/sim/out")
  for(i in allf){
    out <- get(load(paste0("~/hc/sim/out/",i)))
    tmp <- strsplit(i, "=")
    n <- c(n, as.numeric(strsplit(tmp[[1]][2],"_seed")[[1]][1]))
    sl.mse.ns <- rbind(sl.mse.ns, out[[1]][[1]]$diff[1])
    sl.mse.s <- rbind(sl.mse.s, out[[2]][[1]]$diff[1])
    sl.ll.ns <- rbind(sl.ll.ns, out[[3]][[1]]$diff[1])
    sl.ll.s <- rbind(sl.ll.s, out[[4]][[1]]$diff[1])
    cand <- rbind(cand, Reduce(c, lapply(out[[1]][[3]], function(x){x[[3]][1]})))
  }
  algo <- c("SL.glm","SL.gammaLogGLM","SL.gammaIdentityGLM",        
            "SL.logOLS.smear","SL.manningGLM","SL.gengamma",
            "SL.weibull","SL.lognormalsurv",
            "SL.wangZhou","SL.coxph","SL.gilleskie",
            "SL.randomForest","SL.rpart","SL.gbm",
            "SL.rpart.caret1","SL.rf.caret1","SL.gbm.caret1",
            "SL.mean","SL.hal")
  
  out <- data.frame(n, sl.mse.ns,sl.mse.s,sl.ll.ns,sl.ll.s,cand)
  names(out) <- c("n",paste0("sl.",1:4), algo)
  out$truth <- -13796.87
  save(out, file="~/hc/rslt.RData")
}