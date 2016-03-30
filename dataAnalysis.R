################################################################
# Super Learner Analysis for FCS Data
################################################################

#=========================================
# Reading data, loading libraries
#=========================================
# load libraries
library(SuperLearner)
library(caret)
library(foreign)
library(randomForest)
library(rpart)
library(gbm)
library(flexsurv)
library(tmle)
library(moments)

# load SL library functions
dataDir <- "C:/Users/David Benkeser/Dropbox/UW Classes/Consulting/Nita Khandelwal/Methods paper/Data"
codeDir <- "C:/Users/David Benkeser/Dropbox/UW Classes/Consulting/Nita Khandelwal/Methods paper/R code"
# library(snowfall)
# sfInit(parallel=TRUE, cpus=3, type="SOCK")

# dataDir <- "~/Dropbox/UW Classes/Consulting/Nita Khandelwal/Methods paper/Data"
# codeDir <- "~/Dropbox/UW Classes/Consulting/Nita Khandelwal/Methods paper/R code"
# sfExport("dataDir")
# sfExport("codeDir")
source(file.path(codeDir, "superLearnerLibrary.R"))

# read data
dat <- read.dta(file.path(dataDir, "cleanData.dta"),convert.factors=FALSE)
dat <- dat[-which(dat$pid=="HMC0350"),]
#=========================================
# Super Learning
#=========================================

# set up library as combination of learning algorithsm
learners <- c("SL.glm","SL.gammaLogGLM","SL.gammaIdentityGLM",        
              "SL.logOLS.smear","SL.manningGLM","SL.gengamma",
              "SL.weibull","SL.lognormalsurv",
              "SL.wangZhou","SL.gilleskie","SL.coxph",
              "SL.randomForest","SL.rpart","SL.gbm",
              "SL.rpart.caret1","SL.rf.caret1","SL.gbm.caret1"
              )
screens <- c("All","allScreen_noInt","demoScreen","demoScreen_ageInt",
             "demoScreen_raceInt","demoScreen_femaleInt","medScreen",
              "medScreen_sofaInt","medScreen_charlInt"
              )

mySLlibrary <- vector(length(length(learners)*length(screens)),mode="list")
ct <- 0
for(l in 1:length(learners)){ 
  for(s in 1:length(screens)){
    ct <- ct+1
    mySLlibrary[[ct]] <- c(learners[[l]],screens[[s]])
}}

getInd <- function(start){
  c(start,(start+3):(start+5),(start+7):(start+8))
}

mySLlibrary <- mySLlibrary[-c(getInd(145),getInd(136),getInd(127),
                              getInd(118),getInd(109),getInd(100),
                              getInd(82))]

set.seed(128951925)
X <- dat[,which(!(colnames(dat) %in% c("pid","totalcost","directcost")))]
X1 <- X0 <- X
X1$trt <- 1; X0$trt <- 0
X1$ageInt <- X1$age
X1$femaleInt <- X1$female
X1$raceInt <- X1$race
X1$sofaInt <- X1$sofa
X1$score <- X1$score
X0[,grep("Int",names(X0))] <- 0

#====================================
# inference with SuperLearner
#====================================
# total costs
# Y <- dat$totalcost
# fm <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=mySLlibrary,
#                    verbose=FALSE)
# save(fmD, file=file.path(dataDir,"slFit_total_20160321.RData"))

# direct cost
# Y <- dat$directcost
# fmD <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=mySLlibrary,
#                     verbose=FALSE)
# save(fmD, file=file.path(dataDir,"slFit_direct_20160321.RData"))

load(file.path(dataDir,"slFit_total_20160321.RData"))
sl_total <- my.tmle(Y=dat$totalcost,X=X,X0=X0,X1=X1,fm=fm)
coef_total <- fm$coef[fm$coef>0]
fm <- NULL # because R gets fussy with big objects

load(file.path(dataDir,"slFit_direct_20160321.RData"))
sl_direct <- my.tmle(Y=dat$directcost,X=X,X0=X0,X1=X1,fm=fmD)
coef_direct <- fmD$coef[fmD$coef>0]
fmD <- NULL # because R gets fussy with big objects

#===================================
# inference for candidate methods
#===================================
# #total cost
# Y <- dat$totalcost
# allOut <- lapply(mySLlibrary, function(x){
#   fm <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=list(x),
#                      verbose=TRUE,cvControl=list(V=2L))
#   my.tmle(Y=Y,X=X,X0=X0,X1=X1,fm=fm)
# })
# save(allOut,file="inf_AllMethods_total_20160321.RData")
# 
# # direct cost
# Y <- dat$directcost
# allOut <- lapply(mySLlibrary, function(x){
#   fm <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=list(x),
#                      verbose=TRUE,cvControl=list(V=2L))
#   my.tmle(Y=Y,X=X,X0=X0,X1=X1,fm=fm)
# })
# save(allOut,file="inf_AllMethods_direct_20160321.RData")

#==========
# plotting 
#==========
pdf(file.path(dataDir,"inferenceAllMethods.pdf"),height=12,width=4.5)
load(file.path(dataDir, "inf_AllMethods_total_20160321.RData"))
plotInferenceResults(allOut=allOut, SLlibrary=mySLlibrary, xl=c(-50,10),
                     SLrslt=sl_total$diff,t1=60,lmar=11.1)
load(file.path(dataDir, "inf_AllMethods_direct_20160321.RData"))
plotInferenceResults(allOut=allOut, SLlibrary=mySLlibrary, xl=c(-8,2),
                     SLrslt=sl_direct$diff,t1=10,lmar=11.1)
dev.off()


#===============================================
# CV-Super Learner
#===============================================
# set.seed(128951925)
# # total cost
# 
# # cl <- makeCluster(2, type = "PSOCK") # can use different types here
# # clusterSetRNGStream(cl, iseed = 2343)
# fm <- CV.SuperLearner(Y = Y, X = X, V=20,
#                       SL.library = mySLlibrary, verbose=TRUE,
#                       method = "method.NNLS",
#                       parallel="seq")
# 
# save(fm, file=file.path(dataDir,"slCVFit_total_20160321.RData"))
# 
# # direct cost
# Y <- dat$directcost
# fmD <- CV.SuperLearner(Y=Y,X=X,V=20,family=gaussian(), SL.library=mySLlibrary,
#                     verbose=TRUE)
# save(fmD, file=file.path(dataDir,"slCVFit_direct_20160321.RData"))
# #
load(file.path(dataDir,"slCVFit_total_20160321.RData"))
load(file.path(dataDir,"slCVFit_direct_20160321.RData"))

pdf(file.path(dataDir,"cvNEW.pdf"),heigh=15,width=6)
plot(fm)
plot(fmD)
dev.off()





