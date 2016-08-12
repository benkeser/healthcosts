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
codeDir <- "C:/Users/David Benkeser/Dropbox/UW Classes/Consulting/Nita Khandelwal/Methods paper/healthCosts"
outDir <-  "C:/Users/David Benkeser/Dropbox/UW Classes/Consulting/Nita Khandelwal/Methods paper/Manuscript"

# library(snowfall)
# sfInit(parallel=TRUE, cpus=3, type="SOCK")

# dataDir <- "~/Dropbox/UW Classes/Consulting/Nita Khandelwal/Methods paper/Data"
# codeDir <- "~/Dropbox/UW Classes/Consulting/Nita Khandelwal/Methods paper/R code"
# sfExport("dataDir")
# sfExport("codeDir")
source(file.path(codeDir, "healthcost.R"))

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

load(file.path(dataDir,"slFit_total_stratCV_20160801"))
sl_total <- my.tmle(Y=dat$totalcost,X=X,X0=X0,X1=X1,fm=fm)
coef_total <- fm$coef[fm$coef>0]
fm <- NULL # because R gets fussy with big objects

load(file.path(dataDir,"slFit_direct_stratCV_20160801"))
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
load(file.path(dataDir, "inf_AllMethods_total_20160801.RData"))
plotInferenceResults(allOut=allOut, SLlibrary=mySLlibrary, xl=c(-50,10),
                     SLrslt=sl_total$diff,t1=60,lmar=11.1)
load(file.path(dataDir, "inf_AllMethods_direct_20160321.RData"))
plotInferenceResults(allOut=allOut, SLlibrary=mySLlibrary, xl=c(-8,2),
                     SLrslt=sl_direct$diff,t1=10,lmar=11.1)
dev.off()


## new plot
load(file.path(dataDir, "inf_AllMethods_stratCV_total_20160801.RData"))
slTotal <- rbind(totalOut[[1]]$diff, totalOut[[2]]$diff)
candTotal <- Reduce(rbind,lapply(totalOut[[3]],function(x){x$diff}))
K <- length(candTotal[,1]) + 2
ordCand <- order(-candTotal[,1])

slDirect <- rbind(directOut[[1]]$diff, directOut[[2]]$diff)
candDirect <- Reduce(rbind,lapply(directOut[[3]],function(x){x$diff}))
K <- length(candDirect[,1]) + 2
ordCandD <- order(-candDirect[,1])


pdf(file.path(outDir,"allATE.pdf"),height=5.25,width=4)
layout(t(c(1,2)))
par(mar=c(1.1,0.5,0.5,0.5),oma=c(2,1,1,0),mgp=c(2.2,0.6,0))
## total
plot(x=0,y=0,pch="",xlab="",bty="n",
     yaxt="n",xlim=c(min(candTotal[,2]),max(candTotal[,3])),ylab="",
     ylim=c(0,K-2))
for(k in (K-2):1){
  points(x=candTotal[ordCand[k],1],y=k+2,pch=19,col="gray50")
  segments(x0=candTotal[ordCand[k],2],x1=candTotal[ordCand[k],3],
           y0=k+2,lwd=1.5,col="gray50")
}
abline(v=0,lty=3,lwd=1.5)
points(x=slTotal[1,1],y=0,pch=23,cex=1.2,bg=1)
segments(x0=slTotal[1,2],x1=slTotal[1,3],
         y0=0,lwd=1.5)
## direct 
plot(x=0,y=0,pch="",xlab="",bty="n",
     yaxt="n",xlim=c(min(candDirect[,2]),max(candDirect[,3])),ylab="",
     ylim=c(0,K-2))
for(k in (K-2):1){
  points(x=candDirect[ordCandD[k],1],y=k+2,pch=19,col="gray50")
  segments(x0=candDirect[ordCandD[k],2],x1=candDirect[ordCandD[k],3],
           y0=k+2,lwd=1.5,col="gray50")
}
abline(v=0,lty=3,lwd=1.5)
points(x=candDirect[1,1],y=0,pch=23,cex=1.2,bg=1)
segments(x0=candDirect[1,2],x1=candDirect[1,3],
         y0=0,lwd=1.5)
mtext(side=1,outer=TRUE,text=expression(ATE[n]^"*"*paste(" (95% CI)")),
      line=1)
mtext(side=2,outer=TRUE,text="Candidate Methods",line=0)
mtext(side=2,outer=TRUE,text="SL",line=0,at=0.09)
mtext(side=3,outer=TRUE,text="Total ICU costs",line=0,at=0.25)
mtext(side=3,outer=TRUE,text="Direct-variable ICU costs",line=0,at=0.75)

dev.off()






load(file.path(dataDir, "inf_AllMethods_stratCV_total_20160801.RData"))
allEstDirect <- 


#==========
# find 3 largest and 3 smallest estimates
#==========
load(file.path(dataDir, "inf_AllMethods_total_20160321.RData"))
ateEst <- unlist(lapply(allOut, function(x){ x$diff[1] }))
bigThree <- which(rank(ateEst) %in% 1:3)
smallThree <- which(rank(ateEst) %in% max(rank(ateEst)):(max(rank(ateEst))-3))
allOut[bigThree]
allOut[smallThree]
load(file.path(dataDir, "inf_AllMethods_direct_20160321.RData"))
ateEst <- unlist(lapply(allOut, function(x){ x$diff[1] }))
bigThree <- which(rank(ateEst) %in% min(rank(ateEst)):(min(rank(ateEst)+3)))
smallThree <- which(rank(ateEst) %in% max(rank(ateEst)):(max(rank(ateEst))-3))

mySLlibrary[bigThree]
allOut[bigThree]

mySLlibrary[smallThree]
allOut[smallThree]




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
load(file.path(dataDir,"slCVFit_total_20160801.RData"))
load(file.path(dataDir,"slCVFit_direct_20160801.RData"))
## default formatting
pdf(file.path(dataDir,"cvNEW_20160801.pdf"),height=10,width=4.5)
plot.my.CV.SuperLearner(fm)
plot.my.CV.SuperLearner(fmD)
dev.off()

## nicer formatting
pdf(file.path(dataDir,"cvNEW_20160801.pdf"),height=10,width=4.5)
plotCVResult(fm, SLlibrary=mySLlibrary, t1=3.5e9,lmar=10.1)
plotCVResult(fmD, SLlibrary=mySLlibrary, t1=9.5e7, lmar=10.1)
dev.off()





