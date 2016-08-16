#===============================================================
# Super Learner Analysis for FCS Data
#===============================================================

#=========================================
# Reading data, loading libraries
#=========================================
# Install and load packages
pkgs <- c("RCurl","SuperLearner","glmnet","gbm","rpart","randomForest","survival","moments","flexsurv",
          "sandwich","quantreg","caret")
for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
  require(pkg)
}

# load Super Learner functions from github
eval(parse(text=getURL("https://raw.githubusercontent.com/benkeser/healthcosts/master/SuperLearnerWrappers.R")))

# load analysis specific functions from github
eval(parse(text=getURL("https://raw.githubusercontent.com/benkeser/healthcosts/master/Mock%20FCS%20Analysis/AnalysisFunctions.R")))

# load the simulated data from github
healthdata <- getURL("https://raw.githubusercontent.com/benkeser/healthcosts/master/Mock%20FCS%20Analysis/healthdata.csv")
dat <- read.csv(textConnection(healthdata),header=TRUE)

# create interaction variables
dat$sofaInt <- dat$trt*dat$sofa
dat$scoreInt <- dat$trt*dat$score
dat$femaleInt <- dat$trt*dat$female
dat$raceInt <- dat$trt*dat$race
dat$ageInt <- dat$trt*dat$age

# make data frames to be used later by hc.tmle()
X <- dat[,which(!(colnames(dat) %in% c("pid","totalcost","directcost")))]
X1 <- X0 <- X
X1$trt <- 1; X0$trt <- 0
X1$ageInt <- X1$age
X1$femaleInt <- X1$female
X1$raceInt <- X1$race
X1$sofaInt <- X1$sofa
X1$score <- X1$score
X0[,grep("Int",names(X0))] <- 0

#=========================================
# Set up Super Learner Library
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

## add in HAL
mySLlibrary <- append(mySLlibrary, list(c("SL.hal","allScreen_noInt")))
mySLlibrary <- append(mySLlibrary,list(c("SL.hal","demoScreen")))
mySLlibrary <- append(mySLlibrary,list(c("SL.hal","medScreen")))

## add in unadjusted
mySLlibrary <- append(mySLlibrary,list(c("SL.mean","All")))

## to speed up analysis, you could truncate the library
nalgos <- 10
mySLlibrary <- mySLlibrary[1:nalgos]

# make methods using bounded log-likelihood loss
makeBoundedMethod(upperBound=max(dat$totalcost), lowerBound=min(dat$totalcost), name="methodTotal")
makeBoundedMethod(upperBound=max(dat$directcost), lowerBound=min(dat$direct), name="methodDirect")

#==================================================
# inference with SuperLearner + candidate methods
#==================================================
set.seed(1289525)

#==================
# Total ICU costs
#==================
# generate sample splits based on ordered total ICU costs
ordTC <- order(-dat$totalcost)
v <- rep(c(1:10,10:1),10)[1:nrow(dat)]
folds <- split(ordTC, factor(v))

# fit Super Learner for total costs
Y <- dat$totalcost
fm <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=mySLlibrary,
                   verbose=FALSE,method="methodTotal",
                   cvControl=list(V=10L, stratify=FALSE, shuffle=FALSE, validRows=folds))
totalOut <- hc.tmle(Y=Y,X=X,X0=X0,X1=X1,fm=fm,trt="trt")


#=================================
# Direct-variable ICU costs
#=================================
# generate sample splits based on ordered direct-variable ICU costs
ordDC <- order(-dat$directcost)
v <- rep(c(1:10,10:1),10)[1:nrow(dat)]
folds <- split(ordDC, factor(v))

Y <- dat$directcost
fmD <- SuperLearner(Y=Y,X=X,family=gaussian(), SL.library=mySLlibrary,
                    verbose=FALSE,method="methodDirect")
directOut <- hc.tmle(Y=Y,X=X,X0=X0,X1=X1,fm=fmD)


#=================================
# Plot the results
#=================================
# reduce the total cost results
slTotal <- rbind(totalOut[[1]]$diff, totalOut[[2]]$diff)
candTotal <- Reduce(rbind,lapply(totalOut[[3]],function(x){x$diff}))
K <- length(candTotal[,1]) + 2
ordCand <- order(-candTotal[,1])
# reduce the direct-variable cost results
slDirect <- rbind(directOut[[1]]$diff, directOut[[2]]$diff)
candDirect <- Reduce(rbind,lapply(directOut[[3]],function(x){x$diff}))
K <- length(candDirect[,1]) + 2
ordCandD <- order(-candDirect[,1])

# plotting
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


#===============================================
# CV-Super Learner
#===============================================
set.seed(128951925)
# total cost
Y <- dat$totalcost
fm <- CV.SuperLearner(Y = Y, X = X, V=10,
                      SL.library = mySLlibrary, verbose=TRUE,
                      method = "methodTotal")

# cross-validated risk plot
plot.my.CV.SuperLearner(fm)

# direct cost
Y <- dat$directcost
fmD <- CV.SuperLearner(Y=Y,X=X,V=20,family=gaussian(), SL.library=mySLlibrary,
                    verbose=TRUE, method="methodDirect")

# cross-validated risk plot
plot.my.CV.SuperLearner(fmD)






