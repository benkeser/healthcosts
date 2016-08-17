# Code used to plot simulation results
# out
load("rslt.RData") # result of executing sce.sh


# bias result
out$truth <- -13796.87 # computed numerically

biasRslt <- by(out, out$n, function(x){
  abs(colMeans(x[,2:24] - x$truth)/x$truth[1])
})

# variance result
varRslt <- by(out, out$n, function(x){
  means <- matrix(rep(colMeans(x[,2:24]),nrow(x)),nrow=nrow(x),byrow=TRUE)
  colMeans((x[,2:24] - means)^2)/10000
})

# mse results
mseRslt <- by(out, out$n, function(x){
  colMeans((x[,2:24] - x$truth)^2)/10000
})

algoLabels <- c("SL-MSE","SL-MSE2","SL-LLIK","SL-LLIK2",
                "Normal Id. GLM","Gamma Log GLM","Gamma Id. GLM",
                "Normal Id. Log GLM", "Adaptive GLM","AFT Gen. Gamma",
                "AFT Weibull", "AFT Log-Norm.","Quantile Reg.",
                "Cox PH", "Adaptive Hazard","Random Forest",
                "Regression Tree","GBM","CV-Reg. Tree", "CV-Random Forest",
                "CV-GBM","Mean","HAL")

### make results plot
makePlotRow <- function(rslt, 
                        xAdd=0.5, # what to add to x-labels to get them to line up
                        mult=100, # what to multiply rslts by (100 for bias 0.0001 or something small for var/mse)
                        ylim1,ylim2, # passed to plot (for ylim)
                        labMult1,labMult2 # to adjust location of labels
                        ){
## bias plots
par(mar=c(7.1, 1, 1, 2.5), mgp=c(2.1,0.5,0))
plot(0,0,pch="",ylab="", xlab="",
     xaxt="n",bty="n", xlim=c(1,23), ylim=ylim1)
# order of algorithms - SL's
ordA <- order(-rslt[[1]][-(1:4)])
usrRange <- diff(par()$usr[1:2])
text(y=par()$usr[1]-labMult1*usrRange,x=1:(23-4), labels=algoLabels[-(1:4)][ordA], srt=-45,
     pos=4,xpd=TRUE)
# add points for algorithms - SL's
points(x=(1:(23-4))+xAdd, y=rslt[[1]][-(1:4)][ordA]*mult)

# order of SL algorithms
ordSL <- order(-rslt[[1]][(1:4)])
text(y=par()$xpd[1]-labMult1*usrRange,x=20:23, labels=algoLabels[(1:4)][ordSL], srt=-45,
     pos=4,xpd=TRUE)
# add points for algorithms - SL's
points(x=(20:23)+xAdd, y=rslt[[1]][(1:4)][ordSL]*mult,pch=23,bg=1)


## for n=2000
par(mar=c(7.1, 1, 1, 2.5), mgp=c(2.1,0.5,0))
plot(0,0,pch="",ylab="", xlab="",
     xaxt="n",bty="n", xlim=c(1,23),ylim=ylim2)
# order of algorithms - SL's
usrRange <- diff(par()$usr[1:2])
ordA <- order(-rslt[[6]][-(1:4)])
text(y=par()$xpd[1]-labMult2*usrRange,x=1:(23-4), labels=algoLabels[-(1:4)][ordA], srt=-45,
     pos=4,xpd=TRUE)
# add points for algorithms - SL's
points(x=(1:(23-4))+xAdd, y=rslt[[6]][-(1:4)][ordA]*mult)

# order of SL algorithms
ordSL <- order(-rslt[[6]][(1:4)])
text(y=par()$xpd[1]-labMult2*usrRange,x=20:23, labels=algoLabels[(1:4)][ordSL], srt=-45,
     pos=4,xpd=TRUE)
# add points for algorithms - SL's
points(x=(20:23)+xAdd, y=rslt[[6]][(1:4)][ordSL]*mult,pch=23,bg=1)
}


# make the plot
layout(matrix(1:6,ncol=2,byrow=TRUE))
par(oma=c(0,4,2,0))
makePlotRow(biasRslt,ylim1=c(0,15),ylim2=c(0,10),mult=100,labMult1 = 0.05, labMult2=0.05)
makePlotRow(varRslt,ylim1=c(0,3),ylim2=c(0,0.3),mult=0.0005,labMult1 = 0.0025, labMult2 = 0.0025)
makePlotRow(mseRslt,ylim1=c(0,5),ylim2=c(0,0.6),mult=0.001,labMult1 = 0.0025, labMult2 = 0.0025)
mtext(outer=TRUE,side=2, at=0.885, "|Bias| (% of truth)",cex=0.75, line=1.5)
mtext(outer=TRUE,side=2, at=0.55, "Variance (x10e-7)",cex=0.75, line=1.5)
mtext(outer=TRUE,side=2, at=0.22, "Mean-squared error (x10e-7)",cex=0.75, line=1.5)
mtext(outer=TRUE,side=3, at=0.25, "n=100")
mtext(outer=TRUE,side=3, at=0.75, "n=2,000")











