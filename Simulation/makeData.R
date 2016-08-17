makeFCSData <- function(n){
  trt <- rbinom(n,1,0.5)
  age <- runif(n,18,90)
  female <- rbinom(n,1,0.33)
  sofa <- round(runif(n,3,19))
  race <- rbinom(n,1,0.25)
  score <- round(runif(n,0,6))
  ageInt <- age*trt
  femaleInt <- female*trt
  raceInt <- race*trt
  sofaInt <- sofa*trt
  scoreInt <- score*trt
  
  totalcost <- rep(NA,n)

  # first group 
  pg1 <- plogis(-2 + female + femaleInt)
  g1 <- rbinom(n,1,pg1)
  
  # second group 
  pg2 <- plogis(-2 + age/200 + age*score/1000)
  g2 <- rbinom(n,1,pg2)
  
  # assign g1=0, g2=0 costs
  ind <- g1==0 & g2==0
  totalcost[ind] <- 35000 + exp(9.5 + -trt[ind] + female[ind] - age[ind]/100 + rnorm(sum(ind),0,0.9))

  # assign g1=1, g2=0 costs
  ind <- g1==1 & g2==0
  totalcost[ind] <- 35000 + exp(9.5 - trt[ind] + sofa[ind]/20 + rnorm(sum(ind),0,0.9))
  
  # assign g1=0, g2=1 costs
  ind <- g1==0 & g2==1
  totalcost[ind] <- 5000 + exp(9.5 - trt[ind] + age[ind]*score[ind]/500 + rnorm(sum(ind),0,1))
  
  # assign g1=1, g2=1 costs
  ind <- g1==1 & g2==1
  totalcost[ind] <- runif(sum(ind),5000,40000)
  
  return(data.frame(pid=1:n,trt=trt,
                    age=age,female=female,sofa=sofa,totalcost=totalcost,
                    race=race,score=score,ageInt=ageInt,femaleInt=femaleInt,
                    raceInt=raceInt,sofaInt=sofaInt,scoreInt=scoreInt))
}