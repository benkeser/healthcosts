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

