#' @title Tridge Estimator
#' @description generate a tridge estimator for objective function
#' @param Test.case three different function distributions
#' @param X matrix with num.obs x num.par dimension
#' @param y vector, num.obs elements
#' @return tridge estimator vector
#' @details First compute the stationary point for beta by calling optim_ObFn, then compute the ridge estimaor for m equally spaced tuning 
#'         parameters and find the one minimize the objective function.
#' @examples 
#' Test.case="gaussian"
#' X=matrix(c(1:12),3,4)
#' y=c(2,4,5)
#' @rdname TridgeEst
#' @export 
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats rnorm coef
#' 
#'   
TridgeEst<-function(Test.case,X,y){
  qnorm <- function(q, v){
    (sum(abs(v) ^ q)) ^ (1 / q)
  }
  
  if(Test.case=="gaussian") {
    TREX.c.vector <- 2 
  } else {
    TREX.c.vector <- 1
  }
  num.obs <- nrow(X)
  num.par <- ncol(X)
  class.num <- 1
  
  ##### signal-to-noise ratio for gaussian test case 
  if(Test.case == "gaussian") {
    stn <- 10 # signal-to-noise ratio
  }else{
    stn <- NA
  }
  ###compute T-ridge estimators
  par0 <-rep(0,num.par)
  par1 <- optim_ObLs(par0,X,y,Test.case)
  
  GlmTrex.estimators <- matrix(nrow=num.par, ncol=length(TREX.c.vector))
  for(trex.e in 1:length(TREX.c.vector)) {
    TREX.c <- TREX.c.vector[trex.e]
    ###find the stationary point
    GlmTrex.estimation <- optim_ObFn(par1,X,y,Test.case,TREX.c)
    GlmTrex.estimators[, trex.e] <- GlmTrex.estimation
    
  }
  
  length.tuning <- 1000
  edr.lambda <- qnorm(2,GradientLs(GlmTrex.estimators,X,y,Test.case)) / (2 * qnorm(2, GlmTrex.estimators))
  r.max <- edr.lambda + 0.1
  r.min <- max(0.05, (edr.lambda - 0.1))
  if(r.min == Inf || r.max == Inf ) {
    tuning.parameters <- seq(1e+10, 1e+11, length.out=length.tuning)
  } else {
    tuning.parameters <- seq(r.min, r.max, length.out=length.tuning)
  }
  
  init <- cv.glmnet(X, y, nfolds=10, alpha=0, family=Test.case)
  init.estimation <- glmnet(X, y, alpha=0, family=Test.case, lambda=init$lambda.min)
  init.vector <- as.vector(coef(init.estimation))[-1]
  
  Ridge.estimators <- matrix(nrow=num.par, ncol=length.tuning)
  cost <- rep(0, length.tuning)
  for(tune in 1:length.tuning) {
    r <- tuning.parameters[tune]
    if(Test.case=="gaussian"){
      estimation <- glmnet(X, y, alpha=0, family=Test.case, lambda=r)
      Ridge.estimators[, tune] <- as.vector(stats::coef(estimation)[-1])
    }else {
      estimation <- optim_Ridge(init.vector,X,y,Test.case,r)
      Ridge.estimators[, tune] <- estimation
    }
    
    ####compute ObjectFunc value to find minimum
    estimator.r <- as.vector(Ridge.estimators[, tune])
    cost[tune] <- ObjectiveFunction(estimator.r,Test.case,y,X,TREX.c)
  }
  estimator <- Ridge.estimators[, which.min(cost)]
  
  return(estimator)
}