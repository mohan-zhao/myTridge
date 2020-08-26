#' @title Tridge Estimator
#' @description generate a tridge estimator for objective function
#' @param Test.case three different function distributions
#' @param num.obs num of rows of the matrix X
#' @param num.par num of volomns of the matrix y
#' @param k magnitude of the mutual correlations
#' @return tridge estimator vector
#' @details DETAILS
#' @examples TridgeEst(Test.case="gaussian",num.obs=100,num.par=300,k=0.0)
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname TridgeEst
#' @export 
TridgeEst<-function(Test.case,num.obs,num.par,k){
  qnorm <- function(q, v){
    (sum(abs(v) ^ q)) ^ (1 / q)
  }
  
  if(Test.case=="gaussian") {
    TREX.c.vector <- 2 
  } else {
    TREX.c.vector <- 1
  }
  
  class.num <- 1
  
  ##### Required package
  
  library(glmnet) # Generalized linear model
  library(MASS) # Multivariate normal
  
  ##### signal-to-noise ratio for gaussian test case 
  if(Test.case == "gaussian") {
    stn <- 10 # signal-to-noise ratio
  } else {
    stn <- NA
  }
  
  mu <- rep(0, num.par)
  beta <- rnorm(num.par)
  
 
    
    cat("  -> generate test case... ")
    cat("  -> Normalize each column of the design matrix... ")
    
    
    ret <- genDataList(num.obs, mu, num.par,
                       k, beta, stn,Test.case)
    
    X <- ret$normData
    
    beta <-ret$beta
    y <- as.vector(ret$y)
    
    
   
    cat("  -> Perform the K-fold cross validation pipeline... ")
    if(num.obs < num.par) {
      cv <- cv.glmnet(X, y, nfolds=10, alpha=0, family=Test.case)
      
      ##Fit a generalized linear model via penalized maximum likelihood
      
      cv.estimation <- glmnet(X, y, alpha=0, family=Test.case, lambda=cv$lambda.min)
      
      CV.estimator <- as.vector(coef(cv.estimation))[-1]
    } else {
      cv.estimation <- glmnet(X, y, alpha=0, family=Test.case, lambda=0)
      CV.estimator <- as.vector(coef(cv.estimation))[-1]
    }
    
    
    cat("  -> compute T-ridge estimators... ")
    par0 <-rep(0,num.par)
    ###does type=1 make a difference in this optim?
    par1 <- optim_ObLs(par0,X,y,Test.case)
    
    GlmTrex.estimators <- matrix(nrow=num.par, ncol=length(TREX.c.vector))
    for(trex.e in 1:length(TREX.c.vector)) {
      TREX.c <- TREX.c.vector[trex.e]
      ###find the stationary pointc
      GlmTrex.estimation <- optim_ObFn(par1,X,y,Test.case,TREX.c)
      GlmTrex.estimators[, trex.e] <- GlmTrex.estimation
      
    }
    
    length.tuning <- 100
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
      estimation <- glmnet(X, y, alpha=0, family=Test.case, lambda=r)
      ###did not change this following line to improve the speed
      Ridge.estimators[, tune] <- as.vector(coef(estimation)[-1])
      
      ####compute ObjectFunc value to find minimum
      estimator.r <- as.vector(Ridge.estimators[, tune])
      cost[tune] <- ObjectiveFunction(estimator.r,Test.case,y,X,TREX.c)
    }
    estimator <- Ridge.estimators[, which.min(cost)]
    cat("done\n")
    return(estimator)
}
