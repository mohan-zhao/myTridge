#' @title Tridge Estimator
#' @description generate a tridge estimator for objective function
#' @param X input matrix, of dimension nobs x nvars; each row is an observation
#' vector.
#' @param y response variable, num.obs elements
#' @param family Response type, should be a character string representing one of the built-in 
#' families \code{family="gaussian"}, \code{family="poisson"} and \code{family="binomial"}
#' @param nlambda The number of tuning parameters - default is 1000.
#' @param lambda.min.limit the smallest possible value for lambda
#' @param c a constant such that lambda.min = max(lambda.min.limit, (lambda.max - 2*c)) 
#' @param nfolds number of folds - default is 10
#' @param alpha the elastic-net mixing parameter α, default is 0 the ridge penalty. Alpha=1 is the lasso penalty
#' @param stn signal-to-noise ratio. Only useful for the \code{"gaussian"} family such that entries of the noise vector ε are sampled i.i.d. from N (0, stn,
#' default is 10. 
#' @return tridge estimator vector
#' @details First compute the stationary point for beta by calling optim_ObFn, then compute the ridge estimaor for m equally spaced tuning 
#'         parameters and find the one minimize the objective function.
#' @examples 
#' X = matrix(rnorm(20 * 100), 20, 100)
#' y = rnorm(20)
#' family = 'gaussian'
#' nlambda=1000
#' lambda.min.limit=0.05
#' c=0.1
#' nfolds=10
#' alpha=0
#' @rdname TridgeEst
#' @export 
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats rnorm coef
#' 
#'   
TridgeEst<-function(X,y,family=c("gaussian","binomial","poisson"),nlambda=1000,lambda.min.limit=0.05,c=0.1,nfolds=10,alpha=0,stn=10){
  family=match.arg(family)
  if(family=="gaussian") {
    TREX.c.vector <- 2 
  } else {
    TREX.c.vector <- 1
  }
  num.obs <- nrow(X)
  num.par <- ncol(X)
  
  ##### signal-to-noise ratio for gaussian family case 
  if(family == "gaussian") {
    stn <- stn # signal-to-noise ratio
  }else{
    stn <- NA
  }
  ###compute T-ridge estimators
  par0 <-rep(0,num.par)
  par1 <- optim_ObLs(par0,X,y,family)
  
  GlmTrex.estimators <- matrix(nrow=num.par, ncol=length(TREX.c.vector))
  for(trex.e in 1:length(TREX.c.vector)) {
    TREX.c <- TREX.c.vector[trex.e]
    ###find the stationary point
    GlmTrex.estimation <- optim_ObFn(par1,X,y,family,TREX.c)
    GlmTrex.estimators[, trex.e] <- GlmTrex.estimation
    
  }
  
  edr.lambda <- qnorm(2,GradientLs(GlmTrex.estimators,X,y,family)) / (2 * qnorm(2, GlmTrex.estimators))
  r.max <- edr.lambda + c
  r.min <- max(lambda.min.limit, (edr.lambda - c))
  if(r.min == Inf || r.max == Inf ) {
    tuning.parameters <- seq(1e+10, 1e+11, length.out=nlambda)
  } else {
    tuning.parameters <- seq(r.min, r.max, length.out=nlambda)
  }
  
  init <- cv.glmnet(X, y, nfolds=nfolds, alpha=alpha, family=family)
  init.estimation <- glmnet(X, y, alpha=alpha, family=family, lambda=init$lambda.min)
  init.vector <- as.vector(coef(init.estimation))[-1]
  
  Ridge.estimators <- matrix(nrow=num.par, ncol=nlambda)
  cost <- rep(0, nlambda)
  for(tune in 1:nlambda) {
    r <- tuning.parameters[tune]
    if(family=="gaussian"){
      estimation <- glmnet(X, y, alpha=alpha, family=family, lambda=r)
      Ridge.estimators[, tune] <- as.vector(stats::coef(estimation)[-1])
    }else {
      estimation <- optim_Ridge(init.vector,X,y,family,r)
      Ridge.estimators[, tune] <- estimation
    }
    
    ####compute ObjectFunc value to find minimum
    estimator.r <- as.vector(Ridge.estimators[, tune])
    cost[tune] <- ObjectiveFunction(estimator.r,family,y,X,TREX.c)
  }
  estimator <- Ridge.estimators[, which.min(cost)]
  
  return(estimator)
}