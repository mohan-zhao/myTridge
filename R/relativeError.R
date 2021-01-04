#' @title relative error table
#' @description compare errors of estimator for tridge with MLE or 10 fold CV
#' @param family Response type, should be a character string representing one of the built-in 
#' families \code{family="gaussian"}, \code{family="poisson"} and \code{family="binomial"}
#' @param num.obs num of rows of the matrix X
#' @param num.par num of volomns of the matrix y
#' @param k magnitude of the mutual correlations
#' @param num.runs num of runs to estimate and compute error
#' @param stn signal-to-noise ratio. Only useful for the \code{"gaussian"} family such that 
#' entries of the noise vector ε are sampled i.i.d. from N (0, stn, default is 10. 
#' @param nfolds number of folds for comparsion between cross-validation estimator 
#' error and t-ridge estimator error,default is 10
#' @param alpha the elastic-net mixing parameter α, default is 0 the ridge penalty. Alpha=1 is the lasso penalty
#' @return a table showing the errors for different methods
#' @details call TridgeEst to computer the t-ridge estimator then compare the relative errors
#' @examples 
#' family="gaussian"
#' num.obs=100
#' num.par=300
#' k=0.0
#' num.runs=5
#' @rdname relativeError
#' @export 
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom Matrix mean
#' @importFrom stats rnorm coef
#' @importFrom htmlTable htmlTable
#' @importFrom pander pander

relativeError <-function(family=c("gaussian","binomial","poisson"),num.obs,num.par,k,num.runs,stn=10,nfolds=10,alpha=0){
  family=match.arg(family)
  if(family=="gaussian") {
    TREX.c.vector <- 2 
  } else {
    TREX.c.vector <- 1
  }
  
  ##### signal-to-noise ratio for gaussian family case 
  if(family == "gaussian") {
    stn <- stn # signal-to-noise ratio
  } else {
    stn <- NA
  }
  
  # Initialize absolute error
  errors.TREX <- matrix(nrow=num.runs, ncol=length(TREX.c.vector))
  relative.errors.TREX <- matrix(nrow=num.runs, ncol=length(TREX.c.vector))
  errors.CV <- rep(NA, num.runs)
  relative.errors.CV <- rep(NA, num.runs)
  
  
  beta.error.TREX <- matrix(nrow=num.runs, ncol=length(TREX.c.vector))
  relative.beta.error.TREX <- matrix(nrow=num.runs, ncol=length(TREX.c.vector))
  beta.error.CV <- rep(NA, num.runs)
  relative.beta.error.CV <- rep(NA, num.runs)
  
  for (run in 1:num.runs){
    cat(sprintf("run %d out of %d runs\n", run, num.runs))
    mu <- rep(0, num.par)
    beta <- stats::rnorm(num.par)
    
    ret <- genDataList(num.obs, mu, num.par,
                       k, beta, stn,family)
    
    X <- ret$normData
    
    beta <-ret$beta
    y <- as.vector(ret$y)
    
    
    
    
    cat("  -> Perform the K-fold cross validation pipeline... ")
    if(num.obs < num.par) {
      cv <- glmnet::cv.glmnet(X, y, nfolds=nfolds, alpha=alpha, family=family)
      
      ##Fit a generalized linear model via penalized maximum likelihood
      
      cv.estimation <- glmnet(X, y, alpha=alpha, family=family, lambda=cv$lambda.min)
      
      CV.estimator <- as.vector(stats::coef(cv.estimation))[-1]
    } else {
      cv.estimation <- glmnet(X, y, alpha=alpha, family=family, lambda=0)
      CV.estimator <- as.vector(stats::coef(cv.estimation))[-1]
    }
    
    
    cat("  -> compute T-ridge estimators... ")
    estimator <- TridgeEst(X,y,family)
    
    ##### Compute the prediction errors 
    cat("  -> compute errors... ")
    for(trex.err in 1:length(TREX.c.vector)) {
      errors.TREX[run, trex.err] <- qnorm(2, abs(X %*% (estimator - beta))) / sqrt(num.obs)
    }
    errors.CV[run] <-qnorm(2, abs(X %*% (CV.estimator - beta))) / sqrt(num.obs) 
    
    #Relative Prediction error
    relative.errors.TREX[run, ] <- errors.TREX[run, ] / qnorm(2, (X %*% beta)) * sqrt(num.obs)
    relative.errors.CV[run] <- errors.CV[run] / qnorm(2, (X %*% beta)) * sqrt(num.obs)
  
    # beta error
    for(trex.b.err in 1:length(TREX.c.vector)) {
      beta.error.TREX[run, trex.b.err] <- qnorm(2, (estimator - beta))
    }
    beta.error.CV[run] <- qnorm(2, (CV.estimator - beta))
    
    #relative beta error
    relative.beta.error.TREX[run, ] <- beta.error.TREX[run, ] / qnorm(2, beta)
    relative.beta.error.CV[run] <- beta.error.CV[run] / qnorm(2, beta)
    #cat("done\n")
  }
  ##### Output results
  output.Data <- matrix(nrow = (length(TREX.c.vector) + 1), ncol = 4)
  for(case in 1:length(TREX.c.vector)){
    output.Data[case, ] <- c(Matrix::mean(errors.TREX[, case]), Matrix::mean(beta.error.TREX[, case]), Matrix::mean(relative.errors.TREX[, case]), Matrix::mean(relative.beta.error.TREX[, case]))
  }
  output.Data[(case + 1), ] <- c(Matrix::mean(errors.CV), Matrix::mean(beta.error.CV), Matrix::mean(relative.errors.CV), Matrix::mean(relative.beta.error.CV))
  names <- rep("", length(TREX.c.vector))
  for(name in 1:length(TREX.c.vector)){
    names[name] <- paste0("T-ridge", TREX.c.vector[name])
  }
  
  if(num.obs < num.par) {
    rownames(output.Data) <- c(names, "K-fold CV")
  } else {
    rownames(output.Data) <- c(names, "MLE")
  }
  
  colnames(output.Data)<- c("Xbeta", "Beta", "r.Xbeta", "r.Beta")
  #cat("Simulation Results : \n")
  
  output <- 
    matrix(c(round(output.Data[1, 3], 2),
             round(output.Data[2, 3], 2),
             round(output.Data[1, 4], 2),
             round(output.Data[2, 4], 2)), 
           ncol=4, byrow = TRUE)
  
  
  if(num.obs < num.par) {
    H <- c("T-ridge", "K-fold CV",
           "T-ridge", "K-fold CV")
  } else {
    H <- c("T-ridge(sd)", "MLE(sd)",
           "T-ridge(sd)", "MLE(sd)")
  }
  
  htmlTable::htmlTable(output,
                       header =  H,
                       rnames = "Mean relative errors",
                       rgroup = paste0("Case : ", family),
                       n.rgroup = c(1),
                       cgroup = c(pander::pander("$\\frac{||X\\hat{\\beta}_{T-ridge} - X\\beta^{*}||_{2}}{||X\\beta^{*}||_{2}}$"), pander::pander("$\\frac{||\\hat{\\beta}_{T-ridge} - \\beta^{*}||_{2}}{||\\beta^{*}||_{2}}$")),
                       n.cgroup = c(2,2), 
                       caption=paste0("(n,p,k,K)=(", num.obs, ",", num.par, ",", k,",",nfolds, ")"))
}
