library(myTridge)
library(rsvd)
library(ggplot2)
library(glmnet)

num.obs <- 500
num.par <- 5000
family <- "gaussian"

data<-myTridge::genDataList(n = num.obs,
                            mu = rep(0, num.par),
                            p = num.par,
                            rho = 0.2,
                            beta = rnorm(num.par, mean = 0, sd = 1),
                            SNR = 10,
                            family = family)
X <- data$normData
y <- data$y
X_new <- cbind(rep(1,num.obs),X)
r <- 0.5 # tuning parameter


tt <- bench::mark(

  svd = {

    ##### the following commands should only be calculated once for all tuning parameters ####
    mysvd <- base::svd(X_new) # uses the base implementation of svd
    R <- mysvd$u %*% diag(mysvd$d)
    RTR <- crossprod(R) # R^TR
    IN <- diag(num.obs) # NxN identity matrix
    RTY <- crossprod(R,y) # R^TY
    #############################################################################################

    mysvd$v %*% (solve(RTR + r*IN) %*% RTY) # ridge estimator
  },

  rsvd = {

    ##### the following commands should only be calculated once for all tuning parameters ####
    mysvd <- rsvd::rsvd(X_new)   #uses the rsvd package
    R <- mysvd$u %*% diag(mysvd$d)
    RTR <- crossprod(R) # R^TR
    IN <- diag(num.obs) # NxN identity matrix
    RTY <- crossprod(R,y) # R^TY
    #############################################################################################

    mysvd$v %*% (solve(RTR + r*IN) %*% RTY) # ridge estimator
  },

  woodbury = {

    ##### the following commands should only be calculated once for all tuning parameters ####
    XXT <- tcrossprod(X_new)
    XTY <- crossprod(X_new,y)
    IN <- diag(num.obs) # NxN identity matrix
    XXTY <- XXT %*% y
    #############################################################################################

    XTY/r - crossprod(X_new, solve(IN + XXT/r) %*% XXTY) / r^2   # ridge estimator
  },


  glmnet = {
    coef(glmnet(X, y, alpha = 0, lambda = r))
  },

  check = FALSE,
  relative = FALSE
)

tt
ggplot2::autoplot(tt)





plot(mysvd$v %*% (solve(RTR + r*IN) %*% RTY), coef(glmnet(X, y, alpha = 0, lambda = r)))
abline(a=0,b=1, col="red",lwd=3)

fit <- glmnet(X, y, alpha = 0, lambda = r)

library(magrittr)
mysvd$v %*% (solve(RTR + r*IN) %*% RTY) %>% dim
fit$beta %>% dim


coef(fit) %>% dim









