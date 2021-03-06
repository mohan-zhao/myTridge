# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Generate a design matrix and normalize each column of the design matrix
#' @param n number of rows of X
#' @param mu a zero vector of length p where p is number of cols of X
#' @param p number of columns of X
#' @param rho where Σuv = rho to the power of |u−v|
#' @export
mvrnormArma <- function(n, mu, p, rho) {
    .Call(`_myTridge_mvrnormArma`, n, mu, p, rho)
}

#' return a poisson distributed vector
#' @param lambda a vector which is the mean of result vector
#' @export
rpois_rcpp <- function(lambda) {
    .Call(`_myTridge_rpois_rcpp`, lambda)
}

#' return a binomial distributed vector
#' @param n number of rows of X
#' @param size 1
#' @param prob means of output vector
#' @export
cpprbinom <- function(n, size, prob) {
    .Call(`_myTridge_cpprbinom`, n, size, prob)
}

#' Generated X from a p-dimensional normal distribution with mean 0p and covariance matrix 
#' and y = Xβ + ε 
#' β is sampled i.i.d. from N (0, 1) and then projected onto the row space of X to ensure identifiability
#'  
#' @param n number of rows of X
#' @param mu a zero vector of length p where p is number of cols of X
#' @param p number of cols of X
#' @param rho where Σuv = rho to the power of |u−v|
#' @param beta β is sampled i.i.d. from N (0, 1) with length p
#' @param SNR signal-to-noise ratio for, 10 for gaussian family case and NaN for poisson and binomial cases
#' @param family three family cases: gaussian, poisson and binomial
#' @export
genDataList <- function(n, mu, p, rho, beta, SNR, family) {
    .Call(`_myTridge_genDataList`, n, mu, p, rho, beta, SNR, family)
}

#' bfunction
#' @param X high dimensional matrix
#' @param theta unknown regression vector
#' @param family three family cases: gaussian, poisson and binomial
#' @export
bFunction <- function(X, theta, family) {
    .Call(`_myTridge_bFunction`, X, theta, family)
}

#' mean function
#' @param X high dimensional matrix
#' @param theta unknown regression vector
#' @param family three family cases: gaussian, poisson and binomial
#' @export
MeanFunction <- function(X, theta, family) {
    .Call(`_myTridge_MeanFunction`, X, theta, family)
}

#' mean prime function
#' @param X high dimensional matrix
#' @param theta unknown regression vector
#' @param family three family cases: gaussian, poisson and binomial
#' @export
MeanPrime <- function(X, theta, family) {
    .Call(`_myTridge_MeanPrime`, X, theta, family)
}

#' objective function
#' @param theta unknown regression vector
#' @param family three family cases: gaussian, poisson and binomial
#' @param y where y = Xβ + ε
#' @param X high dimensional matrix
#' @param trex_c TREX parameters, 2 for gaussian case, 1 for poisson and binomial case
#' @export
ObjectiveFunction <- function(theta, family, y, X, trex_c) {
    .Call(`_myTridge_ObjectiveFunction`, theta, family, y, X, trex_c)
}

#' Gradient of objective function
#' @param theta unknown regression vector
#' @param family three family cases: gaussian, poisson and binomial
#' @param y where y = Xβ + ε
#' @param X high dimensional matrix
#' @param trex_c TREX parameters, 2 for gaussian case, 1 for poisson and binomial case
#' @export
Gradient <- function(theta, family, y, X, trex_c) {
    .Call(`_myTridge_Gradient`, theta, family, y, X, trex_c)
}

#' GradientLs function
#' @param theta unknown regression vector
#' @param X high dimensional matrix
#' @param y where y = Xβ + ε
#' @param family three family cases: gaussian, poisson and binomial
#' @export
GradientLs <- function(theta, X, y, family) {
    .Call(`_myTridge_GradientLs`, theta, X, y, family)
}

#' find optimized solution for ObjectLs function
#' @param theta vector store our optimized result
#' @param X high dimensional matrix
#' @param y where y = Xβ + ε
#' @param family three family cases: gaussian, poisson and binomial
#' @export
optim_ObLs <- function(theta, X, y, family) {
    .Call(`_myTridge_optim_ObLs`, theta, X, y, family)
}

#' find optimized solution for ObjectRidge function
#' @param theta vector store our optimized result
#' @param X high dimensional matrix
#' @param y where y = Xβ + ε
#' @param family three family cases: gaussian, poisson and binomial
#' @param r tuning.parameter
#' @export
optim_Ridge <- function(theta, X, y, family, r) {
    .Call(`_myTridge_optim_Ridge`, theta, X, y, family, r)
}

#' find optimized solution for Objective function
#' @param theta vector store our optimized result
#' @param X high dimensional matrix
#' @param y where y = Xβ + ε
#' @param family three family cases: gaussian, poisson and binomial
#' @param trex_c TREX parameters, 2 for gaussian case, 1 for poisson and binomial case
#' @export
optim_ObFn <- function(theta, X, y, family, trex_c) {
    .Call(`_myTridge_optim_ObFn`, theta, X, y, family, trex_c)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_myTridge_RcppExport_registerCCallable', PACKAGE = 'myTridge')
})
