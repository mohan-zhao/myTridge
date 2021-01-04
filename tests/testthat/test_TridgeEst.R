context("run TridgeEst with example input ")
### We use matrix X and vector y in data/ to test three cases following:
###"gaussian","poisson"and "binomial"
### X_gau and y_gau are data for "gaussian"
### X_bin and y_bin are data for "binomial"
### X_poi and y_poi are data for "poisson"
estimator1 <- try(TridgeEst(X_gau,y_gau,"gaussian"),
           silent = TRUE)
estimator2 <- try(TridgeEst(X_poi,y_poi,"poisson"),
                  silent = TRUE)
estimator3 <- try(TridgeEst(X_bin,y_bin,"binomial"),
                  silent = TRUE)

test_that("no error regarding the output type in TridgeEst for the input", {
  
  expect_is(estimator1, "numeric")
  
})

test_that("no error regarding the poisson distribution input type in TridgeEst", {
  
  expect_is(estimator2, "numeric")
  
})

test_that("no error regarding the bernoulli distribution input type in TridgeEst", {
  
  expect_is(estimator3, "numeric")
  
})