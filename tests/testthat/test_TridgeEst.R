context("run TridgeEst with example input ")


estimator <- try(TridgeEst("gaussian",100,300,0.0),
           silent = TRUE)

test_that("no error regarding the output type in TridgeEst for the input", {
  
  expect_is(estimator, "numeric")
  
})
