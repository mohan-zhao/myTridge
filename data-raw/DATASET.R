#' R Source code file for creating simulated dataset to be included in the myTridge package
ret_gau<-genDataList(10, rep(0, 30), 30, 0.,
            rnorm(30, mean = 0, sd = 1), 10,"gaussian")
ret_bin<-genDataList(10, rep(0, 30), 30, 0.,
                     rnorm(30, mean = 0, sd = 1), NA,"binomial")
ret_poi<-genDataList(10, rep(0, 30), 30, 0.,
                     rnorm(30, mean = 0, sd = 1), NA,"poisson")
X_gau<-ret_gau$normData
y_gau<-ret_gau$y
Test.case_gau<-ret_gau$Test.case

X_poi<-ret_poi$normData
y_poi<-ret_poi$y
Test.case_poi<-ret_poi$Test.case

X_bin<-ret_bin$normData
y_bin<-ret_bin$y
Test.case_bin<-ret_bin$Test.case
usethis::use_data(X_gau,y_gau,Test.case_gau,X_bin,y_bin,Test.case_bin,X_poi,y_poi,Test.case_poi,overwrite = TRUE)

