
# myTridge

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/mohan-zhao/myTridge.svg?branch=master)](https://travis-ci.com/mohan-zhao/myTridge)
<!-- badges: end -->

The goal of myTridge is to show that ridge esti- mators can be modified
such that tuning parameters can be avoided altogether and apply the
t-ridge estimator to generalized linear models.

## Installation

You can install the released version of myTridge from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("myTridge")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mohan-zhao/myTridge")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(myTridge)
## Here "gaussian" is from three most common cases for the distribution F: Gaussian,Poisson, and Bernoulli.
### high dimensional X :100x300
### k is the magnitude of the mutual correlations for covariance matrix 
TridgeEst(Test.case="gaussian",num.obs = 100,num.par=300,k=0.0)
#>   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators... done
#>   [1]  0.690035964 -0.735042942  0.297691693  0.443587630 -0.369141257
#>   [6] -0.737564216 -0.385692207  0.863998901 -0.245515549  0.151605346
#>  [11] -0.396629029 -0.220858742  0.457065712  1.124022019  0.834713349
#>  [16]  0.145411948 -0.274856132  0.380566804 -0.549235792  0.623252555
#>  [21] -0.822272743 -0.401258708 -0.324819608 -1.916548884 -0.646720531
#>  [26]  0.068863521  0.027989566  0.599081531 -0.917360441 -0.320391011
#>  [31]  0.293179417 -0.421439534  0.054240841 -0.520918007  0.074772726
#>  [36] -0.731551674  1.043859671  0.163425989 -0.666303700  0.129440122
#>  [41] -0.448898848 -0.108372621  0.086802146  1.171042031 -0.019753743
#>  [46] -0.167941323 -0.267249528  0.340197574 -0.123125994 -0.359988944
#>  [51]  0.125291073 -0.501605016  0.305606866  0.551898346 -0.469650316
#>  [56]  0.748042120 -0.816574048 -0.048428136  1.470729171  0.262930515
#>  [61]  0.716996430  0.902288931 -0.509366945 -0.258410986 -0.663663609
#>  [66]  0.378576124 -0.775629343 -0.073592841  0.551962567  0.966722602
#>  [71]  1.209542225  0.492101809  0.251395239 -0.532684370 -0.389499085
#>  [76] -0.058319453  0.508999134 -1.272418561  0.262332291 -0.015514266
#>  [81] -0.286169428  1.030500537 -0.806459388  0.335926988  0.096705875
#>  [86] -0.058610819  0.031031928  1.417746119  0.732245423 -0.042964990
#>  [91] -0.119828301 -0.096012327 -0.205273189 -0.275219135  0.940670682
#>  [96]  0.875547344  0.116303125 -0.341236652 -0.921388040  0.835479534
#> [101] -0.045008766  0.896353669 -0.844714592 -0.231741832  0.486760636
#> [106] -0.258153891  0.232841841 -2.090712544  0.143864903  0.770356103
#> [111] -0.645100921 -0.617638501 -0.982585539 -0.149323600 -0.024007926
#> [116] -0.571309057  0.115164731 -1.272950846 -0.872621407 -0.988665058
#> [121]  0.019513333 -0.213759578 -0.576914635  0.282064260  0.129507323
#> [126]  1.191895138 -0.413436151  0.055331883 -0.532812189 -0.547459164
#> [131]  0.639155548 -0.153101299 -0.392734655  0.109563029 -0.229778834
#> [136]  0.002047006  0.934238261  0.207373229  0.616869173 -0.441762354
#> [141] -1.240922010 -0.187696495  0.203428956 -0.993104198  0.210366748
#> [146]  0.256915994  0.238013405  0.866505884  0.974332108  0.194848894
#> [151]  0.770381054  0.352014835  0.025620949 -1.374250834 -0.761639295
#> [156] -0.037417006  0.197315968  1.273662628 -0.159969498 -0.360174949
#> [161]  0.454746032  0.416630800 -0.359303964 -0.249542999  0.799230740
#> [166]  0.359239691  0.172938850 -0.656151408 -0.825030486 -0.024782614
#> [171]  0.348719158 -0.566793653  1.234402219  0.415830461  0.019219896
#> [176]  0.145809253  1.016200745  0.293256608  0.884274974  0.125157356
#> [181]  0.031966378 -0.097090638  0.031177652  0.741369868 -0.927277295
#> [186]  0.296712353  0.080396558 -0.373639466 -0.201248973 -0.316498701
#> [191] -0.938481751 -0.403320038  0.130907511  0.354099346 -0.230049409
#> [196] -0.325322848  0.510551010 -0.227214253  0.924806980  1.199732209
#> [201]  0.568346638  0.647967153  0.765285749  0.433252028 -0.188265274
#> [206]  0.240638591  0.285588726 -0.359813917 -0.039668233  0.524655053
#> [211] -0.992759331 -0.589826441  0.588411403  0.333186307  0.733573039
#> [216]  0.507615105  0.849323660 -0.210402300  0.529367700 -0.190522341
#> [221] -0.162707483  0.402014305  0.147643334 -0.004038497  0.434498164
#> [226] -0.209372309 -0.476974792 -0.229640652  0.172348556 -0.038057492
#> [231] -0.021157027  0.346603396  0.373778241 -0.259015309 -0.243731564
#> [236] -0.277458173 -0.221486254  0.425700328  0.967402514 -1.160562424
#> [241]  0.950125111 -1.148529160 -1.663338995 -0.202800344 -1.572679505
#> [246]  0.410512125 -0.571055182 -0.314578797  0.060559031 -0.567507634
#> [251] -1.295907001 -0.027757150  1.510706070 -0.805479895  0.337827958
#> [256]  1.030865713  0.367478899  0.906752814 -0.061927860 -0.513815176
#> [261]  0.466485574 -0.348973652 -0.403340422  0.365841349  0.344512961
#> [266]  0.415893282 -0.467959698  0.397353155 -0.748982168  0.026777250
#> [271]  0.623736991 -0.021729907  1.019829499  0.675806180 -0.611207080
#> [276]  0.552254433  0.126695635  0.886452243  0.121170734 -0.814712569
#> [281] -1.338359377  0.519995590 -0.328600811 -0.778198886 -0.296592944
#> [286]  0.106375214 -0.576027164  0.297774935  0.013807069  1.263915834
#> [291] -0.161067293 -0.331291055 -0.029690840  0.546050630  0.777839687
#> [296]  0.217812397  0.458988265  0.028600146 -0.285071058  0.032160500
```

``` r
## example of relativeError() function
###shows t-ridge outperforms K-fold cross-validated ridge
relativeError("gaussian",100,300,0.0,5)
#> run 1 out of 5 runs
#>   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators...   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators... done
#> done
#>   -> compute errors... done
#> run 2 out of 5 runs
#>   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators...   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators... done
#> done
#>   -> compute errors... done
#> run 3 out of 5 runs
#>   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators...   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators... done
#> done
#>   -> compute errors... done
#> run 4 out of 5 runs
#>   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators...   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators... done
#> done
#>   -> compute errors... done
#> run 5 out of 5 runs
#>   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators...   -> generate test case...   -> Normalize each column of the design matrix...   -> Perform the K-fold cross validation pipeline...   -> compute T-ridge estimators... done
#> done
#>   -> compute errors... done
#> Simulation Results :
```

<table class="gmisc_table" style="border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;">

<thead>

<tr>

<td colspan="6" style="text-align: left;">

(n,p,k)=(100,300,0)

</td>

</tr>

<tr>

<th style="border-top: 2px solid grey;">

</th>

<th colspan="2" style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">

\(\frac{||X\hat{\beta}_{T-ridge} - X\beta^{*}||_{2}}{||X\beta^{*}||_{2}}\)

</th>

<th style="border-top: 2px solid grey;; border-bottom: hidden;">

 

</th>

<th colspan="2" style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">

\(\frac{||\hat{\beta}_{T-ridge} - \beta^{*}||_{2}}{||\beta^{*}||_{2}}\)

</th>

</tr>

<tr>

<th style="border-bottom: 1px solid grey;">

</th>

<th style="font-weight: 900; border-bottom: 1px solid grey; text-align: center;">

T-ridge

</th>

<th style="font-weight: 900; border-bottom: 1px solid grey; text-align: center;">

10-fold
CV

</th>

<th style="font-weight: 900; border-bottom: 1px solid grey; text-align: center;" colspan="1">

 

</th>

<th style="font-weight: 900; border-bottom: 1px solid grey; text-align: center;">

T-ridge

</th>

<th style="font-weight: 900; border-bottom: 1px solid grey; text-align: center;">

10-fold CV

</th>

</tr>

</thead>

<tbody>

<tr>

<td colspan="6" style="font-weight: 900;">

Case : gaussian

</td>

</tr>

<tr>

<td style="border-bottom: 2px solid grey; text-align: left;">

  Mean relative
errors

</td>

<td style="border-bottom: 2px solid grey; text-align: center;">

2.97

</td>

<td style="border-bottom: 2px solid grey; text-align: center;">

0.52

</td>

<td style="border-bottom: 2px solid grey; text-align: center;" colspan="1">

 

</td>

<td style="border-bottom: 2px solid grey; text-align: center;">

4.88

</td>

<td style="border-bottom: 2px solid grey; text-align: center;">

0.6

</td>

</tr>

</tbody>

</table>

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.
