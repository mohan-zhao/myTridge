
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
devtools::load_all()
#> Loading myTridge
## Here "gaussian" is from three most common cases for the distribution F: Gaussian,Poisson, and Binomial.
### X and y are matrix and vector under gaussian distribution in data folder
### dimension of X :10x30
### generate a tridge estimator
TridgeEst(Test.case="gaussian",X=X_gau,y=y_gau)
#>  [1]  0.1670192214  0.3444109960 -0.2455607011  0.5691582110  0.0759831119
#>  [6] -0.0008287938 -0.5805231789 -0.6276065836 -0.5171332311  0.1106345505
#> [11] -0.3879062193 -0.7348493696 -0.3388737660  0.5701026405 -0.3441166568
#> [16]  0.2201076189 -0.3350770107 -0.0272578159  0.1429890980  0.1534147531
#> [21] -0.0970010686  0.1400095960  0.4609051226  0.6895894773 -0.3621226935
#> [26]  0.3114799498 -0.3796407771 -0.0615521447 -0.9272345753  0.0806730175
```

``` r

### You can also use genDataList() function in this package to generate a high dimensional matrix X and y randomly then generate a tridge estimator
list<-genDataList(n=100,mu=rep(0,300),p=300, rho=0.,beta=rnorm(300, mean = 0, sd = 1),SNR=10,Test_case="gaussian")
X<-list$normData
y<-list$y
myTridge::TridgeEst(Test.case="gaussian",X=X,y=y)
#>   [1]  0.576012722 -0.744677982 -0.003237774  0.314723457  0.486794088
#>   [6]  0.533943504  0.422719077  0.695709574 -0.612922014 -0.815627168
#>  [11]  0.104739312  0.314199629 -0.329489190  0.325055790 -1.429502718
#>  [16] -0.177188198  0.330998573  0.206948107 -0.219405370  0.265212940
#>  [21] -0.309497061  0.876825309  0.190012681 -0.489564503  0.216992114
#>  [26]  0.620985763  0.747688002  0.387654298 -1.032794253 -0.614994118
#>  [31]  0.534199036 -0.548567308 -1.142429346  0.300904887 -1.061461837
#>  [36]  0.186355902  0.796830821  1.596525250  0.783028626 -0.001740044
#>  [41] -0.053303209  0.020430387  0.144793525  0.079529394  0.088629693
#>  [46]  0.009736881  0.021114218 -0.216470979 -0.730928937 -0.080116591
#>  [51] -1.068705586 -1.312076834  0.209476797 -0.439167733  0.007805762
#>  [56] -0.015294628 -0.911819550 -0.168596459  0.251996208  0.496830102
#>  [61]  0.944055773 -0.429233556  0.403977135 -0.756034064 -0.138673969
#>  [66]  0.952799324 -0.621176717  0.548180036 -0.604107926  0.611087969
#>  [71] -0.254998860 -0.416232896 -0.228883449 -1.143219770 -0.432702518
#>  [76] -0.755050439 -0.754033204  0.380784128  0.081175733 -0.314234774
#>  [81] -1.726630086 -0.223869161  0.481001177 -1.304199985 -0.528314844
#>  [86]  0.256588940  0.845633790  1.083770279 -0.205115412 -0.514725241
#>  [91] -1.375075357  0.672086019  0.303919936  0.316265532 -0.886621630
#>  [96]  0.184174082  0.218639728  0.971098954 -0.006169950  0.538154562
#> [101]  1.170427628 -0.426632522  0.740211833  0.399286912  0.164004948
#> [106]  0.289730808  0.276239981 -0.327793178 -0.030021769  0.363196251
#> [111] -0.409804203 -0.206649509 -0.415812404 -0.139690964 -0.527531043
#> [116]  0.162223334 -0.711432933  0.008568535  1.290440231  0.223765244
#> [121]  0.061077078 -0.131400347 -0.121610015 -0.473245491 -0.398074668
#> [126]  0.248255246 -0.308682685  0.586146188 -0.401160771  0.154957823
#> [131] -0.678766391 -0.520585145  0.359913794 -0.477614017  0.361207562
#> [136] -0.187228256 -0.543993322  0.078082054  0.684090869  0.527885438
#> [141]  1.554904002 -0.172040679 -0.980922272 -0.317298182 -0.086237010
#> [146]  0.167496379  0.211292035  0.044855279 -0.397480626  0.981848048
#> [151] -0.068418930 -0.708156788  0.368691100 -0.366791561 -0.558349823
#> [156] -0.564283093  0.219914037 -0.954450231 -0.567185257  0.257229024
#> [161] -0.799251895 -1.044031505  0.257458675 -0.507092916 -0.516180707
#> [166]  1.426569219 -0.432416876 -0.543623984 -0.997331913  0.156743034
#> [171] -0.247722556 -0.401074433 -0.061858437  0.022701439 -0.476402272
#> [176]  1.287921961  0.194267175  0.139509355  0.256876421  0.434872838
#> [181]  0.156279385  0.940569629 -0.198296737 -0.424071200  0.188822840
#> [186] -0.142691232  0.350679993  0.627540077  0.137306455  0.012844045
#> [191] -0.068795608  0.206269747 -0.650628008 -0.614788469  0.781084744
#> [196]  0.017874476  0.482126253  1.650049608  0.039368966  0.620670608
#> [201]  0.385601864  0.495323631 -0.353899103 -0.395054397 -0.326535235
#> [206] -0.532130525 -0.526594183  0.024960554 -0.441802015 -0.306437759
#> [211] -0.351423977 -0.324682051 -1.172589439  0.691002706 -0.563554682
#> [216]  0.536676962  0.483106709  0.294192995 -0.038881263 -0.661874040
#> [221]  0.071993526 -0.212988650 -0.679634438  1.059574708 -0.392719270
#> [226]  0.760987429 -0.194760459 -0.459523062  0.221099462 -0.794682562
#> [231] -0.020276859 -0.558501839 -0.221578063  0.358186899 -0.435474184
#> [236]  0.561082406  1.271356124  0.791703213 -0.947768116 -0.211621608
#> [241]  0.060550830 -0.792272378  1.007582411  0.599979708  0.067444698
#> [246] -0.247349549 -0.438481747 -0.309550695  0.352479724  0.047080985
#> [251] -0.377153957  0.471305900  0.277745680 -0.339552065  0.219106744
#> [256] -0.663353156 -0.298303353  0.611630466 -0.384929094 -0.134089841
#> [261]  0.176041485 -0.570232350  0.067568023 -0.202293271  0.416515979
#> [266]  1.084802093 -0.425791621 -0.942891311 -0.639342305 -0.271459839
#> [271]  0.241625222  0.679036717  0.046645273 -0.046773930  0.911983883
#> [276] -0.858137667  0.210624778  0.316007129  0.834621777 -0.445543069
#> [281]  0.709130602  0.175083997 -0.878371080 -0.035340068 -0.441733379
#> [286] -0.173560071  0.839937948 -0.271852664 -0.122126964  0.333531561
#> [291]  0.401923837  0.504334541 -0.155621707  0.399904414  0.342230510
#> [296]  0.112696208  0.777230567 -0.185138584 -0.028697167 -0.819482326
```

    #> Simulation Results :

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

0.31

</td>

<td style="border-bottom: 2px solid grey; text-align: center;">

0.52

</td>

<td style="border-bottom: 2px solid grey; text-align: center;" colspan="1">

 

</td>

<td style="border-bottom: 2px solid grey; text-align: center;">

0.37

</td>

<td style="border-bottom: 2px solid grey; text-align: center;">

0.6

</td>

</tr>

</tbody>

</table>

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.
