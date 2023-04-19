
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geessbin

<!-- badges: start -->
<!-- badges: end -->

The geessbin function analyzes small-sample clustered or longitudinal
data using modified generalized estimating equations (GEE) with
bias-adjusted covariance estimator. This function assumes binary outcome
and uses the logit link function. This function provides any combination
of three GEE methods (conventional and two modified GEE methods) and 12
covariance estimators (unadjusted and 11 bias-adjusted estimators).

## Installation

You can install the released version of geessbin from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("geessbin")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rtishii/geessbin")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(geessbin)
data(wheeze)

# analysis of PGEE method with Morel et al. covariance estimator
res <- geessbin(formula = Wheeze ~ City + factor(Age), data = wheeze, id = ID,
                corstr = "ar1", repeated = Age, beta.method = "PGEE",
                SE.method = "MB")

print(res)
#> Call:
#> geessbin(formula = Wheeze ~ City + factor(Age), data = wheeze, 
#>     id = ID, corstr = "ar1", repeated = Age, beta.method = "PGEE", 
#>     SE.method = "MB")
#> 
#> Correlation Structure:  ar1 
#> Estimation Method for Regression Coefficients:  PGEE 
#> Estimation Method for Standard Errors:  MB 
#> 
#> Number of observations:  64 
#> Number of clusters:  16 
#> Maximum cluster size:  4 
#> 
#> Coefficients:
#>   (Intercept)          City factor(Age)10 factor(Age)11 factor(Age)12 
#>        -0.546         0.226        -0.237        -0.511        -0.525 
#> 
#> Estimated Scale Parameter:  1.04
#> Number of Iterations:  8 
#> 
#> Working Correlation:
#>         9    10    11     12
#> 9  1.0000 0.411 0.169 0.0693
#> 10 0.4107 1.000 0.411 0.1687
#> 11 0.1687 0.411 1.000 0.4107
#> 12 0.0693 0.169 0.411 1.0000
#> 
#> Convergence:  1 ( Success )

# hypothesis tests for regression coefficients
summary(res)
#> Call:
#> geessbin(formula = Wheeze ~ City + factor(Age), data = wheeze, 
#>     id = ID, corstr = "ar1", repeated = Age, beta.method = "PGEE", 
#>     SE.method = "MB")
#> 
#> Correlation Structure:  ar1 
#> Estimation Method for Regression Coefficients:  PGEE 
#> Estimation Method for Standard Errors:  MB 
#> 
#> Coefficients:
#>               Estimate Std.err      Z P.value
#> (Intercept)     -0.546   0.782 -0.699   0.625
#> City             0.226   0.827  0.273   0.769
#> factor(Age)10   -0.237   0.752 -0.315   0.759
#> factor(Age)11   -0.511   0.870 -0.588   0.671
#> factor(Age)12   -0.525   0.986 -0.532   0.692
#> 
#> Odds Ratios with 95 % Confidence Intervals :
#>               Odds Ratio Lower Limit Upper Limit
#> City               1.254      0.2477        6.35
#> factor(Age)10      0.789      0.1809        3.44
#> factor(Age)11      0.600      0.1090        3.30
#> factor(Age)12      0.592      0.0857        4.09
#> 
#> Estimated Scale Parameter:  1.04
#> Number of Iterations:  8 
#> 
#> Working Correlation:
#>         9    10    11     12
#> 9  1.0000 0.411 0.169 0.0693
#> 10 0.4107 1.000 0.411 0.1687
#> 11 0.1687 0.411 1.000 0.4107
#> 12 0.0693 0.169 0.411 1.0000
```
