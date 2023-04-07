
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geessbin

<!-- badges: start -->
<!-- badges: end -->

The geessbin function provides results of modified generalized
estimating equations with bias-adjusted covariance estimators for
longitudinal or clustered data with binary outcomes.

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
res <- geessbin(formula = Wheeze ~ City + Time, id = ID,
                repeated = Time, corstr = "ar1", data = wheeze,
                beta.method = "PGEE", SE.method = "MB")

# hypothesis tests for regression coefficients
summary(res)
#> Call:
#> geessbin(formula = Wheeze ~ City + Time, data = wheeze, id = ID, 
#>     corstr = "ar1", repeated = Time, beta.method = "PGEE", SE.method = "MB")
#> 
#> Correlation Structure:  ar1 
#> Estimation Method for Regression Coefficients:  PGEE 
#> Estimation Method for Standard Errors:  MB 
#> 
#> Coefficients:
#>             Estimate Std.err      Z P.value
#> (Intercept)    1.097   3.333  0.329   0.756
#> City           0.228   0.750  0.304   0.762
#> Time          -0.190   0.307 -0.619   0.659
#> 
#> Odds Ratios with 95 % Confidence Intervals :
#>      Odds Ratio Lower Limit Upper Limit
#> City      1.256       0.289        5.46
#> Time      0.827       0.453        1.51
#> 
#> Estimated Scale Parameter:  1.01
#> Number of Iterations:  7 
#> 
#> Working Correlation:
#>         9   10   11     12
#> 9  1.0000 0.40 0.16 0.0642
#> 10 0.4003 1.00 0.40 0.1603
#> 11 0.1603 0.40 1.00 0.4003
#> 12 0.0642 0.16 0.40 1.0000
```
