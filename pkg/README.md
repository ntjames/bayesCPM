
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesCPM

<!-- badges: start -->

<!-- badges: end -->

This package includes functions to fit a Bayesian Cumulative Probability
Model (CPM) using the R interface to Stan.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ntjames/bayes_cpm/pkg")
```

## Example

This is a basic example:

``` r
library(bayesCPM)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

## make example data
set.seed(1567)
n <- 100
x1 <- rnorm(n)
y <- 0.9*x1 + rnorm(n)
dat <- data.frame(y=ordered(y),x1) # outcome must be ordered factor

## sample from Bayes CPM model with probit link
fit <- bayes_cpm(y~x1, data=dat, link="probit", refresh=1000)
#> 
#> SAMPLING FOR MODEL 'bayes_cpm_mod' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.00018 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.8 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.99746 seconds (Warm-up)
#> Chain 1:                1.51528 seconds (Sampling)
#> Chain 1:                3.51274 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bayes_cpm_mod' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 9.4e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.94 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 2.03663 seconds (Warm-up)
#> Chain 2:                1.7976 seconds (Sampling)
#> Chain 2:                3.83422 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'bayes_cpm_mod' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000187 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 1.87 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 2.09917 seconds (Warm-up)
#> Chain 3:                1.56386 seconds (Sampling)
#> Chain 3:                3.66303 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'bayes_cpm_mod' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 9.5e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.95 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 2.13702 seconds (Warm-up)
#> Chain 4:                1.64027 seconds (Sampling)
#> Chain 4:                3.77729 seconds (Total)
#> Chain 4:
```

Using the model fit we can get a summary of the posterior conditional
mean,

``` r
## posterior conditional mean when x1=1
fit_mn <- getMean(fit, newdata=data.frame(x1=c(1)))
fit_mn
#> # A tibble: 1 x 13
#> # Groups:   ndrow [1]
#>   ndrow mean_mn med_mn sd_mn mn_q2.5 mn_q5 mn_q10 mn_q25 mn_q75 mn_q90 mn_q95
#>   <int>   <dbl>  <dbl> <dbl>   <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#> 1     1   0.970  0.974 0.136   0.695 0.734  0.792  0.881   1.06   1.14   1.19
#> # … with 2 more variables: mn_q97.5 <dbl>, x1 <dbl>
```

a posterior conditional quantile,

``` r
## posterior conditional 80th quantile when x1=0.5
fit_q80 <- getQuantile(fit, newdata=data.frame(x1=c(0.5)),q=0.8)
fit_q80
#> # A tibble: 1 x 13
#> # Groups:   ndrow [1]
#>   ndrow mean_qtile med_qtile sd_qtile qtile_q2.5 qtile_q5 qtile_q10 qtile_q25
#>   <int>      <dbl>     <dbl>    <dbl>      <dbl>    <dbl>     <dbl>     <dbl>
#> 1     1       1.47      1.41    0.149       1.26     1.30      1.33      1.36
#> # … with 5 more variables: qtile_q75 <dbl>, qtile_q90 <dbl>, qtile_q95 <dbl>,
#> #   qtile_q97.5 <dbl>, x1 <dbl>
```

or plot the median of the posterior conditional CDF. The true CDF is
shown for reference.

``` r
## get posterior conditional CDF when x1=1
fit_cdf <- getCDF(fit, newdata=data.frame(x1=c(1))) 

plot(fit_cdf$yval, fit_cdf$med_cdf, type="s", xlab="y", ylab="Conditional CDF")
curve(pnorm(x,mean=0.9*1),-3,3, add=TRUE, col="red")
legend("topleft",legend=c("Bayes CPM estimate","True CDF"),lty=c(1,1),col=c(1,2))
```

<img src="man/figures/README-cdf-1.png" width="100%" />

<!--
You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!
-->
