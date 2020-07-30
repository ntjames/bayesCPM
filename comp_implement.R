#### Compare CPM ordinal models ####

# call this once to distribute MCMC chains across cpu cores:
options(mc.cores=parallel::detectCores())

## make example data
set.seed(1567)
n <- 50
x1 <- rnorm(n)
y <- 0.9*x1 + rlogis(n)
dat <- data.frame(y=y,x1)
dat_ord <- data.frame(y=ordered(y),x1)

### bayesCPM 
library(bayesCPM)
library(rstan)
library(dplyr) # needed for bayesCPM::mkStanDat()

## convert data to Stan format (link=1 is logistic, link=2 is probit)
ptm <- proc.time()
dat_stan <- mkStanDat(dat_ord, outcome="y", preds = c("x1"), link=1)

## sample from Bayes CPM model with probit link
bayescpm1 <- bayes_cpm(dat_stan)
proc.time() - ptm
# user  system elapsed 
# 10.588   0.396   4.381 

traceplot(bayescpm1, pars=c("b[1]"))
stan_plot(bayescpm1, pars=c("b[1]"))

stan_diag(bayescpm1)

summary(bayescpm1, pars=c("b[1]"),use_cache=FALSE)$summary %>%
  as_tibble() %>% mutate(par=c("b[1]"))

fit_mn <- getMean(bayescpm1, dat_stan, newdata=data.frame(x1=c(1)))
fit_mn



### bayesCPM w/ QR decomp?


### goodrich w/o QR decomp
# for goodrich mods see https://github.com/bgoodri/covid19_Bayesian_RCT
## ord_mod_goodrich.stan, same as ordered_logit.stan from above

# goodrich w/ QR decomp (blrm?)
## lrmqr.stan from above
library(rmsb)

rm(ptm)

ptm <- proc.time()
blrm1 <- blrm(y~x1)
proc.time() - ptm

# user  system elapsed 
# 14.408   0.748   5.927


blrm1
plot(blrm1)
m<-Mean(blrm1)
m(predict(blrm1,newdata=data.frame(x1=c(1))))

check_divergences(blrm1$rstan)

sbc

### betancourt w/o QR decomp
## ord_mod_betancourt.stan


### betancourt w/ QR decomp


### OpenBUGS - see /home/nathan/Dropbox/njames/stat/stat_book_exercises/congdon_bayes_for_cat_data

