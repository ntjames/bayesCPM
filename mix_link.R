## mixture link from Lang 1999

#### in base R

library(extRemes) # for extreme value dists
# https://cran.r-project.org/web/packages/extRemes/extRemes.pdf

# logistic
x<-seq(-4,4,length=75)
logis_y<-plogis(x)
plot(x,logis_y,"l")

# extreme value maximum (Gumbel maximum)
hist(revd(1000))
evmax_y<-pevd(x)
plot(x,evmax_y,"l")

# extreme value minimum
hist(-revd(1000))
pevd_min<-function(x,...){
  return(pevd(-x,lower.tail=FALSE,...))
}
evmin_y<-pevd_min(x)
plot(x,evmin_y,"l")


# mixture cdf
pmix <- function(q, lambda){
  m1 <- exp(-exp(3.5*lambda+2))
  m3 <- exp(-exp(-3.5*lambda+2))
  m2 <- 1-m1-m3
  
return( m1*pevd_min(q)+m2*plogis(q)+m3*pevd(q) )
}


pmix_log <- function(q, lambda){
  log_m1 <- -exp(3.5*lambda+2)
  log_m3 <- -exp(-3.5*lambda+2)
  log_m2 <- log1p(-(exp(log_m1)+exp(log_m3)))
  
  logA <- log_m1 + pevd_min(q,log.p=TRUE)
  logB <- log_m2 + plogis(q,log.p=TRUE)
  logC <- log_m3 + pevd(q,log.p=TRUE)
  
  #log(A+B+C)=log(A)+log(1+exp(log(B)-log(A))+exp(log(C)-log(A)))
  # https://en.wikipedia.org/wiki/List_of_logarithmic_identities
  return( logA+log1p(exp(logB-logA)+exp(logC-logA)) )
}

log(pmix(0.235,-1))
pmix_log(0.235,-1)


#mixture inverse cdf; smallest x s.t. F(x)<p
qmix0 <- function(p,lambda,lim=c(-10000,10000)){
  f <- function(y) pmix(y,lambda)-p
  return( uniroot(f,lower=lim[1],upper=lim[2])$root )
}

qmix<-function(pv,...){
  sapply(pv, function(p) qmix0(p,...) )
}

qmix(0.3,-1)
qmix(c(0.3,0.4),-1)
qmix(c(2.4e-135,0.00001,0.9998),-1)

# fix this issue
qmix(0.99999999999999994,-1)
qmix(0.99999999999999995,-1)


# grid approx
pseq<-seq(1e-16,1-1e-16,length=1e4)
qseq<-qmix(pseq,-1)
max(qseq[which(pseq<0.61)])
max(qseq[which(pseq<0.523)])

# exact
qmix(0.61,-1)
qmix(0.523,-1)


mix_y_neg2<-pmix(x,-2)
plot(x,mix_y_neg2,"l")
abline(h=0.223)
abline(v=qmix(0.223,-2))

mix_y_neg1<-pmix(x,-1)
plot(x,mix_y_neg1,"l")
abline(h=0.6)
abline(v=qmix(0.6,-1))

mix_y_neg05<-pmix(x,-0.5)
plot(x,mix_y_0,"l")
abline(h=0.6)
abline(v=qmix(0.6,-0.5))

mix_y_0<-pmix(x,0)
plot(x,mix_y_0,"l")

mix_y_1<-pmix(x,1)
plot(x,mix_y_1,"l")

mix_y_2<-pmix(x,2)
plot(x,mix_y_2,"l")

# can also do mixtures with 'mistr' package
# see https://cran.r-project.org/web/packages/mistr/vignettes/mistr-introduction.pdf
# library(mistr)


#### in Stan
libs <- c("rstan", "rstanarm", "rms", "dplyr", "stringr", "readr")
invisible(lapply(libs, library, character.only = TRUE))

dir<-getwd()
source(file.path(dir,"cpm_functions.r"))

# make data
set.seed(930)
n <- 150
x1 <- c(rep(0,75), rep(1,75))
y <- rnorm(n) + 3*x1
dat1 <- data.frame(y=ordered(y),x1)

mod_data_1 <- mkStanDat(dat1, outcome="y", preds = c("x1"), link=2)


# individual links
ord_mod_file1 <- read_file(file.path(getwd(),"ordinal_model_1.stan"))
ord_mod1 <- stan_model(model_code = ord_mod_file1)

# call this once to distribute MCMC chains across cpu cores:
options(mc.cores=parallel::detectCores())

if(0){
#probit link 
fit_ord_1 <- sampling(ord_mod1, data=mod_data_1, seed=6472, 
               iter=4000, warmup=2000, chains=2,
               control = list(adapt_delta = 0.95))

summary(fit_ord_1,pars=c("b[1]"))$summary %>% round(3)

plot(fit_ord_1,pars=c("b[1]"))
traceplot(fit_ord_1,pars=c("b[1]"))
traceplot(fit_ord_1,pars=c("pi[1]"))
traceplot(fit_ord_1,pars=c("cutpoints[15]"))
}


# mixture link

# mixture link model
# NOTE: currently VERY slow (hours!!) because need to numerically solve for
# inverse mixture CDF
# also still have convergence issues
ord_mod_file1b <- read_file(file.path(getwd(),"ordinal_model_1b.stan"))
ord_mod1b <- stan_model(model_code = ord_mod_file1b)

# NOTE: link specification in data doesn't do anything
ptm <- proc.time()
fit_ord_1b <- sampling(ord_mod1b, data=mod_data_1, seed=7723, 
                     iter=6000, warmup=4000, chains=2,
                     control = list(adapt_delta = 0.9))
proc.time() - ptm

summary(fit_ord_1b,pars=c("b[1]"))$summary %>% round(3)
summary(fit_ord_1b,pars=c("lambda"))$summary %>% round(3)

plot(fit_ord_1b,pars=c("b[1]"))
traceplot(fit_ord_1b,pars=c("b[1]"))

traceplot(fit_ord_1b,pars=c("pi[2]"))
traceplot(fit_ord_1b,pars=c("cutpoints[15]"))
traceplot(fit_ord_1b,pars=c("lambda"))






if(0){
#scratch

link_choice <- function(lnk,x){
switch(lnk,
       pevd_min(x),
       plogis(x),
       pevd(x))
}

# qmix <- function(p,xseq,lambda,lim=c(-1000,1000)){
#   f <- function(xseq) pmix(xseq,lambda)-p
# return( uniroot(f,lower=lim[1],upper=lim[2])$root )
# }

# p<-0.2
# f <- function(x) pmix(x,0)-p
# uniroot(f,lower=-100,upper=100)
}
