library(tidyverse)
library(fitdistrplus)
library(rms)
library(MASS)

# See section 3.3 of https://rdrr.io/cran/fitdistrplus/f/inst/doc/paper2JSS.pdf

set.seed(223)


# make data
n<-100
grp<-sort(rep(c(0,1),n/2))
y<-rnorm(n,85)

dat0<-cbind(y,grp)

# detection limits
ll1 <- 85.5
ll2 <- 84.5

# scenario 1 - single detection limit

dat1<-dat0 %>% as_tibble() %>% 
  mutate(lim=ifelse(grp==0,ll2,ll2),left=pmax(lim,y),right=left) %>% 
  mutate(left=ifelse(lim==left,NA,left)) %>% as.data.frame()


# compare uncensored and censored outcome models
coef(lrm(y~1,data=dat1))
coef(lrm(right~1,data=dat1))

plot(ecdf(dat1[,c("y")]),verticals=TRUE,pch="",xlim=c(84.1,88.3))
abline(v=ll2)

plotdistcens(dat1[,c("left","right")])
abline(v=ll2)

#same after lower detection limit

plotdistcens(dat1[,c("left","right")], NPMLE=FALSE)


# scenario 2 - two detection limits
dat2<-dat0 %>% as_tibble() %>% 
  mutate(lim=ifelse(grp==0,ll1,ll2),left=pmax(lim,y),right=left) %>% 
  mutate(left=ifelse(lim==left,NA,left)) %>% as.data.frame()

plot(ecdf(dat2[,c("y")]),verticals=TRUE,pch="",xlim=c(84,88.3))
abline(v=c(ll1,ll2))

plotdistcens(dat2[,c("left","right")])
abline(v=c(ll1,ll2))

# same after max(detection limit)


plotdistcens(dat2[,c("left","right")], NPMLE=FALSE)

## add covariate
## need to update
if (0){

# make data
n<-200
x<-rbinom(n,1,0.4)
beta1 <- 2.5
#beta2 <- -1.5
grp<-sort(rep(c(0,1),n/2))
#y<-rnorm(n,85)+beta1*x+beta2*grp
y<-rnorm(n,85)+beta1*x

dat0<-cbind(y,x,grp)

# scenario 1 - single detection limit
ll1 <- 87.5
ll2 <- 84.5

dat1<-dat0 %>% as_tibble() %>% 
  mutate(lim=ifelse(grp==0,ll2,ll2),left=pmax(lim,y),right=left) %>% 
  mutate(left=ifelse(lim==left,NA,left)) %>% as.data.frame()


table(is.na(dat1[,"left"]))

# compare uncensored and censored outcome models
orm(y~x,data=dat1)

orm(right~x,data=dat1)

lrm(right~x,data=dat1)

polr(factor(y)~x,data=dat1)


#fitdist(y,"norm")
fit_x0<-fitdistcens(dat1[x==0,c("left","right")],"logis")
fit_x1<-fitdistcens(dat1[x==1,c("left","right")],"logis")
plot(fit_x0)
plot(fit_x1)

fit_x0
fit_x1

coef(fit_x0)

plotdistcens(dat1[,c("left","right")])
plotdistcens(dat1[,c("left","right")], NPMLE=FALSE)


# scenario 2 - two detection limits
dat2<-dat0 %>% as_tibble() %>% 
  mutate(lim=ifelse(grp==0,ll1,ll2),left=pmax(lim,y),right=left) %>% 
  mutate(left=ifelse(lim==left,NA,left)) %>% as.data.frame()

orm(y~x,data=dat2)

orm(right~x,data=dat2)

lrm(right~x,data=dat2)

# meta-analysis approach
orm(right~x,subset=c(grp==1),data=dat2)
orm(right~x,subset=c(grp==0),data=dat2)



#fitdist(y,"norm")
fit2_x0<-fitdistcens(dat2[x==0,c("left","right")],"logis")
fit2_x1<-fitdistcens(dat2[x==1,c("left","right")],"logis")
plot(fit2_x0)
plot(fit2_x1)

fit2_x0
fit2_x1

coef(fit2_x1) - coef(fit2_x0)

plotdistcens(dat1[,c("left","right")])
plotdistcens(dat1[,c("left","right")], NPMLE=FALSE)
}

