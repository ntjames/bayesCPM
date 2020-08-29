rm(list=ls())

# code for partially specified Polya tree
# see Hanson et al. (2005) Bayesian Nonparametric Modeling & Data Analysis: An Introduction
# in Handbook of Statistics, Vol. 25, pp. 255-257
# Muller et al. (2015) Bayesian Nonparametric Data Analysis, Chp 3 

# function to get conditional Y probs at level m using 
# Beta(alp,alp) where alp=gam*m^2
condY<-function(gam,m){
  alp<-gam*m^2
  y<-rbeta(2^(m-1),alp,alp)
  c(rbind(y,1-y)) 
}

# partition probabilities 
parprob<-function(gam,m){
  cprobs<-matrix(NA,nrow=m,ncol=2^m)
  for (i in 1:m){
    cprobs[i,]<-rep(condY(gam=gam,i),each=2^(m-i))
  }
  apply(cprobs,2,prod)
}

# partition intervals for centering distribution F0
# F0 is inverse cdf function of centering dist
# i is index of partition at level m
# m is partition level
# ... is params for dist
PI<-function(F0,i,m,...){
  sort(c(F0((i-1)/(2^m),...), F0(i/(2^m),...)))
}


# midpoints of cut() function
# from https://www.r-bloggers.com/finding-the-midpoint-when-creating-intervals/
midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
  upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
  return(round(lower+(upper-lower)/2, dp))
}


if(0){
c(PI(qnorm,1,1),PI(qnorm,2,1)) # level m=1 partition
c(PI(qnorm,1,2),PI(qnorm,2,2),PI(qnorm,3,2),PI(qnorm,4,2)) # level m=2 partition
  
PI(qnorm,1:2^1,1) # level m=1 partition
PI(qnorm,1:2^2,2) # level m=2 partition
PI(qnorm,1:2^3,3) # level m=3 partition
}


#F0 inverse cdf function of centering dist
# gg gamma for calculating Beta/Dirichlet params - current uses gam*m^2 for partition level m
# mm number of levels for partition (default is mm=10 for partition of 2^10=1024 at lowest lev)
# br number of breaks to summarize draw
# ... params for centering dist
draw_Polya<-function(F0,gg=5,mm=10,br=50,...){
  br0<-unique(PI(F0=F0,i=1:2^mm,m=mm,...))
  ypr<-parprob(gam=gg,m=mm)
  mrg0<-cbind(i0=br0[1:(2^mm)],i1=br0[2:(2^mm+1)],ypr)
  mrg<-data.frame(mrg0[-c(1,2^mm),]) # trim intervals at ends w/ -Inf and Inf

  ctlev<-cut(seq(mrg[1,"i0"],mrg[2^(mm)-2,"i1"],length=br),br)
  lvs0<-paste0(gsub("\\(||\\]","",levels(ctlev)),sep=",",collapse="")
  lvs1<-substr(lvs0,1,nchar(lvs0)-1)
  lvs<-unique(as.numeric(unlist(strsplit(lvs1, split=","))))
  
  ct <- factor(findInterval(mrg$i1,lvs),levels=1:br)
  prplt<-data.frame(cbind(mp=unique(midpoints(ctlev)),pr=by(mrg$ypr,ct,sum)))
  
  prplt[is.na(prplt[,"pr"]),"pr"]<-0
  
 prplt
}


# example std. Normal centering dist
nn<-draw_Polya(F0=qnorm,gg=5,mm=10,br=100)
plot(nn$mp,nn$pr,type="l")

# example exp(1) centering dist
ee<-draw_Polya(F0=qexp,gg=25,mm=10,br=120)
plot(ee$mp,ee$pr,type="l")

# example gamma(9,2) centering dist
gg<-draw_Polya(F0=qgamma,gg=50,mm=10,br=80,shape=9,rate=2)
plot(gg$mp,gg$pr,type="l")


# reproduce Fig 3.1 from Muller (2015)
library(dplyr)
library(magrittr)
library(ggplot2)

set.seed(2673)
out_c1000 <- out_c100 <- out_c5 <- data.frame()
for (n in 1:10){
  c5<-draw_Polya(F0=qnorm,gg=5,mm=10,br=100)
  c100<-draw_Polya(F0=qnorm,gg=100,mm=10,br=100)  
  c1000<-draw_Polya(F0=qnorm,gg=1000,mm=10,br=100) 
  c1000$draw<-c100$draw<-c5$draw<-n
  out_c5<-rbind(out_c5,c5)
  out_c100<-rbind(out_c100,c100)
  out_c1000<-rbind(out_c1000,c1000)
}

out_c5 %<>% mutate(param=5)
out_c100 %<>% mutate(param=100)
out_c1000 %<>% mutate(param=1000)
fctlabs<-c("c=5","c=100","c=1000")

pdraws <-bind_rows(out_c5,out_c100,out_c1000) %>% 
  mutate(param=factor(param,labels=fctlabs))

pdrawsmn <- pdraws %>% group_by(param,mp) %>% summarize(pr=mean(pr)) %>% mutate(draw=-1)

scale<-16
pdraws %>% ggplot(aes(x=mp,y=pr*scale))+
  geom_line(aes(group=draw),alpha=0.2)+
  geom_line(data=pdrawsmn)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),col="red",alpha=0.6)+
  facet_grid(~param)+xlab("")+ylab("")



# scratch below
if(0){
mm<-15
br0<-unique(PI(F0=qnorm,i=1:2^mm,m=mm))
ypr<-parprob(gam=100,m=mm)
mrg0<-cbind(i0=br0[1:(2^mm)],i1=br0[2:(2^mm+1)],ypr)
mrg<-data.frame(mrg0[-c(1,2^mm),]) # trim intervals at ends w/ -Inf and Inf

hbr<-50
ctlev<-cut(seq(mrg[1,"i0"],mrg[2^(mm)-2,"i1"],length=hbr),hbr)
lvs0<-paste0(gsub("\\(||\\]","",levels(ctlev)),sep=",",collapse="")
lvs1<-substr(lvs0,1,nchar(lvs0)-1)
lvs<-unique(as.numeric(unlist(strsplit(lvs1, split=","))))

ct <- factor(findInterval(mrg$i1,lvs),levels=1:hbr)
prplt<-data.frame(cbind(mp=unique(midpoints(ctlev)),pr=by(mrg$ypr,ct,sum)))

plot(prplt$mp,prplt$pr,type="l")


condY(gam=1,1)  # Y level m=1 
condY(gam=1,2)  # Y level m=2
condY(gam=1,3)  # Y level m=3

sample(c(0,1),3,replace=TRUE) # binary seq

rep(condY(gam=1,1),each=2^(3-1))
rep(condY(gam=1,2),each=2^(3-2))
rep(condY(gam=1,3),each=2^(3-3))

condY<-function(gam,m){
  i=1
  y<-vector("numeric",m)
  while(i<=m){
    # print(i) 
    #  print(rbeta(1,alp,alp))
    alp=gam*i^2
    if (i==1) { 
      y[i]<-rbeta(i,alp,alp)}
    else {
      y<-c(y,rbeta(2^(m-1),alp,alp))
    }
    i=i+1
  }
  c(rbind(y,1-y)) 
}


br1<-br0[2:(2^mm-1)]
br<-br0
br[1] <- br[2]-sd(br1)
br[(2^mm+1)] <- br[(2^mm)]+sd(br1)

xbr<-(br[1:(2^mm)]+br[2:(2^mm+1)])/2

plot(xbr,ypr,type="l")

mpbr<-data.frame(cbind(xbr,ypr))

breaks<-cut(xbr,5)

by(mpbr$ypr,breaks,sum)

xseq<-data.frame(xbr=seq(xbr[1],xbr[2^mm],length=10),ypr=0)
mpbre<-rbind(mpbr,xseq)
mpbre$yprnew<-0
mpbre[order(mpbre$xbr),]
}