library(MCMCpack) # for dirichlet dist
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# cutpoint functions
funa <- function(J) 2/max(1, J-5)
funb <- function(J) 1/J            # 0.042
func <-function(J) 1/(J + 2)
fund <-function(J) 1/sqrt(J)
fune <-function(J) 1/(max(3, J))
funf <-function(J) 1/min(J, 20)   # not bad esp late 0.026
fung <-function(J) 1/min(J, 30)  # not as good
funh <-function(J) 1 / (2 + (J/3))  # winner below
#funi <-function(J) 1/(0.8 + 0.35 * max(J, 3))
funi <-function(J) 1/(0.8 + 0.35 * J)

funj <- function(g) 1
funk <- function(g) 1/2

#induced prior on cutpoints

set.seed(958)
#n <- 10
#cats<-50

if(0){

cnc <- funb(50)
ddraws <- rdirichlet(1e4,rep(cnc,50))
cdraws <- apply(ddraws,1,cumsum) # cdf

ctpt <- apply(cdraws,1,qlogis) # cutpoints

# ctpoints based on MLE when beta_x is 0
mlectpt <- qlogis(cumsum(rep(1/50,50)))

ctpt_qt <- apply(ctpt, 2, quantile, probs=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
if(mnplt) ctpt_mn <- apply(ctpt, 2, mean, na.rm=TRUE)

dt <- as_tibble(t(ctpt_qt)) %>% add_rownames() %>% 
  rename(cat=rowname) %>% mutate(cat=as.numeric(cat))

ggplot(dt,aes(x=cat,y=`50%`))+
  geom_ribbon(aes(ymin=`25%`,ymax=`75%`),fill="blue", alpha=0.5)+
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),fill="blue", alpha=0.25)+
  geom_line(col="blue") +
  coord_cartesian(ylim=c(-15,15))+
  ylab(expression(gamma[j]))+xlab("j")
}

plot_induced_full_gg <- function(conc_fun, cats, link=qlogis,
                              mleplt=TRUE, n=1e4,
                              prbs=c(0.025,0.25,0.5,0.75,0.975),
                              yl=c(-15,15)){
  cnc <- conc_fun(cats)
  ddraws <- rdirichlet(n,rep(cnc,cats))
  cdraws <- apply(ddraws,1,cumsum) # cdf
  
  ctpt <- apply(cdraws,1,link)[,1:(cats-1)] # cutpoints
  
  # ctpoints based on MLE when beta_x is 0
  mlectpt <- link(cumsum(rep(1/cats,cats)))[1:(cats-1)]
  
  ctpt_qt <- apply(ctpt, 2, quantile, probs=prbs, na.rm=TRUE)
  
  dt <- as_tibble(t(ctpt_qt)) %>% add_rownames() %>% 
    rename(cpt=rowname) %>% mutate(cpt=as.numeric(cpt))
  
  V<-paste0(prbs*100,"%")
  
  plt<-ggplot(dt,aes(x=cpt,y=get(V[3])))+
    geom_ribbon(aes(ymin=get(V[2]),ymax=get(V[4])),fill="blue", alpha=0.5)+
    geom_ribbon(aes(ymin=get(V[1]),ymax=get(V[5])),fill="blue", alpha=0.25)+
    geom_line(col="blue") +
    coord_cartesian(ylim=yl)+
    ylab(expression(gamma[j]))+xlab("j")
  
  if (mleplt) plt<-plt+geom_line(data=data.frame(cpt=1:(cats-1),mle=mlectpt),
                            aes(x=cpt,y=mle),col="red",linetype="dashed")

plt
}

plot_induced_full_gg(funb,cats=100)


plot_induced_full_gg_dat <- function(conc_fun, cats, link=qlogis,
                                 mleplt=TRUE, n=1e4,
                                 prbs=c(0.025,0.25,0.5,0.75,0.975),
                                 yl=c(-15,15)){
  cnc <- conc_fun(cats)
  ddraws <- rdirichlet(n,rep(cnc,cats))
  cdraws <- apply(ddraws,1,cumsum) # cdf
  
  ctpt <- apply(cdraws,1,link)[,1:(cats-1)] # cutpoints
  
  ctpt_qt <- apply(ctpt, 2, quantile, probs=prbs, na.rm=TRUE)
  
  dt <- as_tibble(t(ctpt_qt)) %>% add_rownames() %>% 
    rename(cpt=rowname) %>% 
    mutate(cpt=as.numeric(cpt),ncats=cats,fun=deparse(conc_fun)[2])
  
  if(mleplt){
  # ctpoints based on MLE when beta_x is 0
  mlectpt <- link(cumsum(rep(1/cats,cats)))[1:(cats-1)]
  dt <- dt %>% mutate(mlectpt=mlectpt)
  }
  
dt
}

#plot_induced_full_gg_dat(funb,cats=100,n=10)

plotmle<-TRUE
lnk <- qnorm
a1<-plot_induced_full_gg_dat(funb,cats=25,link=lnk,mleplt=plotmle) # 1/ncats
a2<-plot_induced_full_gg_dat(funi,cats=25,link=lnk,mleplt=plotmle)
a3<-plot_induced_full_gg_dat(funh,cats=25,link=lnk,mleplt=plotmle)
a4<-plot_induced_full_gg_dat(funk,cats=25,link=lnk,mleplt=plotmle) # 1/2
a5<-plot_induced_full_gg_dat(funj,cats=25,link=lnk,mleplt=plotmle) # 1

b1<-plot_induced_full_gg_dat(funb,cats=50,link=lnk,mleplt=plotmle) # 1/ncats
b2<-plot_induced_full_gg_dat(funi,cats=50,link=lnk,mleplt=plotmle)
b3<-plot_induced_full_gg_dat(funh,cats=50,link=lnk,mleplt=plotmle)
b4<-plot_induced_full_gg_dat(funk,cats=50,link=lnk,mleplt=plotmle) # 1/2
b5<-plot_induced_full_gg_dat(funj,cats=50,link=lnk,mleplt=plotmle) # 1

c1<-plot_induced_full_gg_dat(funb,cats=100,link=lnk,mleplt=plotmle) # 1/ncats
c2<-plot_induced_full_gg_dat(funi,cats=100,link=lnk,mleplt=plotmle)
c3<-plot_induced_full_gg_dat(funh,cats=100,link=lnk,mleplt=plotmle)
c4<-plot_induced_full_gg_dat(funk,cats=100,link=lnk,mleplt=plotmle) # 1/2
c5<-plot_induced_full_gg_dat(funj,cats=100,link=lnk,mleplt=plotmle) # 1

ppr_fig_dir<-file.path('/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper/fig')

pltw<-12
plth<-9
atxtsz<-9
fctsiz<-13

pltdt<-bind_rows(a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,c1,c2,c3,c4,c5) %>% 
  mutate(conc=factor(fun,levels=c('1/J','1/(2 + (J/3))','1/(0.8 + 0.35 * J)','1/2','1'))) 

pltdt %>% 
  ggplot(aes(x=cpt,y=`50%`)) + 
  geom_line(aes(color='median')) +
  geom_ribbon(aes(ymin=`25%`,ymax=`75%`,fill='50% int'))+
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill='95% int'))+
  scale_color_manual(labels = c('median'),
                     values = c('blue')) +
  scale_fill_manual(labels= c('50% int','95% int'),
                    values = alpha('blue',c(0.4,0.2))) +
  facet_wrap(~ncats+conc, labeller= label_both, scales = "free_x", nrow=4, ncol=5) +
  labs(color='', fill='') +
  coord_cartesian(ylim=c(-10,10))+
  ylab(expression(gamma[j]))+xlab("j")+
  theme(axis.title.x=element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        legend.text = element_text(size=fctsiz))

ggsave(file.path(ppr_fig_dir,'probit_induced.png'),width=pltw,height=plth)


pltdt %>% pivot_longer(cols=c(`50%`,mlectpt)) %>% 
  ggplot(aes(x=cpt,y=value)) + 
  geom_ribbon(aes(ymin=`25%`,ymax=`75%`,fill='50% int'))+
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill='95% int'))+
  geom_line(aes(color=name, linetype = name)) +
  scale_color_manual(name='',labels=c('median','MLE'), values = c('blue','red')) +
  scale_linetype_manual(name='',labels=c('median','MLE'),values = c('solid','dashed')) +
  scale_fill_manual(labels= c('50% int','95% int'),
                    values = alpha('blue',c(0.4,0.2))) +
  facet_wrap(~ncats+conc, labeller= label_both, scales = "free_x", nrow=4, ncol=5) +
  labs(color='', fill='') +
  coord_cartesian(ylim=c(-10,10))+
  ylab(expression(gamma[j]))+xlab("j")+
  theme(axis.title.x=element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        legend.text = element_text(size=fctsiz))

ggsave(file.path(ppr_fig_dir,'probit_induced_mle.png'),width=pltw,height=plth)

# can try this for better legend
# https://stackoverflow.com/questions/28714492/legend-with-geom-line-and-geom-ribbon










pltdt %>% 
  ggplot(aes(x=cpt,y=`50%`)) + 
  geom_line(color='blue') +
  geom_ribbon(aes(ymin=`25%`,ymax=`75%`),fill="blue", alpha=0.4)+
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),fill="blue", alpha=0.2)+
  coord_cartesian(ylim=c(-10,10))+
  ylab(expression(gamma[j]))+xlab("j")+
  facet_wrap(~ncats+conc, labeller= label_both, scales = "free_x", nrow=4, ncol=4)




# logit
a1<-plot_induced_full_gg(funb,cats=25,mleplt=plotmle) # 1/ncats
a2<-plot_induced_full_gg(funi,cats=25,mleplt=plotmle) # 1/(0.8 + 0.35 * max(g, 3))
a3<-plot_induced_full_gg(funk,cats=25,mleplt=plotmle) # 1/2
a4<-plot_induced_full_gg(funj,cats=25,mleplt=plotmle) # 1

b1<-plot_induced_full_gg(funb,cats=50,mleplt=plotmle) # 1/ncats
b2<-plot_induced_full_gg(funi,cats=50,mleplt=plotmle) # 1/(0.8 + 0.35 * max(g, 3))
b3<-plot_induced_full_gg(funk,cats=50,mleplt=plotmle) # 1/2
b4<-plot_induced_full_gg(funj,cats=50,mleplt=plotmle) # 1

c1<-plot_induced_full_gg(funb,cats=100,mleplt=plotmle) # 1/ncats
c2<-plot_induced_full_gg(funi,cats=100,mleplt=plotmle) # 1/(0.8 + 0.35 * max(g, 3))
c3<-plot_induced_full_gg(funk,cats=100,mleplt=plotmle) # 1/2
c4<-plot_induced_full_gg(funj,cats=100,mleplt=plotmle) # 1


plot_grid(a1,a2,a3,a4,
          b1,b2,b3,b4,
          c1,c2,c3,c4,
          ncol=4,nrow=3)




# probit
lnk <- qnorm
yrng <- c(-10,10)
d1<-plot_induced_full_gg(funb,link=lnk,cats=25,yl=yrng,mleplt=plotmle) # 1/ncats
d2<-plot_induced_full_gg(funi,link=lnk,cats=25,yl=yrng,mleplt=plotmle) # 1/(0.8 + 0.35 * max(g, 3))
d3<-plot_induced_full_gg(funk,link=lnk,cats=25,yl=yrng,mleplt=plotmle) # 1/2
d4<-plot_induced_full_gg(funj,link=lnk,cats=25,yl=yrng,mleplt=plotmle) # 1

e1<-plot_induced_full_gg(funb,link=lnk,cats=50,yl=yrng,mleplt=plotmle) # 1/ncats
e2<-plot_induced_full_gg(funi,link=lnk,cats=50,yl=yrng,mleplt=plotmle) # 1/(0.8 + 0.35 * max(g, 3))
e3<-plot_induced_full_gg(funk,link=lnk,cats=50,yl=yrng,mleplt=plotmle) # 1/2
e4<-plot_induced_full_gg(funj,link=lnk,cats=50,yl=yrng,mleplt=plotmle) # 1

f1<-plot_induced_full_gg(funb,link=lnk,cats=100,yl=yrng,mleplt=plotmle) # 1/ncats
f2<-plot_induced_full_gg(funi,link=lnk,cats=100,yl=yrng,mleplt=plotmle) # 1/(0.8 + 0.35 * max(g, 3))
f3<-plot_induced_full_gg(funk,link=lnk,cats=100,yl=yrng,mleplt=plotmle) # 1/2
f4<-plot_induced_full_gg(funj,link=lnk,cats=100,yl=yrng,mleplt=plotmle) # 1


plot_grid(d1,d2,d3,d4,
          e1,e2,e3,e4,
          f1,f2,f3,f4,
          ncol=4,nrow=3)





library(rms)
x<-runif(100)
y1<-rlogis(100)+0*x
orm_logit1<-orm(y1~x)

par(mfrow=c(1,1))
plot_induced_full(funb,cats=100,mleplt=T)
lines(1:99,-coef(orm_logit1)[1:99],col="blue")

library(bayesCPM)
library(dplyr)
dat <- data.frame(y=ordered(y1),x)
dat_stan  <- mkStanDat(dat, outcome="y", preds = c("x"), link=1)

## sample from Bayes CPM model with probit link
options(mc.cores=parallel::detectCores())
bcpm_logit1 <- bayes_cpm(dat_stan)

data.frame(bcpm_logit1)



y2<-rlogis(100)+2.5*x
orm_logit2<-orm(y2~x)
plot_induced_full(funb,cats=100,mleplt=T)
lines(1:99,-coef(orm_logit2)[1:99],col="blue")



# induced prior is centered around theoretical MLE value when beta_x = 0 
# very small values of alpha provide little info about location of cutpoints

# large alpha values are more concentrated around theoretical MLE value when beta_x = 0
# and will have large amount of shrinkage to these values



# scratch below
if(0){
  plot_induced<-function(conc_fun, cats, n=3000, brks=30){
    cnc <- conc_fun(cats)
    ddraws<-rdirichlet(n,rep(cnc,cats))
    cdraws<-apply(ddraws,1,cumsum) # cdf
    
    # NB: cumsum() can return values very slightly greater than 1 
    # due to numerical imprecision
    # these values will return NaN when passed to qlogis()
    # fix by setting max to 1?
    
    ctpt<-apply(cdraws,1,qlogis) # cutpoints
    
    # ctpoints based on MLE 
    mlectpt<-qlogis(cumsum(rep(1/cats,cats)))
    
    par(mfrow=c(1,5))
    qnt<-quantile(1:cats,probs=c(0.05,0.25,0.5,0.75,0.95)) # change to probs 5, 25, 50, 75, 95
    
    hist(ctpt[,qnt[1]],xlim=c(-50,5),breaks=brks)
    abline(v=mlectpt[qnt[1]],col="red")
    
    hist(ctpt[,qnt[2]],xlim=c(-25,15),breaks=brks)
    abline(v=mlectpt[qnt[2]],col="red")
    
    hist(ctpt[,qnt[3]],xlim=c(-20,20),breaks=brks)
    abline(v=mlectpt[qnt[3]],col="red")
    
    hist(ctpt[,qnt[4]],xlim=c(-15,25),breaks=brks)
    abline(v=mlectpt[qnt[4]],col="red")
    
    hist(ctpt[,qnt[5]],xlim=c(-5,50),breaks=brks)
    abline(v=mlectpt[qnt[5]],col="red")
  }
  
  
  plot_induced(funj,cats=5)
  
  plot_induced(funj,cats=100)
  plot_induced(funj,cats=250)
  plot_induced(funj,cats=500)
  
  plot_induced(funb,cats=100)
  plot_induced(funb,cats=250)
  plot_induced(funb,cats=500)
  
  plot_induced(funi,cats=100)
  plot_induced(funi,cats=250)
  plot_induced(funi,cats=500)
  
  
  plot_induced(funj,cats=100)
  
  
  
  plot_induced_full <- function(conc_fun, cats, link=qlogis,
                                mleplt=TRUE, mnplt=FALSE, n=1e4){
    cnc <- conc_fun(cats)
    ddraws<-rdirichlet(n,rep(cnc,cats))
    cdraws<-apply(ddraws,1,cumsum) # cdf
    
    ctpt<-apply(cdraws,1,link) # cutpoints
    
    # ctpoints based on MLE when beta_x is 0
    mlectpt<-link(cumsum(rep(1/cats,cats)))
    
    ctpt_qt <- apply(ctpt, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=TRUE)
    if(mnplt) ctpt_mn <- apply(ctpt, 2, mean, na.rm=TRUE)
    
    if(mleplt) {
      plot(1:cats,mlectpt,col="red",ylim=c(-18,18),type="l",
           xlab="j",ylab="gamma_j")
    } else {
      plot(1:cats,mlectpt,col="red",ylim=c(-18,18),type="n",
           xlab="j",ylab="gamma_j")
    }
    lines(1:cats,ctpt_qt[1,])
    lines(1:cats,ctpt_qt[2,])
    lines(1:cats,ctpt_qt[3,])
    if(mnplt) lines(1:cats,ctpt_mn,lty="dashed")
  }
  
  # logit
  par(mfrow=c(3,3))
  plot_induced_full(funb,cats=25,mleplt=F) # 1/ncats
  plot_induced_full(funi,cats=25,mleplt=F) # 1/(0.8 + 0.35 * max(g, 3))
  #plot_induced_full(fund,cats=25,mleplt=F) # 1/sqrt(ncats)
  plot_induced_full(funj,cats=25,mleplt=F) # 1
  
  plot_induced_full(funb,cats=50,mleplt=F) # 1/ncats
  plot_induced_full(funi,cats=50,mleplt=F) # 1/(0.8 + 0.35 * max(g, 3))
  #plot_induced_full(fund,cats=50,mleplt=F) # 1/sqrt(ncats)
  plot_induced_full(funj,cats=50,mleplt=F) # 1
  
  plot_induced_full(funb,cats=100,mleplt=F) # 1/ncats
  plot_induced_full(funi,cats=100,mleplt=F) # 1/(0.8 + 0.35 * max(g, 3))
  #plot_induced_full(fund,cats=100,mleplt=F) # 1/sqrt(ncats)
  plot_induced_full(funj,cats=100,mleplt=F) # 1
  
  
  # probit
  plot_induced_full(funb,link=qnorm, cats=100, mleplt=T) # 1/ncats
  plot_induced_full(funi,link=qnorm, cats=100, mleplt=T) # 1/(0.8 + 0.35 * max(g, 3)
  plot_induced_full(fund,link=qnorm, cats=100, mleplt=T) # 1/sqrt(ncats)
  plot_induced_full(funj,link=qnorm, cats=100, mleplt=T) # 1
  
  
  
}

if(0){
#for what conc does ddirichlet(adraw,rep(conc,nn)) * abs(Jac()) = 1 ?
# i.e. for what conc does
# log(ddirichlet(adraw,rep(conc,nn))) + log(abs(Jac())) = 0 ?

nn<-5
adraw<-rdirichlet(1,rep(0.5,nn))

ddirichlet(adraw,rep(0.5,nn))

# make cutpoints from dirichlet draws

cuts<-qlogis(cumsum(adraw)[1:(nn-1)]) # cutpoints

# make Jacobian function
Jac <- function(ctpt, printJ=TRUE){
  ln<-length(ctpt)
  g <- dlogis(ctpt)
  A0 <- diag(-dlogis(ctpt))
  
  A <- rbind(A0,rep(0,nrow(A0)))
  B <- rbind(rep(0,nrow(A0)),-A0)
  C <- A+B

  J <- cbind(rep(1,ln+1),C)
  if(printJ) print(J)
  
  det(J)
}

term <- function(alpha,nn=5,draws=500){
  # adraw<-rdirichlet(5,rep(alpha,nn))
  # cuts<-qlogis(cumsum(adraw)[1:(nn-1)])
  # 
  # ddirichlet(adraw,rep(alpha,nn))*abs(Jac(cuts,F))
  
  adraw<-rdirichlet(draws,rep(alpha,nn))
  
  cdraw<-apply(adraw,1,cumsum)
  
  cuts<-apply(cdraw,1,qlogis)[,1:(nn-1)]
  
  #cuts<-qlogis(cumsum(adraw)[1:(nn-1)])
  
  ddir<-ddirichlet(adraw,matrix(alpha,nrow=nrow(adraw),ncol=ncol(adraw)))
  
  absJ<-abs(apply(cuts,1,Jac,printJ=F))
  
  ddir*absJ
  
}


nn<-5
alpha<-0.4
draws<-10
adraw<-rdirichlet(draws,rep(alpha,nn))

cuts<-apply(adraw,1,qlogis)

#cuts<-qlogis(cumsum(adraw)[1:(nn-1)])

ddir<-ddirichlet(adraw,matrix(alpha,nrow=nrow(adraw),ncol=ncol(adraw)))

absJ<-abs(apply(cuts,2,Jac,printJ=F))

ddir*absJ

ddirichlet(adraw,rep(alpha,nn))*abs(Jac(cuts,F))

mean(term(1/5,draws=1000))
mean(term(0.3921569,draws=1000))

1/(0.8 + 0.35 * max(5, 3))

mean(term(1/50,nn=50,draws=1000))
mean(term(0.05464481,nn=50,draws=1000))

1/(0.8 + 0.35 * max(50, 3))

}


if(0){
# plot rank vs cdf value
par(mfrow=c(1,1))
for(i in 1:ncol(cdraws)){
  if(i==1){
    plot(1:cats,cdraws[,i],type="s")
  } else {
    lines(1:cats,cdraws[,i],type="s")
  }
  
}


curve(dlogis,-20,20,n=1000)
abline(v=ctpt[3,],col='green')
abline(v=ctpt[4,],col='red')
}



