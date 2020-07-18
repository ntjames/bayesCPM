
# check Dirichlet density
library(plotly)
library(tibble)
library(MCMCpack)
library(gtools)

# from Wikipedia
# Values of the concentration parameter above 1 prefer variates that are dense, evenly distributed distributions, i.e. all the values within a single sample are similar to each other. Values of the concentration parameter below 1 prefer sparse distributions, i.e. most of the values within a single sample will be close to 0, and the vast majority of the mass will be concentrated in a few of the values. 

# alpha (conc) = 1, uniform over discrete dists with 4 categories
d1<-rdirichlet(20, c(1,1,1,1))

# alpha (conc) = 15, discrete dists with 4 balanced category probabilities more likely
d2<-rdirichlet(20, c(15,15,15,15))

# alpha (conc) = 0.5, discrete dists with 4 unbalanced category probabilities more likely
d3<-rdirichlet(20, c(0.5,0.5,0.5,0.5))


a<-rdirichlet(5,rep(1,7)) 
b<-rdirichlet(5,rep(10,7)) 
c<-rdirichlet(5,rep(0.1,7))

plot(a[1,], type="h",ylim=c(0,1))
plot(a[2,], type="h",ylim=c(0,1))
plot(a[3,], type="h",ylim=c(0,1))

cumsum(a[2,])
plot(ecdf(cumsum(a[2,])))
plot(ecdf(a[2,]))


plot(b[1,], type="h",ylim=c(0,1))
plot(b[2,], type="h",ylim=c(0,1))
plot(b[3,], type="h",ylim=c(0,1))


plot(c[1,], type="h",ylim=c(0,1))
plot(c[2,], type="h",ylim=c(0,1))
plot(c[3,], type="h",ylim=c(0,1))

mat<-rdirichlet(10,c(1,1,1))

s<-seq(0,1,by=0.05)
sg<-expand.grid(s,s)
sgm<-sg[rowSums(sg)<=1,]

sdiff <- 1-rowSums(sgm)

mat<-cbind(sgm,sdiff)

z<-ddirichlet(mat,c(2,2,2))

dat<-cbind(x=mat[,1],y=mat[,3],z) %>% as_tibble()

plot_ly(dat, x=~x, y=~y, z= ~z, type="scatter3d", mode="markers", color=~z, size=1)




ddirichlet(c(0.3,0.7),c(1,1))
ddirichlet(c(.1,.2,.6,.1), c(1,1,1,1))

alp<-c(1,1,1,1)
thet<-c(.1,.2,.6,.1)

(gamma(sum(alp))/prod(sapply(alp,gamma)))*prod(thet^(alp-1))








# prior
library(tidyverse)
ncats<-50
conc<-1
alp_pr <- rep(conc,ncats)
ndraws<-12
pr<-rdirichlet(ndraws,alp_pr) %>% 
  as_tibble() %>% 
  mutate(draw=1:ndraws) %>% 
  pivot_longer(cols=paste0("V",1:ncats)) %>% 
  mutate(cat=gsub("V","",name))

ggplot(pr, aes(x=cat,y=value))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~draw)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# post
alp_pst <- alp_pr+1
alp_pst <- rep(0,length(alp_pr))+1
pst<-rdirichlet(ndraws,alp_pst) %>% 
  as_tibble() %>% 
  mutate(draw=1:ndraws) %>% 
  pivot_longer(cols=paste0("V",1:ncats)) %>% 
  mutate(cat=gsub("V","",name))

ggplot(pst, aes(x=cat,y=value))+
  geom_bar(stat="identity",width=0.7)+
  facet_wrap(~draw)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

post<-function(x,n,c){
  unif<-(x+1)/(n+c)
  jeffreys<-(x+1/2)/(n+c/2)
  ref<-(x+1/2)/(n+1)
  obj<-(x+1/c)/(n+1)
  data.frame(unif,jeffreys,ref,obj)
}

post(c(2,1),3,1000)

post(1,1000,1000)

pst_wide <- pst %>% pivot_wider(id_cols=draw)

cutpts <- as.numeric(pst_wide[1,2:(ncats+1)]) %>% cumsum() %>% qlogis()

cdf_est <- plogis(cutpts)
