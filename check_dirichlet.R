
# check Dirichlet density
library(plotly)
library(MCMCpack)

rdirichlet(10,c(1,1,1,1,1))
rdirichlet(10,c(10,10,10,10,10))
rdirichlet(10,c(0.1,0.1,0.1,0.1,0.1))

mat<-rdirichlet(10,c(1,1,1))

s<-seq(0,1,by=0.05)
sg<-expand.grid(s,s)
sgm<-sg[rowSums(sg)<=1,]

sdiff <- 1-rowSums(sgm)

mat<-cbind(sgm,sdiff)

z<-ddirichlet(mat,c(2,2,2))

dat<-cbind(x=mat[,1],y=mat[,3],z) %>% as_tibble()

plot_ly(dat, x=~x, y=~y, z= ~z, type="scatter3d", mode="markers", color=~z, size=1)
