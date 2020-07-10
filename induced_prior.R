library(MCMCpack)

# cutpoint functions
funa <- function(g) 2/max(1, g-5)
funb <- function(g) 1/g            # 0.042
func <-function(g) 1/(g + 2)
fund <-function(g) 1/sqrt(g)
fune <-function(g) 1/(max(3, g))
funf <-function(g) 1/min(g, 20)   # not bad esp late 0.026
fung <-function(g) 1/min(g, 30)  # not as good
funh <-function(g) 1 / (2 + (g/3))  # winner below
funi <-function(g) 1/(0.8 + 0.35 * max(g, 3))

funj <-function(g) 1

#induced prior on cutpoints


set.seed(958)
#n <- 10
#cats<-50

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


fun0 <- function(g) 1
plot_induced(fun0,cats=5)

plot_induced(fun0,cats=100)
plot_induced(fun0,cats=250)
plot_induced(fun0,cats=500)

plot_induced(funb,cats=100)
plot_induced(funb,cats=250)
plot_induced(funb,cats=500)

plot_induced(funi,cats=100)
plot_induced(funi,cats=250)
plot_induced(funi,cats=500)


plot_induced(funj,cats=100)


nn<-5
adraw<-rdirichlet(1,rep(0.5,nn))

log(ddirichlet(adraw,rep(0.5,nn)))

# make cutpoints from dirichlet draws

qlogis(cumsum(adraw)[1:(nn-1)]) # cutpoints


# make Jacobian function
Jac <- function(ctpt){
  g <- dlogis(ctpt)
  
}

#for what conc does ddirichlet(adraw,rep(conc,nn)) * abs(Jac()) = 1 ?
# i.e. for what conc does
# log(ddirichlet(adraw,rep(conc,nn))) + log(abs(Jac())) = 0 ?


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



