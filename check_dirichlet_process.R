
### check Dirichlet process ###

## manual ##

# draw from prior (stick-breaking process)
set.seed(2743)
nn <- 1000 # finite approx to DP
alph <- 100
base_dis <- function(...) rnorm(...)

# is there a distribution that produces just the weights from stick-breaking process?? 
# if so it would be similar to draws from dirichlet, right??

samp_prior <- function(nn,alph){
V <- rbeta(nn,1,alph)
U <- base_dis(nn)
w <- vector("numeric",length(V))

w[1]<-V[1]
V[nn]<-1
for(i in 2:length(w)){
  w[i]<-V[i]*prod(1-V[1:(i-1)])
}

dat <- cbind(U,w)
dat <- dat[order(dat[,1]),]
dat <- cbind(dat,cw=cumsum(dat[,2]))

return(dat)
}

pr_draw <- samp_prior(nn,alph)

plot(pr_draw[,1], pr_draw[,2], type="h") # prior prob mass function
plot(pr_draw[,1], pr_draw[,3], type="s") # prior cdf

# draw from marginal
nd <- 50 # number of draws

# method 1 (sample F, then draw X1,..X20 from F)

# i) redraw new F for each X_i ?? think should use ii
# marg1i <- vector("numeric",nd)
# for(i in 1:nd){
#   pr_draw <- samp_prior(nn,alph)
#   marg1i[i] <- sample(pr_draw[,1],1,replace=TRUE,prob=pr_draw[,2])
# }

# ii) sample all X_i from same F ??
marg1ii <- sample(pr_draw[,1],nd,replace=TRUE,prob=pr_draw[,2])

# method 2 (Chinese restaurant process)
marg2 <- vector("numeric",nd)
marg2[1] <- base_dis(1)

for(i in 2:nd){
  emp <- sample(marg2[1:(i-1)],1) # empirical dist of X_1 to X_{i-1}
  base <- base_dis(1) # base measure
  p_emp <- (i-1)/(i+alph-1)
  
  marg2[i]<- ifelse(runif(1)<=p_emp,emp,base)
}

#qqplot(marg1i,marg2)
#qqplot(marg1ii,marg2)

#hist(marg1i, freq=FALSE, breaks=50)
hist(marg1ii, freq=FALSE, breaks=100)
hist(marg2, freq=FALSE, breaks=100)


#plot(ecdf(marg1i))
plot(ecdf(marg1ii))
plot(ecdf(marg2))


# samples
x <- rnorm(100,4,0.5)


# draw from posterior
#! check/collapse rows with same U in dat?
samp_post <- function(x,nn,alph){
  samp_n<-length(x)
  p_samp<-samp_n/(samp_n+alph)
  
  # posterior = comb. of empirical x and prior
  U <- replicate(nn, ifelse(runif(1)<=p_samp,sample(x,1),base_dis(1)))
  V <- rbeta(nn,1,alph+samp_n)
  w <- vector("numeric",length(V))
  
  w[1]<-V[1]
  V[nn]<-1
  for(i in 2:length(w)){
    w[i]<-V[i]*prod(1-V[1:(i-1)])
  }

  dat <- cbind(U,w)
  dat <- dat[order(dat[,1]),]
  dat <- cbind(dat,cw=cumsum(dat[,2]))
  
  return(dat)
}

#functions to calc mean of DP samples
#nn <- 15
#sp1<-samp_post(x1,nn,alph)
#sp2<-samp_post(x1,nn,alph)

draws_mean_grid<-function(draws, grid_seq=seq(-4,6,length=20)){
  gr<-data.frame(U=grid_seq, gridU=1)
  mrg<-merge(gr,draws,all=TRUE)
  
  mrg$w0<-mrg$w
  mrg[is.na(mrg$w0),"w0"] <- 0
  mrg$mcw<-cumsum(mrg$w0)
  subset(mrg,gridU==1,c("U","mcw"))
}

#draws_mean_grid(sp1)

draws_mean_calc<-function(draws_reps, grid=seq(-4,6,length=20)){
  df<-do.call(data.frame,apply(draws_reps,3,function(x) draws_mean_grid(x,grid)))
  
  mcw_nms<-names(df)[grep("mcw",names(df))]
  df[,c("U",mcw_nms)]
  
  data.frame(U=df$U,mn_cw=rowMeans(df[,mcw_nms]))
}

#pdl<-replicate(10,samp_post(x1,nn,alph))
#draws_mean_calc(pdl,grid=seq(-10,10,length=100))



plot_combo<-function(x, prdraws=5, postdraws=10){

# draws from prior
pr_draw <- replicate(prdraws,samp_prior(nn,alph))

# draws from posterior
post_draw <- replicate(postdraws,samp_post(x,nn,alph))

xmin<-floor(min(pr_draw[,1,],post_draw[,1,],x))*1.05 
xmax<-ceiling(max(pr_draw[,1,],post_draw[,1,],x))*1.05 

pr_draw_mean<-draws_mean_calc(pr_draw,grid=seq(xmin,xmax,length=nn*3))
post_draw_mean<-draws_mean_calc(post_draw,grid=seq(xmin,xmax,length=nn*3))

# blank plot
plot(x=0,y=0.5, typ="n",xlim=c(xmin,xmax),ylim=c(0,1),
     xlab="",ylab="cdf")

# plot draws from prior
for (i in 1:prdraws){
  lines(pr_draw[,1,i], pr_draw[,3,i], type="s",
        col=rgb(red=1, green=0, blue=0, alpha=0.3))
}
lines(pr_draw_mean[,1],pr_draw_mean[,2], type="s",
      col=rgb(red=1, green=0, blue=0, alpha=1),lwd=2,lty=2)

# plot draws from posterior
for (i in 1:postdraws){
  lines(post_draw[,1,i], post_draw[,3,i], type="s",
        col=rgb(red=0, green=0, blue=1, alpha=0.3))
}
lines(post_draw_mean[,1],post_draw_mean[,2], type="s",
      col=rgb(red=0, green=0, blue=1, alpha=1),lwd=2,lty=2)

# plot sample
samp_ecdf <- ecdf(x)
lines(sort(x), samp_ecdf(sort(x)), type="s", col="green", lwd=2)

legend("topleft",legend=c("prior draws", "prior mean",
                          "posterior draws", "posterior mean",
                          "likelihood"),
       lty=c(1,2,1,2,1), lwd=c(1,2,1,2,2),
       col=c("red","red","blue","blue","green"))
}

nn <- 1000 # finite approx to DP
alph <- 5

# continuous
x0 <- rnorm(15,3,2)
plot_combo(x0)

# discrete
x1 <- rnorm(50,3,0.5)
x1 <- as.numeric(cut(x1,7))-2
plot_combo(x1)

# continuous below threshold, categorical above
x2 <- rnorm(50,2,0.5)
x2[x2>2]<-as.numeric(cut(x2[x2>2],4))+1
plot_combo(x2)

# lower limit
x3 <- rnorm(60,2.5,1.5)
x3[x3<2]<-2
plot_combo(x3)

# skewed
x4<-rgamma(40,3,0.75)
plot_combo(x4)

# bimodal
x5<-rep(NA,75)
for(i in 1:length(x5)){
  if(runif(1)<0.6){
    x5[i]<-rnorm(1,0.9,1)
  } else {
    x5[i]<-rnorm(1,4,0.2)
  }
}
hist(x5)
plot_combo(x5)

# try with uniform
nn <- 1000 # finite approx to DP
alph <- 5 # concentration param
base_dis <-function(...) runif(...)

pr_draw_unif <- samp_prior(nn,alph)

plot(pr_draw_unif[,1], pr_draw_unif[,2], type="h") # prior prob mass function

plot(pr_draw_unif[,2], type="h")

# samples from discrete uniform
x6<-sample(seq(0,1,by=1/64),50,replace=TRUE)
plot_combo(x6, postdraws=25)

# samples from uniform
x7 <- rbeta(58,1,1)
alph <- 10 # concentration param
plot_combo(x7, prdraws=5, postdraws=15)

## with dirichletprocess ##
#library(dirichletprocess)


