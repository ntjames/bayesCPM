require(rms)
# set.seed(1) # for blrm_plots_conc_*_1
#set.seed(77283) # for blrm_plots_conc_*_2
set.seed(9149) # for blrm_plots_conc_*_3
x <- rnorm(800)
y <- rnorm(800) + x

#stanCompile() # run once
options(stancompiled='~/R/stan')

# use all cores
options(mc.cores=parallel::detectCores())

blrmdir<-file.path('~','Dropbox','njames','school','PhD','orise_ra','bayes_cpm','blrm_ex')

f <- orm(y ~ x)
concs <- c(0.0001, 0.001, .005, .01, .02, 0.03, 0.04, .05, .07, 0.1, 0.5, 1)
for(conc in concs) {
g <- blrm(y ~ x, priorsd=100, conc=conc)
pdf(file.path(blrmdir,paste0('blrm_plots_conc',conc,'_3.pdf'))) 
plot(coef(f), coef(g)); abline(a=0,b=1,col='red')
d <- abs(coef(f) - coef(g))
mad <- round(mean(d), 4)
d90 <- round(quantile(d, 0.9), 4)
title(paste('conc:', conc, ' mad:', mad, ' d90:', d90))
dev.off()
}
# Best of MAD 0.005    D90 0.01    default=0.0036 1/k=0.001



if(0){
set.seed(3478)
n<-50
x1 <- rlogis(n)
x2 <- rbinom(n,1,0.8)
x3 <- rlogis(n)
x4 <- rexp(n)
x5 <- rbinom(n,1,0.3)

y0 <- rlogis(n) + 0*x1
y1 <- rlogis(n) + x1
y2 <- rlogis(n) + x1 + 3.7*x2 - 2.4*x3 + 0.5*x4 + 0.8*x5

f0 <- orm(y0 ~ x1)
g0a <- blrm(y0 ~ x1, priorsd=100, conc = 1/n)
g0b <- blrm(y0 ~ x1, priorsd=100, conc = 1/(0.8 + 0.35 * max(n, 3)))

plot(coef(f0), -qlogis(cumsum(rep(1/50,50)))); abline(a=0,b=1,col='red')

plot(coef(f0), coef(g0a)); abline(a=0,b=1,col='red')
points(coef(f0), coef(g0b), col='blue')

mean(abs(coef(f0) - coef(g0a)))
mean(abs(coef(f0) - coef(g0b)))



f1 <- orm(y1 ~ x1)
g1a <- blrm(y1 ~ x1, priorsd=100, conc = 1/n)
g1b <- blrm(y1 ~ x1, priorsd=100, conc = 1/(0.8 + 0.35 * max(n, 3)))

plot(coef(f1), coef(g1a)); abline(a=0,b=1,col='red')
points(coef(f1), coef(g1b), col='blue')

mean(abs(coef(f1) - coef(g1a)))
mean(abs(coef(f1) - coef(g1b)))

mad <- round(mean(d), 4)
d90 <- round(quantile(d, 0.9), 4)


f2 <- orm(y2 ~ x1+x2+x3+x4+x5)
g2a <- blrm(y2 ~ x1+x2+x3+x4+x5, priorsd=100, conc = 1/n)
g2b <- blrm(y2 ~ x1+x2+x3+x4+x5, priorsd=100, conc = 1/(0.8 + 0.35 * max(n, 3)))
g2c <- blrm(y2 ~ x1+x2+x3+x4+x5, priorsd=100, conc = 2.85/n)

mean(abs(coef(f2) - coef(g2a)))
mean(abs(coef(f2) - coef(g2b)))
mean(abs(coef(f2) - coef(g2c)))



plot(coef(f2), coef(g2a)); abline(a=0,b=1,col='red')
points(coef(f2), coef(g2b), col='blue')
points(coef(f2), coef(g2c), col='green')

plot(coef(g2a),coef(g2b)); abline(a=0,b=1,col='red')
  
}

foo <- function(g) 1/(0.8 + 0.35 * max(g, 3))
sapply(1:25,foo) * 1:25

tail(sapply(1:10000,foo) * 1:10000)



# try case when beta_x=0, i.e. x has no effect
set.seed(78593) # for blrm_plots_conc_*_4

# probit link
x <- rnorm(800)

y <- rnorm(800) + 0*x
f <- orm(y ~ x, family="probit")
plot(-coef(f)[1:799],qnorm(cumsum(rep(1/800,800)))[1:799]); abline(a=0,b=1,col='red')


y <- rnorm(800) + 1*x
sd(rnorm(1e5)+rnorm(1e5))  # var = 1 + 1^2
f <- orm(y ~ x, family="probit")
plot(-coef(f)[1:799],qnorm(cumsum(rep(1/800,800)),sd=sqrt(2))[1:799]); abline(a=0,b=1,col='red')

y <- rnorm(800) + 1.8*x
sd(rnorm(1e5)+1.8*rnorm(1e5)) # var = 1 + 1.8^2
f <- orm(y ~ x, family="probit")
plot(-coef(f)[1:799],qnorm(cumsum(rep(1/800,800)),sd=sqrt(4.24))[1:799]); abline(a=0,b=1,col='red')


y <- rnorm(800) + 2.3*x
var(rnorm(1e5)+2.3*rnorm(1e5)) # var = 1 + 2.3^2
f <- orm(y ~ x, family="probit")
plot(-coef(f)[1:799],qnorm(cumsum(rep(1/800,800)),sd=sqrt(6.29))[1:799]); abline(a=0,b=1,col='red')

x2 <- rbinom(800,1,0.8)
y <- rnorm(800) + 2.3*x + x2
var(rnorm(1e5)+2.3*rnorm(1e5)+rbinom(1e5,1,0.8)) # var = 1 + 2.3^2 + 0.8*0.2
f <- orm(y ~ x+x2, family="probit")
plot(-coef(f)[1:799],qnorm(cumsum(rep(1/800,800)))[1:799]); abline(a=0,b=1,col='red')
plot(1:799,-coef(f)[1:799])


x3<-rlogis(800)
x4<-rexp(800)
y <- rnorm(800) + 2.3*x + x2 + 1.7*x3-3.1*x4
f <- orm(y ~ x+x2+x3+x4, family="probit")
plot(1:799,-coef(f)[1:799])


# logit link
x <- rlogis(800)
y <- rlogis(800) + 0*x
f <- orm(y ~ x)
plot(-coef(f)[1:799],qlogis(cumsum(rep(1/800,800)))[1:799]); abline(a=0,b=1,col='red')


x <- rlogis(800)
y <- rlogis(800) + 1*x
# sum of logistics is NOT logistic
sd(rlogis(1e5)+1*rlogis(1e5)) # var = pi^2/3*1^2 + pi^2/3*1^2
f <- orm(y ~ x)
plot(-coef(f)[1:799],qlogis(cumsum(rep(1/800,800)),scale=1)[1:799]); abline(a=0,b=1,col='red')

# cloglog
library(evd)
x <- rgumbel(800)
y <- rgumbel(800) + 0*x
f <- orm(y ~ x, family="cloglog")
plot(-coef(f)[1:799],qgumbel(cumsum(rep(1/800,800))[1:799])); abline(a=0,b=1,col='red')


# try case when beta_x=0, i.e. x has no effect
set.seed(78593) # for blrm_plots_conc_*_4
x <- rlogis(800)
y <- rlogis(800) + 0*x
f <- orm(y ~ x)
concs <- c(0.0001, 0.001, .005, .01, .02, 0.03, 0.04, .05, .07, 0.1, 0.5, 1, 10)
for(conc in concs) {
  g <- blrm(y ~ x, priorsd=100, conc=conc)
  pdf(file.path(blrmdir,paste0('blrm_plots_conc',conc,'_4.pdf'))) 
  plot(coef(f), coef(g)); abline(a=0,b=1,col='red')
  d <- abs(coef(f) - coef(g))
  mad <- round(mean(d), 4)
  d90 <- round(quantile(d, 0.9), 4)
  title(paste('conc:', conc, ' mad:', mad, ' d90:', d90))
  dev.off()
}

 
doit <- FALSE
if(doit) {

require(rms)
#stanSet()
set.seed(1)
n  <- 100
x  <- rnorm(n)
x2 <- rnorm(n)
x  <- x  - mean(x)
x2 <- x2 - mean(x2)
Y <- x + 2 * x2 + rnorm(n)

y <- round(3*Y)
i <- rep(1 : 20, 5)
f <- lrm(y ~ x + x2)
b <- blrm(y ~ x + x2, priorsd=1000)
bc <- blrm(y ~ x + x2 + cluster(i), priorsd=1000)
cbind(lrm=coef(f), blrm=coef(b, 'mean'), 'cluster blrm'=coef(bc, 'mean'))

d <- blrm(y ~ x + x2, priorsd=1000, standata=TRUE)
saveRDS(d, '~/tmp/d.rds')
}


if(doit) {
#set.seed(1) # for blrm_plots_conc_*_1.pdf and blrm_plots_conc_*_4.pdf (use stat='mode')
#set.seed(34678) # for blrm_plots_conc_*_2.pdf
set.seed(8891) # for blrm_plots_conc_*_3.pdf
n  <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 3 * x1 + 4 * x2 + 2 * rnorm(n)
#lrm(cut2(y, g=20) ~ pol(x1,2) + pol(x2,2))

do <- function(fun, eq=TRUE) {
  m <- 0
  gs <- c(2:4, seq(5, 100, by=5))
  #gs <- c(2:4, seq(5, 100, by=10))
  #gs<- c(10,20) # for test
  nm <- gsub('\\.','',make.names(enquote(fun))[2]) # make fun ok for filename
  pdf(file.path(blrmdir,paste0('blrm_plots_conc_',nm,'_3.pdf')),
      width=12, height=12) 
  par(mfrow=c(5,5), mar=c(3, 2.5, 1, 1))
  for(g in gs) {
    u <- if(eq) cut2(y, g=g)
         else cut2(y, if(g==2) 0 else  seq(-10, 10, length=g-1))
    if(length(unique(u)) != g) next
    f <- lrm (u ~ pol(x1,2) + x2)
    b <- blrm(u ~ pol(x1,2) + x2, conc=fun(g), priorsd=10000)
    u <- coef(f)[1 : (g-1)] 
    #v <- coef(b,stat="mode")[1 : (g-1)] # use coef(f,"mode") instead of default mean
    v <- coef(b)[1 : (g-1)]
    plot(u, v)
    abline(a=0, b=1, col='red')
    d <- abs(u - v)
    w <- paste("ncat:", g,
               "maxdif:", round(max(d), 3),
               "avgdif:", round(mean(d), 3))
    title(w, cex.main=0.8)
    m <- m + mean(d)
  }
  cat('Avg mean abs diff:', m / length(gs), '\n')
  pt <- paste('Avg mean abs diff:', m / length(gs))
  mtext(pt, side=1, at=50)
  dev.off()
}

do(function(g) 2/max(1, g-5))
do(function(g) 1/g)            # 0.042
do(function(g) 1/(g + 2))
do(function(g) 1/sqrt(g))
do(function(g) 1/(max(3, g)))
do(function(g) 1/min(g, 20))   # not bad esp late 0.026
do(function(g) 1/min(g, 30))   # not as good
do(function(g) 1 / (2 + (g/3)))  # winner below
do(function(g) 1/(0.8 + 0.35 * max(g, 3)), eq=FALSE)  # good w/n=1000
do(function(g) 1/(0.8 + 0.35 * max(g, 3) + 0), eq=TRUE)   # good w/n=1000
do(function(g) 1)   # try uniform Dirichlet
do(function(g) 10)   # try 'strong' Dirichlet





funa <- function(g) 2/max(1, g-5)
funb <-function(g) 1/g            # 0.042
func <-function(g) 1/(g + 2)
fund <-function(g) 1/sqrt(g)
fune <-function(g) 1/(max(3, g))
funf <-function(g) 1/min(g, 20)   # not bad esp late 0.026
fung <-function(g) 1/min(g, 30)  # not as good
funh <-function(g) 1 / (2 + (g/3))  # winner below
funi <-function(g) 1/(0.8 + 0.35 * max(g, 3))

allfuns<-function(g){
  c(funa(g),
  funb(g),
  func(g),
  fund(g),
  fune(g),
  funf(g),
  fung(g),
  funh(g),
  funi(g))
}

sapply(1:7,allfuns)

allfuns(1000)


set.seed(1)
n  <- 200
x1  <- rnorm(n)
x2 <- rnorm(n)
y <- 3 * x1 + 4 * x2 + 2 * rnorm(n)

z <- expand.grid(g=c(2:5, seq(10, 90, by=5)),
                 conc=c(0.01, 0.02, .03, .04, seq(0.05, 1, by=0.05)))
z <- data.frame(z, mad=NA)
for(i in 1 : nrow(z)) {
  g    <- z$g[i]
  conc <- z$conc[i]
  cat('conc=', conc, '\n')
  u <- cut2(y, g=g)
  f <- lrm(u ~ x1 + x2)
  b <- blrm(u ~ x1 + x2, conc=conc, priorsd=10000)
  u <- coef(f)[1 : (g-1)]
  v <- coef(b)[1 : (g-1)]
  d <- abs(u - v)
  z$mad[i] <- mean(d)
}

saveRDS(z, 'blrm-po.rds') 

ggplot(z, aes(x=g, y=conc, size=log(mad))) + geom_point()

gs <- unique(z$g)
bestconc <- numeric(length(gs))
i <- 0
for(k in gs) {
  i <- i + 1
  u <- subset(z, g==k)
  bestconc[i] <- with(u, conc[which.min(mad)])
}
plot(gs, bestconc)
plot(gs, 1/bestconc)
abline(lsfit(gs, 1/bestconc))
summary(lm(1/bestconc ~ gs))
summary(lm(1/bestconc ~ pol(gs, 2)))

## 1 / (2 + k/3)

1/bestconc[gs==2]
1/bestconc[gs==3]

plot(gs, 1/bestconc)
abline(lsfit(gs[gs > 2], 1/bestconc[gs > 2]))
lm(1/bestconc ~ gs, subset=gs > 2)
}
