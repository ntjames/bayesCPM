################# code for sims and applications in orm paper ################
library(rms)

################################### FUNCTIONS ################

#### function to estimate conditional mean and its standard error for orm model
mean.orm <- function(mod, new.data, se=TRUE){
  if(is.null(mod$yunique)) {
    stop("Need to set x=TRUE and y=TRUE for orm") 
  } else{
    order.y <- mod$yunique
    n.alpha <- length(order.y)-1
    xb <- as.matrix(new.data)%*%matrix(coef(mod)[colnames(new.data)])
    m.alpha <- mod$coef[1:n.alpha]
    lb <- t(outer(m.alpha, xb, "+")[,,1])
    m.s <- mod$trans$cumprob(lb)
    m.f <- t(apply(m.s, 1, FUN=function(x) c(1,x[1:n.alpha]) - c(x[1:n.alpha], 0)))
    m.mean <- apply(m.f, 1, FUN=function(x) sum(x*order.y))
    
    if(se){
      if(mod$family=="logistic") mod$trans$deriv <- function(x) exp(-x)/(1+exp(-x))^2
     
      dmean.dalpha <- t(apply(mod$trans$deriv(lb), 1, FUN=function(x) x*(order.y[2:length(order.y)] - order.y[1:n.alpha])))
      dmean.dbeta <- apply(dmean.dalpha, 1, sum)*as.matrix(new.data)
      dmean.dtheta <- cbind(dmean.dalpha, dmean.dbeta)   
      mean.var <-diag(dmean.dtheta%*%solve(mod$info.matrix)%*%t(dmean.dtheta))
      mean.se <- sqrt(mean.var)   
      result <- cbind(m.mean, mean.se)
      ci <- t(apply(result, 1, FUN=function(x) c(x[1]- qnorm(0.975)*x[2], x[1]+ qnorm(0.975)*x[2])))
      result <- cbind(result, ci)
      colnames(result) <- c("est", "se", "lb", "ub")
    } else{
      result <- matrix(m.mean)
      colnames(result) <- c("est")
    }
    
    
    return(result)
    
    
  } 
  
}


#### function to estimate conditional quantiles and confidence intervals for orm model
quantile.orm <- function(mod, new.data, probs=0.5, se=TRUE){
  
  quantile <- matrix(NA, nrow=dim(new.data)[1], ncol=length(probs))
  order.y <- mod$yunique
  #n.alpha <- length(order.y)-1
  xb <- as.matrix(new.data)%*%matrix(coef(mod)[colnames(new.data)])
  alpha <- mod$coef[1:(length(unique(order.y))-1)]
  lb <- t(outer(alpha, xb, "+")[,,1])
  m.cdf <- 1- mod$trans$cumprob(lb)
  m.cdf <- cbind(0, m.cdf, 1)
  for(i in 1: length(probs)){
    try({
      index.1 <- apply(m.cdf, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
      index.2 <- apply(m.cdf, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
      
      index.y1 <- ifelse(index.1>length(order.y), Inf, order.y[index.1])
      index.y2 <- ifelse(index.2>length(order.y), Inf, order.y[index.2])
      
      index.y1.cdf <- ifelse(index.1==0, 0, m.cdf[cbind(1:dim(new.data)[1], index.1)])
      index.y2.cdf <- ifelse(index.2>length(order.y), 1, m.cdf[cbind(1:dim(new.data)[1], index.2)])
      
      
      quantile[,i] <- ifelse(index.1==index.2, index.y1, 
                             (index.y2-index.y1)/(index.y2.cdf - index.y1.cdf)*(probs[i]-index.y1.cdf) + index.y1) 
      quantile[, i] <- ifelse(is.infinite(quantile[,i]), max(order.y), quantile[, i])
    })
    
  }
  result <- quantile
  
  if(se){
    if(mod$family=="logistic") mod$trans$deriv <- function(x) exp(-x)/(1+exp(-x))^2
    
    quantile.lb <- quantile.ub <- matrix(NA, nrow=dim(new.data)[1], ncol=length(probs))
    lb.se <- matrix(NA, ncol=dim(lb)[2], nrow=dim(new.data)[1])
    var <- as.matrix(solve(mod$info.matrix))
    
    for(i in 1:dim(lb)[2]){
      var.i <- var[c(i, which(names(coef(mod)) %in% colnames(new.data))), 
                   c(i, which(names(coef(mod)) %in% colnames(new.data)))]
      
      dcdf.dtheta <- cbind(-mod$trans$deriv(lb[,i]),  
                           -mod$trans$deriv(lb[,i])*as.matrix(new.data) )
      dlb.dtheta <- as.matrix(cbind(1, new.data))
      lb.se[,i] <- sqrt(diag(dlb.dtheta%*%var.i%*% t(dlb.dtheta)))
    }
    
    ci.lb <- sapply(1:dim(lb)[2], FUN=function(i) { 1- mod$trans$cumprob(lb[, i] +qnorm(0.975)*lb.se[, i])})
    ci.ub <- sapply(1:dim(lb)[2], FUN=function(i) { 1- mod$trans$cumprob(lb[, i] -qnorm(0.975)*lb.se[, i])})
    ci.lb <- matrix(ci.lb, nrow=dim(new.data)[1])
    ci.ub <- matrix(ci.ub, nrow=dim(new.data)[1])
    
    ci.lb <- cbind(0, ci.lb, 1)
    ci.ub <- cbind(0, ci.ub, 1)
    
    for(i in 1: length(probs)){
      try({
        index.1 <- apply(ci.lb, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
        index.2 <- apply(ci.lb, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
 
        index.y1 <- ifelse(index.1>length(order.y), Inf, order.y[index.1])
        index.y2 <- ifelse(index.2>length(order.y),Inf,order.y[index.2])
        
        index.y1.cdf <- ifelse(index.1==0, 0, ci.lb[cbind(1:dim(new.data)[1], index.1)])
        
        index.y2.cdf <- ifelse(index.2>length(order.y), 1, ci.lb[cbind(1:dim(new.data)[1], index.2)])
        
        
        quantile.lb[,i] <- ifelse(index.1==index.2, index.y1, 
                                  (index.y2-index.y1)/(index.y2.cdf - index.y1.cdf)*(probs[i]-index.y1.cdf) + index.y1) 
        quantile.lb[, i] <- ifelse(is.infinite(quantile.lb[,i]), max(order.y), quantile.lb[, i])
        
        index.1 <- apply(ci.ub, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
        index.2 <- apply(ci.ub, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
        
        index.y1 <- ifelse(index.1>length(order.y), Inf, order.y[index.1])
        index.y2 <- ifelse(index.2>length(order.y),Inf,order.y[index.2])
        
        index.y1.cdf <- ifelse(index.1==0, 0, ci.ub[cbind(1:dim(new.data)[1], index.1)])
        
        index.y2.cdf <- ifelse(index.2>length(order.y), 1, ci.ub[cbind(1:dim(new.data)[1], index.2)])
        
        
        quantile.ub[,i] <- ifelse(index.1==index.2, index.y1, 
                                  (index.y2-index.y1)/(index.y2.cdf - index.y1.cdf)*(probs[i]-index.y1.cdf) + index.y1) 
        quantile.ub[, i] <- ifelse(is.infinite(quantile.ub[,i]), max(order.y), quantile.ub[, i])
        
        
      })
      
    }
    
    result <- list(quantile=quantile,
                   lb=quantile.ub,
                   ub=quantile.lb)
    
    
    
  }
  
  
  
  return(result)
  
}


#### function to estimate conditional CDF and its standard error for orm models
cdf.orm <- function(mod, new.data, at.y=0,se=TRUE){
  if(is.null(mod$yunique)) {
    stop("Need to set x=TRUE and y=TRUE for orm") 
  } else{
    order.y <- mod$yunique
    xb <- as.matrix(new.data)%*%matrix(coef(mod)[colnames(new.data)])
    
    index <- sapply(at.y, FUN=function(x) {if(x<min(order.y)[1]) result <- Inf 
    else if (x==min(order.y)[1]) result <- 1
    else if(x >= max(order.y)[1]) result <- -Inf
    else which(order.y>=x)[1]-1})
    
    m.alpha <- mod$coef[index]
    m.alpha <- ifelse(is.infinite(index), index, m.alpha)
    if(length(at.y)==1){
      lb <- as.matrix(outer(m.alpha, xb, "+")[,,1])
    } else lb <- t(outer(m.alpha, xb, "+")[,,1])
    m.cdf <- 1- mod$trans$cumprob(lb)
    
    
    if(se){
      if(mod$family=="logistic") mod$trans$deriv <- function(x) exp(-x)/(1+exp(-x))^2
      cdf.se <- matrix(NA, ncol=length(at.y), nrow=dim(new.data)[1])
      lb.se <- matrix(NA, ncol=length(at.y), nrow=dim(new.data)[1])
      
      var <- as.matrix(solve(mod$info.matrix))
      
      for(i in 1:length(at.y)) {
        
        var.i <- var[c(index[i], which(names(coef(mod)) %in% colnames(new.data))), 
                     c(index[i], which(names(coef(mod)) %in% colnames(new.data)))]
        dcdf.dtheta <- cbind(-mod$trans$deriv(lb[,i]),  
                             -mod$trans$deriv(lb[,i])*as.matrix(new.data) )
        dlb.dtheta <- as.matrix(cbind(1, new.data))
        cdf.se[,i] <- sqrt(diag(dcdf.dtheta %*% var.i%*% t(dcdf.dtheta)))
        lb.se[, i] <- sqrt(diag(dlb.dtheta%*%var.i%*% t(dlb.dtheta)))
        
      }
      ci.lb <- sapply(1:length(at.y), FUN=function(i) { 1- mod$trans$cumprob(lb[, i] +qnorm(0.975)*lb.se[, i])})
      ci.ub <- sapply(1:length(at.y), FUN=function(i) { 1- mod$trans$cumprob(lb[, i] -qnorm(0.975)*lb.se[, i])})
      
   
      result <- list(est=m.cdf,
                     se=cdf.se,
                     lb=ci.lb,
                     ub=ci.ub)
    } else{
      result <- list(est=m.cdf)
    }
    
    
    return(result)
    
  } 
  
}


############################ SIMULATIONS ##############

###########Simulation 1: orm with properly specified link functions under  ####################

################### i) normal error  ###################
generate.data.1 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  log.y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=log.y, z1=z1, z2=z2)
  return(data)
}


######### Sim 1 (i) regression coefficients ################
sim1_coeff.fun <- function(sim=100,seed=1, n=50,
                       p=0.5, alpha=0, beta=c(1,-0.5), sigma=1,
                       log.y=c(-1, -0.33, 0.5, 1.33, 2), se=TRUE){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  alpha.y <- matrix(NA, ncol=length(log.y), nrow=sim)
  alpha.y.se <- matrix(NA, ncol=length(log.y), nrow=sim)
  beta.est <- matrix(NA, ncol=length(beta), nrow=sim)
  beta.se <- matrix(NA, ncol=length(beta), nrow=sim)
  
  for(i in 1:sim){
    try({
      data <- generate.data.1(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta, sigma=sigma)
      order.y <- data$y[order(data$y)]
      m.orm <- orm(y~z1+z2, data=data, family="probit")
      beta.est[i,] <- m.orm$coef[c(n, n+1)]
      beta.se[i, 1] <- sqrt(m.orm$var["z1", "z1"])
      beta.se[i, 2] <- sqrt(m.orm$var["z2", "z2"])
      alpha.y[i,] <- -sapply(log.y, FUN=function(x) {if(x<=min(data$y)) result <- Inf 
      else if(x > max(data$y)) result <- -Inf
      else m.orm$coef[which(order.y>=x)[1]-1]})
      
      if(se){
        var <- as.matrix(solve(m.orm$info.matrix))
        
        alpha.y.se[i,] <- sqrt(sapply(log.y, FUN=function(x) {if(x<=min(data$y)) result <- Inf 
        else if(x > max(data$y)) result <- Inf
        else var[which(order.y>=x)[1]-1, which(order.y>=x)[1]-1]}) )
      }
    })
  }
  
  return(list(beta.est=beta.est,
              beta.se=beta.se,
              alpha.y=alpha.y,
              alpha.y.se=alpha.y.se))
  
}

seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1
log.y <- c(-1, -0.33, 0.5, 1.33, 2)
sim=10^4
#sim=100

sim1_coeff.n25 <- sim1_coeff.fun(n=25, sim=sim)
sim1_coeff.n50 <- sim1_coeff.fun(n=50, sim=sim)
sim1_coeff.n100 <- sim1_coeff.fun(n=100, sim=sim)
sim1_coeff.n200 <- sim1_coeff.fun(n=200, sim=sim)
sim1_coeff.n500 <- sim1_coeff.fun(n=500, sim=sim)
sim1_coeff.n1000 <- sim1_coeff.fun(n=1000, sim=sim)


sim1_coeff.summary <- function(result, beta=c(1, -0.5),alpha.y=c(-1, -0.33, 0.5, 1.33, 2)){
  
  mean.est <- c(apply(result$beta.est, 2, FUN=function(x) mean(x, na.rm=TRUE)),
                apply(ifelse(is.infinite(result$alpha.y), NA, result$alpha.y), 2, FUN=function(x) mean(x, na.rm=TRUE)))
  
  mean.se <- c(apply(result$beta.se, 2, FUN=function(x) mean(x, na.rm=TRUE)),
               apply(ifelse(is.infinite(result$alpha.y.se), NA, result$alpha.y.se), 2, FUN=function(x) mean(x, na.rm=TRUE)))
  
  emp.se <- c(apply(result$beta.est, 2, FUN=function(x) sd(x, na.rm=TRUE)),
              apply(ifelse(is.infinite(result$alpha.y), NA, result$alpha.y), 2, FUN=function(x) sd(x, na.rm=TRUE)))
  
  mse <- c(apply((result$beta.est-matrix(beta, nrow=dim(result$beta.est)[1], ncol=length(beta),byrow=TRUE))^2, 2, FUN=function(x) mean(x, na.rm=TRUE)),
           apply((ifelse(is.infinite(result$alpha.y), NA, result$alpha.y)-matrix(alpha.y, nrow=dim(result$alpha.y)[1], ncol=length(alpha.y), byrow=TRUE))^2, 2 ,FUN=function(x) mean(x, na.rm=TRUE)))
  
  
  beta.lb <-result$beta.est+qnorm(0.025, 0, 1)*result$beta.se 
  beta.ub <- result$beta.est-qnorm(0.025, 0, 1)*result$beta.se
  
  alpha.lb <- cbind(result$alpha.y+qnorm(0.025, 0, 1)*result$alpha.y.se)
  alpha.ub <- cbind(result$alpha.y-qnorm(0.025, 0, 1)*result$alpha.y.se)
  alpha.lb <- ifelse(is.infinite(alpha.lb), NA, alpha.lb)
  alpha.ub <- ifelse(is.infinite(alpha.ub), NA, alpha.ub)
  
  alpha.coverage <- rep(NA, length(alpha.y))
  for(i in 1:length(alpha.y)){
    
    alpha.coverage[i] <-mean(apply(cbind(alpha.lb[,i], alpha.ub[,i]), 1, FUN=function(x) alpha.y[i]>=x[1] & alpha.y[i]<=x[2]), na.rm=TRUE)
    
  }
  
  beta.coverage <- rep(NA, length(beta))
  for(i in 1:length(beta)){
    
    beta.coverage[i] <-mean(apply(cbind(beta.lb[,i], beta.ub[,i]), 1, FUN=function(x) beta[i]>=x[1] & beta[i]<=x[2]), na.rm=TRUE)
    
  }
  
  coverage <- c(beta.coverage, alpha.coverage)
  
  result <- cbind(format(c(beta, alpha.y), nsmall=1),
                  format(round(mean.est, 2), nsmall=2),
                  format(round(mean.se, 3) , nsmall=3),
                  format(round(emp.se, 3) , nsmall=3),
                  format(round(mse, 4) , nsmall=3),
                  format(round(coverage, 3) , nsmall=3))
  return(result)
  
}

########    Table S.2; Figure 6  (i) error ~ Normal ################
sim1_coeff.probit.n25 <- sim1_coeff.summary(sim1_coeff.n25)
sim1_coeff.probit.n50 <- sim1_coeff.summary(sim1_coeff.n50)
sim1_coeff.probit.n100 <- sim1_coeff.summary(sim1_coeff.n100)
sim1_coeff.probit.n200 <- sim1_coeff.summary(sim1_coeff.n200)
sim1_coeff.probit.n500 <- sim1_coeff.summary(sim1_coeff.n500)
sim1_coeff.probit.n1000 <- sim1_coeff.summary(sim1_coeff.n1000)

######### Sim 1 (i) regression coefficients compared with Box-Cox transformation ######### ################

box.cox.trans <- function(x,lambda){
  if(lambda==0){
    x <- log(x)
  }else{
    x <- (x^lambda-1)/lambda
  }
  
  return(x)
}

generate.data.1 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  log.y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=exp(log.y), z1=z1, z2=z2)
  return(data)
}

sim1_boxcox.fun <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1, -0.5), sigma=1,
                         log.y=c(-1, -0.33, 0.5, 1.33, 2), se=TRUE){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  t.y <- matrix(NA, ncol=length(log.y), nrow=sim)
  
  beta.est <- matrix(NA, ncol=length(beta), nrow=sim)
  beta.se <- matrix(NA, ncol=length(beta), nrow=sim)
  beta.confint <- matrix(NA, ncol=4, nrow=sim)
  for(i in 1:sim){
    try({
      data <- generate.data.1(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta, sigma=sigma)
      
      b.m <- boxcox(y~z1+z2, data=data, plotit=F, interp = TRUE)
      lambda <-  b.m$x[which(b.m$y==max(b.m$y))][1]
      data$y.trans <- box.cox.trans(data$y, lambda)
      mod <- ols(y.trans~z1+z2, data=data)
      
      beta.est[i, ] <- mod$coef[c("z1", "z2")]
      beta.se[i, 1] <- sqrt(mod$var["z1", "z1"])
      beta.se[i, 2] <- sqrt(mod$var["z2", "z2"])
      beta.confint[i,c(1,2)] <- confint(mod)["z1",]
      beta.confint[i,c(3,4)] <- confint(mod)["z2",]
      
      t.y[i,] <- box.cox.trans(exp(log.y), lambda)      
      
    })
  }
  
  return(list(beta.est=beta.est,
              beta.se=beta.se,
              beta.confint=beta.confint,
              t.y=t.y))
  
}


seed=1
p <- 0.5
alpha <- 0
beta <- c(1,-0.5)
sigma <- 1
log.y=c(-1, -0.33, 0.5, 1.33, 2)

sim=10^4
sim1_boxcox.n25 <- sim1_boxcox.fun(n=25, sim=sim)
sim1_boxcox.n50 <- sim1_boxcox.fun(n=50, sim=sim)
sim1_boxcox.n100 <- sim1_boxcox.fun(n=100, sim=sim)
sim1_boxcox.n200 <- sim1_boxcox.fun(n=200, sim=sim)
sim1_boxcox.n500 <- sim1_boxcox.fun(n=500, sim=sim)
sim1_boxcox.n1000 <- sim1_boxcox.fun(n=1000, sim=sim)

###### relative efficiency: orm vs. box-cox
RE.summary <- function(result1, result2,beta=c(1, -0.5)){
  
  MSE1 <- apply(result1$beta.est -matrix(beta, nrow=dim(result1$beta.est)[1], ncol=length(beta),byrow=TRUE), 2, FUN=function(x) mean(x^2, na.rm=TRUE))
  
  MSE2 <- apply(result2$beta.est -matrix(beta, nrow=dim(result2$beta.est)[1], ncol=length(beta),byrow=TRUE), 2, FUN=function(x) mean(x^2, na.rm=TRUE))
  
  re <- MSE2/MSE1
  
  return(re)              
  
}

sim1_beta.re.n25 <- RE.summary(sim1_coeff.n25, sim1_boxcox.n25)
sim1_beta.re.n50 <- RE.summary(sim1_coeff.n50, sim1_boxcox.n50)
sim1_beta.re.n100 <- RE.summary(sim1_coeff.n100, sim1_boxcox.n100)
sim1_beta.re.n200 <- RE.summary(sim1_coeff.n200, sim1_boxcox.n200)
sim1_beta.re.n500 <- RE.summary(sim1_coeff.n500, sim1_boxcox.n500)
sim1_beta.re.n1000 <- RE.summary(sim1_coeff.n1000, sim1_boxcox.n1000)

##### Figure 7: relative  efficiency (a) 

sim1_re.probit.beta.1 <- c(sim1_beta.re.n25[1],
                           sim1_beta.re.n50[1],
                           sim1_beta.re.n100[1],
                           sim1_beta.re.n200[1],
                           sim1_beta.re.n500[1],
                           sim1_beta.re.n1000[1])
sim1_re.probit.beta.2 <- c(sim1_beta.re.n25[2],
                           sim1_beta.re.n50[2],
                           sim1_beta.re.n100[2],
                           sim1_beta.re.n200[2],
                           sim1_beta.re.n500[2],
                           sim1_beta.re.n1000[2])



####### Sim 1 (i) cdf and mean

generate.data.1 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  log.y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=exp(log.y), z1=z1, z2=z2)
  return(data)
}

sim1_cdfmean.fun <- function(sim=100,seed=1, n=50,
                       p=0.5, alpha=0, beta=c(1, -0.5), sigma=1,
                       log.y=c(-1, -0.33, 0.5, 1.33, 2), se=TRUE){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]  
  control.est <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.est <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  control.est.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.est.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  
  control.se <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.se <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  control.se.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.se.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  
  control.lb <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  control.ub <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.lb <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.ub <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  
  control.lb.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  control.ub.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.lb.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.ub.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  
  for(i in 1:sim){
    try({
      data <- generate.data.1(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta, sigma=sigma)
      #order.y <- data$y[order(data$y)]
      m.orm <- orm(y~z1+z2, data=data, family="probit", x=TRUE, y=TRUE)
      cdf.result <- cdf.orm(m.orm, new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1)), at.y=exp(log.y),se=se)
      mean.result <- data.frame(mean.orm(m.orm, new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1)), se=se))
      
      control.est[i, ] <- c(mean.result$est[1], cdf.result$est[1,])
      trt.est[i, ] <- c(mean.result$est[2], cdf.result$est[2,])
      control.est.1[i, ] <- c(mean.result$est[3], cdf.result$est[3,])
      trt.est.1[i, ] <- c(mean.result$est[4], cdf.result$est[4,])
      
      if(se){
        control.se[i, ] <- c(mean.result$se[1], cdf.result$se[1,])
        trt.se[i, ] <- c(mean.result$se[2], cdf.result$se[2,])
        
        control.se.1[i, ] <- c(mean.result$se[3], cdf.result$se[3,])
        trt.se.1[i, ] <- c(mean.result$se[4], cdf.result$se[4,])
        
        control.lb[i, ] <- c(mean.result$lb[1], cdf.result$lb[1,])
        control.ub[i, ] <- c(mean.result$ub[1], cdf.result$ub[1,])
        control.lb.1[i, ] <- c(mean.result$lb[3], cdf.result$lb[3,])
        control.ub.1[i, ] <- c(mean.result$ub[3], cdf.result$ub[3,])
        
        trt.lb[i, ] <- c(mean.result$lb[2], cdf.result$lb[2,])
        trt.ub[i, ] <- c(mean.result$ub[2], cdf.result$ub[2,])
        trt.lb.1[i, ] <- c(mean.result$lb[4], cdf.result$lb[4,])
        trt.ub.1[i, ] <- c(mean.result$ub[4], cdf.result$ub[4,])
      }
      
      
    })
  }
  
  return(list(control.est=control.est,
              trt.est=trt.est,
              control.se=control.se,
              trt.se=trt.se,
              control.lb=control.lb,
              trt.lb=trt.lb,
              control.ub=control.ub,
              trt.ub=trt.ub,
              control.est.1=control.est.1,
              trt.est.1=trt.est.1,
              control.se.1=control.se.1,
              trt.se.1=trt.se.1,
              control.lb.1=control.lb.1,
              trt.lb.1=trt.lb.1,
              control.ub.1=control.ub.1,
              trt.ub.1=trt.ub.1))
  
}



seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1
log.y=c(-1, -0.33, 0.5, 1.33, 2)

sim=10^4
sim1_cdfmean.n25 <- sim1_cdfmean.fun(n=25, sim=sim)
sim1_cdfmean.n50 <- sim1_cdfmean.fun(n=50, sim=sim)
sim1_cdfmean.n100 <- sim1_cdfmean.fun(n=100, sim=sim)
sim1_cdfmean.n200 <- sim1_cdfmean.fun(n=200, sim=sim)
sim1_cdfmean.n500 <- sim1_cdfmean.fun(n=500, sim=sim)
sim1_cdfmean.n1000 <- sim1_cdfmean.fun(n=1000, sim=sim)

sim1_cdfmean.summary <- function(result, beta=c(1, -0.5),alpha=0,sigma=1,alpha.y=c(-1, -0.33, 0.5, 1.33, 2), family=c("probit", "loglog")){
  if(family[1]=="probit"){
    true.cdf.01 <- pnorm(log.y, alpha+c(0,1)%*%beta, sigma)
    true.cdf.11 <- pnorm(log.y, alpha+c(1,1)%*%beta)
    
    true.mean.01 <- exp(alpha+c(0,1)%*%beta+sigma^2/2)
    true.mean.11 <- exp(alpha+c(1,1)%*%beta+sigma^2/2)
    true.mean.00 <- exp(alpha+c(0,0)%*%beta+sigma^2/2)
    true.mean.10 <- exp(alpha+c(1,0)%*%beta+sigma^2/2)
    true.mean <- c(true.mean.00, true.mean.10, true.mean.01, true.mean.11)
  } else if (family[1]=="loglog"){
    true.cdf.01 <- pGumbelMin(log.y, alpha+c(0,1)%*%beta, sigma)
    true.cdf.11 <- pGumbelMin(log.y, alpha+c(1,1)%*%beta, sigma)
    true.mean.01 <- mean(exp(rGumbelMin(10^6, alpha+c(0,1)%*%beta, sigma)))
    true.mean.11 <- mean(exp(rGumbelMin(10^6, alpha+c(1,1)%*%beta, sigma)))
    true.mean.00 <- mean(exp(rGumbelMin(10^6, alpha+c(0,0)%*%beta, sigma)))
    true.mean.10 <- mean(exp(rGumbelMin(10^6, alpha+c(1,0)%*%beta, sigma)))
    true.mean <- c(true.mean.00, true.mean.10, true.mean.01, true.mean.11)
    
  }
  
  
  mean.est <- c(mean(result$control.est[,1], na.rm=TRUE),
                mean(result$trt.est[,1], na.rm=TRUE),
                mean(result$control.est.1[,1], na.rm=TRUE),
                mean(result$trt.est.1[,1], na.rm=TRUE))
  
  mean.emp.se <-c(sd(result$control.est[,1], na.rm=TRUE),
                  sd(result$trt.est[,1], na.rm=TRUE),
                  sd(result$control.est.1[,1], na.rm=TRUE),
                  sd(result$trt.est.1[,1], na.rm=TRUE))
  
  mean.est.se <- c(mean(result$control.se[,1], na.rm=TRUE),
                   mean(result$trt.se[,1], na.rm=TRUE),
                   mean(result$control.se.1[,1], na.rm=TRUE),
                   mean(result$trt.se.1[,1], na.rm=TRUE))
  
  mean.coverage <- c(mean(true.mean[1] >=result$control.lb[, 1] & result$control.ub[,1]>=true.mean[1], na.rm=TRUE),
                     mean(true.mean[2] >=result$trt.lb[, 1] & result$trt.ub[,1]>=true.mean[2], na.rm=TRUE),
                     mean(true.mean[3] >=result$control.lb.1[, 1] & result$control.ub.1[,1]>=true.mean[3], na.rm=TRUE),
                     mean(true.mean[4] >=result$trt.lb.1[, 1] & result$trt.ub.1[,1]>=true.mean[4], na.rm=TRUE) )
  
 
  
  
  cdf.01.est <- apply(result$control.est.1[,-1], 2, FUN=function(x) mean(x, na.rm=TRUE))
  cdf.11.est <- apply(result$trt.est.1[,-1], 2, FUN=function(x) mean(x, na.rm=TRUE))
  
  cdf.01.est.se <- apply(result$control.est.1[,-1], 2, FUN=function(x) sd(x, na.rm=TRUE))
  cdf.11.est.se <- apply(result$trt.est.1[,-1], 2, FUN=function(x) sd(x, na.rm=TRUE))
  
  cdf.01.emp.se <- apply(result$control.se.1[,-1], 2, FUN=function(x) mean(x, na.rm=TRUE))
  cdf.11.emp.se <- apply(result$trt.se.1[,-1], 2, FUN=function(x) mean(x, na.rm=TRUE))
  
  cdf.01.coverage <- sapply(1:(length(log.y)), FUN=function(i) mean( true.cdf.01[i]>=result$control.lb.1[,i+1]&true.cdf.01[i]<=result$control.ub.1[, i+1], na.rm=TRUE))
  
  cdf.11.coverage <- sapply(1:(length(log.y)), FUN=function(i) mean( true.cdf.11[i]>=result$trt.lb.1[,i+1]&true.cdf.11[i]<=result$trt.ub.1[, i+1], na.rm=TRUE))
  
  
  mean.out <- cbind(format(round(true.mean, 2), nsmall=2),
                    format(round(mean.est, 3), nsmall=3),
                    format(round(mean.est.se, 3), nsmall=3),
                    format(round(mean.emp.se, 3), nsmall=3),
                    format(round(mean.coverage, 3), nsmall=3))
  rownames(mean.out) <- c("mean|00", "mean|10", "mean|01", "mean|11")
  
 
  
  cdf.01.out <- cbind(format(round(true.cdf.01, 4),4),
                      format(round(cdf.01.est, 4),4),
                      format(round(cdf.01.est.se, 4),4),
                      format(round(cdf.01.emp.se, 4),4),
                      format(round(cdf.01.coverage, 3),3)
  )
  
  cdf.11.out <- cbind(format(round(true.cdf.11, 4),4),
                      format(round(cdf.11.est, 4),4),
                      format(round(cdf.11.est.se, 4),4),
                      format(round(cdf.11.emp.se, 4),4),
                      format(round(cdf.11.coverage, 3),3)
  )
  result <- list(mean.out=mean.out,
                 cdf.01.out=cdf.01.out,
                 cdf.11.out=cdf.11.out)
  
  return(result)
  
}

sim1_cdfmean.probit.n25 <- sim1_cdfmean.summary(sim1_cdfmean.n25)
sim1_cdfmean.probit.n50 <- sim1_cdfmean.summary(sim1_cdfmean.n50)
sim1_cdfmean.probit.n100 <- sim1_cdfmean.summary(sim1_cdfmean.n100)
sim1_cdfmean.probit.n200 <- sim1_cdfmean.summary(sim1_cdfmean.n200)
sim1_cdfmean.probit.n500 <- sim1_cdfmean.summary(sim1_cdfmean.n500)
sim1_cdfmean.probit.n1000 <- sim1_cdfmean.summary(sim1_cdfmean.n1000)

###### Figure 8 (i) error ~ Normal, Table S.3 (i) error ~ Normal
sim1_cdf.probit.out <- rbind(rep("", 5),
                        rbind(sim1_cdfmean.probit.n25$cdf.11.out),
                        rep("", 5),
                        rbind(sim1_cdfmean.probit.n50$cdf.11.out),
                        rep("", 5),
                        rbind(sim1_cdfmean.probit.n100$cdf.11.out),
                        rep("", 5),
                        rbind(sim1_cdfmean.probit.n200$cdf.11.out),
                        rep("", 5),
                        rbind(sim1_cdfmean.probit.n500$cdf.11.out),
                        rep("", 5),
                        rbind(sim1_cdfmean.probit.n1000$cdf.11.out))

###### Figure S.1 (i) error ~ Normal, Table S.4 (i) error ~ Normal
sim1_mean.probit.out <- rbind(rep("", 5),
                              sim1_cdfmean.probit.n25$mean.out,
                              rep("", 5),
                              sim1_cdfmean.probit.n50$mean.out,
                              rep("", 5),
                              sim1_cdfmean.probit.n100$mean.out,
                              rep("", 5),
                              sim1_cdfmean.probit.n200$mean.out,
                              rep("", 5),
                              sim1_cdfmean.probit.n500$mean.out,
                              rep("", 5),
                              sim1_cdfmean.probit.n1000$mean.out)


########### Sim 1 (i) quantiles
generate.data.1 <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  log.y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(y=exp(log.y), z1=z1, z2=z2)
  return(data)
}


sim1_quantile.fun <- function(sim=100,seed=1, n=50,
                        p=0.5, alpha=0, beta=c(1, -0.5), sigma=1,
                        probs=c(0.1, 0.25, 0.5, 0.75, 0.9), se=TRUE){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  est.00 <- est.10 <- est.01<- est.11 <- matrix(NA, ncol=length(probs), nrow=sim)
  lb.00 <- lb.10 <- lb.01<- lb.11 <- matrix(NA, ncol=length(probs), nrow=sim)
  ub.00 <- ub.10 <- ub.01<- ub.11 <- matrix(NA, ncol=length(probs), nrow=sim)
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.1(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta, sigma=sigma)
      m.orm <- orm(y~z1+z2, data=data, family="probit", x=TRUE, y=TRUE)
      result <- quantile.orm(mod=m.orm, new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1)), probs=probs, se=TRUE)
      est.00[i, ] <- result$quantile[1,]
      est.10[i, ] <- result$quantile[2,]
      est.01[i, ] <- result$quantile[3,]
      est.11[i, ] <- result$quantile[4,]
      
      lb.00[i, ] <- result$lb[1,]
      lb.10[i, ] <- result$lb[2,]
      lb.01[i, ] <- result$lb[3,]
      lb.11[i, ] <- result$lb[4,]
      
      ub.00[i, ] <- result$ub[1,]
      ub.10[i, ] <- result$ub[2,]
      ub.01[i, ] <- result$ub[3,]
      ub.11[i, ] <- result$ub[4,]
      
      
    })
  }
  
  result <- list(result.00=list(est=est.00,
                                lb=lb.00,
                                ub=ub.00),
                 result.10=list(est=est.10,
                                lb=lb.10,
                                ub=ub.10),
                 result.01=list(est=est.01,
                                lb=lb.01,
                                ub=ub.01),
                 result.11=list(est=est.11,
                                lb=lb.11,
                                ub=ub.11))
  return(result)
}

seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1
probs=c(0.1, 0.25, 0.5, 0.75, 0.9)
se=TRUE

sim=10^4
sim1_quantile.n25 <- sim1_quantile.fun(n=25, sim=sim)
sim1_quantile.n50 <- sim1_quantile.fun(n=50, sim=sim)
sim1_quantile.n100 <- sim1_quantile.fun(n=100, sim=sim)
sim1_quantile.n200 <- sim1_quantile.fun(n=200, sim=sim)
sim1_quantile.n500 <- sim1_quantile.fun(n=500, sim=sim)
sim1_quantile.n1000 <- sim1_quantile.fun(n=1000, sim=sim)

sim1_quantile.summary <- function(result, true.quantile){
  
  true <- true.quantile
  est <- apply(result$est, 2, FUN=function(x) mean(x, na.rm=TRUE))
  emp.se <- apply(result$est, 2, FUN=function(x) sd(x, na.rm=TRUE))
  
  mse <- apply(result$est -matrix(true.quantile, ncol=length(true.quantile), nrow=dim(result$est)[1], byrow=T ),
               2, FUN=function(x) mean(x^2, na.rm=TRUE))
  
  cover <- sapply(1:length(true.quantile), FUN=function(i) true.quantile[i]>=result$lb[,i] & true.quantile[i]<=result$ub[,i]) 
  coverage <- apply(cover, 2, FUN=function(x) mean(x, na.rm=TRUE))
  
  
  result <- cbind(format(round(true, 3), nsmall=3),
                  format(round(est, 3), nsmall=3),
                  format(round(emp.se, 4), nsmall=4),
                  format(round(coverage, 3), nsmall=3)
  )
  
  
  
}

new.data <- rbind(c(0,1), c(1, 1))
true.mean <- beta %*% new.data

sim1_quantile.probit.n25 <- sim1_quantile.summary(result=sim1_quantile.n25$result.11, 
                                                  true.quantile=exp(qnorm(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1))) 

sim1_quantile.probit.n50 <- sim1_quantile.summary(result=sim1_quantile.n50$result.11, 
                                                  true.quantile=exp(qnorm(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))

sim1_quantile.probit.n100 <- sim1_quantile.summary(result=sim1_quantile.n100$result.11, 
                                                       true.quantile=exp(qnorm(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))

sim1_quantile.probit.n200 <- sim1_quantile.summary(result=sim1_quantile.n200$result.11, 
                                                       true.quantile=exp(qnorm(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))

sim1_quantile.probit.n500 <- sim1_quantile.summary(result=sim1_quantile.n500$result.11, 
                                                       true.quantile=exp(qnorm(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))

sim1_quantile.probit.n1000 <- sim1_quantile.summary(result=sim1_quantile.n1000$result.11, 
                                                       true.quantile=exp(qnorm(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))


###### Figure S.2 (i) error ~ Normal, Table S.5 (i) error ~ Normal

sim1_quantile.probit.out <- rbind(rep("", 4),
                                  rbind(sim1_quantile.probit.n25),
                                  rep("", 4),
                                  rbind(sim1_quantile.probit.n50),
                                  rep("", 4),
                                  rbind(sim1_quantile.probit.n100),
                                  rep("", 4),
                                  rbind(sim1_quantile.probit.n200),
                                  rep("", 4),
                                  rbind(sim1_quantile.probit.n500),
                                  rep("", 4),
                                  rbind(sim1_quantile.probit.n1000)
                                  )

               



##ii) extreme value error distribution

##### Sim 1 (ii) regression coefficients ###########

rGumbelMin<- function(n, mu=0, sigma=1){
  u <- runif(n, min=0, max=1)
  x <- mu + sigma*log(-log(1-u))
  return(x)
}



dGumbelMin <- function(x, mu=0, sigma=1){
  df <- 1/sigma *exp((x-mu)/sigma)*exp(-exp((x-mu)/sigma))
  return(df)
}

pGumbelMin <- function(x, mu=0, sigma=1){
  pf <- 1-exp(-exp((x-mu)/sigma))
  return(pf)
}

qGumbelMin <- function(p, mu=0, sigma=1){
  x <- mu + sigma*(-log(-log(1-p)))
  
}

generate.data.2 <- function(seed, n, p=0.5,alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  
  y <- rGumbelMin(n, mu=alpha+beta[1]*z1+beta[2]*z2, sigma=sigma)
  
  data <- data.frame(y, z1=z1, z2=z2)
  return(data) 
}




sim2_coeff.fun <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5), sigma=1,
                         log.y=c(-1, -0.33, 0.5, 1.33, 2), se=TRUE){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  alpha.y <- matrix(NA, ncol=length(log.y), nrow=sim)
  alpha.y.se <- matrix(NA, ncol=length(log.y), nrow=sim)
  beta.est <- matrix(NA, ncol=length(beta), nrow=sim)
  beta.se <- matrix(NA, ncol=length(beta), nrow=sim)
  
  for(i in 1:sim){
    try({
      data <- generate.data.2(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta, sigma=sigma)
      order.y <- data$y[order(data$y)]
      m.orm <- orm(y~z1+z2, data=data, family=loglog)
      beta.est[i,] <- m.orm$coef[c(n, n+1)]
      beta.se[i, 1] <- sqrt(m.orm$var["z1", "z1"])
      beta.se[i, 2] <- sqrt(m.orm$var["z2", "z2"])
      alpha.y[i,] <- -sapply(log.y, FUN=function(x) {if(x<=min(data$y)) result <- Inf 
      else if(x > max(data$y)) result <- -Inf
      else m.orm$coef[which(order.y>=x)[1]-1]})
      
      
      
      if(se){
        var <- as.matrix(solve(m.orm$info.matrix))
        
        alpha.y.se[i,] <- sqrt(sapply(log.y, FUN=function(x) {if(x<=min(data$y)) result <- Inf 
        else if(x > max(data$y)) result <- Inf
        else var[which(order.y>=x)[1]-1, which(order.y>=x)[1]-1]}) )
      }
    })
  }
  
  return(list(beta.est=beta.est,
              beta.se=beta.se,
              alpha.y=alpha.y,
              alpha.y.se=alpha.y.se))
  
}
seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1
log.y=c(-1, -0.33, 0.5, 1.33, 2)
sim=10^4
#sim=100

sim2_coeff.n25 <- sim2_coeff.fun(n=25, sim=sim)
sim2_coeff.n50 <- sim2_coeff.fun(n=50, sim=sim)
sim2_coeff.n100 <- sim2_coeff.fun(n=100, sim=sim)
sim2_coeff.n200 <- sim2_coeff.fun(n=200, sim=sim)
sim2_coeff.n500 <- sim2_coeff.fun(n=500, sim=sim)
sim2_coeff.n1000 <- sim2_coeff.fun(n=1000, sim=sim)

sim2_coeff.summary <- function(result, beta=c(1, -0.5),alpha.y=c(-1, -0.33, 0.5, 1.33, 2)){
  
  mean.est <- c(apply(result$beta.est, 2, FUN=function(x) mean(x, na.rm=TRUE)),
                apply(ifelse(is.infinite(result$alpha.y), NA, result$alpha.y), 2, FUN=function(x) mean(x, na.rm=TRUE)))
  
  mean.se <- c(apply(result$beta.se, 2, FUN=function(x) mean(x, na.rm=TRUE)),
               apply(ifelse(is.infinite(result$alpha.y.se), NA, result$alpha.y.se), 2, FUN=function(x) mean(x, na.rm=TRUE)))
  
  emp.se <- c(apply(result$beta.est, 2, FUN=function(x) sd(x, na.rm=TRUE)),
              apply(ifelse(is.infinite(result$alpha.y), NA, result$alpha.y), 2, FUN=function(x) sd(x, na.rm=TRUE)))
  
  mse <- c(apply((result$beta.est-matrix(beta, nrow=dim(result$beta.est)[1], ncol=length(beta),byrow=TRUE))^2, 2, FUN=function(x) mean(x, na.rm=TRUE)),
           apply((ifelse(is.infinite(result$alpha.y), NA, result$alpha.y)-matrix(alpha.y, nrow=dim(result$alpha.y)[1], ncol=length(alpha.y), byrow=TRUE))^2, 2 ,FUN=function(x) mean(x, na.rm=TRUE)))
  
  
  beta.lb <-result$beta.est+qnorm(0.025, 0, 1)*result$beta.se 
  beta.ub <- result$beta.est-qnorm(0.025, 0, 1)*result$beta.se
  
  alpha.lb <- cbind(result$alpha.y+qnorm(0.025, 0, 1)*result$alpha.y.se)
  alpha.ub <- cbind(result$alpha.y-qnorm(0.025, 0, 1)*result$alpha.y.se)
  alpha.lb <- ifelse(is.infinite(alpha.lb), NA, alpha.lb)
  alpha.ub <- ifelse(is.infinite(alpha.ub), NA, alpha.ub)
  
  alpha.coverage <- rep(NA, length(alpha.y))
  for(i in 1:length(alpha.y)){
    
    alpha.coverage[i] <-mean(apply(cbind(alpha.lb[,i], alpha.ub[,i]), 1, FUN=function(x) alpha.y[i]>=x[1] & alpha.y[i]<=x[2]), na.rm=TRUE)
    
  }
  
  beta.coverage <- rep(NA, length(beta))
  for(i in 1:length(beta)){
    
    beta.coverage[i] <-mean(apply(cbind(beta.lb[,i], beta.ub[,i]), 1, FUN=function(x) beta[i]>=x[1] & beta[i]<=x[2]), na.rm=TRUE)
    
  }
  
  coverage <- c(beta.coverage, alpha.coverage)
  
  result <- cbind(format(c(beta, alpha.y), nsmall=1),
                  format(round(mean.est, 2), nsmall=2),
                  format(round(mean.se, 3) , nsmall=3),
                  format(round(emp.se, 3) , nsmall=3),
                  format(round(mse, 4) , nsmall=3),
                  format(round(coverage, 3) , nsmall=3))
  return(result)
  
}

########  Table S.2; Figure 6  (ii) error ~ extreme Type I  ################
sim2_coeff.loglog.n25 <- sim2_coeff.summary(sim2_coeff.n25)
sim2_coeff.loglog.n50 <- sim2_coeff.summary(sim2_coeff.n50)
sim2_coeff.loglog.n100 <- sim2_coeff.summary(sim2_coeff.n100)
sim2_coeff.loglog.n200 <- sim2_coeff.summary(sim2_coeff.n200)
sim2_coeff.loglog.n500 <- sim2_coeff.summary(sim2_coeff.n500)
sim2_coeff.loglog.n1000 <- sim2_coeff.summary(sim2_coeff.n1000)



######## sim 1 (ii) regression coefficients compared with Cox proportional hazards model

sim2_cox.fun <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5), sigma=1,
                         log.y=c(-1, -0.33, 0.5, 1.33, 2)){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  alpha.y <- matrix(NA, ncol=length(log.y), nrow=sim)
  alpha.y.se <- matrix(NA, ncol=length(log.y), nrow=sim)
  beta.est <- matrix(NA, ncol=length(beta), nrow=sim)
  beta.se <- matrix(NA, ncol=length(beta), nrow=sim)
  
  for(i in 1:sim){
    try({
      data <- generate.data.2(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta, sigma=sigma)
      order.y <- data$y[order(data$y)]
      index <- sapply(log.y, FUN=function(x) {if(x<=min(order.y)[1]) result <- Inf 
      else if(x > max(order.y)[1]) result <- -Inf
      else which(order.y>=x)[1]-1})
      m.cox <- coxph(Surv(y, rep(1, n))~z1+z2, data=data)
      beta.est[i,] <- -m.cox$coef
      beta.se[i, 1] <- sqrt(m.cox$var[1, 1])
      beta.se[i, 2] <- sqrt(m.cox$var[2, 2])
      chz <- basehaz(m.cox, centered=FALSE)[,1]
      est.F <- 1-exp(-chz)
      
      alpha.y[i,] <- -qGumbelMin(est.F[index])
      
      
      
    })
  }
  
  return(list(beta.est=beta.est,
              beta.se=beta.se,
              alpha.y=alpha.y,
              alpha.y.se=alpha.y.se))
  
}


seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1
log.y=c(-1, -0.33, 0.5, 1.33, 2)

sim=10^4
sim2_cox.n25 <- sim2_cox.fun(n=25, sim=sim)
sim2_cox.n50 <- sim2_cox.fun(n=50, sim=sim)
sim2_cox.n100 <- sim2_cox.fun(n=100, sim=sim)
sim2_cox.n200 <- sim2_cox.fun(n=200, sim=sim)
sim2_cox.n500 <- sim2_cox.fun(n=500, sim=sim)
sim2_cox.n1000 <- sim2_cox.fun(n=1000, sim=sim)

RE.summary <- function(result1, result2,beta=c(1, -0.5)){
  
  MSE1 <- apply(result1$beta.est -matrix(beta, nrow=dim(result1$beta.est)[1], ncol=length(beta),byrow=TRUE), 2, FUN=function(x) mean(x^2, na.rm=TRUE))
  
  MSE2 <- apply(result2$beta.est -matrix(beta, nrow=dim(result2$beta.est)[1], ncol=length(beta),byrow=TRUE), 2, FUN=function(x) mean(x^2, na.rm=TRUE))
  
  re <- MSE2/MSE1
  
  return(re)              
  
}

sim2_beta.re.n25 <- RE.summary(sim2_coeff.n25, sim2_cox.n25)
sim2_beta.re.n50 <- RE.summary(sim2_coeff.n50, sim2_cox.n50)
sim2_beta.re.n100 <- RE.summary(sim2_coeff.n100, sim2_cox.n100)
sim2_beta.re.n200 <- RE.summary(sim2_coeff.n200, sim2_cox.n200)
sim2_beta.re.n500 <- RE.summary(sim1_coeff.n500, sim2_cox.n500)
sim2_beta.re.n1000 <- RE.summary(sim1_coeff.n1000, sim2_cox.n1000)

##### Figure 7: relative  efficiency (b) 

sim2_re.loglog.beta.1 <- c(sim2_beta.re.n25[1],
                           sim2_beta.re.n50[1],
                           sim2_beta.re.n100[1],
                           sim2_beta.re.n200[1],
                           sim2_beta.re.n500[1],
                           sim2_beta.re.n1000[1])
sim2_re.loglog.beta.2 <- c(sim2_beta.re.n25[2],
                           sim2_beta.re.n50[2],
                           sim2_beta.re.n100[2],
                           sim2_beta.re.n200[2],
                           sim2_beta.re.n500[2],
                           sim2_beta.re.n1000[2])


############ Sim 1 (ii): cdf and mean ##########

generate.data.2 <- function(seed, n, p=0.5,alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  
  log.y <- rGumbelMin(n, mu=alpha+beta[1]*z1+beta[2]*z2, sigma=sigma)
  
  data <- data.frame(y=exp(log.y), z1=z1, z2=z2)
  return(data) 
}

sim2_cdfmean.fun <- function(sim=100,seed=1, n=50,
                       p=0.5, alpha=0, beta=c(1, -0.5), sigma=1,
                       log.y=c(-1, -0.33, 0.5, 1.33, 2), se=TRUE){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]  
  control.est <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.est <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  control.est.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.est.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  
  control.se <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.se <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  control.se.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.se.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  
  control.lb <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  control.ub <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.lb <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.ub <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  
  control.lb.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  control.ub.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.lb.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  trt.ub.1 <- matrix(NA, ncol=length(log.y)+1, nrow=sim)
  

  for(i in 1:sim){
    try({
      data <- generate.data.2(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta, sigma=sigma)
      
      m.orm <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      cdf.result <- cdf.orm(m.orm, new.data=data.frame(z1=c( 0, 1), z2=c( 1, 1)), at.y=exp(log.y),se=se)
      mean.result <- data.frame(mean.orm(m.orm, new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1)), se=se))
      
      control.est[i, ] <- c(mean.result$est[1])
      trt.est[i, ] <- c(mean.result$est[2])
      control.est.1[i, ] <- c(mean.result$est[3], cdf.result$est[1,])
      trt.est.1[i, ] <- c(mean.result$est[4], cdf.result$est[2,])

      if(se){
        control.se[i, ] <- c(mean.result$se[1])
        trt.se[i, ] <- c(mean.result$se[2])
        
        control.se.1[i, ] <- c(mean.result$se[3], cdf.result$se[1,])
        trt.se.1[i, ] <- c(mean.result$se[4], cdf.result$se[2,])
        
        control.lb[i, ] <- c(mean.result$lb[1])
        control.ub[i, ] <- c(mean.result$ub[1])
        control.lb.1[i, ] <- c(mean.result$lb[3], cdf.result$lb[1,])
        control.ub.1[i, ] <- c(mean.result$ub[3], cdf.result$ub[1,])
        
        trt.lb[i, ] <- c(mean.result$lb[2])
        trt.ub[i, ] <- c(mean.result$ub[2])
        trt.lb.1[i, ] <- c(mean.result$lb[4], cdf.result$lb[2,])
        trt.ub.1[i, ] <- c(mean.result$ub[4], cdf.result$ub[2,])
      }
      
      
    })
  }
  
  return(list(control.est=control.est,
              trt.est=trt.est,
              control.se=control.se,
              trt.se=trt.se,
              control.lb=control.lb,
              trt.lb=trt.lb,
              control.ub=control.ub,
              trt.ub=trt.ub,
              control.est.1=control.est.1,
              trt.est.1=trt.est.1,
              control.se.1=control.se.1,
              trt.se.1=trt.se.1,
              control.lb.1=control.lb.1,
              trt.lb.1=trt.lb.1,
              control.ub.1=control.ub.1,
              trt.ub.1=trt.ub.1))
  
}



seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1
log.y=c(-1, -0.33, 0.5, 1.33, 2)

sim=10^4
sim2_cdfmean.n25 <- sim2_cdfmean.fun(n=25, sim=sim)
sim2_cdfmean.n50 <- sim2_cdfmean.fun(n=50, sim=sim)
sim2_cdfmean.n100 <- sim2_cdfmean.fun(n=100, sim=sim)
sim2_cdfmean.n200 <- sim2_cdfmean.fun(n=200, sim=sim)
sim2_cdfmean.n500 <- sim2_cdfmean.fun(n=500, sim=sim)
sim2_cdfmean.n1000 <- sim2_cdfmean.fun(n=1000, sim=sim)


sim2_cdfmean.summary <- function(result, beta=c(1, -0.5),alpha=0,sigma=1,alpha.y=c(-1, -0.33, 0.5, 1.33, 2), family=c("probit", "loglog")){
  if(family[1]=="probit"){
    true.cdf.01 <- pnorm(log.y, alpha+c(0,1)%*%beta, sigma)
    true.cdf.11 <- pnorm(log.y, alpha+c(1,1)%*%beta)
    
    true.mean.01 <- exp(alpha+c(0,1)%*%beta+sigma^2/2)
    true.mean.11 <- exp(alpha+c(1,1)%*%beta+sigma^2/2)
    true.mean.00 <- exp(alpha+c(0,0)%*%beta+sigma^2/2)
    true.mean.10 <- exp(alpha+c(1,0)%*%beta+sigma^2/2)
    true.mean <- c(true.mean.00, true.mean.10, true.mean.01, true.mean.11)
  } else if (family[1]=="loglog"){
    true.cdf.01 <- pGumbelMin(log.y, alpha+c(0,1)%*%beta, sigma)
    true.cdf.11 <- pGumbelMin(log.y, alpha+c(1,1)%*%beta, sigma)
    true.mean.01 <- mean(exp(rGumbelMin(10^6, alpha+c(0,1)%*%beta, sigma)))
    true.mean.11 <- mean(exp(rGumbelMin(10^6, alpha+c(1,1)%*%beta, sigma)))
    true.mean.00 <- mean(exp(rGumbelMin(10^6, alpha+c(0,0)%*%beta, sigma)))
    true.mean.10 <- mean(exp(rGumbelMin(10^6, alpha+c(1,0)%*%beta, sigma)))
    true.mean <- c(true.mean.00, true.mean.10, true.mean.01, true.mean.11)
    
  }
  
  
  mean.est <- c(mean(result$control.est[,1], na.rm=TRUE),
                mean(result$trt.est[,1], na.rm=TRUE),
                mean(result$control.est.1[,1], na.rm=TRUE),
                mean(result$trt.est.1[,1], na.rm=TRUE))
  
  mean.emp.se <-c(sd(result$control.est[,1], na.rm=TRUE),
                  sd(result$trt.est[,1], na.rm=TRUE),
                  sd(result$control.est.1[,1], na.rm=TRUE),
                  sd(result$trt.est.1[,1], na.rm=TRUE))
  
  mean.est.se <- c(mean(result$control.se[,1], na.rm=TRUE),
                   mean(result$trt.se[,1], na.rm=TRUE),
                   mean(result$control.se.1[,1], na.rm=TRUE),
                   mean(result$trt.se.1[,1], na.rm=TRUE))
  
  mean.coverage <- c(mean(true.mean[1] >=result$control.lb[, 1] & result$control.ub[,1]>=true.mean[1], na.rm=TRUE),
                     mean(true.mean[2] >=result$trt.lb[, 1] & result$trt.ub[,1]>=true.mean[2], na.rm=TRUE),
                     mean(true.mean[3] >=result$control.lb.1[, 1] & result$control.ub.1[,1]>=true.mean[3], na.rm=TRUE),
                     mean(true.mean[4] >=result$trt.lb.1[, 1] & result$trt.ub.1[,1]>=true.mean[4], na.rm=TRUE) )
  
  
  
  
  cdf.01.est <- apply(result$control.est.1[,-1], 2, FUN=function(x) mean(x, na.rm=TRUE))
  cdf.11.est <- apply(result$trt.est.1[,-1], 2, FUN=function(x) mean(x, na.rm=TRUE))
  
  cdf.01.est.se <- apply(result$control.est.1[,-1], 2, FUN=function(x) sd(x, na.rm=TRUE))
  cdf.11.est.se <- apply(result$trt.est.1[,-1], 2, FUN=function(x) sd(x, na.rm=TRUE))
  
  cdf.01.emp.se <- apply(result$control.se.1[,-1], 2, FUN=function(x) mean(x, na.rm=TRUE))
  cdf.11.emp.se <- apply(result$trt.se.1[,-1], 2, FUN=function(x) mean(x, na.rm=TRUE))
  
  cdf.01.coverage <- sapply(1:(length(log.y)), FUN=function(i) mean( true.cdf.01[i]>=result$control.lb.1[,i+1]&true.cdf.01[i]<=result$control.ub.1[, i+1], na.rm=TRUE))
  
  cdf.11.coverage <- sapply(1:(length(log.y)), FUN=function(i) mean( true.cdf.11[i]>=result$trt.lb.1[,i+1]&true.cdf.11[i]<=result$trt.ub.1[, i+1], na.rm=TRUE))
  
  
  mean.out <- cbind(format(round(true.mean, 2), nsmall=2),
                    format(round(mean.est, 3), nsmall=3),
                    format(round(mean.est.se, 3), nsmall=3),
                    format(round(mean.emp.se, 3), nsmall=3),
                    format(round(mean.coverage, 3), nsmall=3))
  rownames(mean.out) <- c("mean|00", "mean|10", "mean|01", "mean|11")
  
  
  
  cdf.01.out <- cbind(format(round(true.cdf.01, 4),4),
                      format(round(cdf.01.est, 4),4),
                      format(round(cdf.01.est.se, 4),4),
                      format(round(cdf.01.emp.se, 4),4),
                      format(round(cdf.01.coverage, 3),3)
  )
  
  cdf.11.out <- cbind(format(round(true.cdf.11, 4),4),
                      format(round(cdf.11.est, 4),4),
                      format(round(cdf.11.est.se, 4),4),
                      format(round(cdf.11.emp.se, 4),4),
                      format(round(cdf.11.coverage, 3),3)
  )
  result <- list(mean.out=mean.out,
                 cdf.01.out=cdf.01.out,
                 cdf.11.out=cdf.11.out)
  
  return(result)
  
}

sim2_cdfmean.loglog.n25 <- sim1_cdfmean.summary(sim2_cdfmean.n25)
sim2_cdfmean.loglog.n50 <- sim1_cdfmean.summary(sim2_cdfmean.n50)
sim2_cdfmean.loglog.n100 <- sim1_cdfmean.summary(sim2_cdfmean.n100)
sim2_cdfmean.loglog.n200 <- sim1_cdfmean.summary(sim2_cdfmean.n200)
sim2_cdfmean.loglog.n500 <- sim1_cdfmean.summary(sim2_cdfmean.n500)
sim2_cdfmean.loglog.n1000 <- sim1_cdfmean.summary(sim2_cdfmean.n1000)

###### Figure 8 (ii) error ~ extreme value (type I), Table S.3 (ii) error ~ extreme value (type I)
sim2_cdf.loglog.out <- rbind(rep("", 5),
                             rbind(sim2_cdfmean.loglog.n25$cdf.11.out),
                             rep("", 5),
                             rbind(sim2_cdfmean.loglog.n50$cdf.11.out),
                             rep("", 5),
                             rbind(sim2_cdfmean.loglog.n100$cdf.11.out),
                             rep("", 5),
                             rbind(sim2_cdfmean.loglog.n200$cdf.11.out),
                             rep("", 5),
                             rbind(sim2_cdfmean.loglog.n500$cdf.11.out),
                             rep("", 5),
                             rbind(sim2_cdfmean.loglog.n1000$cdf.11.out))


###### Figure S.1 (ii) error ~ extreme value (type I), Table S.4 (ii) error ~ extreme value (type I)
sim2_mean.loglog.out <- rbind(rep("", 5),
                              sim2_cdfmean.loglog.n25$mean.out,
                              rep("", 5),
                              sim2_cdfmean.loglog.n50$mean.out,
                              rep("", 5),
                              sim2_cdfmean.loglog.n100$mean.out,
                              rep("", 5),
                              sim2_cdfmean.loglog.n200$mean.out,
                              rep("", 5),
                              sim2_cdfmean.loglog.n500$mean.out,
                              rep("", 5),
                              sim2_cdfmean.loglog.n1000$mean.out)

############# Sim 1(ii): quantiles  ############
generate.data.2 <- function(seed, n, p=0.5,alpha=0, beta=c(1, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  
  log.y <- rGumbelMin(n, mu=alpha+beta[1]*z1+beta[2]*z2, sigma=sigma)
  
  data <- data.frame(y=exp(log.y), z1=z1, z2=z2)
  return(data) 
}

sim2_quantile.fun <- function(sim=100,seed=1, n=50,
                        p=0.5, alpha=0, beta=c(1, -0.5), sigma=1,
                        probs=c(0.1, 0.25, 0.5, 0.75, 0.9), se=TRUE){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  est.00 <- est.10 <- est.01<- est.11 <- matrix(NA, ncol=length(probs), nrow=sim)
  lb.00 <- lb.10 <- lb.01<- lb.11 <- matrix(NA, ncol=length(probs), nrow=sim)
  ub.00 <- ub.10 <- ub.01<- ub.11 <- matrix(NA, ncol=length(probs), nrow=sim)
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.2(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta, sigma=sigma)
      m.orm <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      result <- quantile.orm(mod=m.orm, new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1)), probs=probs, se=TRUE)
      est.00[i, ] <- result$quantile[1,]
      est.10[i, ] <- result$quantile[2,]
      est.01[i, ] <- result$quantile[3,]
      est.11[i, ] <- result$quantile[4,]
      
      lb.00[i, ] <- result$lb[1,]
      lb.10[i, ] <- result$lb[2,]
      lb.01[i, ] <- result$lb[3,]
      lb.11[i, ] <- result$lb[4,]
      
      ub.00[i, ] <- result$ub[1,]
      ub.10[i, ] <- result$ub[2,]
      ub.01[i, ] <- result$ub[3,]
      ub.11[i, ] <- result$ub[4,]
      
      
    })
  }
  
  result <- list(result.00=list(est=est.00,
                                lb=lb.00,
                                ub=ub.00),
                 result.10=list(est=est.10,
                                lb=lb.10,
                                ub=ub.10),
                 result.01=list(est=est.01,
                                lb=lb.01,
                                ub=ub.01),
                 result.11=list(est=est.11,
                                lb=lb.11,
                                ub=ub.11))
  return(result)
}

seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1
probs=c(0.1, 0.25, 0.5, 0.75, 0.9)
se=TRUE
sim=10^4
#n=100
sim2_quantile.n25 <- sim2_quantile.fun(n=25, sim=sim)
sim2_quantile.n50 <- sim2_quantile.fun(n=50, sim=sim)
sim2_quantile.n100 <- sim2_quantile.fun(n=100, sim=sim)
sim2_quantile.n200 <- sim2_quantile.fun(n=200, sim=sim)
sim2_quantile.n500 <- sim2_quantile.fun(n=500, sim=sim)
sim2_quantile.n1000 <- sim2_quantile.fun(n=1000, sim=sim)

sim2_quantile.summary <- function(result, true.quantile){
  
  true <- true.quantile
  est <- apply(result$est, 2, FUN=function(x) mean(x, na.rm=TRUE))
  emp.se <- apply(result$est, 2, FUN=function(x) sd(x, na.rm=TRUE))
  
  mse <- apply(result$est -matrix(true.quantile, ncol=length(true.quantile), nrow=dim(result$est)[1], byrow=T ),
               2, FUN=function(x) mean(x^2, na.rm=TRUE))
  
  cover <- sapply(1:length(true.quantile), FUN=function(i) true.quantile[i]>=result$lb[,i] & true.quantile[i]<=result$ub[,i]) 
  coverage <- apply(cover, 2, FUN=function(x) mean(x, na.rm=TRUE))
  
  
  result <- cbind(format(round(true, 3), nsmall=3),
                  format(round(est, 3), nsmall=3),
                  format(round(emp.se, 4), nsmall=4),
                  format(round(coverage, 3), nsmall=3)
  )
  
  
  
}

new.data <- rbind(c(0,1), c(1, 1))
true.mean <- beta %*% new.data
qGumbelMin <- function(p, mean=0, sd=1){
  x <- mu + sigma*(log(-log(1-p)))
  
}

sim2_quantile.loglog.n25 <- sim2_quantile.summary(result=sim2_quantile.n25.$result.11, 
                                                  true.quantile=exp(GumbelMin(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1))) 

sim2_quantile.loglog.n50 <- sim2_quantile.summary(result=sim2_quantile.n50$result.11, 
                                                  true.quantile=exp(GumbelMin(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))

sim2_quantile.loglog.n100 <- sim2_quantile.summary(result=sim2_quantile.n100$result.11, 
                                                   true.quantile=exp(GumbelMin(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))

sim2_quantile.loglog.n200 <- sim2_quantile.summary(result=sim2_quantile.n200$result.11, 
                                                   true.quantile=exp(GumbelMin(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))

sim2_quantile.loglog.n500 <- sim2_quantile.summary(result=sim2_quantile.n500$result.11, 
                                                   true.quantile=exp(GumbelMin(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))

sim2_quantile.loglog.n1000 <- sim2_quantile.summary(result=sim2_quantile.n1000.$result.11, 
                                                    true.quantile=exp(GumbelMin(c(0.1, 0.25, 0.5, 0.75, 0.9), mean=true.mean[2], sd=1)))


###### Figure S.2 (ii) error ~ extreme value (Type I), Table S.5 (ii) error ~ extreme value (Type I)

sim2_quantile.loglog.out <- rbind(rep("", 4),
                                  rbind(sim2_quantile.loglog.n25),
                                  rep("", 4),
                                  rbind(sim2_quantile.loglog.n50),
                                  rep("", 4),
                                  rbind(sim2_quantile.loglog.n100),
                                  rep("", 4),
                                  rbind(sim2_quantile.loglog.n200),
                                  rep("", 4),
                                  rbind(sim2_quantile.loglog.n500),
                                  rep("", 4),
                                  rbind(sim2_quantile.loglog.n1000)
)


######################### sim 2: with misspecified link functions   ##############
rGumbelMin<- function(n, mu=0, sigma=1){
  u <- runif(n, min=0, max=1)
  x <- mu + sigma*log(-log(1-u))
  return(x)
}

rGumbelMax<- function(n, mu=0, sigma=1){
  u <- runif(n, min=0, max=1)
  x <- mu + sigma*(-log(-log(u)))
  return(x)
}


generate.data.a <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, 1)
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}


generate.data.b <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  y <-  rlogis(n, alpha+beta[1]*z1 + beta[2]*z2, 1)
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}


generate.data.c <- function(seed, n, p=0.5,alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  error <-rGumbelMin(n, mu=0, sigma=1)-digamma(1)
  y <- alpha+beta[1]*z1+beta[2]*z2 +error
  
  data <- data.frame(y, z1=z1, z2=z2)
  return(data) 
}




generate.data.d <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #error <- rgamma(n, 1, 1)-1
  error <- rGumbelMax(n, 0, 1)+digamma(1)
  y <- alpha+beta[1]*z1+beta[2]*z2 +error
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}


generate.data.e <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  y <- alpha+beta[1]*z1 + beta[2]*z2 + rt(n, df=5)
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}

generate.data.f <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #y.star <- rnorm(n, alpha+beta*z, sigma)
  #y <- exp(y.star)
  y <- alpha+beta[1]*z1 + beta[2]*z2 + runif(n, -5, 5)
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}


generate.data.g <- function(seed, n, p=0.5,alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  error <-(rbeta(n, 5, 2) -5/(5+2))*sqrt((5+2)^2*(5+2+1)/5/2)
  y <- alpha+beta[1]*z1+beta[2]*z2 +error
  
  data <- data.frame(y, z1=z1, z2=z2)
  return(data) 
}

generate.data.h <- function(seed=1, n=200, p=0.5, alpha=0, beta=c(1, -0.5)){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  #error <- rgamma(n, 1, 1)-1
  error <-(rbeta(n, 2, 5) -2/(5+2))*sqrt((5+2)^2*(5+2+1)/5/2)
  y <- alpha+beta[1]*z1+beta[2]*z2 +error
  data <- data.frame(y=y, z1=z1, z2=z2)
  return(data)
}

##### Sim 2: means ##########

sim3_mean.fun.a <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=4, nrow=sim)
  probit.10 <- matrix(NA, ncol=4, nrow=sim)
  probit.01 <- matrix(NA, ncol=4, nrow=sim)
  probit.11 <- matrix(NA, ncol=4, nrow=sim)
  colnames(probit.00) <- c("est", "se", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=4, nrow=sim)
  logit.10 <- matrix(NA, ncol=4, nrow=sim)
  logit.01 <- matrix(NA, ncol=4, nrow=sim)
  logit.11 <- matrix(NA, ncol=4, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=4, nrow=sim)
  loglog.10 <- matrix(NA, ncol=4, nrow=sim)
  loglog.01 <- matrix(NA, ncol=4, nrow=sim)
  loglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  ols.00 <- matrix(NA, ncol=4, nrow=sim)
  ols.10 <- matrix(NA, ncol=4, nrow=sim)
  ols.01 <- matrix(NA, ncol=4, nrow=sim)
  ols.11 <- matrix(NA, ncol=4, nrow=sim)
  
  for(i in 1:sim){
    try({
      data <- generate.data.a(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      ols.mod <- ols(y~ z1+z2, data=data)
      
      mean.probit <- mean.orm(orm.probit, new.data, se=TRUE)
      mean.logit <- mean.orm(orm.logit, new.data, se=TRUE)
      mean.loglog <- mean.orm(orm.loglog, new.data, se=TRUE)
      mean.cloglog <- mean.orm(orm.cloglog, new.data, se=TRUE)
      
      mean.ols <- predict(ols.mod, new.data, se.fit=TRUE, conf.int=0.95)
      mean.ols <- cbind(mean.ols$linear.predictors, mean.ols$se.fit, mean.ols$lower, mean.ols$upper)
      
      probit.00[i,] <- mean.probit[1,]
      probit.10[i,] <- mean.probit[2,]
      probit.01[i,] <- mean.probit[3,]
      probit.11[i,] <- mean.probit[4,]
      
      logit.00[i,] <- mean.logit[1,]
      logit.10[i,] <- mean.logit[2,]
      logit.01[i,] <- mean.logit[3,]
      logit.11[i,] <- mean.logit[4,]
      
      loglog.00[i,] <- mean.loglog[1,]
      loglog.10[i,] <- mean.loglog[2,]
      loglog.01[i,] <- mean.loglog[3,]
      loglog.11[i,] <- mean.loglog[4,]
      
      cloglog.00[i,] <- mean.cloglog[1,]
      cloglog.10[i,] <- mean.cloglog[2,]
      cloglog.01[i,] <- mean.cloglog[3,]
      cloglog.11[i,] <- mean.cloglog[4,]
      
      ols.00[i,] <- mean.ols[1,]
      ols.10[i,] <- mean.ols[2,]
      ols.01[i,] <- mean.ols[3,]
      ols.11[i,] <- mean.ols[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              ols.00=ols.00,
              ols.10=ols.10,
              ols.01=ols.01,
              ols.11=ols.11
  ))
  
}


sim3_mean.fun.b <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=4, nrow=sim)
  probit.10 <- matrix(NA, ncol=4, nrow=sim)
  probit.01 <- matrix(NA, ncol=4, nrow=sim)
  probit.11 <- matrix(NA, ncol=4, nrow=sim)
  colnames(probit.00) <- c("est", "se", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=4, nrow=sim)
  logit.10 <- matrix(NA, ncol=4, nrow=sim)
  logit.01 <- matrix(NA, ncol=4, nrow=sim)
  logit.11 <- matrix(NA, ncol=4, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=4, nrow=sim)
  loglog.10 <- matrix(NA, ncol=4, nrow=sim)
  loglog.01 <- matrix(NA, ncol=4, nrow=sim)
  loglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  ols.00 <- matrix(NA, ncol=4, nrow=sim)
  ols.10 <- matrix(NA, ncol=4, nrow=sim)
  ols.01 <- matrix(NA, ncol=4, nrow=sim)
  ols.11 <- matrix(NA, ncol=4, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.b(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      ols.mod <- ols(y~ z1+z2, data=data)
      
      mean.probit <- mean.orm(orm.probit, new.data, se=TRUE)
      mean.logit <- mean.orm(orm.logit, new.data, se=TRUE)
      mean.loglog <- mean.orm(orm.loglog, new.data, se=TRUE)
      mean.cloglog <- mean.orm(orm.cloglog, new.data, se=TRUE)
      
      mean.ols <- predict(ols.mod, new.data, se.fit=TRUE, conf.int=0.95)
      mean.ols <- cbind(mean.ols$linear.predictors, mean.ols$se.fit, mean.ols$lower, mean.ols$upper)
      
      probit.00[i,] <- mean.probit[1,]
      probit.10[i,] <- mean.probit[2,]
      probit.01[i,] <- mean.probit[3,]
      probit.11[i,] <- mean.probit[4,]
      
      logit.00[i,] <- mean.logit[1,]
      logit.10[i,] <- mean.logit[2,]
      logit.01[i,] <- mean.logit[3,]
      logit.11[i,] <- mean.logit[4,]
      
      loglog.00[i,] <- mean.loglog[1,]
      loglog.10[i,] <- mean.loglog[2,]
      loglog.01[i,] <- mean.loglog[3,]
      loglog.11[i,] <- mean.loglog[4,]
      
      cloglog.00[i,] <- mean.cloglog[1,]
      cloglog.10[i,] <- mean.cloglog[2,]
      cloglog.01[i,] <- mean.cloglog[3,]
      cloglog.11[i,] <- mean.cloglog[4,]
      
      ols.00[i,] <- mean.ols[1,]
      ols.10[i,] <- mean.ols[2,]
      ols.01[i,] <- mean.ols[3,]
      ols.11[i,] <- mean.ols[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              ols.00=ols.00,
              ols.10=ols.10,
              ols.01=ols.01,
              ols.11=ols.11
  ))
  
}

sim3_mean.fun.c <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=4, nrow=sim)
  probit.10 <- matrix(NA, ncol=4, nrow=sim)
  probit.01 <- matrix(NA, ncol=4, nrow=sim)
  probit.11 <- matrix(NA, ncol=4, nrow=sim)
  colnames(probit.00) <- c("est", "se", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=4, nrow=sim)
  logit.10 <- matrix(NA, ncol=4, nrow=sim)
  logit.01 <- matrix(NA, ncol=4, nrow=sim)
  logit.11 <- matrix(NA, ncol=4, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=4, nrow=sim)
  loglog.10 <- matrix(NA, ncol=4, nrow=sim)
  loglog.01 <- matrix(NA, ncol=4, nrow=sim)
  loglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  ols.00 <- matrix(NA, ncol=4, nrow=sim)
  ols.10 <- matrix(NA, ncol=4, nrow=sim)
  ols.01 <- matrix(NA, ncol=4, nrow=sim)
  ols.11 <- matrix(NA, ncol=4, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.c(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      ols.mod <- ols(y~ z1+z2, data=data)
      
      mean.probit <- mean.orm(orm.probit, new.data, se=TRUE)
      mean.logit <- mean.orm(orm.logit, new.data, se=TRUE)
      mean.loglog <- mean.orm(orm.loglog, new.data, se=TRUE)
      mean.cloglog <- mean.orm(orm.cloglog, new.data, se=TRUE)
      
      mean.ols <- predict(ols.mod, new.data, se.fit=TRUE, conf.int=0.95)
      mean.ols <- cbind(mean.ols$linear.predictors, mean.ols$se.fit, mean.ols$lower, mean.ols$upper)
      
      probit.00[i,] <- mean.probit[1,]
      probit.10[i,] <- mean.probit[2,]
      probit.01[i,] <- mean.probit[3,]
      probit.11[i,] <- mean.probit[4,]
      
      logit.00[i,] <- mean.logit[1,]
      logit.10[i,] <- mean.logit[2,]
      logit.01[i,] <- mean.logit[3,]
      logit.11[i,] <- mean.logit[4,]
      
      loglog.00[i,] <- mean.loglog[1,]
      loglog.10[i,] <- mean.loglog[2,]
      loglog.01[i,] <- mean.loglog[3,]
      loglog.11[i,] <- mean.loglog[4,]
      
      cloglog.00[i,] <- mean.cloglog[1,]
      cloglog.10[i,] <- mean.cloglog[2,]
      cloglog.01[i,] <- mean.cloglog[3,]
      cloglog.11[i,] <- mean.cloglog[4,]
      
      ols.00[i,] <- mean.ols[1,]
      ols.10[i,] <- mean.ols[2,]
      ols.01[i,] <- mean.ols[3,]
      ols.11[i,] <- mean.ols[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              ols.00=ols.00,
              ols.10=ols.10,
              ols.01=ols.01,
              ols.11=ols.11
  ))
  
}


sim3_mean.fun.d <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=4, nrow=sim)
  probit.10 <- matrix(NA, ncol=4, nrow=sim)
  probit.01 <- matrix(NA, ncol=4, nrow=sim)
  probit.11 <- matrix(NA, ncol=4, nrow=sim)
  colnames(probit.00) <- c("est", "se", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=4, nrow=sim)
  logit.10 <- matrix(NA, ncol=4, nrow=sim)
  logit.01 <- matrix(NA, ncol=4, nrow=sim)
  logit.11 <- matrix(NA, ncol=4, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=4, nrow=sim)
  loglog.10 <- matrix(NA, ncol=4, nrow=sim)
  loglog.01 <- matrix(NA, ncol=4, nrow=sim)
  loglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  ols.00 <- matrix(NA, ncol=4, nrow=sim)
  ols.10 <- matrix(NA, ncol=4, nrow=sim)
  ols.01 <- matrix(NA, ncol=4, nrow=sim)
  ols.11 <- matrix(NA, ncol=4, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.d(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      ols.mod <- ols(y~ z1+z2, data=data)
      
      mean.probit <- mean.orm(orm.probit, new.data, se=TRUE)
      mean.logit <- mean.orm(orm.logit, new.data, se=TRUE)
      mean.loglog <- mean.orm(orm.loglog, new.data, se=TRUE)
      mean.cloglog <- mean.orm(orm.cloglog, new.data, se=TRUE)
      
      mean.ols <- predict(ols.mod, new.data, se.fit=TRUE, conf.int=0.95)
      mean.ols <- cbind(mean.ols$linear.predictors, mean.ols$se.fit, mean.ols$lower, mean.ols$upper)
      
      probit.00[i,] <- mean.probit[1,]
      probit.10[i,] <- mean.probit[2,]
      probit.01[i,] <- mean.probit[3,]
      probit.11[i,] <- mean.probit[4,]
      
      logit.00[i,] <- mean.logit[1,]
      logit.10[i,] <- mean.logit[2,]
      logit.01[i,] <- mean.logit[3,]
      logit.11[i,] <- mean.logit[4,]
      
      loglog.00[i,] <- mean.loglog[1,]
      loglog.10[i,] <- mean.loglog[2,]
      loglog.01[i,] <- mean.loglog[3,]
      loglog.11[i,] <- mean.loglog[4,]
      
      cloglog.00[i,] <- mean.cloglog[1,]
      cloglog.10[i,] <- mean.cloglog[2,]
      cloglog.01[i,] <- mean.cloglog[3,]
      cloglog.11[i,] <- mean.cloglog[4,]
      
      ols.00[i,] <- mean.ols[1,]
      ols.10[i,] <- mean.ols[2,]
      ols.01[i,] <- mean.ols[3,]
      ols.11[i,] <- mean.ols[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              ols.00=ols.00,
              ols.10=ols.10,
              ols.01=ols.01,
              ols.11=ols.11
  ))
  
}


sim3_mean.fun.e <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=4, nrow=sim)
  probit.10 <- matrix(NA, ncol=4, nrow=sim)
  probit.01 <- matrix(NA, ncol=4, nrow=sim)
  probit.11 <- matrix(NA, ncol=4, nrow=sim)
  colnames(probit.00) <- c("est", "se", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=4, nrow=sim)
  logit.10 <- matrix(NA, ncol=4, nrow=sim)
  logit.01 <- matrix(NA, ncol=4, nrow=sim)
  logit.11 <- matrix(NA, ncol=4, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=4, nrow=sim)
  loglog.10 <- matrix(NA, ncol=4, nrow=sim)
  loglog.01 <- matrix(NA, ncol=4, nrow=sim)
  loglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  ols.00 <- matrix(NA, ncol=4, nrow=sim)
  ols.10 <- matrix(NA, ncol=4, nrow=sim)
  ols.01 <- matrix(NA, ncol=4, nrow=sim)
  ols.11 <- matrix(NA, ncol=4, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.e(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      ols.mod <- ols(y~ z1+z2, data=data)
      
      mean.probit <- mean.orm(orm.probit, new.data, se=TRUE)
      mean.logit <- mean.orm(orm.logit, new.data, se=TRUE)
      mean.loglog <- mean.orm(orm.loglog, new.data, se=TRUE)
      mean.cloglog <- mean.orm(orm.cloglog, new.data, se=TRUE)
      
      mean.ols <- predict(ols.mod, new.data, se.fit=TRUE, conf.int=0.95)
      mean.ols <- cbind(mean.ols$linear.predictors, mean.ols$se.fit, mean.ols$lower, mean.ols$upper)
      
      probit.00[i,] <- mean.probit[1,]
      probit.10[i,] <- mean.probit[2,]
      probit.01[i,] <- mean.probit[3,]
      probit.11[i,] <- mean.probit[4,]
      
      logit.00[i,] <- mean.logit[1,]
      logit.10[i,] <- mean.logit[2,]
      logit.01[i,] <- mean.logit[3,]
      logit.11[i,] <- mean.logit[4,]
      
      loglog.00[i,] <- mean.loglog[1,]
      loglog.10[i,] <- mean.loglog[2,]
      loglog.01[i,] <- mean.loglog[3,]
      loglog.11[i,] <- mean.loglog[4,]
      
      cloglog.00[i,] <- mean.cloglog[1,]
      cloglog.10[i,] <- mean.cloglog[2,]
      cloglog.01[i,] <- mean.cloglog[3,]
      cloglog.11[i,] <- mean.cloglog[4,]
      
      ols.00[i,] <- mean.ols[1,]
      ols.10[i,] <- mean.ols[2,]
      ols.01[i,] <- mean.ols[3,]
      ols.11[i,] <- mean.ols[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              ols.00=ols.00,
              ols.10=ols.10,
              ols.01=ols.01,
              ols.11=ols.11
  ))
  
}


sim3_mean.fun.f <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=4, nrow=sim)
  probit.10 <- matrix(NA, ncol=4, nrow=sim)
  probit.01 <- matrix(NA, ncol=4, nrow=sim)
  probit.11 <- matrix(NA, ncol=4, nrow=sim)
  colnames(probit.00) <- c("est", "se", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=4, nrow=sim)
  logit.10 <- matrix(NA, ncol=4, nrow=sim)
  logit.01 <- matrix(NA, ncol=4, nrow=sim)
  logit.11 <- matrix(NA, ncol=4, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=4, nrow=sim)
  loglog.10 <- matrix(NA, ncol=4, nrow=sim)
  loglog.01 <- matrix(NA, ncol=4, nrow=sim)
  loglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  ols.00 <- matrix(NA, ncol=4, nrow=sim)
  ols.10 <- matrix(NA, ncol=4, nrow=sim)
  ols.01 <- matrix(NA, ncol=4, nrow=sim)
  ols.11 <- matrix(NA, ncol=4, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.f(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      ols.mod <- ols(y~ z1+z2, data=data)
      
      mean.probit <- mean.orm(orm.probit, new.data, se=TRUE)
      mean.logit <- mean.orm(orm.logit, new.data, se=TRUE)
      mean.loglog <- mean.orm(orm.loglog, new.data, se=TRUE)
      mean.cloglog <- mean.orm(orm.cloglog, new.data, se=TRUE)
      
      mean.ols <- predict(ols.mod, new.data, se.fit=TRUE, conf.int=0.95)
      mean.ols <- cbind(mean.ols$linear.predictors, mean.ols$se.fit, mean.ols$lower, mean.ols$upper)
      
      probit.00[i,] <- mean.probit[1,]
      probit.10[i,] <- mean.probit[2,]
      probit.01[i,] <- mean.probit[3,]
      probit.11[i,] <- mean.probit[4,]
      
      logit.00[i,] <- mean.logit[1,]
      logit.10[i,] <- mean.logit[2,]
      logit.01[i,] <- mean.logit[3,]
      logit.11[i,] <- mean.logit[4,]
      
      loglog.00[i,] <- mean.loglog[1,]
      loglog.10[i,] <- mean.loglog[2,]
      loglog.01[i,] <- mean.loglog[3,]
      loglog.11[i,] <- mean.loglog[4,]
      
      cloglog.00[i,] <- mean.cloglog[1,]
      cloglog.10[i,] <- mean.cloglog[2,]
      cloglog.01[i,] <- mean.cloglog[3,]
      cloglog.11[i,] <- mean.cloglog[4,]
      
      ols.00[i,] <- mean.ols[1,]
      ols.10[i,] <- mean.ols[2,]
      ols.01[i,] <- mean.ols[3,]
      ols.11[i,] <- mean.ols[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              ols.00=ols.00,
              ols.10=ols.10,
              ols.01=ols.01,
              ols.11=ols.11
  ))
  
}

sim3_mean.fun.g <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=4, nrow=sim)
  probit.10 <- matrix(NA, ncol=4, nrow=sim)
  probit.01 <- matrix(NA, ncol=4, nrow=sim)
  probit.11 <- matrix(NA, ncol=4, nrow=sim)
  colnames(probit.00) <- c("est", "se", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=4, nrow=sim)
  logit.10 <- matrix(NA, ncol=4, nrow=sim)
  logit.01 <- matrix(NA, ncol=4, nrow=sim)
  logit.11 <- matrix(NA, ncol=4, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=4, nrow=sim)
  loglog.10 <- matrix(NA, ncol=4, nrow=sim)
  loglog.01 <- matrix(NA, ncol=4, nrow=sim)
  loglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  ols.00 <- matrix(NA, ncol=4, nrow=sim)
  ols.10 <- matrix(NA, ncol=4, nrow=sim)
  ols.01 <- matrix(NA, ncol=4, nrow=sim)
  ols.11 <- matrix(NA, ncol=4, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.g(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      ols.mod <- ols(y~ z1+z2, data=data)
      
      mean.probit <- mean.orm(orm.probit, new.data, se=TRUE)
      mean.logit <- mean.orm(orm.logit, new.data, se=TRUE)
      mean.loglog <- mean.orm(orm.loglog, new.data, se=TRUE)
      mean.cloglog <- mean.orm(orm.cloglog, new.data, se=TRUE)
      
      mean.ols <- predict(ols.mod, new.data, se.fit=TRUE, conf.int=0.95)
      mean.ols <- cbind(mean.ols$linear.predictors, mean.ols$se.fit, mean.ols$lower, mean.ols$upper)
      
      probit.00[i,] <- mean.probit[1,]
      probit.10[i,] <- mean.probit[2,]
      probit.01[i,] <- mean.probit[3,]
      probit.11[i,] <- mean.probit[4,]
      
      logit.00[i,] <- mean.logit[1,]
      logit.10[i,] <- mean.logit[2,]
      logit.01[i,] <- mean.logit[3,]
      logit.11[i,] <- mean.logit[4,]
      
      loglog.00[i,] <- mean.loglog[1,]
      loglog.10[i,] <- mean.loglog[2,]
      loglog.01[i,] <- mean.loglog[3,]
      loglog.11[i,] <- mean.loglog[4,]
      
      cloglog.00[i,] <- mean.cloglog[1,]
      cloglog.10[i,] <- mean.cloglog[2,]
      cloglog.01[i,] <- mean.cloglog[3,]
      cloglog.11[i,] <- mean.cloglog[4,]
      
      ols.00[i,] <- mean.ols[1,]
      ols.10[i,] <- mean.ols[2,]
      ols.01[i,] <- mean.ols[3,]
      ols.11[i,] <- mean.ols[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              ols.00=ols.00,
              ols.10=ols.10,
              ols.01=ols.01,
              ols.11=ols.11
  ))
  
}

sim3_mean.fun.h <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=4, nrow=sim)
  probit.10 <- matrix(NA, ncol=4, nrow=sim)
  probit.01 <- matrix(NA, ncol=4, nrow=sim)
  probit.11 <- matrix(NA, ncol=4, nrow=sim)
  colnames(probit.00) <- c("est", "se", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=4, nrow=sim)
  logit.10 <- matrix(NA, ncol=4, nrow=sim)
  logit.01 <- matrix(NA, ncol=4, nrow=sim)
  logit.11 <- matrix(NA, ncol=4, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=4, nrow=sim)
  loglog.10 <- matrix(NA, ncol=4, nrow=sim)
  loglog.01 <- matrix(NA, ncol=4, nrow=sim)
  loglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=4, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=4, nrow=sim)
  
  ols.00 <- matrix(NA, ncol=4, nrow=sim)
  ols.10 <- matrix(NA, ncol=4, nrow=sim)
  ols.01 <- matrix(NA, ncol=4, nrow=sim)
  ols.11 <- matrix(NA, ncol=4, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.h(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      ols.mod <- ols(y~ z1+z2, data=data)
      
      mean.probit <- mean.orm(orm.probit, new.data, se=TRUE)
      mean.logit <- mean.orm(orm.logit, new.data, se=TRUE)
      mean.loglog <- mean.orm(orm.loglog, new.data, se=TRUE)
      mean.cloglog <- mean.orm(orm.cloglog, new.data, se=TRUE)
      
      mean.ols <- predict(ols.mod, new.data, se.fit=TRUE, conf.int=0.95)
      mean.ols <- cbind(mean.ols$linear.predictors, mean.ols$se.fit, mean.ols$lower, mean.ols$upper)
      
      probit.00[i,] <- mean.probit[1,]
      probit.10[i,] <- mean.probit[2,]
      probit.01[i,] <- mean.probit[3,]
      probit.11[i,] <- mean.probit[4,]
      
      logit.00[i,] <- mean.logit[1,]
      logit.10[i,] <- mean.logit[2,]
      logit.01[i,] <- mean.logit[3,]
      logit.11[i,] <- mean.logit[4,]
      
      loglog.00[i,] <- mean.loglog[1,]
      loglog.10[i,] <- mean.loglog[2,]
      loglog.01[i,] <- mean.loglog[3,]
      loglog.11[i,] <- mean.loglog[4,]
      
      cloglog.00[i,] <- mean.cloglog[1,]
      cloglog.10[i,] <- mean.cloglog[2,]
      cloglog.01[i,] <- mean.cloglog[3,]
      cloglog.11[i,] <- mean.cloglog[4,]
      
      ols.00[i,] <- mean.ols[1,]
      ols.10[i,] <- mean.ols[2,]
      ols.01[i,] <- mean.ols[3,]
      ols.11[i,] <- mean.ols[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              ols.00=ols.00,
              ols.10=ols.10,
              ols.01=ols.01,
              ols.11=ols.11
  ))
  
}

seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1

sim=10^4

sim3_mean.a.n50 <-sim3_mean.fun.a(sim=sim, n=50) 
sim3_mean.b.n50 <-sim3_mean.fun.b(sim=sim, n=50) 
sim3_mean.c.n50 <-sim3_mean.fun.c(sim=sim, n=50) 
sim3_mean.d.n50 <-sim3_mean.fun.d(sim=sim, n=50) 
sim3_mean.e.n50 <-sim3_mean.fun.e(sim=sim, n=50) 
sim3_mean.f.n50 <-sim3_mean.fun.f(sim=sim, n=50) 
sim3_mean.g.n50 <-sim3_mean.fun.g(sim=sim, n=50) 
sim3_mean.h.n50 <-sim3_mean.fun.h(sim=sim, n=50) 

sim3_mean.a.n100 <-sim3_mean.fun.a(sim=sim, n=50) 
sim3_mean.b.n100 <-sim3_mean.fun.b(sim=sim, n=50) 
sim3_mean.c.n100 <-sim3_mean.fun.c(sim=sim, n=50) 
sim3_mean.d.n100 <-sim3_mean.fun.d(sim=sim, n=50) 
sim3_mean.e.n100 <-sim3_mean.fun.e(sim=sim, n=50) 
sim3_mean.f.n100 <-sim3_mean.fun.f(sim=sim, n=50) 
sim3_mean.g.n100 <-sim3_mean.fun.g(sim=sim, n=50) 
sim3_mean.h.n100 <-sim3_mean.fun.h(sim=sim, n=50) 


sim3_mean.a.n50 <-sim3_mean.fun.a(sim=sim, n=100) 
sim3_mean.b.n50 <-sim3_mean.fun.b(sim=sim, n=100) 
sim3_mean.c.n50 <-sim3_mean.fun.c(sim=sim, n=100) 
sim3_mean.d.n50 <-sim3_mean.fun.d(sim=sim, n=100) 
sim3_mean.e.n50 <-sim3_mean.fun.e(sim=sim, n=100) 
sim3_mean.f.n50 <-sim3_mean.fun.f(sim=sim, n=100) 
sim3_mean.g.n50 <-sim3_mean.fun.g(sim=sim, n=100) 
sim3_mean.h.n50 <-sim3_mean.fun.h(sim=sim, n=100) 

sim3_mean.a.n200 <-sim3_mean.fun.a(sim=sim, n=200) 
sim3_mean.b.n200 <-sim3_mean.fun.b(sim=sim, n=200) 
sim3_mean.c.n200 <-sim3_mean.fun.c(sim=sim, n=200) 
sim3_mean.d.n200 <-sim3_mean.fun.d(sim=sim, n=200) 
sim3_mean.e.n200 <-sim3_mean.fun.e(sim=sim, n=200) 
sim3_mean.f.n200 <-sim3_mean.fun.f(sim=sim, n=200) 
sim3_mean.g.n200 <-sim3_mean.fun.g(sim=sim, n=200) 
sim3_mean.h.n200 <-sim3_mean.fun.h(sim=sim, n=200)

new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))
beta <- c(1, -0.5)
true=as.matrix(new.data)%*%beta

sim3_mean.summary <- function(true=as.matrix(new.data)%*%beta, result){
  
  result.00.est <- c(mean(result$probit.00[,1], na.rm=TRUE),
                     mean(result$logit.00[,1], na.rm=TRUE),
                     mean(result$loglog.00[,1], na.rm=TRUE),
                     mean(result$cloglog.00[,1], na.rm=TRUE))
  result.00.coverage <- c(mean(result$probit.00[,3]<=true[1,1]&result$probit.00[,4]>=true[1,1], na.rm=TRUE),
                          mean(result$logit.00[,3]<=true[1,1]&result$logit.00[,4]>=true[1,1], na.rm=TRUE),
                          mean(result$loglog.00[,3]<=true[1,1]&result$loglog.00[,4]>=true[1,1], na.rm=TRUE),
                          mean(result$cloglog.00[,3]<=true[1,1]&result$cloglog.00[,4]>=true[1,1], na.rm=TRUE))
  
  result.00.mse <- c(mean((result$probit.00[,1]-true[1,1])^2, na.rm=TRUE),
                     mean( (result$logit.00[,1]-true[1,1])^2, na.rm=TRUE),
                     mean((result$loglog.00[,1]-true[1,1])^2, na.rm=TRUE),
                     mean((result$cloglog.00[,1]-true[1,1])^2, na.rm=TRUE),
                     mean((result$ols.00[,1]-true[1,1])^2, na.rm=TRUE)  )
  
  result.00.re <- (result.00.mse[5]/result.00.mse)[-5]
  
  
  result.00 <- cbind(format(round(result.00.est,2), nsmall=2), 
                     format(round(result.00.coverage, 3), nsmall=3),
                     format(round(result.00.re,3), nsmall=3))
  result.00.out <- c(format(round(true[1,1],2),nsmall=0),
                     result.00[1,],
                     result.00[2,],
                     result.00[3,],
                     result.00[4,])
  
  ########  10
  result.10.est <- c(mean(result$probit.10[,1], na.rm=TRUE),
                     mean(result$logit.10[,1], na.rm=TRUE),
                     mean(result$loglog.10[,1], na.rm=TRUE),
                     mean(result$cloglog.10[,1], na.rm=TRUE))
  result.10.coverage <- c(mean(result$probit.10[,3]<=true[2,1]&result$probit.10[,4]>=true[2,1], na.rm=TRUE),
                          mean(result$logit.10[,3]<=true[2,1]&result$logit.10[,4]>=true[2,1], na.rm=TRUE),
                          mean(result$loglog.10[,3]<=true[2,1]&result$loglog.10[,4]>=true[2,1], na.rm=TRUE),
                          mean(result$cloglog.10[,3]<=true[2,1]&result$cloglog.10[,4]>=true[2,1], na.rm=TRUE))
  
  result.10.mse <- c(mean((result$probit.10[,1]-true[2,1])^2, na.rm=TRUE),
                     mean( (result$logit.10[,1]-true[2,1])^2, na.rm=TRUE),
                     mean((result$loglog.10[,1]-true[2,1])^2, na.rm=TRUE),
                     mean((result$cloglog.10[,1]-true[2,1])^2, na.rm=TRUE),
                     mean((result$ols.10[,1]-true[2,1])^2, na.rm=TRUE)  )
  
  result.10.re <- (result.10.mse[5]/result.10.mse)[-5]
  
  
  result.10 <- cbind(format(round(result.10.est,2), nsmall=2), 
                     format(round(result.10.coverage, 3), nsmall=3),
                     format(round(result.10.re,3), nsmall=3))
  result.10.out <- c(format(round(true[2,1],2),nsmall=0),
                     result.10[1,],
                     result.10[2,],
                     result.10[3,],
                     result.10[4,])
  
  ###### 01
  result.01.est <- c(mean(result$probit.01[,1], na.rm=TRUE),
                     mean(result$logit.01[,1], na.rm=TRUE),
                     mean(result$loglog.01[,1], na.rm=TRUE),
                     mean(result$cloglog.01[,1], na.rm=TRUE))
  result.01.coverage <- c(mean(result$probit.01[,3]<=true[3,1]&result$probit.01[,4]>=true[3,1], na.rm=TRUE),
                          mean(result$logit.01[,3]<=true[3,1]&result$logit.01[,4]>=true[3,1], na.rm=TRUE),
                          mean(result$loglog.01[,3]<=true[3,1]&result$loglog.01[,4]>=true[3,1], na.rm=TRUE),
                          mean(result$cloglog.01[,3]<=true[3,1]&result$cloglog.01[,4]>=true[3,1], na.rm=TRUE))
  
  result.01.mse <- c(mean((result$probit.01[,1]-true[3,1])^2, na.rm=TRUE),
                     mean( (result$logit.01[,1]-true[3,1])^2, na.rm=TRUE),
                     mean((result$loglog.01[,1]-true[3,1])^2, na.rm=TRUE),
                     mean((result$cloglog.01[,1]-true[3,1])^2, na.rm=TRUE),
                     mean((result$ols.01[,1]-true[3,1])^2, na.rm=TRUE)  )
  
  result.01.re <- (result.01.mse[5]/result.01.mse)[-5]
  
  
  result.01 <- cbind(format(round(result.01.est,2), nsmall=2), 
                     format(round(result.01.coverage, 3), nsmall=3),
                     format(round(result.01.re,3), nsmall=3))
  result.01.out <- c(format(round(true[3,1],2),nsmall=0),
                     result.01[1,],
                     result.01[2,],
                     result.01[3,],
                     result.01[4,])
  ###### 11
  
  result.11.est <- c(mean(result$probit.11[,1], na.rm=TRUE),
                     mean(result$logit.11[,1], na.rm=TRUE),
                     mean(result$loglog.11[,1], na.rm=TRUE),
                     mean(result$cloglog.11[,1], na.rm=TRUE))
  result.11.coverage <- c(mean(result$probit.11[,3]<=true[4,1]&result$probit.11[,4]>=true[4,1], na.rm=TRUE),
                          mean(result$logit.11[,3]<=true[4,1]&result$logit.11[,4]>=true[4,1], na.rm=TRUE),
                          mean(result$loglog.11[,3]<=true[4,1]&result$loglog.11[,4]>=true[4,1], na.rm=TRUE),
                          mean(result$cloglog.11[,3]<=true[4,1]&result$cloglog.11[,4]>=true[4,1], na.rm=TRUE))
  
  result.11.mse <- c(mean((result$probit.11[,1]-true[4,1])^2, na.rm=TRUE),
                     mean( (result$logit.11[,1]-true[4,1])^2, na.rm=TRUE),
                     mean((result$loglog.11[,1]-true[4,1])^2, na.rm=TRUE),
                     mean((result$cloglog.11[,1]-true[4,1])^2, na.rm=TRUE),
                     mean((result$ols.11[,1]-true[4,1])^2, na.rm=TRUE)  )
  
  result.11.re <- (result.11.mse[5]/result.11.mse)[-5]
  
  
  result.11 <- cbind(format(round(result.11.est,2), nsmall=2), 
                     format(round(result.11.coverage, 3), nsmall=3),
                     format(round(result.11.re,3), nsmall=3))
  result.11.out <- c(format(round(true[4,1],2),nsmall=0),
                     result.11[1,],
                     result.11[2,],
                     result.11[3,],
                     result.11[4,])
  
  result.out <- rbind(result.00.out,
                      result.10.out,
                      result.01.out,
                      result.11.out)
  
  return(result.out)
  
}

beta=c(1, -0.5)
new.data=as.matrix(data.frame(z1=c(0, 1, 0, 1), z2=c(0, 0, 1, 1)))
result.sim3_mean.a.n100 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.a.n100)
result.sim3_mean.b.n100 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.b.n100)
result.sim3_mean.c.n100 <- sim3_mean.summary(true=digamma(1)+as.matrix(new.data)%*%beta, result=sim3_mean.c.n100)
result.sim3_mean.d.n100 <- sim3_mean.summary(true=-digamma(1)+as.matrix(new.data)%*%beta, result=sim3_mean.d.n100)

result.sim3_mean.e.n100 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.e.n100)
result.sim3_mean.f.n100 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.f.n100)
result.sim3_mean.g.n100 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.g.n100)
result.sim3_mean.h.n100 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.h.n100)

###### Figure 9, Tables S.10 and S.11

sim3_mean.out.n100 <- rbind(rep("", 13),
                            result.sim3_mean.a.n100,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.b.n100,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.c.n100,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.d.n100,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.e.n100,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.f.n100,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.g.n100,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.h.n100)

###### Tables S.6 and S.7 ####################

result.sim3_mean.a.n50 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.a.n50)
result.sim3_mean.b.n50 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.b.n50)
result.sim3_mean.c.n50 <- sim3_mean.summary(true=digamma(1)+as.matrix(new.data)%*%beta, result=sim3_mean.c.n50)
result.sim3_mean.d.n50 <- sim3_mean.summary(true=-digamma(1)+as.matrix(new.data)%*%beta, result=sim3_mean.d.n50)

result.sim3_mean.e.n50 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.e.n50)
result.sim3_mean.f.n50 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.f.n50)
result.sim3_mean.g.n50 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.g.n50)
result.sim3_mean.h.n50 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.h.n50)


sim3_mean.out.n50 <- rbind(rep("", 13),
                           result.sim3_mean.a.n50,
                           rep("", 13),
                           rep("", 13),
                           result.sim3_mean.b.n50,
                           rep("", 13),
                           rep("", 13),
                           result.sim3_mean.c.n50,
                           rep("", 13),
                           rep("", 13),
                           result.sim3_mean.d.n50,
                           rep("", 13),
                           rep("", 13),
                           result.sim3_mean.e.n50,
                           rep("", 13),
                           rep("", 13),
                           result.sim3_mean.f.n50,
                           rep("", 13),
                           rep("", 13),
                           result.sim3_mean.g.n50,
                           rep("", 13),
                           rep("", 13),
                           result.sim3_mean.h.n50)


############# Tables S.14 and S.15 ################

result.sim3_mean.a.n200 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.a.n200)
result.sim3_mean.b.n200 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.b.n200)
result.sim3_mean.c.n200 <- sim3_mean.summary(true=digamma(1)+as.matrix(new.data)%*%beta, result=sim3_mean.c.n200)
result.sim3_mean.d.n200 <- sim3_mean.summary(true=-digamma(1)+as.matrix(new.data)%*%beta, result=sim3_mean.d.n200)

result.sim3_mean.e.n200 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.e.n200)
result.sim3_mean.f.n200 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.f.n200)
result.sim3_mean.g.n200 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.g.n200)
result.sim3_mean.h.n200 <- sim3_mean.summary(true=as.matrix(new.data)%*%beta, result=sim3_mean.h.n200)

sim3_mean.out.n200 <- rbind(rep("", 13),
                            result.sim3_mean.a.n200,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.b.n200,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.c.n200,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.d.n200,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.e.n200,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.f.n200,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.g.n200,
                            rep("", 13),
                            rep("", 13),
                            result.sim3_mean.h.n200)




######### Sim 2: medians #######
library(quantreg)

sim4_median.fun.a <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=3, nrow=sim)
  probit.10 <- matrix(NA, ncol=3, nrow=sim)
  probit.01 <- matrix(NA, ncol=3, nrow=sim)
  probit.11 <- matrix(NA, ncol=3, nrow=sim)
  colnames(probit.00) <- c("est", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=3, nrow=sim)
  logit.10 <- matrix(NA, ncol=3, nrow=sim)
  logit.01 <- matrix(NA, ncol=3, nrow=sim)
  logit.11 <- matrix(NA, ncol=3, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=3, nrow=sim)
  loglog.10 <- matrix(NA, ncol=3, nrow=sim)
  loglog.01 <- matrix(NA, ncol=3, nrow=sim)
  loglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  rq.00 <- matrix(NA, ncol=3, nrow=sim)
  rq.10 <- matrix(NA, ncol=3, nrow=sim)
  rq.01 <- matrix(NA, ncol=3, nrow=sim)
  rq.11 <- matrix(NA, ncol=3, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.a(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      rq.mod <- rq(y~ z1+z2, data=data)
      
      quantile.probit <- quantile.orm(orm.probit, new.data,probs=0.5, se=TRUE)
      quantile.probit <- cbind(quantile.probit$quantile, quantile.probit$lb, quantile.probit$ub)
      
      quantile.logit <- quantile.orm(orm.logit, new.data, probs=0.5, se=TRUE)
      quantile.logit <- cbind(quantile.logit$quantile, quantile.logit$lb, quantile.logit$ub)
      
      quantile.loglog <- quantile.orm(orm.loglog, new.data,probs=0.5, se=TRUE)
      quantile.loglog <- cbind(quantile.loglog$quantile, quantile.loglog$lb, quantile.loglog$ub)
      
      quantile.cloglog <- quantile.orm(orm.cloglog, new.data, probs=0.5,se=TRUE)
      quantile.cloglog <- cbind(quantile.cloglog$quantile, quantile.cloglog$lb, quantile.cloglog$ub)
      
      
      quantile.rq <- predict(rq.mod, new.data,  interval =  "confidence", conf.int=0.95)
      
      probit.00[i,] <- quantile.probit[1,]
      probit.10[i,] <- quantile.probit[2,]
      probit.01[i,] <- quantile.probit[3,]
      probit.11[i,] <- quantile.probit[4,]
      
      logit.00[i,] <- quantile.logit[1,]
      logit.10[i,] <- quantile.logit[2,]
      logit.01[i,] <- quantile.logit[3,]
      logit.11[i,] <- quantile.logit[4,]
      
      loglog.00[i,] <- quantile.loglog[1,]
      loglog.10[i,] <- quantile.loglog[2,]
      loglog.01[i,] <- quantile.loglog[3,]
      loglog.11[i,] <- quantile.loglog[4,]
      
      cloglog.00[i,] <- quantile.cloglog[1,]
      cloglog.10[i,] <- quantile.cloglog[2,]
      cloglog.01[i,] <- quantile.cloglog[3,]
      cloglog.11[i,] <- quantile.cloglog[4,]
      
      rq.00[i,] <- quantile.rq[1,]
      rq.10[i,] <- quantile.rq[2,]
      rq.01[i,] <- quantile.rq[3,]
      rq.11[i,] <- quantile.rq[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              rq.00=rq.00,
              rq.10=rq.10,
              rq.01=rq.01,
              rq.11=rq.11
  ))
  
}



sim4_median.fun.b <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=3, nrow=sim)
  probit.10 <- matrix(NA, ncol=3, nrow=sim)
  probit.01 <- matrix(NA, ncol=3, nrow=sim)
  probit.11 <- matrix(NA, ncol=3, nrow=sim)
  colnames(probit.00) <- c("est", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=3, nrow=sim)
  logit.10 <- matrix(NA, ncol=3, nrow=sim)
  logit.01 <- matrix(NA, ncol=3, nrow=sim)
  logit.11 <- matrix(NA, ncol=3, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=3, nrow=sim)
  loglog.10 <- matrix(NA, ncol=3, nrow=sim)
  loglog.01 <- matrix(NA, ncol=3, nrow=sim)
  loglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  rq.00 <- matrix(NA, ncol=3, nrow=sim)
  rq.10 <- matrix(NA, ncol=3, nrow=sim)
  rq.01 <- matrix(NA, ncol=3, nrow=sim)
  rq.11 <- matrix(NA, ncol=3, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.b(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      rq.mod <- rq(y~ z1+z2, data=data)
      
      quantile.probit <- quantile.orm(orm.probit, new.data,probs=0.5, se=TRUE)
      quantile.probit <- cbind(quantile.probit$quantile, quantile.probit$lb, quantile.probit$ub)
      
      quantile.logit <- quantile.orm(orm.logit, new.data, probs=0.5, se=TRUE)
      quantile.logit <- cbind(quantile.logit$quantile, quantile.logit$lb, quantile.logit$ub)
      
      quantile.loglog <- quantile.orm(orm.loglog, new.data,probs=0.5, se=TRUE)
      quantile.loglog <- cbind(quantile.loglog$quantile, quantile.loglog$lb, quantile.loglog$ub)
      
      quantile.cloglog <- quantile.orm(orm.cloglog, new.data, probs=0.5,se=TRUE)
      quantile.cloglog <- cbind(quantile.cloglog$quantile, quantile.cloglog$lb, quantile.cloglog$ub)
      
      
      quantile.rq <- predict(rq.mod, new.data,  interval =  "confidence", conf.int=0.95)
      
      probit.00[i,] <- quantile.probit[1,]
      probit.10[i,] <- quantile.probit[2,]
      probit.01[i,] <- quantile.probit[3,]
      probit.11[i,] <- quantile.probit[4,]
      
      logit.00[i,] <- quantile.logit[1,]
      logit.10[i,] <- quantile.logit[2,]
      logit.01[i,] <- quantile.logit[3,]
      logit.11[i,] <- quantile.logit[4,]
      
      loglog.00[i,] <- quantile.loglog[1,]
      loglog.10[i,] <- quantile.loglog[2,]
      loglog.01[i,] <- quantile.loglog[3,]
      loglog.11[i,] <- quantile.loglog[4,]
      
      cloglog.00[i,] <- quantile.cloglog[1,]
      cloglog.10[i,] <- quantile.cloglog[2,]
      cloglog.01[i,] <- quantile.cloglog[3,]
      cloglog.11[i,] <- quantile.cloglog[4,]
      
      rq.00[i,] <- quantile.rq[1,]
      rq.10[i,] <- quantile.rq[2,]
      rq.01[i,] <- quantile.rq[3,]
      rq.11[i,] <- quantile.rq[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              rq.00=rq.00,
              rq.10=rq.10,
              rq.01=rq.01,
              rq.11=rq.11
  ))
  
}




sim4_median.fun.c <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=3, nrow=sim)
  probit.10 <- matrix(NA, ncol=3, nrow=sim)
  probit.01 <- matrix(NA, ncol=3, nrow=sim)
  probit.11 <- matrix(NA, ncol=3, nrow=sim)
  colnames(probit.00) <- c("est", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=3, nrow=sim)
  logit.10 <- matrix(NA, ncol=3, nrow=sim)
  logit.01 <- matrix(NA, ncol=3, nrow=sim)
  logit.11 <- matrix(NA, ncol=3, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=3, nrow=sim)
  loglog.10 <- matrix(NA, ncol=3, nrow=sim)
  loglog.01 <- matrix(NA, ncol=3, nrow=sim)
  loglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  rq.00 <- matrix(NA, ncol=3, nrow=sim)
  rq.10 <- matrix(NA, ncol=3, nrow=sim)
  rq.01 <- matrix(NA, ncol=3, nrow=sim)
  rq.11 <- matrix(NA, ncol=3, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.c(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      rq.mod <- rq(y~ z1+z2, data=data)
      
      quantile.probit <- quantile.orm(orm.probit, new.data,probs=0.5, se=TRUE)
      quantile.probit <- cbind(quantile.probit$quantile, quantile.probit$lb, quantile.probit$ub)
      
      quantile.logit <- quantile.orm(orm.logit, new.data, probs=0.5, se=TRUE)
      quantile.logit <- cbind(quantile.logit$quantile, quantile.logit$lb, quantile.logit$ub)
      
      quantile.loglog <- quantile.orm(orm.loglog, new.data,probs=0.5, se=TRUE)
      quantile.loglog <- cbind(quantile.loglog$quantile, quantile.loglog$lb, quantile.loglog$ub)
      
      quantile.cloglog <- quantile.orm(orm.cloglog, new.data, probs=0.5,se=TRUE)
      quantile.cloglog <- cbind(quantile.cloglog$quantile, quantile.cloglog$lb, quantile.cloglog$ub)
      
      
      quantile.rq <- predict(rq.mod, new.data,  interval =  "confidence", conf.int=0.95)
      
      probit.00[i,] <- quantile.probit[1,]
      probit.10[i,] <- quantile.probit[2,]
      probit.01[i,] <- quantile.probit[3,]
      probit.11[i,] <- quantile.probit[4,]
      
      logit.00[i,] <- quantile.logit[1,]
      logit.10[i,] <- quantile.logit[2,]
      logit.01[i,] <- quantile.logit[3,]
      logit.11[i,] <- quantile.logit[4,]
      
      loglog.00[i,] <- quantile.loglog[1,]
      loglog.10[i,] <- quantile.loglog[2,]
      loglog.01[i,] <- quantile.loglog[3,]
      loglog.11[i,] <- quantile.loglog[4,]
      
      cloglog.00[i,] <- quantile.cloglog[1,]
      cloglog.10[i,] <- quantile.cloglog[2,]
      cloglog.01[i,] <- quantile.cloglog[3,]
      cloglog.11[i,] <- quantile.cloglog[4,]
      
      rq.00[i,] <- quantile.rq[1,]
      rq.10[i,] <- quantile.rq[2,]
      rq.01[i,] <- quantile.rq[3,]
      rq.11[i,] <- quantile.rq[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              rq.00=rq.00,
              rq.10=rq.10,
              rq.01=rq.01,
              rq.11=rq.11
  ))
  
}




sim4_median.fun.d <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=3, nrow=sim)
  probit.10 <- matrix(NA, ncol=3, nrow=sim)
  probit.01 <- matrix(NA, ncol=3, nrow=sim)
  probit.11 <- matrix(NA, ncol=3, nrow=sim)
  colnames(probit.00) <- c("est", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=3, nrow=sim)
  logit.10 <- matrix(NA, ncol=3, nrow=sim)
  logit.01 <- matrix(NA, ncol=3, nrow=sim)
  logit.11 <- matrix(NA, ncol=3, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=3, nrow=sim)
  loglog.10 <- matrix(NA, ncol=3, nrow=sim)
  loglog.01 <- matrix(NA, ncol=3, nrow=sim)
  loglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  rq.00 <- matrix(NA, ncol=3, nrow=sim)
  rq.10 <- matrix(NA, ncol=3, nrow=sim)
  rq.01 <- matrix(NA, ncol=3, nrow=sim)
  rq.11 <- matrix(NA, ncol=3, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.d(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      rq.mod <- rq(y~ z1+z2, data=data)
      
      quantile.probit <- quantile.orm(orm.probit, new.data,probs=0.5, se=TRUE)
      quantile.probit <- cbind(quantile.probit$quantile, quantile.probit$lb, quantile.probit$ub)
      
      quantile.logit <- quantile.orm(orm.logit, new.data, probs=0.5, se=TRUE)
      quantile.logit <- cbind(quantile.logit$quantile, quantile.logit$lb, quantile.logit$ub)
      
      quantile.loglog <- quantile.orm(orm.loglog, new.data,probs=0.5, se=TRUE)
      quantile.loglog <- cbind(quantile.loglog$quantile, quantile.loglog$lb, quantile.loglog$ub)
      
      quantile.cloglog <- quantile.orm(orm.cloglog, new.data, probs=0.5,se=TRUE)
      quantile.cloglog <- cbind(quantile.cloglog$quantile, quantile.cloglog$lb, quantile.cloglog$ub)
      
      
      quantile.rq <- predict(rq.mod, new.data,  interval =  "confidence", conf.int=0.95)
      
      probit.00[i,] <- quantile.probit[1,]
      probit.10[i,] <- quantile.probit[2,]
      probit.01[i,] <- quantile.probit[3,]
      probit.11[i,] <- quantile.probit[4,]
      
      logit.00[i,] <- quantile.logit[1,]
      logit.10[i,] <- quantile.logit[2,]
      logit.01[i,] <- quantile.logit[3,]
      logit.11[i,] <- quantile.logit[4,]
      
      loglog.00[i,] <- quantile.loglog[1,]
      loglog.10[i,] <- quantile.loglog[2,]
      loglog.01[i,] <- quantile.loglog[3,]
      loglog.11[i,] <- quantile.loglog[4,]
      
      cloglog.00[i,] <- quantile.cloglog[1,]
      cloglog.10[i,] <- quantile.cloglog[2,]
      cloglog.01[i,] <- quantile.cloglog[3,]
      cloglog.11[i,] <- quantile.cloglog[4,]
      
      rq.00[i,] <- quantile.rq[1,]
      rq.10[i,] <- quantile.rq[2,]
      rq.01[i,] <- quantile.rq[3,]
      rq.11[i,] <- quantile.rq[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              rq.00=rq.00,
              rq.10=rq.10,
              rq.01=rq.01,
              rq.11=rq.11
  ))
  
}


sim4_median.fun.e <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=3, nrow=sim)
  probit.10 <- matrix(NA, ncol=3, nrow=sim)
  probit.01 <- matrix(NA, ncol=3, nrow=sim)
  probit.11 <- matrix(NA, ncol=3, nrow=sim)
  colnames(probit.00) <- c("est", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=3, nrow=sim)
  logit.10 <- matrix(NA, ncol=3, nrow=sim)
  logit.01 <- matrix(NA, ncol=3, nrow=sim)
  logit.11 <- matrix(NA, ncol=3, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=3, nrow=sim)
  loglog.10 <- matrix(NA, ncol=3, nrow=sim)
  loglog.01 <- matrix(NA, ncol=3, nrow=sim)
  loglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  rq.00 <- matrix(NA, ncol=3, nrow=sim)
  rq.10 <- matrix(NA, ncol=3, nrow=sim)
  rq.01 <- matrix(NA, ncol=3, nrow=sim)
  rq.11 <- matrix(NA, ncol=3, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.e(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      rq.mod <- rq(y~ z1+z2, data=data)
      
      quantile.probit <- quantile.orm(orm.probit, new.data,probs=0.5, se=TRUE)
      quantile.probit <- cbind(quantile.probit$quantile, quantile.probit$lb, quantile.probit$ub)
      
      quantile.logit <- quantile.orm(orm.logit, new.data, probs=0.5, se=TRUE)
      quantile.logit <- cbind(quantile.logit$quantile, quantile.logit$lb, quantile.logit$ub)
      
      quantile.loglog <- quantile.orm(orm.loglog, new.data,probs=0.5, se=TRUE)
      quantile.loglog <- cbind(quantile.loglog$quantile, quantile.loglog$lb, quantile.loglog$ub)
      
      quantile.cloglog <- quantile.orm(orm.cloglog, new.data, probs=0.5,se=TRUE)
      quantile.cloglog <- cbind(quantile.cloglog$quantile, quantile.cloglog$lb, quantile.cloglog$ub)
      
      
      quantile.rq <- predict(rq.mod, new.data,  interval =  "confidence", conf.int=0.95)
      
      probit.00[i,] <- quantile.probit[1,]
      probit.10[i,] <- quantile.probit[2,]
      probit.01[i,] <- quantile.probit[3,]
      probit.11[i,] <- quantile.probit[4,]
      
      logit.00[i,] <- quantile.logit[1,]
      logit.10[i,] <- quantile.logit[2,]
      logit.01[i,] <- quantile.logit[3,]
      logit.11[i,] <- quantile.logit[4,]
      
      loglog.00[i,] <- quantile.loglog[1,]
      loglog.10[i,] <- quantile.loglog[2,]
      loglog.01[i,] <- quantile.loglog[3,]
      loglog.11[i,] <- quantile.loglog[4,]
      
      cloglog.00[i,] <- quantile.cloglog[1,]
      cloglog.10[i,] <- quantile.cloglog[2,]
      cloglog.01[i,] <- quantile.cloglog[3,]
      cloglog.11[i,] <- quantile.cloglog[4,]
      
      rq.00[i,] <- quantile.rq[1,]
      rq.10[i,] <- quantile.rq[2,]
      rq.01[i,] <- quantile.rq[3,]
      rq.11[i,] <- quantile.rq[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              rq.00=rq.00,
              rq.10=rq.10,
              rq.01=rq.01,
              rq.11=rq.11
  ))
  
}



sim4_median.fun.f <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=3, nrow=sim)
  probit.10 <- matrix(NA, ncol=3, nrow=sim)
  probit.01 <- matrix(NA, ncol=3, nrow=sim)
  probit.11 <- matrix(NA, ncol=3, nrow=sim)
  colnames(probit.00) <- c("est", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=3, nrow=sim)
  logit.10 <- matrix(NA, ncol=3, nrow=sim)
  logit.01 <- matrix(NA, ncol=3, nrow=sim)
  logit.11 <- matrix(NA, ncol=3, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=3, nrow=sim)
  loglog.10 <- matrix(NA, ncol=3, nrow=sim)
  loglog.01 <- matrix(NA, ncol=3, nrow=sim)
  loglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  rq.00 <- matrix(NA, ncol=3, nrow=sim)
  rq.10 <- matrix(NA, ncol=3, nrow=sim)
  rq.01 <- matrix(NA, ncol=3, nrow=sim)
  rq.11 <- matrix(NA, ncol=3, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.f(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      rq.mod <- rq(y~ z1+z2, data=data)
      
      quantile.probit <- quantile.orm(orm.probit, new.data,probs=0.5, se=TRUE)
      quantile.probit <- cbind(quantile.probit$quantile, quantile.probit$lb, quantile.probit$ub)
      
      quantile.logit <- quantile.orm(orm.logit, new.data, probs=0.5, se=TRUE)
      quantile.logit <- cbind(quantile.logit$quantile, quantile.logit$lb, quantile.logit$ub)
      
      quantile.loglog <- quantile.orm(orm.loglog, new.data,probs=0.5, se=TRUE)
      quantile.loglog <- cbind(quantile.loglog$quantile, quantile.loglog$lb, quantile.loglog$ub)
      
      quantile.cloglog <- quantile.orm(orm.cloglog, new.data, probs=0.5,se=TRUE)
      quantile.cloglog <- cbind(quantile.cloglog$quantile, quantile.cloglog$lb, quantile.cloglog$ub)
      
      
      quantile.rq <- predict(rq.mod, new.data,  interval =  "confidence", conf.int=0.95)
      
      probit.00[i,] <- quantile.probit[1,]
      probit.10[i,] <- quantile.probit[2,]
      probit.01[i,] <- quantile.probit[3,]
      probit.11[i,] <- quantile.probit[4,]
      
      logit.00[i,] <- quantile.logit[1,]
      logit.10[i,] <- quantile.logit[2,]
      logit.01[i,] <- quantile.logit[3,]
      logit.11[i,] <- quantile.logit[4,]
      
      loglog.00[i,] <- quantile.loglog[1,]
      loglog.10[i,] <- quantile.loglog[2,]
      loglog.01[i,] <- quantile.loglog[3,]
      loglog.11[i,] <- quantile.loglog[4,]
      
      cloglog.00[i,] <- quantile.cloglog[1,]
      cloglog.10[i,] <- quantile.cloglog[2,]
      cloglog.01[i,] <- quantile.cloglog[3,]
      cloglog.11[i,] <- quantile.cloglog[4,]
      
      rq.00[i,] <- quantile.rq[1,]
      rq.10[i,] <- quantile.rq[2,]
      rq.01[i,] <- quantile.rq[3,]
      rq.11[i,] <- quantile.rq[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              rq.00=rq.00,
              rq.10=rq.10,
              rq.01=rq.01,
              rq.11=rq.11
  ))
  
}


sim4_median.fun.g <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=3, nrow=sim)
  probit.10 <- matrix(NA, ncol=3, nrow=sim)
  probit.01 <- matrix(NA, ncol=3, nrow=sim)
  probit.11 <- matrix(NA, ncol=3, nrow=sim)
  colnames(probit.00) <- c("est", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=3, nrow=sim)
  logit.10 <- matrix(NA, ncol=3, nrow=sim)
  logit.01 <- matrix(NA, ncol=3, nrow=sim)
  logit.11 <- matrix(NA, ncol=3, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=3, nrow=sim)
  loglog.10 <- matrix(NA, ncol=3, nrow=sim)
  loglog.01 <- matrix(NA, ncol=3, nrow=sim)
  loglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  rq.00 <- matrix(NA, ncol=3, nrow=sim)
  rq.10 <- matrix(NA, ncol=3, nrow=sim)
  rq.01 <- matrix(NA, ncol=3, nrow=sim)
  rq.11 <- matrix(NA, ncol=3, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.g(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      rq.mod <- rq(y~ z1+z2, data=data)
      
      quantile.probit <- quantile.orm(orm.probit, new.data,probs=0.5, se=TRUE)
      quantile.probit <- cbind(quantile.probit$quantile, quantile.probit$lb, quantile.probit$ub)
      
      quantile.logit <- quantile.orm(orm.logit, new.data, probs=0.5, se=TRUE)
      quantile.logit <- cbind(quantile.logit$quantile, quantile.logit$lb, quantile.logit$ub)
      
      quantile.loglog <- quantile.orm(orm.loglog, new.data,probs=0.5, se=TRUE)
      quantile.loglog <- cbind(quantile.loglog$quantile, quantile.loglog$lb, quantile.loglog$ub)
      
      quantile.cloglog <- quantile.orm(orm.cloglog, new.data, probs=0.5,se=TRUE)
      quantile.cloglog <- cbind(quantile.cloglog$quantile, quantile.cloglog$lb, quantile.cloglog$ub)
      
      
      quantile.rq <- predict(rq.mod, new.data,  interval =  "confidence", conf.int=0.95)
      
      probit.00[i,] <- quantile.probit[1,]
      probit.10[i,] <- quantile.probit[2,]
      probit.01[i,] <- quantile.probit[3,]
      probit.11[i,] <- quantile.probit[4,]
      
      logit.00[i,] <- quantile.logit[1,]
      logit.10[i,] <- quantile.logit[2,]
      logit.01[i,] <- quantile.logit[3,]
      logit.11[i,] <- quantile.logit[4,]
      
      loglog.00[i,] <- quantile.loglog[1,]
      loglog.10[i,] <- quantile.loglog[2,]
      loglog.01[i,] <- quantile.loglog[3,]
      loglog.11[i,] <- quantile.loglog[4,]
      
      cloglog.00[i,] <- quantile.cloglog[1,]
      cloglog.10[i,] <- quantile.cloglog[2,]
      cloglog.01[i,] <- quantile.cloglog[3,]
      cloglog.11[i,] <- quantile.cloglog[4,]
      
      rq.00[i,] <- quantile.rq[1,]
      rq.10[i,] <- quantile.rq[2,]
      rq.01[i,] <- quantile.rq[3,]
      rq.11[i,] <- quantile.rq[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              rq.00=rq.00,
              rq.10=rq.10,
              rq.01=rq.01,
              rq.11=rq.11
  ))
  
}


sim4_median.fun.h <- function(sim=100,seed=1, n=50,
                         p=0.5, alpha=0, beta=c(1,-0.5),new.data=data.frame(z1=c(0,1, 0, 1), z2=c(0, 0, 1, 1))){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  
  probit.00 <- matrix(NA, ncol=3, nrow=sim)
  probit.10 <- matrix(NA, ncol=3, nrow=sim)
  probit.01 <- matrix(NA, ncol=3, nrow=sim)
  probit.11 <- matrix(NA, ncol=3, nrow=sim)
  colnames(probit.00) <- c("est", "lb", "ub")
  
  logit.00 <- matrix(NA, ncol=3, nrow=sim)
  logit.10 <- matrix(NA, ncol=3, nrow=sim)
  logit.01 <- matrix(NA, ncol=3, nrow=sim)
  logit.11 <- matrix(NA, ncol=3, nrow=sim)
  
  loglog.00 <- matrix(NA, ncol=3, nrow=sim)
  loglog.10 <- matrix(NA, ncol=3, nrow=sim)
  loglog.01 <- matrix(NA, ncol=3, nrow=sim)
  loglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  cloglog.00 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.10 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.01 <- matrix(NA, ncol=3, nrow=sim)
  cloglog.11 <- matrix(NA, ncol=3, nrow=sim)
  
  rq.00 <- matrix(NA, ncol=3, nrow=sim)
  rq.10 <- matrix(NA, ncol=3, nrow=sim)
  rq.01 <- matrix(NA, ncol=3, nrow=sim)
  rq.11 <- matrix(NA, ncol=3, nrow=sim)
  
  
  
  
  for(i in 1:sim){
    try({
      data <- generate.data.h(seed=seeds[i], n=n,p=p, alpha=alpha, beta=beta)
      orm.probit <- orm(y~z1+z2, data=data, family=probit, x=TRUE, y=TRUE)
      orm.logit <- orm(y~z1+z2, data=data, family=logistic, x=TRUE, y=TRUE)
      orm.loglog <- orm(y~z1+z2, data=data, family=loglog, x=TRUE, y=TRUE)
      orm.cloglog <- orm(y~z1+z2, data=data, family=cloglog, x=TRUE, y=TRUE)
      rq.mod <- rq(y~ z1+z2, data=data)
      
      quantile.probit <- quantile.orm(orm.probit, new.data,probs=0.5, se=TRUE)
      quantile.probit <- cbind(quantile.probit$quantile, quantile.probit$lb, quantile.probit$ub)
      
      quantile.logit <- quantile.orm(orm.logit, new.data, probs=0.5, se=TRUE)
      quantile.logit <- cbind(quantile.logit$quantile, quantile.logit$lb, quantile.logit$ub)
      
      quantile.loglog <- quantile.orm(orm.loglog, new.data,probs=0.5, se=TRUE)
      quantile.loglog <- cbind(quantile.loglog$quantile, quantile.loglog$lb, quantile.loglog$ub)
      
      quantile.cloglog <- quantile.orm(orm.cloglog, new.data, probs=0.5,se=TRUE)
      quantile.cloglog <- cbind(quantile.cloglog$quantile, quantile.cloglog$lb, quantile.cloglog$ub)
      
      
      quantile.rq <- predict(rq.mod, new.data,  interval =  "confidence", conf.int=0.95)
      
      probit.00[i,] <- quantile.probit[1,]
      probit.10[i,] <- quantile.probit[2,]
      probit.01[i,] <- quantile.probit[3,]
      probit.11[i,] <- quantile.probit[4,]
      
      logit.00[i,] <- quantile.logit[1,]
      logit.10[i,] <- quantile.logit[2,]
      logit.01[i,] <- quantile.logit[3,]
      logit.11[i,] <- quantile.logit[4,]
      
      loglog.00[i,] <- quantile.loglog[1,]
      loglog.10[i,] <- quantile.loglog[2,]
      loglog.01[i,] <- quantile.loglog[3,]
      loglog.11[i,] <- quantile.loglog[4,]
      
      cloglog.00[i,] <- quantile.cloglog[1,]
      cloglog.10[i,] <- quantile.cloglog[2,]
      cloglog.01[i,] <- quantile.cloglog[3,]
      cloglog.11[i,] <- quantile.cloglog[4,]
      
      rq.00[i,] <- quantile.rq[1,]
      rq.10[i,] <- quantile.rq[2,]
      rq.01[i,] <- quantile.rq[3,]
      rq.11[i,] <- quantile.rq[4,]
      
    })
  }
  
  return(list(probit.00=probit.00,
              probit.10=probit.10,
              probit.01=probit.01,
              probit.11=probit.11,
              logit.00=logit.00,
              logit.10=logit.10,
              logit.01=logit.01,
              logit.11=logit.11,
              loglog.00=loglog.00,
              loglog.10=loglog.10,
              loglog.01=loglog.01,
              loglog.11=loglog.11,
              cloglog.00=cloglog.00,
              cloglog.10=cloglog.10,
              cloglog.01=cloglog.01,
              cloglog.11=cloglog.11,
              rq.00=rq.00,
              rq.10=rq.10,
              rq.01=rq.01,
              rq.11=rq.11
  ))
  
}


seed=1
p <- 0.5
alpha <- 0
beta <- c(1, -0.5)
sigma <- 1

sim=10^4
#sim=1

sim4_median.a.n50 <-sim4_median.fun.a(sim=sim, n=50) 
sim4_median.b.n50 <-sim4_median.fun.b(sim=sim, n=50) 
sim4_median.c.n50 <-sim4_median.fun.c(sim=sim, n=50) 
sim4_median.d.n50 <-sim4_median.fun.d(sim=sim, n=50) 
sim4_median.e.n50 <-sim4_median.fun.e(sim=sim, n=50) 
sim4_median.f.n50 <-sim4_median.fun.f(sim=sim, n=50) 
sim4_median.g.n50 <-sim4_median.fun.g(sim=sim, n=50) 
sim4_median.h.n50 <-sim4_median.fun.h(sim=sim, n=50) 

sim4_median.a.n100 <-sim4_median.fun.a(sim=sim, n=100) 
sim4_median.b.n100 <-sim4_median.fun.b(sim=sim, n=100) 
sim4_median.c.n100 <-sim4_median.fun.c(sim=sim, n=100) 
sim4_median.d.n100 <-sim4_median.fun.d(sim=sim, n=100) 
sim4_median.e.n100 <-sim4_median.fun.e(sim=sim, n=100) 
sim4_median.f.n100 <-sim4_median.fun.f(sim=sim, n=100) 
sim4_median.g.n100 <-sim4_median.fun.g(sim=sim, n=100) 
sim4_median.h.n100 <-sim4_median.fun.h(sim=sim, n=100) 


sim4_median.a.n200 <-sim4_median.fun.a(sim=sim, n=200) 
sim4_median.b.n200 <-sim4_median.fun.b(sim=sim, n=200) 
sim4_median.c.n200 <-sim4_median.fun.c(sim=sim, n=200) 
sim4_median.d.n200 <-sim4_median.fun.d(sim=sim, n=200) 
sim4_median.e.n200 <-sim4_median.fun.e(sim=sim, n=200) 
sim4_median.f.n200 <-sim4_median.fun.f(sim=sim, n=200) 
sim4_median.g.n200 <-sim4_median.fun.g(sim=sim, n=200) 
sim4_median.h.n200 <-sim4_median.fun.h(sim=sim, n=200) 


sim4_median.summary <- function(true=as.matrix(new.data)%*%beta, result){
  
  result.00.est <- c(mean(result$probit.00[,1], na.rm=TRUE),
                     mean(result$logit.00[,1], na.rm=TRUE),
                     mean(result$loglog.00[,1], na.rm=TRUE),
                     mean(result$cloglog.00[,1], na.rm=TRUE))
  result.00.coverage <- c(mean(result$probit.00[,2]<=true[1,1]&result$probit.00[,3]>=true[1,1], na.rm=TRUE),
                          mean(result$logit.00[,2]<=true[1,1]&result$logit.00[,3]>=true[1,1], na.rm=TRUE),
                          mean(result$loglog.00[,2]<=true[1,1]&result$loglog.00[,3]>=true[1,1], na.rm=TRUE),
                          mean(result$cloglog.00[,2]<=true[1,1]&result$cloglog.00[,3]>=true[1,1], na.rm=TRUE))
  
  result.00.mse <- c(mean((result$probit.00[,1]-true[1,1])^2, na.rm=TRUE),
                     mean( (result$logit.00[,1]-true[1,1])^2, na.rm=TRUE),
                     mean((result$loglog.00[,1]-true[1,1])^2, na.rm=TRUE),
                     mean((result$cloglog.00[,1]-true[1,1])^2, na.rm=TRUE),
                     mean((result$rq.00[,1]-true[1,1])^2, na.rm=TRUE)  )
  
  result.00.re <- (result.00.mse[5]/result.00.mse)[-5]
  
  
  result.00 <- cbind(format(round(result.00.est,2), nsmall=2), 
                     format(round(result.00.coverage, 3), nsmall=3),
                     format(round(result.00.re,3), nsmall=3))
  result.00.out <- c(format(round(true[1,1],2),nsmall=0),
                     result.00[1,],
                     result.00[2,],
                     result.00[3,],
                     result.00[4,])
  
  ########  10
  result.10.est <- c(mean(result$probit.10[,1], na.rm=TRUE),
                     mean(result$logit.10[,1], na.rm=TRUE),
                     mean(result$loglog.10[,1], na.rm=TRUE),
                     mean(result$cloglog.10[,1], na.rm=TRUE))
  result.10.coverage <- c(mean(result$probit.10[,2]<=true[2,1]&result$probit.10[,3]>=true[2,1], na.rm=TRUE),
                          mean(result$logit.10[,2]<=true[2,1]&result$logit.10[,3]>=true[2,1], na.rm=TRUE),
                          mean(result$loglog.10[,2]<=true[2,1]&result$loglog.10[,3]>=true[2,1], na.rm=TRUE),
                          mean(result$cloglog.10[,2]<=true[2,1]&result$cloglog.10[,3]>=true[2,1], na.rm=TRUE))
  
  result.10.mse <- c(mean((result$probit.10[,1]-true[2,1])^2, na.rm=TRUE),
                     mean( (result$logit.10[,1]-true[2,1])^2, na.rm=TRUE),
                     mean((result$loglog.10[,1]-true[2,1])^2, na.rm=TRUE),
                     mean((result$cloglog.10[,1]-true[2,1])^2, na.rm=TRUE),
                     mean((result$rq.10[,1]-true[2,1])^2, na.rm=TRUE)  )
  
  result.10.re <- (result.10.mse[5]/result.10.mse)[-5]
  
  
  result.10 <- cbind(format(round(result.10.est,2), nsmall=2), 
                     format(round(result.10.coverage, 3), nsmall=3),
                     format(round(result.10.re,3), nsmall=3))
  result.10.out <- c(format(round(true[2,1],2),nsmall=0),
                     result.10[1,],
                     result.10[2,],
                     result.10[3,],
                     result.10[4,])
  
  ###### 01
  result.01.est <- c(mean(result$probit.01[,1], na.rm=TRUE),
                     mean(result$logit.01[,1], na.rm=TRUE),
                     mean(result$loglog.01[,1], na.rm=TRUE),
                     mean(result$cloglog.01[,1], na.rm=TRUE))
  result.01.coverage <- c(mean(result$probit.01[,2]<=true[3,1]&result$probit.01[,3]>=true[3,1], na.rm=TRUE),
                          mean(result$logit.01[,2]<=true[3,1]&result$logit.01[,3]>=true[3,1], na.rm=TRUE),
                          mean(result$loglog.01[,2]<=true[3,1]&result$loglog.01[,3]>=true[3,1], na.rm=TRUE),
                          mean(result$cloglog.01[,2]<=true[3,1]&result$cloglog.01[,3]>=true[3,1], na.rm=TRUE))
  
  result.01.mse <- c(mean((result$probit.01[,1]-true[3,1])^2, na.rm=TRUE),
                     mean( (result$logit.01[,1]-true[3,1])^2, na.rm=TRUE),
                     mean((result$loglog.01[,1]-true[3,1])^2, na.rm=TRUE),
                     mean((result$cloglog.01[,1]-true[3,1])^2, na.rm=TRUE),
                     mean((result$rq.01[,1]-true[3,1])^2, na.rm=TRUE)  )
  
  result.01.re <- (result.01.mse[5]/result.01.mse)[-5]
  
  
  result.01 <- cbind(format(round(result.01.est,2), nsmall=2), 
                     format(round(result.01.coverage, 3), nsmall=3),
                     format(round(result.01.re,3), nsmall=3))
  result.01.out <- c(format(round(true[3,1],2),nsmall=0),
                     result.01[1,],
                     result.01[2,],
                     result.01[3,],
                     result.01[4,])
  ###### 11
  
  result.11.est <- c(mean(result$probit.11[,1], na.rm=TRUE),
                     mean(result$logit.11[,1], na.rm=TRUE),
                     mean(result$loglog.11[,1], na.rm=TRUE),
                     mean(result$cloglog.11[,1], na.rm=TRUE))
  result.11.coverage <- c(mean(result$probit.11[,2]<=true[4,1]&result$probit.11[,3]>=true[4,1], na.rm=TRUE),
                          mean(result$logit.11[,2]<=true[4,1]&result$logit.11[,3]>=true[4,1], na.rm=TRUE),
                          mean(result$loglog.11[,2]<=true[4,1]&result$loglog.11[,3]>=true[4,1], na.rm=TRUE),
                          mean(result$cloglog.11[,2]<=true[4,1]&result$cloglog.11[,3]>=true[4,1], na.rm=TRUE))
  
  result.11.mse <- c(mean((result$probit.11[,1]-true[4,1])^2, na.rm=TRUE),
                     mean( (result$logit.11[,1]-true[4,1])^2, na.rm=TRUE),
                     mean((result$loglog.11[,1]-true[4,1])^2, na.rm=TRUE),
                     mean((result$cloglog.11[,1]-true[4,1])^2, na.rm=TRUE),
                     mean((result$rq.11[,1]-true[4,1])^2, na.rm=TRUE)  )
  
  result.11.re <- (result.11.mse[5]/result.11.mse)[-5]
  
  
  result.11 <- cbind(format(round(result.11.est,2), nsmall=2), 
                     format(round(result.11.coverage, 3), nsmall=3),
                     format(round(result.11.re,3), nsmall=3))
  result.11.out <- c(format(round(true[4,1],2),nsmall=0),
                     result.11[1,],
                     result.11[2,],
                     result.11[3,],
                     result.11[4,])
  
  result.out <- rbind(result.00.out,
                      result.10.out,
                      result.01.out,
                      result.11.out)
  
  return(result.out)
  
}
############### Figure 10; Tables S.13 and S.14 ###################
result.sim4_median.a.n100 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.a.n100)
result.sim4_median.b.n100 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.b.n100)
result.sim4_median.c.n100 <- sim4_median.summary(true=log(-log(0.5))+as.matrix(new.data)%*%beta, result=sim4_median.c.n100)
result.sim4_median.d.n100 <- sim4_median.summary(true=-log(-log(0.5))+as.matrix(new.data)%*%beta, result=sim4_median.d.n100)
result.sim4_median.e.n100 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.e.n100)
result.sim4_median.f.n100 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.f.n100)
result.sim4_median.g.n100 <- sim4_median.summary(true=(qbeta(0.5, 5,2) -5/7)*sqrt((5+2)^2*(5+2+1)/5/2)+as.matrix(new.data)%*%beta, result=sim4_median.g.n100)
result.sim4_median.h.n100 <- sim4_median.summary(true=(qbeta(0.5, 2,5) -2/7)*sqrt((5+2)^2*(5+2+1)/5/2)+as.matrix(new.data)%*%beta, result=sim4_median.h.n100)

sim4_median.out.n100 <- rbind(rep("", 13),
                         result.sim4_median.a.n100,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.b.n100,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.c.n100,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.d.n100,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.e.n100,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.f.n100,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.g.n100,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.h.n100)



############### Tables S.8 and S.9 ###################
result.sim4_median.a.n50 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.a.n50)
result.sim4_median.b.n50 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.b.n50)
result.sim4_median.c.n50 <- sim4_median.summary(true=log(-log(0.5))+as.matrix(new.data)%*%beta, result=sim4_median.c.n50)
result.sim4_median.d.n50 <- sim4_median.summary(true=-log(-log(0.5))+as.matrix(new.data)%*%beta, result=sim4_median.d.n50)
result.sim4_median.e.n50 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.e.n50)
result.sim4_median.f.n50 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.f.n50)
result.sim4_median.g.n50 <- sim4_median.summary(true=(qbeta(0.5, 5,2) -5/7)*sqrt((5+2)^2*(5+2+1)/5/2)+as.matrix(new.data)%*%beta, result=sim4_median.g.n50)
result.sim4_median.h.n50 <- sim4_median.summary(true=(qbeta(0.5, 2,5) -2/7)*sqrt((5+2)^2*(5+2+1)/5/2)+as.matrix(new.data)%*%beta, result=sim4_median.h.n50)

sim4_median.out.n50 <- rbind(rep("", 13),
                         result.sim4_median.a.n50,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.b.n50,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.c.n50,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.d.n50,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.e.n50,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.f.n50,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.g.n50,
                         rep("", 13),
                         rep("", 13),
                         result.sim4_median.h.n50)


############### Tables S.16 and S.17 ###################
result.sim4_median.a.n200 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.a.n200)
result.sim4_median.b.n200 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.b.n200)
result.sim4_median.c.n200 <- sim4_median.summary(true=log(-log(0.5))+as.matrix(new.data)%*%beta, result=sim4_median.c.n200)
result.sim4_median.d.n200 <- sim4_median.summary(true=-log(-log(0.5))+as.matrix(new.data)%*%beta, result=sim4_median.d.n200)
result.sim4_median.e.n200 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.e.n200)
result.sim4_median.f.n200 <- sim4_median.summary(true=as.matrix(new.data)%*%beta, result=sim4_median.f.n200)
result.sim4_median.g.n200 <- sim4_median.summary(true=(qbeta(0.5, 5,2) -5/7)*sqrt((5+2)^2*(5+2+1)/5/2)+as.matrix(new.data)%*%beta, result=sim4_median.g.n200)
result.sim4_median.h.n200 <- sim4_median.summary(true=(qbeta(0.5, 2,5) -2/7)*sqrt((5+2)^2*(5+2+1)/5/2)+as.matrix(new.data)%*%beta, result=sim4_median.h.n200)

sim4_median.out.n200 <- rbind(rep("", 13),
                              result.sim4_median.a.n200,
                              rep("", 13),
                              rep("", 13),
                              result.sim4_median.b.n200,
                              rep("", 13),
                              rep("", 13),
                              result.sim4_median.c.n200,
                              rep("", 13),
                              rep("", 13),
                              result.sim4_median.d.n200,
                              rep("", 13),
                              rep("", 13),
                              result.sim4_median.e.n200,
                              rep("", 13),
                              rep("", 13),
                              result.sim4_median.f.n200,
                              rep("", 13),
                              rep("", 13),
                              result.sim4_median.g.n200,
                              rep("", 13),
                              rep("", 13),
                              result.sim4_median.h.n200)


####### SIMULATIONS IN SUPPLEMENTAL MATERIALS #############
###### Sim 3: orm with an automatic link function selection procedure
library(rms)

rGumbelMin<- function(n, mu=0, sigma=1){
  u <- runif(n, min=0, max=1)
  x <- mu + sigma*log(-log(1-u))
  return(x)
}

rGumbelMax<- function(n, mu=0, sigma=1){
  u <- runif(n, min=0, max=1)
  x <- mu + sigma*(-log(-log(u)))
  return(x)
}


generate.data <- function(seed=1, n=100, 
                          family=c("probit", "logistic", "gumbelMin", "gumbelMax", "cauchy"),
                          z.type=c("continuous", "binary"),
                          p=0.5,sigma.z=1,
                          alpha=0, beta=1, sigma=1){
  set.seed(seed)
  if(z.type[1]=="continuous"){
    z <- rnorm(n, 0, sigma.z)
  } else if (z.type=="binary"){
    z <- sample(c(0,1),size=n, replace=TRUE, prob=c(1-p, p))
  } else stop("Need to specify the covariate Z to be continuous or binary")
  
  if (family[1]=="probit"){
    x <- rnorm(n, alpha+beta*z, sigma)
  } else if (family[1]=="logistic"){
    x <- rlogis(n, alpha+beta*z, sigma)
  } else if (family[1]=="gumbelMin"){
    x <- rGumbelMin(n, mu=alpha+beta*z, sigma=sigma)
  } else if (family[1]=="gumbelMax"){
    x <- rGumbelMax(n, mu=alpha+beta*z, sigma=sigma)
  } else if(family[1]=="cauchy"){
    x <- rcauchy(n, alpha+beta*z, sigma)  
  } else x <- rep(NA, n)
  data <- data.frame(x=x, z=z)
  return(data) 
}

orm.link.sim <- function(sim=100,seed=1, n=50, z.type=c("continuous", "binary"),family=c("probit", "logistic", "gumbelMin","gumbelMax", "cauchy"),
                         p=0.5,sigma.z=1,
                         alpha=0, beta=1, sigma=1){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
  LR.result <- matrix(NA, ncol=4, nrow=sim)
  
  p.result <- matrix(NA, ncol=4, nrow=sim)
  colnames(LR.result) <- c("probit","logit","loglog", "cloglog")
  colnames(p.result) <- c("probit","logit","loglog", "cloglog")
  
  for(i in 1:sim){
    try({
      data <- generate.data(seed=seeds[i], n=n, family=family,
                            z.type=z.type, p=p,sigma.z=sigma.z,
                            alpha=alpha, beta=beta, sigma=sigma)
      
      
      mod.probit <- orm(x ~z, data=data, family=probit)
      mod.logit <- orm(x ~z, data=data)
      
      mod.loglog <- orm(x~z, data=data, family=loglog)
      mod.cloglog <- orm(x~z, data=data, family=cloglog)
      
      
      LR.result[i, ] <- -c(mod.probit$deviance[2], mod.logit$deviance[2],mod.loglog$deviance[2], mod.cloglog$deviance[2])
      p.result[i, ] <- c(anova(mod.probit)["z", "P"],
                         anova(mod.logit)["z", "P"],
                         anova(mod.loglog)["z", "P"],
                         anova(mod.cloglog)["z", "P"])
      

      
    })
  }
  
  return(list(LR.result=LR.result,
              # beta.result=beta.result,
              #se.result=se.result,
              p.result=p.result))
  
}

z.type <- "continuous"

seed <- 1
sim <- 10^4

sigma.z <- 1
alpha <- 0
beta <- 0.1
sigma <- 1


###### true error distribution: normal  ##################

probit.25.null <- orm.link.sim(sim=sim, seed=seed, n=25,family="probit",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0,
                               sigma=sigma)

probit.25.power <- orm.link.sim(sim=sim, seed=seed, n=25,family="probit",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0.1,
                                sigma=sigma)


probit.50.null <- orm.link.sim(sim=sim, seed=seed, n=50,family="probit",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0,
                               sigma=sigma)

probit.50.power <- orm.link.sim(sim=sim, seed=seed, n=50,family="probit",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0.1,
                                sigma=sigma)


probit.100.null <- orm.link.sim(sim=sim, seed=seed, n=100,family="probit",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

probit.100.power <- orm.link.sim(sim=sim, seed=seed, n=100,family="probit",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)


probit.200.null <- orm.link.sim(sim=sim, seed=seed, n=200,family="probit",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

probit.200.power <- orm.link.sim(sim=sim, seed=seed, n=200,family="probit",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)


probit.500.null <- orm.link.sim(sim=sim, seed=seed, n=500,family="probit",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

probit.500.power <- orm.link.sim(sim=sim, seed=seed, n=500,family="probit",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)

probit.1000.null <- orm.link.sim(sim=sim, seed=seed, n=1000,family="probit",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0,
                                 sigma=sigma)
probit.1000.power <- orm.link.sim(sim=sim, seed=seed, n=1000,family="probit",
                                  z.type=z.type,p=p, sigma.z=sigma.z,
                                  alpha=alpha, beta=0.1,
                                  sigma=sigma)


###### true error distribution: logistic ##############
logit.25.null <- orm.link.sim(sim=sim, seed=seed, n=25,family ="logistic",
                              z.type=z.type,p=p, sigma.z=sigma.z,
                              alpha=alpha, beta=0,
                              sigma=sigma)

logit.25.power <- orm.link.sim(sim=sim, seed=seed, n=25,family ="logistic",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0.1,
                               sigma=sigma)


logit.50.null <- orm.link.sim(sim=sim, seed=seed, n=50,family ="logistic",
                              z.type=z.type,p=p, sigma.z=sigma.z,
                              alpha=alpha, beta=0,
                              sigma=sigma)

logit.50.power <- orm.link.sim(sim=sim, seed=seed, n=50,family ="logistic",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0.1,
                               sigma=sigma)


logit.100.null <- orm.link.sim(sim=sim, seed=seed, n=100,family ="logistic",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0,
                               sigma=sigma)

logit.100.power <- orm.link.sim(sim=sim, seed=seed, n=100,family ="logistic",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0.1,
                                sigma=sigma)


logit.200.null <- orm.link.sim(sim=sim, seed=seed, n=200,family ="logistic",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0,
                               sigma=sigma)

logit.200.power <- orm.link.sim(sim=sim, seed=seed, n=200,family ="logistic",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0.1,
                                sigma=sigma)


logit.500.null <- orm.link.sim(sim=sim, seed=seed, n=500,family ="logistic",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0,
                               sigma=sigma)

logit.500.power <- orm.link.sim(sim=sim, seed=seed, n=500,family ="logistic",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0.1,
                                sigma=sigma)


logit.1000.null <- orm.link.sim(sim=sim, seed=seed, n=1000,family ="logistic",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

logit.1000.power <- orm.link.sim(sim=sim, seed=seed, n=1000,family ="logistic",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)


###### true error distribution: extreme value (type I)   ##################

loglog.25.null <- orm.link.sim(sim=sim, seed=seed, n=25,family ="gumbelMin",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0,
                               sigma=sigma)

loglog.25.power <- orm.link.sim(sim=sim, seed=seed, n=25,family ="gumbelMin",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0.1,
                                sigma=sigma)


loglog.50.null <- orm.link.sim(sim=sim, seed=seed, n=50,family ="gumbelMin",
                               z.type=z.type,p=p, sigma.z=sigma.z,
                               alpha=alpha, beta=0,
                               sigma=sigma)

loglog.50.power <- orm.link.sim(sim=sim, seed=seed, n=50,family ="gumbelMin",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0.1,
                                sigma=sigma)


loglog.100.null <- orm.link.sim(sim=sim, seed=seed, n=100,family ="gumbelMin",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

loglog.100.power <- orm.link.sim(sim=sim, seed=seed, n=100,family ="gumbelMin",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)


loglog.200.null <- orm.link.sim(sim=sim, seed=seed, n=200,family ="gumbelMin",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

loglog.200.power <- orm.link.sim(sim=sim, seed=seed, n=200,family ="gumbelMin",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)


loglog.500.null <- orm.link.sim(sim=sim, seed=seed, n=500,family ="gumbelMin",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

loglog.500.power <- orm.link.sim(sim=sim, seed=seed, n=500,family ="gumbelMin",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)


loglog.1000.null <- orm.link.sim(sim=sim, seed=seed, n=1000,family ="gumbelMin",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0,
                                 sigma=sigma)

loglog.1000.power <- orm.link.sim(sim=sim, seed=seed, n=1000,family ="gumbelMin",
                                  z.type=z.type,p=p, sigma.z=sigma.z,
                                  alpha=alpha, beta=0.1,
                                  sigma=sigma)




###### true error distribution: extreme value (type II) ##########

cloglog.25.null <- orm.link.sim(sim=sim, seed=seed, n=25,family = "gumbelMax",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

cloglog.25.power <- orm.link.sim(sim=sim, seed=seed, n=25,family = "gumbelMax",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)


cloglog.50.null <- orm.link.sim(sim=sim, seed=seed, n=50,family = "gumbelMax",
                                z.type=z.type,p=p, sigma.z=sigma.z,
                                alpha=alpha, beta=0,
                                sigma=sigma)

cloglog.50.power <- orm.link.sim(sim=sim, seed=seed, n=50,family = "gumbelMax",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0.1,
                                 sigma=sigma)


cloglog.100.null <- orm.link.sim(sim=sim, seed=seed, n=100,family = "gumbelMax",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0,
                                 sigma=sigma)

cloglog.100.power <- orm.link.sim(sim=sim, seed=seed, n=100,family = "gumbelMax",
                                  z.type=z.type,p=p, sigma.z=sigma.z,
                                  alpha=alpha, beta=0.1,
                                  sigma=sigma)


cloglog.200.null <- orm.link.sim(sim=sim, seed=seed, n=200,family = "gumbelMax",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0,
                                 sigma=sigma)

cloglog.200.power <- orm.link.sim(sim=sim, seed=seed, n=200,family = "gumbelMax",
                                  z.type=z.type,p=p, sigma.z=sigma.z,
                                  alpha=alpha, beta=0.1,
                                  sigma=sigma)


cloglog.500.null <- orm.link.sim(sim=sim, seed=seed, n=500,family = "gumbelMax",
                                 z.type=z.type,p=p, sigma.z=sigma.z,
                                 alpha=alpha, beta=0,
                                 sigma=sigma)

cloglog.500.power <- orm.link.sim(sim=sim, seed=seed, n=500,family = "gumbelMax",
                                  z.type=z.type,p=p, sigma.z=sigma.z,
                                  alpha=alpha, beta=0.1,
                                  sigma=sigma)


cloglog.1000.null <- orm.link.sim(sim=sim, seed=seed, n=1000,family = "gumbelMax",
                                  z.type=z.type,p=p, sigma.z=sigma.z,
                                  alpha=alpha, beta=0,
                                  sigma=sigma)

cloglog.1000.power <- orm.link.sim(sim=sim, seed=seed, n=1000,family = "gumbelMax",
                                   z.type=z.type,p=p, sigma.z=sigma.z,
                                   alpha=alpha, beta=0.1,
                                   sigma=sigma)


####### PRINT OUT TABLE S.1
tab.s1.summary <- function(result){
  
  power <- apply(result$p.result[, c(1,3,4)], 2, FUN=function(x) mean(x<0.05, na.rm=TRUE))
  power <- format(round(power,3), nsmall=3)
  winner <- apply(result$LR.result[, c(1,3,4)], 1, FUN=function(x) which(x==max(x, na.rm=TRUE))[1])
  
  winner.perc <- format(c(mean(winner==1, na.rm=TRUE), 
                          mean(winner==2, na.rm=TRUE),
                          mean(winner==3, na.rm=TRUE))*100, nsmall=2)
  p.winner <- result$p.result[, c(1,3,4)][cbind(1:sim, winner)]
  power.winner <- mean(p.winner<0.05, na.rm=TRUE)
  winner.out <- paste(format(round(power.winner,3), nsmall=3),
                      " (",winner.perc[1], "%, ",
                      winner.perc[2], "%, ",
                      winner.perc[3], "%)",sep="")
  out <- c(winner.out, power)
  names(out)[1] <- "link function selection"

  return(out)
  
}

probit.result <- rbind(tab.s1.summary(probit.25.null),
                      tab.s1.summary(probit.25.power),
                      tab.s1.summary(probit.50.null),
                      tab.s1.summary(probit.50.power),
                      tab.s1.summary(probit.100.null),
                      tab.s1.summary(probit.100.power),
                      tab.s1.summary(probit.200.null),
                      tab.s1.summary(probit.200.power),
                      tab.s1.summary(probit.500.null),
                      tab.s1.summary(probit.500.power),
                      tab.s1.summary(probit.1000.null),
                      tab.s1.summary(probit.1000.power))

logit.result <- rbind(tab.s1.summary(logit.25.null),
                     tab.s1.summary(logit.25.power),
                     tab.s1.summary(logit.50.null),
                     tab.s1.summary(logit.50.power),
                     tab.s1.summary(logit.100.null),
                     tab.s1.summary(logit.100.power),
                     tab.s1.summary(logit.200.null),
                     tab.s1.summary(logit.200.power),
                     tab.s1.summary(logit.500.null),
                     tab.s1.summary(logit.500.power),
                     tab.s1.summary(logit.1000.null),
                     tab.s1.summary(logit.1000.power))


loglog.result <- rbind(tab.1s.summary(loglog.25.null),
                      tab.s1.summary(loglog.25.power),
                      tab.s1.summary(loglog.50.null),
                      tab.s1.summary(loglog.50.power),
                      tab.s1.summary(loglog.100.null),
                      tab.s1.summary(loglog.100.power),
                      tab.s1.summary(loglog.200.null),
                      tab.s1.summary(loglog.200.power),
                      tab.s1.summary(loglog.500.null),
                      tab.s1.summary(loglog.500.power),
                      tab.s1.summary(loglog.1000.null),
                      tab.s1.summary(loglog.1000.power))

cloglog.result <- rbind(tab.1s.summary(cloglog.25.null),
                       tab.s1.summary(cloglog.25.power),
                       tab.s1.summary(cloglog.50.null),
                       tab.s1.summary(cloglog.50.power),
                       tab.s1.summary(cloglog.100.null),
                       tab.s1.summary(cloglog.100.power),
                       tab.s1.summary(cloglog.200.null),
                       tab.s1.summary(cloglog.200.power),
                       tab.s1.summary(cloglog.500.null),
                       tab.s1.summary(cloglog.500.power),
                       tab.s1.summary(cloglog.1000.null),
                       tab.s1.summary(cloglog.1000.power))

tab.s1.out <- rbind(rep("", 4),
                      probit.result,
                      rep("", 4),
                      logit.result,
                      rep("", 4),
                      loglog.result,
                      rep("", 4),
                      cloglog.result
)


rownames(tab.s1.out) <- c("Normal", "~ ~ ~ ~ $n=25$ ~ ~ $H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=50$ ~ ~ $H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$", "~ ~ ~ ~ $n=100$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$", "~ ~ ~ ~ $n=200$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=500$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=1000$ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  $H_1$",
                            "Logistic", "~ ~ ~ ~ $n=25$ ~ ~ $H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=50$ ~ ~ $H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$", "~ ~ ~ ~ $n=100$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$", "~ ~ ~ ~ $n=200$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=500$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=1000$ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  $H_1$",
                            "Extreme Type I", "~ ~ ~ ~ $n=25$ ~ ~ $H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=50$ ~ ~ $H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$", "~ ~ ~ ~ $n=100$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$", "~ ~ ~ ~ $n=200$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=500$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=1000$ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  $H_1$",
                            "Extreme Type II", "~ ~ ~ ~ $n=25$ ~ ~ $H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=50$ ~ ~ $H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$", "~ ~ ~ ~ $n=100$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$", "~ ~ ~ ~ $n=200$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=500$ ~ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ $H_1$","~ ~ ~ ~ $n=1000$ ~$H_0$","~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  $H_1$")

latex(tab.s1.out, where="!htbp", file="",title="True error distribution",caption="link function selection", size="scriptsize")


######### sim 4: orm for outcomes subject to dection limit ###########
generate.data <- function(seed=1, n=100, alpha=3, beta=0.25, sigma=1){
  set.seed(seed)
  z <- rnorm(n, 0, 1)
  y <- rnorm(n, alpha+beta*z, sigma)
  data <- data.frame(z=z, y=y)
  return(data)
}


detection.limit.sim <- function(sim=1000,seed=1, n=100,
                                alpha=3, beta=0.25, sigma=1, cut.off=3){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]  
  pval <- matrix(NA, ncol=5, nrow=sim)
  colnames(pval) <- c("orm (probit)","orm (logit)", "logistic", "ols 1", "ols 2")
  est <- matrix(NA, ncol=5, nrow=sim)
  colnames(est) <- c("orm (probit)","orm (logit)", "logistic", "ols 1", "ols 2")
  est.se <- matrix(NA, ncol=5, nrow=sim)
  colnames(est.se) <- c("orm (probit)","orm (logit)", "logistic", "ols 1", "ols 2")
  for(i in 1:sim){
    try({
      data <- generate.data(seed=seeds[i], n=n, alpha=alpha, beta=beta, sigma=sigma)
      data$y.d <- ifelse(data$y<cut.off, 0, data$y)
      
      m.orm <- orm(y.d~z, data=data, family=probit)
      pval[i,1 ] <- anova(m.orm)["z", "P"]  
      est[i,1] <- m.orm$coef["z"]
      est.se[i,1] <- sqrt(m.orm$var["z", "z"])
      
      m.orm <- orm(y.d~z, data=data)
      pval[i,2 ] <- anova(m.orm)["z", "P"]  
      est[i,2] <- m.orm$coef["z"]
      est.se[i,2] <- sqrt(m.orm$var["z", "z"])
      
      data$y.binary <-  ifelse(data$y<cut.off, 0, 1)
      #m.binary <- lrm(y.binary~z,  data=data)
      m.binary <- glm(y.binary~z, data=data, family="binomial")
      #pval[i,2] <- anova(m.binary)["z", "P"]
      pval[i, 3] <- summary(m.binary)$coef["z", 4]
      est[i, 3] <- summary(m.binary)$coef["z", 1]
      est.se[i, 3] <- summary(m.binary)$coef["z", 2]
      
      data$y.impute <- ifelse(data$y<cut.off, cut.off, data$y)
      m.ols.1 <- ols(y.impute ~z, data)
      pval[i,4] <- anova(m.ols.1)["z", "P"]
      est[i,4] <- m.ols.1$coef["z"]
      est.se[i,4] <- sqrt(m.ols.1$var["z", "z"])
      
      
      m.ols.2 <- ols(y.d ~z, data)
      pval[i,5] <- anova(m.ols.2)["z", "P"]
      est[i,5] <- m.ols.2$coef["z"]
      est.se[i,5] <- sqrt(m.ols.2$var["z", "z"])
      
      
      
    })
  }
  
  return(list(pval=pval,
              est=est,
              est.se=est.se))
  
}

library(rms)
set.seed(1)
alpha=3
beta=0.25
sigma=1
n <- 100
cut.off.p <- c( 0.1, 0.25, 0.5, 0.75, 0.9)
cut.off <- qnorm(cut.off.p, alpha, 1+beta^2)

sim <- 10^4
result.1.null <- detection.limit.sim(sim=sim,seed=1, n=100,
                                     alpha=3, beta=0, sigma=1, 
                                     cut.off=cut.off[1])

result.1.power <- detection.limit.sim(sim=sim,seed=1, n=100,
                                      alpha=3, beta=0.25, sigma=1, 
                                      cut.off=cut.off[1])

result.2.null <- detection.limit.sim(sim=sim,seed=1, n=100,
                                     alpha=3, beta=0, sigma=1, 
                                     cut.off=cut.off[2])

result.2.power <- detection.limit.sim(sim=sim,seed=1, n=100,
                                      alpha=3, beta=0.25, sigma=1, 
                                      cut.off=cut.off[2])


result.3.null <- detection.limit.sim(sim=sim,seed=1, n=100,
                                     alpha=3, beta=0, sigma=1, 
                                     cut.off=cut.off[3])

result.3.power <- detection.limit.sim(sim=sim,seed=1, n=100,
                                      alpha=3, beta=0.25, sigma=1, 
                                      cut.off=cut.off[3])


result.4.null <- detection.limit.sim(sim=sim,seed=1, n=100,
                                     alpha=3, beta=0, sigma=1, 
                                     cut.off=cut.off[4])

result.4.power <- detection.limit.sim(sim=sim,seed=1, n=100,
                                      alpha=3, beta=0.25, sigma=1, 
                                      cut.off=cut.off[4])

result.5.null <- detection.limit.sim(sim=sim,seed=1, n=100,
                                     alpha=3, beta=0, sigma=1, 
                                     cut.off=cut.off[5])

result.5.power <- detection.limit.sim(sim=sim,seed=1, n=100,
                                      alpha=3, beta=0.25, sigma=1, 
                                      cut.off=cut.off[5])


tab.s2.summary <- function(result){
  out <- apply(result$pval, 2, FUN=function(x) mean(x<0.05, na.rm=TRUE)) 
  out <- format(round(out,3), nsmall=3)
  return(out) 
}

tab.s2.out <- rbind(rep("", 5),
                      tab.s2.summary(result.1.null),
                      tab.s2.summary(result.1.power),
                      rep("", 5),
                      tab.s2.summary(result.2.null),
                      tab.s2.summary(result.2.power),
                      rep("", 5),
                      tab.s2.summary(result.3.null),
                      tab.s2.summary(result.3.power),
                      rep("", 5),
                      tab.s2.summary(result.4.null),
                      tab.s2.summary(result.4.power),
                      rep("", 5),
                      tab.s2.summary(result.5.null),
                      tab.s2.summary(result.5.power)                    
)

rownames(tab.s2.out) <- c("$10\\%$ ",
                            "~ ~ ~ ~ ~ ~ $H_0$",
                            "~ ~ ~ ~ ~ ~ $H_1$",
                            "$25\\%$ ",
                            "~ ~ ~ ~ ~ ~ $H_0$",
                            "~ ~ ~ ~ ~ ~ $H_1$",
                            "$50\\%$ ",
                            "~ ~ ~ ~ ~ ~ $H_0$",
                            "~ ~ ~ ~ ~ ~ $H_1$",
                            "$75\\%$",
                            "~ ~ ~ ~ ~ ~ $H_0$",
                            "~ ~ ~ ~ ~ ~ $H_1$",
                            "$90\\%$",
                            "~ ~ ~ ~ ~ ~ $H_0$",
                            "~ ~ ~ ~ ~ ~ $H_1$")

colnames(tab.s2.out) <- c("probit", "logit", "binary outcome", "impute DL", "impute 0" )

###### Table S.26
latex(tab.s2.out, where="!htbp", file="", n.cgroup=c(2, 1, 2),cgroup=c("orm", "logistic regression", "linear regression"),title="below detection limits",caption="detection limits", size="scriptsize")



########################## APPLICAITONS ###############

load("ccasanet_post6.RData")
##### function for calculate probility-scale residuals. 
########This has also been implemented as presid() in "PResiduals" package
presid.orm <- function(object){
  k <- object$non.slopes
  L <- object$linear.predictors
  trans <- object$trans
  cumprob <- if(length(trans)) trans$cumprob else plogis
  if(length(L)==0) stop('you did not use linear.predictors=TRUE for the fit')
  if(length(Y <- object$y) == 0) stop("you did not specify y=TRUE in the fit")
  if(!is.factor(Y)) Y <- factor(Y)
  Y <- unclass(Y) - 1
  cof <- object$coefficients
  if(length(X <- object$x)==0)
    stop("you did not use x=TRUE in the fit")
  
  
  N <- length(Y)
  px <- 1 - cumprob(outer(cof[1:k],
                          as.vector(X %*% cof[- (1:k)]), "+"))
  low.x = rbind(0, px)[cbind(Y + 1L, 1:N)]
  hi.x  = 1 - rbind(px, 1)[cbind(Y + 1L, 1:N)]
  return(low.x - hi.x)
}

#### Or use presid(), see details in the help documents
library(PResiduals)
help(presid)

############### CD4 ############################

###### cumulative probability models with different link functions #####
orm.logit <- orm(post.cd4~ site + male + rcs(age) + init.class + route + 
                   rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
                   preARTAIDS.clinical, data=d, x=TRUE, y=TRUE)

orm.probit <- orm(post.cd4~ site + male + rcs(age) + init.class + route + 
                    rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
                    preARTAIDS.clinical, data=d, x=TRUE, y=TRUE, family=probit)
##### not converge

orm.cloglog <- orm(post.cd4~ site + male + rcs(age) + init.class + route + 
                     rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
                     preARTAIDS.clinical, data=d, x=TRUE, y=TRUE, family=cloglog)


orm.loglog <- orm(post.cd4~ site + male + rcs(age) + init.class + route + 
                    rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
                    preARTAIDS.clinical, data=d, x=TRUE, y=TRUE, family=loglog)### loglog in orm correspoonds to cloglog in conventional cumulative link models

orm.cloglog <- orm.loglog


##### probability-scale residuals
# likelihoods (effectively compare the AIC since number of parameters are the same)
log.L <- c(-orm.probit$deviance[2], -orm.logit$deviance[2], -orm.cloglog$deviance[2])/2
log.L <- c(-orm.probit$deviance[2], -orm.logit$deviance[2], -orm.cloglog$deviance[2])/2

names(log.L) <- c("probit", "logit", "cloglog")
log.L
# 
# probit     logit   cloglog 
# -28796.42 -28709.87 -29185.60 

presid.logit <- presid.orm(orm.logit)

presid.probit <- presid.orm(orm.probit)

presid.cloglog <- presid.orm(orm.cloglog)

##### Box-Cox transformation

library(MASS)

lambda <- boxcox.lambda$x[which(boxcox.lambda$y==max(boxcox.lambda$y))]
box.cox.trans <- function(x,lambda){
  if(lambda==0){
    x <- log(x)
  }else{
    x <- (x^lambda-1)/lambda
  }
  
  return(x)
}


box.cox.trans.inverse <- function(y, lambda){
  if(lambda==0){
    x <- exp(y)
  }else{
    x <- (lambda*y+1)^{1/lambda}
  }
  return(x)
}

d$post.cd4.boxcox <- box.cox.trans(d$post.cd4, lambda)
lm.boxcox <-  ols(post.cd4.boxcox ~ site + male + rcs(age) + init.class + route + 
                    rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
                    preARTAIDS.clinical, data=d, y=TRUE, x=TRUE)

presid.boxcox <- presid(lm.boxcox)


###### Figure 12

postscript("fig13-app_1.eps", height=6, width=8.5,  horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(3,4))
par(mar=c(2.5,2.7,0.5,0.2), mgp = c(1.5, 0.5, 0), oma=c(0.5,3,2,0.5))


alpha.probit <- -orm.probit$coefficients[1:(length(orm.probit$yunique)-1)]
plot(orm.probit$yunique, c(alpha.probit, Inf), type="s", ylab=expression(paste(hat(alpha), "(y)")), xlab="CD4")


alpha.logit <- -orm.logit$coefficients[1:(length(orm.logit$yunique)-1)]
plot(orm.logit$yunique, c(alpha.logit, Inf), type="s", ylab=expression(paste(hat(alpha), "(y)")), xlab="CD4")

alpha.cloglog <- -orm.cloglog$coefficients[1:(length(orm.cloglog$yunique)-1)]
plot(orm.cloglog$yunique, c(alpha.cloglog, Inf), type="s", ylab=expression(paste(hat(alpha), "(y)")), xlab="CD4")

plot(c(min(d$post.cd4):max(d$post.cd4)), box.cox.trans(c(min(d$post.cd4):max(d$post.cd4)), lambda), type="n",
     xlab="CD4", ylab="Estimated transformation")
lines(c(min(d$post.cd4):max(d$post.cd4)), box.cox.trans(c(min(d$post.cd4):max(d$post.cd4)), lambda), lty=1)

order.presid.probit <- presid.probit[order(presid.probit)] 
order.presid.logit <- presid.logit[order(presid.logit)] 
order.presid.cloglog <- presid.cloglog[order(presid.cloglog)] 
order.presid.boxcox <- presid.boxcox[order(presid.boxcox)]
plot( qunif( (1:length(presid.logit))/length(presid.logit), -1, 1),
      order.presid.probit,
      cex=0.1, xlab="Quantiles of Uniform Distribution", ylab="Quantiles of PSR", main="")


#qqline(order.presid.probit, distribution=function(x) qunif(x,-1,1), col=2)
abline(a=0, b=1, col=2)

plot(qunif( (1:length(presid.logit))/length(presid.logit), -1, 1), 
     order.presid.logit, cex=0.1, xlab="Quantiles of Uniform Distribution", ylab="Quantiles of PSR", main="")
abline(a=0, b=1, col=2)
#qqline(order.presid.logit, distribution=function(x) qunif(x,-1,1), col=2)

plot( qunif( (1:length(presid.logit))/length(presid.logit), -1, 1), 
      order.presid.cloglog,cex=0.1, xlab="Quantiles of Uniform Distribution", ylab="Quantiles of PSR", main="")
abline(a=0, b=1, col=2)

#qqline(order.presid.cloglog, distribution=function(x) qunif(x,-1,1), col=2)

plot( qunif( (1:length(presid.logit))/length(presid.logit), -1, 1), 
      order.presid.boxcox,cex=0.1, xlab="Quantiles of Uniform Distribution", ylab="Quantiles of PSR", main="")
abline(a=0, b=1, col=2)

orm.omer.logit <- -orm.logit$coefficients[sapply(d$post.cd4, FUN=function(x) which(orm.logit$yunique==x))] - (predict(orm.logit, type="lp", kint=1) - orm.logit$coefficients[1])
orm.omer.logit <- orm.omer.logit[-which(d$post.cd4==max(d$post.cd4))]

orm.omer.probit <- -orm.probit$coefficients[sapply(d$post.cd4, FUN=function(x) which(orm.probit$yunique==x))] - (predict(orm.probit, type="lp", kint=1) - orm.probit$coefficients[1])
orm.omer.probit <- orm.omer.probit[-which(d$post.cd4==max(d$post.cd4))]

orm.omer.cloglog <- -orm.cloglog$coefficients[sapply(d$post.cd4, FUN=function(x) which(orm.cloglog$yunique==x))] - (predict(orm.cloglog, type="lp", kint=1) - orm.cloglog$coefficients[1])
orm.omer.cloglog <- orm.omer.cloglog[-which(d$post.cd4==max(d$post.cd4))]


order.orm.omer.probit <- orm.omer.probit[order(orm.omer.probit)] 
order.orm.omer.logit <- orm.omer.logit[order(orm.omer.logit)] 
order.orm.omer.cloglog <- orm.omer.cloglog[order(orm.omer.cloglog)] 

plot( qnorm( (1:length(orm.omer.logit))/length(orm.omer.logit)),
      order.orm.omer.probit,
      cex=0.25,, xlab="Quantiles of Normal Distribution", ylab="Quantiles of OMER ", main="")

abline(a=0, b=1, col=2)

plot(qlogis( (1:length(orm.omer.logit))/length(orm.omer.logit)), 
     order.orm.omer.logit, cex=0.25, xlab="Quantiles of Logistic Distribution", ylab="Quantiles of OMER", main="")
abline(a=0, b=1, col=2)


qGumbelMax <- function(p, mu=0, sigma=1){
  x <- mu + -sigma*(log(-log(p)))
  
}

qGumbelMin <- function(p, mu=0, sigma=1){
  x <- mu + sigma*(log(-log(1-p)))
  
}


plot( qGumbelMin((1:length(orm.omer.logit))/length(orm.omer.logit), 0, 1), 
      order.orm.omer.cloglog, cex=0.25, xlab="Quantiles of Extreme Value (I) ", ylab="Quantiles of OMER", main="")
abline(a=0, b=1, col=2)

qqnorm(lm.boxcox$residuals/sd(lm.boxcox$residuals), cex=0.25,  main="", xlab="Quantiles of Normal Distribution", ylab="Quantiles of OMER")
abline(a=0, b=1, col=2)

mtext( "(a)", side=2, outer=TRUE, at=1, cex=1, las=1)
mtext( "(b)", side=2, outer=TRUE, at=2/3, cex=1, las=1)
mtext( "(c)", side=2, outer=TRUE, at=1/3, cex=1, las=1)

mtext("probit", side=3, at=1/8, outer=TRUE, cex=1)
mtext("logit", side=3, at=3/8, outer=TRUE, cex=1)
mtext("cloglog", side=3, at=5/8, outer=TRUE, cex=1)
mtext("Box-Cox", side=3, at=7/8, outer=TRUE, cex=1)

dev.off()

##### Figure 13: illustrate residual-by-predictor plot ####
#### not include nadir cd4
orm.logit.2 <- orm(post.cd4~ site + male  + rcs(age)+init.class + route + 
                     rcs(init.year)+ rcs(log(pre.rna))+
                     preARTAIDS.clinical, data=d, x=TRUE, y=TRUE)
presid.logit.2 <- presid.orm(orm.logit.2)

postscript("fig14-app_2.eps", height=6, width=6,  horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfcol=c(2,2))
par(mar=c(2.5,2.5,0.5,0.5), mgp = c(1.5, 0.5, 0), oma=c(0.5,0.5,2,0.5))
 
plot(sqrt(d$nadir.cd4), presid.logit.2, cex=0.1, type="n",
     xlab="square root of nadir CD4", ylab="PSR", main="")
points(sqrt(d$nadir.cd4), presid.logit.2,pch=1,cex=.5,col=gray(.6))
lines(supsmu(sqrt(d$nadir.cd4), presid.logit.2, bass=1),col=1,lwd=2, cex=1)
abline(h=0,lty=2)

orm.omer.logit.2 <- -orm.logit.2$coefficients[sapply(d$post.cd4, FUN=function(x) which(orm.logit.2$yunique==x))] - (predict(orm.logit.2, type="lp", kint=1) - orm.logit.2$coefficients[1])
orm.omer.logit.2 <- orm.omer.logit.2[-which(d$post.cd4==max(d$post.cd4))]
plot(sqrt(d$nadir.cd4[-which(d$post.cd4==max(d$post.cd4))]), orm.omer.logit.2, cex=0.1, type="n", xlab="square root of nadir CD4", ylab="OMER")
points(sqrt(d$nadir.cd4[-which(d$post.cd4==max(d$post.cd4))]), orm.omer.logit.2,pch=1,cex=.5,col=gray(.6))
lines(supsmu(sqrt(d$nadir.cd4[-which(d$post.cd4==max(d$post.cd4))]), orm.omer.logit.2, bass=1),col=1,lwd=2, cex=1)
abline(h=0,lty=2)

plot(d$nadir.cd4, presid.logit, cex=0.1, type="n",
     xlab="square root of nadir CD4", ylab="PSR", main="")
points(d$nadir.cd4,presid.logit,pch=1,cex=.5,col=gray(.6))
lines(supsmu(d$nadir.cd4, presid.logit, bass=1),col=1,lwd=2, cex=1)
abline(h=0,lty=2)

plot(d$nadir.cd4[-which(d$post.cd4==max(d$post.cd4))], orm.omer.logit, cex=0.1, type="n", xlab="square root of nadir CD4", ylab="OMER ")
points(d$nadir.cd4[-which(d$post.cd4==max(d$post.cd4))], orm.omer.logit,pch=1,cex=.5,col=gray(.6))
lines(supsmu(d$nadir.cd4[-which(d$post.cd4==max(d$post.cd4))], orm.omer.logit, bass=1),col=1,lwd=2, cex=1)
abline(h=0,lty=2)

mtext("not include nadir CD4", side=3, at=1/4, outer=TRUE, cex=1)
mtext("include nadir CD4", side=3, at=3/4, outer=TRUE, cex=1)
dev.off()

##### Figure 14: estimate conditional distribution

##### conditional mean: orm with logit vs Box-Cox
new.age <- (180:820)/10
new.data <- matrix(NA, ncol=dim(orm.logit$x)[2], nrow=length(new.age))

colnames(new.data) <- colnames(orm.logit$x)
for(i in 1:5) new.data[,i] <- rep(orm.logit$x[which(d$site=="brazil")[1], i], length(new.age))

new.data[, "male"] <- 1
new.data[, 7:10] <-  rcspline.eval(new.age, knots=orm.logit$Design$parms$age, inclx=TRUE)

for(i in 11:13) new.data[,i] <- rep(orm.logit$x[which(d$init.class=="NNRTI")[1], i], length(new.age))
for(i in 14:17) new.data[,i] <- rep(orm.logit$x[which(d$route=="homo/bisexual")[1], i], length(new.age))

for(i in 18:21) new.data[,i] <- rep(orm.logit$x[which(d$nadir.cd4==167)[1], i], length(new.age))
for(i in 22:25) new.data[,i] <- rep(orm.logit$x[which(d$init.year==2010)[1], i], length(new.age))
new.data[,26:29] <- matrix(rcspline.eval(log(91728),  knots=orm.logit$Design$parms$pre.rna, inclx=TRUE),nrow=length(new.age), ncol=4, byrow=TRUE)
for(i in 30:31) new.data[,i] <- rep(orm.logit$x[which(d$preARTAIDS.clinical=="not AIDS")[1], i], length(new.age)) 



##### conditional mean as a function of age, fix other covariates at median (for continuous variables) or mode (for discrete variables)
orm.logit <- orm(post.cd4~ site + male + rcs(age) + init.class + route + 
                   rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
                   preARTAIDS.clinical, data=d, x=TRUE, y=TRUE)

est.mean.orm.logit <- mean.orm(orm.logit,new.data,se=TRUE )

new.data.ols <- data.frame(site="brazil", male=1, age=new.age, 
                           init.class="NNRTI", route="homo/bisexual", 
                           nadir.cd4=167, init.year=2010, pre.rna=91728,
                           preARTAIDS.clinical="not AIDS")

est.lm.boxcox.t <- predict(lm.boxcox, new.data.ols, se.fit=TRUE, conf.int=0.95)


est.mean.boxcox <- cbind(box.cox.trans.inverse(est.lm.boxcox.t$linear.predictors, lambda),
                         box.cox.trans.inverse(est.lm.boxcox.t$lower, lambda),
                         box.cox.trans.inverse(est.lm.boxcox.t$upper, lambda))

##### conditional median: orm logit vs median regression

est.median.orm.logit <- quantile.orm(mod=orm.logit,new.data=new.data, probs=0.5,se=TRUE )
library(quantreg)

rq.mod <- rq(post.cd4~ site + male + rcs(age) + init.class + route + 
               rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
               preARTAIDS.clinical, data=d)

X <-  cbind(1, new.data)
pred <-drop(X%*%rq.mod$coef)
V <- summary(rq.mod, cov = TRUE)
df <- V$rdf
tfrac <- qt((1 - 0.95)/2, df)
sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))

est.median.rq <- cbind(pred, pred + tfrac * sdpred %o% 
                         c(1, -1))
colnames(est.median.rq) <- c("fit", "lower", "higher")


est.cdf.orm.logit <- cdf.orm(mod=orm.logit, new.data=new.data, at.y=500)
est.cdf.orm.logit2 <- cdf.orm(mod=orm.logit, new.data=new.data, at.y=350)
#est.cdf.orm.logit3 <- cdf.orm(mod=orm.logit, new.data=new.data, at.y=150)

d$post.cd4.350 <- ifelse(d$post.cd4>350, 1, 0)
d$post.cd4.500 <- ifelse(d$post.cd4>500, 1, 0)
#d$post.cd4.150 <- ifelse(d$post.cd4>150, 1, 0)

##### logistic regression models
lrm.500 <- lrm(post.cd4.500 ~ site + male + rcs(age) + init.class + route + 
                 rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
                 preARTAIDS.clinical, data=d, x=TRUE)
test <- predict(lrm.500, new.data.ols, se.fit=TRUE)
est.cdf.glm.logit <- cbind(plogis(test$linear.predictors), plogis(with(test, linear.predictors + 1.96*cbind(-se.fit,se.fit))))

lrm.350 <- lrm(post.cd4.350 ~ site + male + rcs(age) + init.class + route + 
                 rcs(sqrt(nadir.cd4)) +rcs(init.year)+ rcs(log(pre.rna))+
                 preARTAIDS.clinical, data=d, x=TRUE)
test <- predict(lrm.350, new.data.ols, se.fit=TRUE)
est.cdf.glm.logit2 <- cbind(plogis(test$linear.predictors), plogis(with(test, linear.predictors + 1.96*cbind(-se.fit,se.fit))))




postscript("fig15-app_3.eps", height=4.5, width=8,  horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(2,4))
par(mar=c(2.5,2.5,0.5,0.5), mgp = c(1.5, 0.5, 0), oma=c(0.5,3,2,0.5))
plot(new.age, est.mean.orm.logit[,1], cex= 0.01, ylim=c(280, 410), ylab="CD4", xlab="age", type="n", main="")


polygon(c(new.age, rev(new.age)), c(est.mean.boxcox[,2], rev(est.mean.boxcox[,3])), col =colors()[560], border = NA)
polygon(c(new.age, rev(new.age)), c(est.mean.orm.logit[,3], rev(est.mean.orm.logit[,4])), col = "grey", border = NA)
lines(new.age, est.mean.orm.logit[,1], cex=0.01)
lines(new.age, est.mean.orm.logit[,3], cex=0.01, lty=2)
lines(new.age, est.mean.orm.logit[,4], cex=0.01, lty=2)

lines(new.age, est.mean.boxcox[,1], cex=0.01, col=2)
lines(new.age, est.mean.boxcox[,2], cex=0.01, col=2, lty=2)
lines(new.age, est.mean.boxcox[,3], cex=0.01, col=2, lty=2)

legend("top", c("orm (logit)", "Box-Cox" ), col=c(1,2), lty=c(1,1),bty="n")



plot(new.age, est.median.orm.logit$quantile, cex=0.01, ylim=c(280, 410), ylab="CD4", xlab="age", type="n")


polygon(c(new.age, rev(new.age)), c(est.median.rq[,2], rev(est.median.rq[,3])),col =colors()[560], border = NA)
polygon(c(new.age, rev(new.age)), c(est.median.orm.logit$lb, rev(est.median.orm.logit$ub)), col = 'grey', border = NA)
lines(new.age, est.median.orm.logit$quantile, cex=0.01)
lines(new.age, est.median.orm.logit$lb, cex=0.01, lty=2)
lines(new.age, est.median.orm.logit$ub, cex=0.01, lty=2)

lines(new.age, est.median.rq[,1], cex=0.01, col=2)
lines(new.age, est.median.rq[,2], cex=0.01, col=2, lty=2)
lines(new.age, est.median.rq[,3], cex=0.01, col=2, lty=2)

legend("top", c("orm (logit)", "Median regression"), col=c(1,2), lty=c(1,1),bty="n")

plot(new.age,1-est.cdf.orm.logit2$est[,1], cex=0.01, xlab="age", ylab="Probability", ylim=c(0.25,0.75), typ="n")
polygon(c(new.age, rev(new.age)),  c(est.cdf.glm.logit2[,2], rev(est.cdf.glm.logit2[,3])), col =colors()[560] , border = NA)
polygon(c(new.age, rev(new.age)), c(1-est.cdf.orm.logit2$ub[,1], rev(1-est.cdf.orm.logit2$lb[,1])), col = 'grey', border = NA)


lines(new.age, 1-est.cdf.orm.logit2$est[,1], cex=0.01)
lines(new.age, 1-est.cdf.orm.logit2$ub[,1], cex=0.01, lty=2)
lines(new.age, 1-est.cdf.orm.logit2$lb[,1], cex=0.01, lty=2)

lines(new.age, est.cdf.glm.logit2[,1], cex=0.01, col=2)
lines(new.age, est.cdf.glm.logit2[,2], cex=0.01, col=2, lty=2)
lines(new.age, est.cdf.glm.logit2[,3], cex=0.01, col=2, lty=2)

legend("top", c("orm (logit)", "logistic regression"), col=c(1,2), lty=c(1,1),bty="n")


plot(new.age,1-est.cdf.orm.logit$est[,1], cex=0.01, xlab="age", ylab="Probability", ylim=c(0,0.41), typ="n")
polygon(c(new.age, rev(new.age)),  c(est.cdf.glm.logit[,2], rev(est.cdf.glm.logit[,3])), col = colors()[560], border = NA)
polygon(c(new.age, rev(new.age)), c(1-est.cdf.orm.logit$ub[,1], rev(1-est.cdf.orm.logit$lb[,1])), col = 'grey', border = NA)

lines(new.age, 1-est.cdf.orm.logit$est[,1], cex=0.01)
lines(new.age, 1-est.cdf.orm.logit$ub[,1], cex=0.01, lty=2)
lines(new.age, 1-est.cdf.orm.logit$lb[,1], cex=0.01, lty=2)
lines(new.age, est.cdf.glm.logit[,1], cex=0.01, col=2)
lines(new.age, est.cdf.glm.logit[,2], cex=0.01, col=2, lty=2)
lines(new.age, est.cdf.glm.logit[,3], cex=0.01, col=2, lty=2)

legend("top", c("orm (logit)", "logistic regression"), col=c(1,2), lty=c(1,1),bty="n")

mtext("Mean", side=3, at=1/8, outer=TRUE, cex=1)
mtext("Median", side=3, at=3/8, outer=TRUE, cex=1)
mtext("P[CD4>350|Z]", side=3, at=5/8, outer=TRUE, cex=1)
mtext("P[CD4>500|Z]", side=3, at=7/8, outer=TRUE, cex=1)




####### by treatment class
new.data2 <- new.data[1:4,]

new.data2[,7:10] <- matrix(rcspline.eval(35,  knots=orm.logit$Design$parms$age, inclx=TRUE),nrow=4, ncol=4, byrow=TRUE)

new.data2[1, 11:13] <- orm.logit$x[which(d$init.class=="Boosted PI")[1], 11:13]

new.data2[2, 11:13] <- orm.logit$x[which(d$init.class=="NNRTI")[1], 11:13]
new.data2[3, 11:13] <- orm.logit$x[which(d$init.class=="Unboosted PI")[1], 11:13]
new.data2[4, 11:13] <- orm.logit$x[which(d$init.class=="Other")[1], 11:13]


est.mean.orm.logit <- mean.orm(orm.logit, new.data2,se=TRUE )

new.data.ols <- data.frame(site="brazil", male=1, age=35, 
                           init.class=c("Boosted PI","NNRTI","Unboosted PI","Other"), route="homo/bisexual", 
                           nadir.cd4=167, init.year=2010, pre.rna=91728,
                           preARTAIDS.clinical="not AIDS")



est.lm.boxcox.t <- predict(lm.boxcox, new.data.ols, se.fit=TRUE, conf.int=0.95)


est.mean.boxcox <- cbind(box.cox.trans.inverse(est.lm.boxcox.t$linear.predictors, lambda),
                         box.cox.trans.inverse(est.lm.boxcox.t$lower, lambda),
                         box.cox.trans.inverse(est.lm.boxcox.t$upper, lambda))

X <-  cbind(1, new.data2)
pred <-drop(X%*%rq.mod$coef)
V <- summary(rq.mod, cov = TRUE)
df <- V$rdf
tfrac <- qt((1 - 0.95)/2, df)
sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))

est.median.rq <- cbind(pred, pred + tfrac * sdpred %o% 
                         c(1, -1))
colnames(est.median.rq) <- c("fit", "lower", "higher")



est.median.orm.logit <- quantile.orm(mod=orm.logit,new.data=new.data2, probs=0.5,se=TRUE )


test <- predict(lrm.500, new.data.ols, se.fit=TRUE)
est.cdf.glm.logit <- cbind(plogis(test$linear.predictors), plogis(with(test, linear.predictors + 1.96*cbind(-se.fit,se.fit))))


test <- predict(lrm.350, new.data.ols, se.fit=TRUE)
est.cdf.glm.logit2 <- cbind(plogis(test$linear.predictors), plogis(with(test, linear.predictors + 1.96*cbind(-se.fit,se.fit))))


est.cdf.orm.logit <- cdf.orm(mod=orm.logit, new.data=new.data2, at.y=500)
est.cdf.orm.logit2 <- cdf.orm(mod=orm.logit, new.data=new.data2, at.y=350)




plot(1:4,est.mean.orm.logit[,1], type="n",
     axes=FALSE, ylab="CD4", xlab="", xlim=c(0.5, 4.5), ylim=c(300, 450), main="")
#abline(v=1)
axis(2, at=c(300,350,400, 450), labels=c("300", "350","400", "450"))
segments(1:4,est.mean.orm.logit[,3], 1:4, est.mean.orm.logit[,4])
points(1:4,est.mean.orm.logit[,1], pch="-", cex=1)

segments((1:4)+0.2,est.mean.boxcox[,2], (1:4)+0.2, est.mean.boxcox[,3], col=2)
points((1:4)+0.2,est.mean.boxcox[,1], pch="-", cex=1, col=2)
#segments(hernan.700.sw[,2], cd4+3, hernan.700.sw[,3], cd4+3, col=2)
#points(hernan.700.sw[,1], cd4+3, pch="|", cex=0.4, col=2)
legend("top", c("orm (logit)", "Box-Cox"), col=c(1,2), lty=c(1,1),bty="n")
axis(1, 1:4, labels=c("BPI","NNRTI","UBPI","Other"),
     las=1, tck=-0.01, col=1, cex.axis=0.8)
box()


plot(1:4,est.median.orm.logit$quantile, type="n",
     axes=FALSE, ylab="CD4", xlab="", xlim=c(0.5, 4.5), ylim=c(300, 450), main="")
#abline(v=1)
axis(2, at=c(300,350,400, 450), labels=c("300", "350","400", "450"))
segments(1:4,est.median.orm.logit$lb, 1:4, est.median.orm.logit$ub)
points(1:4,est.median.orm.logit$quantile, pch="-", cex=1)

segments((1:4)+0.2,est.median.rq[,2], (1:4)+0.2, est.median.rq[,3], col=2)
points((1:4)+0.2,est.median.rq[,1], pch="-", cex=1, col=2)
#segments(hernan.700.sw[,2], cd4+3, hernan.700.sw[,3], cd4+3, col=2)
#points(hernan.700.sw[,1], cd4+3, pch="|", cex=0.4, col=2)
legend("top", c("orm (logit)", "Median regression"), col=c(1,2), lty=c(1,1),bty="n")
axis(1, 1:4, labels=c("BPI","NNRTI","UBPI","Other"),
     las=1, tck=-0.01, col=1, cex.axis=0.8)
box()

plot(1:4,1-est.cdf.orm.logit2$est, type="n",
     axes=FALSE, ylab="probability", xlab="", xlim=c(0.5, 4.5), ylim=c(0.3, 0.8), main="")

axis(2, at=c(0.4, 0.5, 0.6,0.7), labels=c("0.4", "0.5","0.6","0.7"))
segments(1:4,1-est.cdf.orm.logit2$lb, 1:4, 1-est.cdf.orm.logit2$ub)
points(1:4,1-est.cdf.orm.logit2$est, pch="-", cex=1)

segments((1:4)+0.2,est.cdf.glm.logit2[,2], (1:4)+0.2, est.cdf.glm.logit2[,3], col=2)
points((1:4)+0.2,est.cdf.glm.logit2[,1], pch="-", cex=1, col=2)

legend("top", c("orm (logit)", "Logistic regression"), col=c(1,2), lty=c(1,1),bty="n")
axis(1, 1:4, labels=c("BPI","NNRTI","UBPI","Other"),
     las=1, tck=-0.01, col=1, cex.axis=0.8)
box()

plot(1:4,1-est.cdf.orm.logit$est, type="n",
     axes=FALSE, ylab="probability", xlab="", xlim=c(0.5, 4.5), ylim=c(0.05, 0.4), main="")
axis(2, at=c( 0.1, 0.2, 0.3, 0.4), labels=c("0.1", "0.2","0.3", "0.4"))
segments(1:4,1-est.cdf.orm.logit$lb, 1:4, 1-est.cdf.orm.logit$ub)
points(1:4,1-est.cdf.orm.logit$est, pch="-", cex=1)
segments((1:4)+0.2,est.cdf.glm.logit[,2], (1:4)+0.2, est.cdf.glm.logit[,3], col=2)
points((1:4)+0.2,est.cdf.glm.logit[,1], pch="-", cex=1, col=2)
legend("top", c("orm (logit)", "Logistic regression"), col=c(1,2), lty=c(1,1),bty="n")
axis(1, 1:4, labels=c("BPI","NNRTI","UBPI","Other"),
     las=1, tck=-0.01, col=1, cex.axis=0.8)
box()

dev.off()


##### viral load as outcome

###########################################################
########## orm with different link functions ##############
###########################################################
orm.logit <- orm(post.rna.d~ site + male + rcs(age,4) + init.class + route + 
                   rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
                   preARTAIDS.clinical, data=d, x=TRUE, y=TRUE)

orm.probit <- orm(post.rna.d~ site + male + rcs(age,4) + init.class + route + 
                    rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
                    preARTAIDS.clinical, data=d, x=TRUE, y=TRUE, family=probit)

orm.cloglog <- orm(post.rna.d~ site + male + rcs(age,4) + init.class + route + 
                     rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
                     preARTAIDS.clinical, data=d, x=TRUE, y=TRUE, family=cloglog)


orm.loglog <- orm(post.rna.d~ site + male + rcs(age,4) + init.class + route + 
                    rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
                    preARTAIDS.clinical, data=d, x=TRUE, y=TRUE, family=loglog,eps=1e-3)


log.L <- c(-orm.probit$deviance[2], -orm.logit$deviance[2], -orm.loglog$deviance[2], -orm.cloglog$deviance[2])/2
### cloglog in orm correspoonds to loglog in conventional cumulative link models
names(log.L) <- c("probit", "logit", "loglog", "cloglog")
log.L
# probit     logit    loglog   cloglog 
# -5213.691 -5185.259 -5245.030 -5162.703 

presid.logit <- presid.orm(orm.logit)

presid.probit <- presid.orm(orm.probit)

presid.cloglog <- presid.orm(orm.cloglog)

presid.loglog <- presid.orm(orm.loglog)

###############################################################
###### Figure S.8 estimated transformation and PSR QQ plot ##############
###############################################################

postscript("fig19-app_4.eps", height=6, width=8,  horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(3,4))
par(mar=c(2.5,2.7,0.5,0.2), mgp = c(1.5, 0.5, 0), oma=c(0.5,3,2,0.5))


alpha.probit <- -orm.probit$coefficients[1:(length(orm.probit$yunique)-1)]


plot(c(0, log(orm.probit$yunique[-1],10)), c(alpha.probit, Inf), type="s",
     ylab=expression(paste(hat(alpha))), xlab="log(viral load)", col=1, axes=FALSE)
abline(v=log(orm.probit$yunique[-1],10)[1], lty=2)
axis(side=2)
axis(side=1, at=c(3,4,5), labels=c("3", "4", "5"))
box()
mtext("undetectable", side=1, at=1, cex=0.6)


alpha.logit <- -orm.logit$coefficients[1:(length(orm.logit$yunique)-1)]

plot(c(0, log(orm.logit$yunique[-1],10)), c(alpha.logit, Inf), type="s",
     ylab=expression(paste(hat(alpha))),  xlab="log(viral load)", col=1, axes=FALSE)
abline(v=log(orm.logit$yunique[-1],10)[1], lty=2)
axis(side=2)
axis(side=1, at=c(3, 4,5), labels=c("3", "4", "5"))
box()
mtext("undetectable", side=1, at=1, cex=0.6)



alpha.loglog <- -orm.loglog$coefficients[1:(length(orm.loglog$yunique)-1)]
plot(c(0, log(orm.loglog$yunique[-1],10)), c(alpha.loglog, Inf), type="s",
     ylab=expression(paste(hat(alpha))),  xlab="log(viral load)", col=1, axes=FALSE)
abline(v=log(orm.loglog$yunique[-1],10)[1], lty=2)
axis(side=2)
axis(side=1, at=c(3, 4,5), labels=c("3", "4", "5"))
box()
mtext("undetectable", side=1, at=1, cex=0.6)


alpha.cloglog <- -orm.cloglog$coefficients[1:(length(orm.cloglog$yunique)-1)]
plot(c(0, log(orm.cloglog$yunique[-1],10)), c(alpha.cloglog, Inf), type="s",
     ylab=expression(paste(hat(alpha))),  xlab="log(viral load)", col=1, axes=FALSE)
abline(v=log(orm.cloglog$yunique[-1],10)[1], lty=2)
axis(side=2)
axis(side=1, at=c(3, 4,5), labels=c("3", "4", "5"))
box()
mtext("undetectable", side=1, at=1, cex=0.6)



order.presid.probit <- presid.probit[order(presid.probit)] 
order.presid.logit <- presid.logit[order(presid.logit)] 
order.presid.loglog <- presid.loglog[order(presid.loglog)]
order.presid.cloglog <- presid.cloglog[order(presid.cloglog)]


plot( qunif( (1:length(presid.probit))/length(presid.probit), -1, 1),
      order.presid.probit,
      cex=0.1, xlab="Quantiles of Uniform Distribution", 
      ylab="Quantiles of PSR", main="",
      col=ifelse( (d$post.rna.d==0)[order(presid.probit)], 4,1), ylim=c(-1,1))

abline(a=0, b=1, col=2)
legend("top", c("undetectable", "detectable"), col=c(4,1), lty=1, bty="n")


plot(qunif( (1:length(presid.logit))/length(presid.logit), -1, 1), 
     order.presid.logit, cex=0.1, 
     xlab="Quantiles of Uniform Distribution", 
     ylab="Quantiles of PSR", main="",
     col=ifelse( (d$post.rna.d==0)[order(presid.logit)], 4,1), ylim=c(-1, 1))
abline(a=0, b=1, col=2)
legend("top", c("undetectable", "detectable"), col=c(4,1), lty=1, bty="n")


plot( qunif( (1:length(presid.loglog))/length(presid.loglog), -1, 1), 
      order.presid.loglog,cex=0.1, xlab="Quantiles of Uniform Distribution", ylab="Quantiles of PSR", main="",
      col=ifelse( (d$post.rna.d==0)[order(presid.loglog)], 4,1), ylim=c(-1,1))
abline(a=0, b=1, col=2)
legend("top", c("undetectable", "detectable"), col=c(4,1), lty=1, bty="n")

plot( qunif( (1:length(presid.cloglog))/length(presid.cloglog), -1, 1), 
      order.presid.cloglog,cex=0.1, 
      xlab="Quantiles of Uniform Distribution",
      ylab="Quantiles of PSR", main="",
      col=ifelse( (d$post.rna.d==0)[order(presid.cloglog)], 4,1), ylim=c(-1,1))
abline(a=0, b=1, col=2)
legend("top", c("undetectable", "detectable"), col=c(4,1), lty=1, bty="n")


mtext( "(a)", side=2, outer=TRUE, at=1, cex=1, las=1)
mtext( "(b)", side=2, outer=TRUE, at=2/3, cex=1, las=1)
mtext( "(c)", side=2, outer=TRUE, at=1/3, cex=1, las=1)


mtext("probit", side=3, at=1/8, outer=TRUE, cex=1)
mtext("logit", side=3, at=3/8, outer=TRUE, cex=1)
mtext("cloglog", side=3, at=5/8, outer=TRUE, cex=1)
mtext("loglog", side=3, at=7/8, outer=TRUE, cex=1)



orm.omer.logit <- -orm.logit$coefficients[sapply(d$post.rna.d, FUN=function(x) which(orm.logit$yunique==x))] - (predict(orm.logit, type="lp", kint=1) - orm.logit$coefficients[1])
orm.omer.logit <- orm.omer.logit[-which(d$post.rna.d==max(d$post.rna.d))]

orm.omer.probit <- -orm.probit$coefficients[sapply(d$post.rna.d, FUN=function(x) which(orm.probit$yunique==x))] - (predict(orm.probit, type="lp", kint=1) - orm.probit$coefficients[1])
orm.omer.probit <- orm.omer.probit[-which(d$post.rna.d==max(d$post.rna.d))]

orm.omer.cloglog <- -orm.cloglog$coefficients[sapply(d$post.rna.d, FUN=function(x) which(orm.cloglog$yunique==x))] - (predict(orm.cloglog, type="lp", kint=1) - orm.cloglog$coefficients[1])
orm.omer.cloglog <- orm.omer.cloglog[-which(d$post.rna.d==max(d$post.rna.d))]

orm.omer.loglog <- -orm.loglog$coefficients[sapply(d$post.rna.d, FUN=function(x) which(orm.loglog$yunique==x))] - (predict(orm.loglog, type="lp", kint=1) - orm.loglog$coefficients[1])
orm.omer.loglog <- orm.omer.loglog[-which(d$post.rna.d==max(d$post.rna.d))]


order.orm.omer.probit <- orm.omer.probit[order(orm.omer.probit)] 
order.orm.omer.logit <- orm.omer.logit[order(orm.omer.logit)] 
order.orm.omer.loglog <- orm.omer.loglog[order(orm.omer.loglog)]
order.orm.omer.cloglog <- orm.omer.cloglog[order(orm.omer.cloglog)] 




plot( qnorm( (1:length(orm.omer.probit))/length(orm.omer.probit)),
      order.orm.omer.probit,
      cex=0.25, xlab="Quantiles of Normal Distribution", 
      ylab="Quantiles of OMER", main="",
      col=ifelse( (d$post.rna.d[-which(d$post.rna.d==max(d$post.rna.d)) ]==0)[order(orm.omer.probit)], 4,1))

abline(a=0, b=1, col=2)
legend("top", c("undetectable", "detectable"), col=c(4,1), lty=1, bty="n")


plot(qlogis( (1:length(orm.omer.logit))/length(orm.omer.logit)), 
     order.orm.omer.logit, cex=0.25, 
     xlab="Quantiles of Logistic Distribution", 
     ylab="Quantiles of OMER", main="",
     col=ifelse( (d$post.rna.d[-which(d$post.rna.d==max(d$post.rna.d)) ]==0)[order(orm.omer.logit)], 4,1))
abline(a=0, b=1, col=2)
legend("top", c("undetectable", "detectable"), col=c(4,1), lty=1, bty="n")



qGumbelMax <- function(p, mu=0, sigma=1){
  x <- mu + -sigma*(log(-log(p)))
  
}


qGumbelMin <- function(p, mu=0, sigma=1){
  x <- mu + sigma*(log(-log(1-p)))
  
}


plot( qGumbelMin((1:length(orm.omer.logit))/length(orm.omer.logit), 0, 1), 
      order.orm.omer.loglog,cex=0.25, xlab="Quantiles of Extreme Value (I) ",
      ylab="Quantiles of OMER", main="",
      col=ifelse( (d$post.rna.d[-which(d$post.rna.d==max(d$post.rna.d)) ]==0)[order(orm.omer.loglog)], 4,1))
#qqline(order.orm.omer.cloglog, distribution=function(x) qGumbelMax(x,), col=3)
abline(a=0, b=1, col=2)
legend("top", c("undetectable", "detectable"), col=c(4,1), lty=1, bty="n")


plot( qGumbelMax((1:length(orm.omer.logit))/length(orm.omer.logit), 0, 1), 
      order.orm.omer.cloglog,cex=0.25,
      xlab="Quantiles of Extreme Value (II) ", 
      ylab="Quantiles of OMER", main="",
      col=ifelse( (d$post.rna.d[-which(d$post.rna.d==max(d$post.rna.d)) ]==0)[order(orm.omer.cloglog)], 4,1))

abline(a=0, b=1, col=2)
legend("top", c("undetectable", "detectable"), col=c(4,1), lty=1, bty="n")


dev.off()

###################################################
####### Figure S.9: PSR by predictor plot #########
##################################################

######### model not include age ##############
orm.logit.2 <- orm(post.rna.d~ site + male  + init.class + route + rcs(log(pre.rna),4) +
                     rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+
                     preARTAIDS.clinical, data=d, x=TRUE, y=TRUE)

presid.logit.2 <- presid.orm(orm.logit.2)

orm.cloglog.2 <- orm(post.rna.d~ site + male  + init.class + route + rcs(log(pre.rna),4) +
                       rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+
                       preARTAIDS.clinical, data=d,family=cloglog, x=TRUE, y=TRUE)

presid.cloglog.2 <- presid.orm(orm.logit.2)


postscript("fig20-app_5.eps", height=6, width=6,  horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfcol=c(2,2))
par(mar=c(2.5,2.5,0.5,0.5), mgp = c(1.5, 0.5, 0), oma=c(0.5,0.5,2,0.5))

plot(d$age, presid.logit.2, cex=0.1, type="n", ylim=c(-1, 1),
     xlab="age", ylab="PSR (logit) ", main="")
points(d$age, presid.logit.2,cex=.5,col=gray(.6), pch=ifelse(d$post.rna.detect==0, 1,3))
lines(supsmu(d$age, presid.logit.2, bass=1),col=1,lwd=2, cex=1)

#lines(lowess(presid.logit.2~log(d$pre.rna,10)),lwd=2)
abline(h=0,lty=2)
legend("bottomleft", c("undetectable", "detectable"), cex=0.8,pch=c(1,3), bty="n", col=gray(0.6))


plot(d$age, presid.cloglog.2, cex=0.1, type="n", ylim=c(-1, 1),
     xlab="age", ylab="PSR (loglog)", main="")
points(d$age, presid.cloglog.2,cex=.5,col=gray(.6), pch=ifelse(d$post.rna.detect==0, 1,3))
lines(supsmu(d$age, presid.cloglog.2, bass=1),col=1,lwd=2, cex=1)

#lines(lowess(presid.logit.2~log(d$pre.rna,10)),lwd=2)
abline(h=0,lty=2)
legend("bottomleft", c("undetectable", "detectable"), cex=0.8,pch=c(1,3), bty="n", col=gray(0.6))



plot(d$age, presid.logit, cex=0.1, type="n", ylim=c(-1, 1),
     xlab="age", ylab="PSR (logit)",
     main="")
points(d$age, presid.logit,cex=0.5,col=gray(.6), pch=ifelse(d$post.rna.detect==0, 1,3) )
lines(supsmu(d$age, presid.logit, bass=1),col=1,lwd=2, cex=1)
#lines(lowess(presid.logit~log(d$pre.rna,10)),lwd=2)
abline(h=0,lty=2)
legend("bottomleft", c("undetectable", "detectable"), cex=0.8,pch=c(1,3), bty="n", col=gray(0.6))


plot(d$age, presid.cloglog, cex=0.1, type="n", ylim=c(-1, 1),
     xlab="age", ylab="PSR (loglog) ", main="")
points(d$age, presid.cloglog,cex=0.5,col=gray(.6), pch=ifelse(d$post.rna.detect==0, 1,3) )
lines(supsmu(d$age, presid.cloglog, bass=1),col=1,lwd=2, cex=1)
#lines(lowess(presid.logit~log(d$pre.rna,10)),lwd=2)
abline(h=0,lty=2)
legend("bottomleft", c("undetectable", "detectable"), cex=0.8,pch=c(1,3), bty="n", col=gray(0.6))



mtext("not include age", side=3, at=1/4, outer=TRUE, cex=1)
mtext("include age", side=3, at=3/4, outer=TRUE, cex=1)

dev.off()

############################################################################
####### Figure 15: performance of orm on estimation conditional quantiles, CDFs #########
#############################################################################
postscript("fig16-app_6.eps", height=4, width=6,  horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(2,3))
par(mar=c(2.5,2.5,0.5,0.5), mgp = c(1.5, 0.5, 0), oma=c(0.5,3,2,0.5))


new.age <- (180:820)/10
new.data <- matrix(NA, ncol=dim(orm.logit$x)[2], nrow=length(new.age))

colnames(new.data) <- colnames(orm.logit$x)

for(i in 1:5) new.data[,i] <- rep(orm.logit$x[which(d$site=="brazil")[1], i], length(new.age))

new.data[, "male"] <- 1
new.data[, 7:9] <-  rcspline.eval(new.age, knots=orm.logit$Design$parms$age, inclx=TRUE)

for(i in 10:12) new.data[,i] <- rep(orm.logit$x[which(d$init.class=="NNRTI")[1], i], length(new.age))
for(i in 13:16) new.data[,i] <- rep(orm.logit$x[which(d$route=="homo/bisexual")[1], i], length(new.age))

for(i in 17:19) new.data[,i] <- rep(orm.logit$x[which(d$nadir.cd4==167)[1], i], length(new.age))
for(i in 20:22) new.data[,i] <- rep(orm.logit$x[which(d$init.year==2010)[1], i], length(new.age))
new.data[,23:25] <- matrix(rcspline.eval(log(91728),  knots=orm.logit$Design$parms$pre.rna, inclx=TRUE),nrow=length(new.age), ncol=3, byrow=TRUE)
for(i in 26:27) new.data[,i] <- rep(orm.logit$x[which(d$preARTAIDS.clinical=="not AIDS")[1], i], length(new.age)) 

new.data.ols <- data.frame(site="brazil", male=1, age=new.age, 
                           init.class="NNRTI", route="homo/bisexual", 
                           nadir.cd4=167, init.year=2010, pre.rna=91728,
                           preARTAIDS.clinical="not AIDS")

est.cdf.orm.logit <- cdf.orm(mod=orm.logit, new.data=new.data, at.y=0)
est.cdf.orm.cloglog <- cdf.orm(mod=orm.cloglog, new.data=new.data, at.y=0)
binary.mod <- lrm(post.rna.detect~ site + male + rcs(age,4) + init.class + route + 
                    rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
                    preARTAIDS.clinical, data=d, x=TRUE, y=TRUE)

test <- predict(binary.mod, new.data.ols, se.fit=TRUE)
est.cdf.binary <- cbind(plogis(test$linear.predictors), plogis(with(test, linear.predictors + 1.96*cbind(-se.fit,se.fit))))

plot(new.age,1-est.cdf.orm.logit$est[,1], cex=0.01, 
     xlab="age", ylab="Probability", ylim=c(0,0.15), typ="n")
polygon(c(new.age, rev(new.age)),  c(est.cdf.binary[,2], rev(est.cdf.binary[,3])), col = colors()[560], border = NA)
polygon(c(new.age, rev(new.age)), c(1-est.cdf.orm.logit$ub[,1], rev(1-est.cdf.orm.logit$lb[,1])), col = colors()[618], border = NA)
polygon(c(new.age, rev(new.age)), c(1-est.cdf.orm.cloglog$ub[,1], rev(1-est.cdf.orm.cloglog$lb[,1])), col = 'grey', border = NA)
#polygon(c(new.age, rev(new.age)),  c(est.cdf.boxcox[,2], rev(est.cdf.boxcox[,3])), col = "grey", border = NA)

lines(new.age, 1-est.cdf.orm.cloglog$est[,1], cex=0.01, col=1)
lines(new.age, 1-est.cdf.orm.cloglog$ub[,1], cex=0.01, lty=2, col=1)
lines(new.age, 1-est.cdf.orm.cloglog$lb[,1], cex=0.01, lty=2, col=1)

lines(new.age, 1-est.cdf.orm.logit$est[,1], cex=0.01, col=4)
lines(new.age, 1-est.cdf.orm.logit$ub[,1], cex=0.01, lty=2, col=4)
lines(new.age, 1-est.cdf.orm.logit$lb[,1], cex=0.01, lty=2, col=4)


lines(new.age, est.cdf.binary[,1], cex=0.01, col=2)
lines(new.age, est.cdf.binary[,2], cex=0.01, col=2, lty=2)
lines(new.age, est.cdf.binary[,3], cex=0.01, col=2, lty=2)

legend("top", c("orm (loglog)","orm (logit)", "logistic regression"), col=c(1,4,2), lty=c(1,1),bty="n")


est.cdf.orm.logit.g1000<- cdf.orm(mod=orm.logit, new.data=new.data, at.y=1000)
est.cdf.orm.cloglog.g1000 <- cdf.orm(mod=orm.cloglog, new.data=new.data, at.y=1000)

post.rna.g1000 <- ifelse(d$post.rna.d>1000, 1, 0)

binary.mod.g1000 <- lrm(post.rna.g1000~ site + male + rcs(age,4) + init.class + route + 
                          rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
                          preARTAIDS.clinical, data=d, x=TRUE, y=TRUE)

test <- predict(binary.mod.g1000, new.data.ols, se.fit=TRUE)

est.cdf.binary.g1000 <- cbind(plogis(test$linear.predictors), plogis(with(test, linear.predictors + 1.96*cbind(-se.fit,se.fit))))

plot(new.age,1-est.cdf.orm.logit.g1000$est[,1], cex=0.01, 
     xlab="age", ylab="Probability", ylim=c(0,0.15), typ="n")
polygon(c(new.age, rev(new.age)),  c(est.cdf.binary.g1000[,2], rev(est.cdf.binary.g1000[,3])), col = colors()[560], border = NA)
polygon(c(new.age, rev(new.age)), c(1-est.cdf.orm.logit.g1000$ub[,1], rev(1-est.cdf.orm.logit.g1000$lb[,1])), col = colors()[618], border = NA)
polygon(c(new.age, rev(new.age)), c(1-est.cdf.orm.cloglog.g1000$ub[,1], rev(1-est.cdf.orm.cloglog.g1000$lb[,1])), col = 'grey', border = NA)
#polygon(c(new.age, rev(new.age)),  c(est.cdf.boxcox[,2], rev(est.cdf.boxcox[,3])), col = "grey", border = NA)

lines(new.age, 1-est.cdf.orm.cloglog.g1000$est[,1], cex=0.01, col=1)
lines(new.age, 1-est.cdf.orm.cloglog.g1000$ub[,1], cex=0.01, lty=2, col=1)
lines(new.age, 1-est.cdf.orm.cloglog.g1000$lb[,1], cex=0.01, lty=2, col=1)

lines(new.age, 1-est.cdf.orm.logit.g1000$est[,1], cex=0.01, col=4)
lines(new.age, 1-est.cdf.orm.logit.g1000$ub[,1], cex=0.01, lty=2, col=4)
lines(new.age, 1-est.cdf.orm.logit.g1000$lb[,1], cex=0.01, lty=2, col=4)


lines(new.age, est.cdf.binary.g1000[,1], cex=0.01, col=2)
lines(new.age, est.cdf.binary.g1000[,2], cex=0.01, col=2, lty=2)
lines(new.age, est.cdf.binary.g1000[,3], cex=0.01, col=2, lty=2)

legend("top", c("orm (loglog)","orm (logit)", "logistic regression"), col=c(1,4,2), lty=c(1,1),bty="n")



est.quantile.orm.logit <- quantile.orm(mod=orm.logit,new.data=new.data, probs=0.95,se=TRUE )


est.quantile.orm.cloglog <- quantile.orm(mod=orm.cloglog,new.data=new.data, probs=0.95,se=TRUE )

plot(new.age, est.quantile.orm.logit$quantile,type="n",
     cex=0.1, col=4, ylim=c(300, 6000), axes=FALSE, log="y", xlab="age", ylab="viral load")

axis(2, at=c(400, 800, 1500, 3000, 6000), labels=c("DL", "800","1500", "3000", "6000"))
axis(1, at=seq(20, 80, by=10), labels=c("20", "30","40","50", "60", "70","80"),
     las=1, col=1, cex.axis=0.8)
polygon(c(new.age, rev(new.age)), c(est.quantile.orm.logit$ub, rev(est.quantile.orm.logit$lb[,1])), col = colors()[618], border = NA)
polygon(c(new.age, rev(new.age)), c(est.quantile.orm.cloglog$ub, rev(est.quantile.orm.cloglog$lb[,1])), col ='grey' , border = NA)

lines(new.age, est.quantile.orm.logit$quantile, cex=0.01, col=4)
lines(new.age, est.quantile.orm.logit$ub, cex=0.01, col=4, lty=2)
lines(new.age, est.quantile.orm.logit$lb, cex=0.01, col=4, lty=2)
lines(new.age, est.quantile.orm.cloglog$quantile, cex=0.01)
lines(new.age, est.quantile.orm.cloglog$lb cex=0.01, lty=2)
lines(new.age, est.quantile.orm.cloglog$ub, cex=0.01, lty=2)

abline(h=400, lty=2, col=1)
legend("top", c("orm (loglog)","orm (logit)"), col=c(1,4), lty=c(1,1), bty="n")
box()


library(quantreg)
# 
rq.mod <- rq(post.rna.d~ site + male + rcs(age,4) + init.class + route + 
               rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
               preARTAIDS.clinical, data=d, tau=0.95)


rq.mod.2 <- rq(tmp~ site + male + rcs(age,4) + init.class + route + 
                 rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
                 preARTAIDS.clinical, data=d, tau=0.95)
X <-  cbind(1, new.data)
pred <-drop(X%*%rq.mod$coef)
V <- summary(rq.mod, cov = TRUE)
df <- V$rdf
tfrac <- qt((1 - 0.95)/2, df)
sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))

est.quantile.rq <- cbind(pred, pred + tfrac * sdpred %o% 
                           c(1, -1))
colnames(est.quantile.rq) <- c("fit", "lower", "higher")

X <-  cbind(1, new.data)
pred <-drop(X%*%rq.mod.2$coef)
V <- summary(rq.mod.2, cov = TRUE)
df <- V$rdf
tfrac <- qt((1 - 0.95)/2, df)
sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))

est.quantile.rq.2 <- cbind(pred, pred + tfrac * sdpred %o% 
                             c(1, -1))
colnames(est.quantile.rq.2) <- c("fit", "lower", "higher")

new.data2 <- new.data[1:4,]

new.data2[,7:9] <- matrix(rcspline.eval(35,  knots=orm.logit$Design$parms$age, inclx=TRUE),nrow=4, ncol=3, byrow=TRUE)

new.data2[1, 10:12] <- orm.logit$x[which(d$init.class=="Boosted PI")[1], 10:12]

new.data2[2, 10:12] <- orm.logit$x[which(d$init.class=="NNRTI")[1], 10:12]
new.data2[3, 10:12] <- orm.logit$x[which(d$init.class=="Unboosted PI")[1], 10:12]
new.data2[4, 10:12] <- orm.logit$x[which(d$init.class=="Other")[1], 10:12]


#est.mean.orm.logit <- mean.orm(orm.logit, new.data2,se=TRUE )

new.data2.ols <- data.frame(site="brazil", male=1, age=35, 
                            init.class=c("Boosted PI","NNRTI","Unboosted PI","Other"), route="homo/bisexual", 
                            nadir.cd4=167, init.year=2010, pre.rna=91728,
                            preARTAIDS.clinical="not AIDS")

est.cdf.orm.logit <- cdf.orm(mod=orm.logit, new.data=new.data2, at.y=0)
est.cdf.orm.cloglog <- cdf.orm(mod=orm.cloglog, new.data=new.data2, at.y=0)


test <- predict(binary.mod, new.data2.ols, se.fit=TRUE)
est.cdf.binary <- cbind(plogis(test$linear.predictors), plogis(with(test, linear.predictors + 1.96*cbind(-se.fit,se.fit))))


plot(1:4,1-est.cdf.orm.cloglog$est, type="n",
     axes=FALSE, ylab="Probability", xlab="", xlim=c(0.5, 4.5), ylim=c(0.02, 0.22), main="")
#abline(v=1)
axis(2, at=c(0, 0.05, 0.1, 0.15, 0.2), labels=c("0", "0.05", "0.10", "0.15", "0.20"))
segments(1:4,1-est.cdf.orm.cloglog$lb, 1:4, 1-est.cdf.orm.cloglog$ub)
points(1:4,1-est.cdf.orm.cloglog$est, pch="-", cex=1)


segments((1:4)+0.15,1-est.cdf.orm.logit$lb, (1:4)+0.15, 1-est.cdf.orm.logit$ub, col=4)
points((1:4)+0.15,1-est.cdf.orm.logit$est, pch="-", cex=1, col=4)

segments((1:4)+0.3, est.cdf.binary[,2], (1:4)+0.3, est.cdf.binary[,3], col=2)
points((1:4)+0.3, est.cdf.binary[,1], pch="-", cex=1, col=2)

legend("top", c("orm (loglog)","orm (logit)", "logistic regression"), col=c(1,4,2), lty=c(1,1),bty="n")
axis(1, 1:4, labels=c("BPI","NNRTI","UBPI","Other"),
     las=1, tck=-0.01, col=1, cex.axis=0.8)
box()

est.cdf.orm.logit.g1000<- cdf.orm(mod=orm.logit, new.data=new.data2, at.y=1000)
est.cdf.orm.cloglog.g1000 <- cdf.orm(mod=orm.cloglog, new.data=new.data2, at.y=1000)

post.rna.g1000 <- ifelse(d$post.rna.d>1000, 1, 0)

binary.mod.g1000 <- lrm(post.rna.g1000~ site + male + rcs(age,4) + init.class + route + 
                          rcs(sqrt(nadir.cd4),4) +rcs(init.year,4)+ rcs(log(pre.rna),4)+
                          preARTAIDS.clinical, data=d, x=TRUE, y=TRUE)

test <- predict(binary.mod.g1000, new.data2.ols, se.fit=TRUE)

est.cdf.binary.g1000 <- cbind(plogis(test$linear.predictors), plogis(with(test, linear.predictors + 1.96*cbind(-se.fit,se.fit))))


plot(1:4,1-est.cdf.orm.cloglog.g1000$est, type="n",
     axes=FALSE, ylab="Probability", xlab="", xlim=c(0.5, 4.5), ylim=c(0.02, 0.22), main="")
#abline(v=1)
axis(2, at=c(0, 0.05, 0.1, 0.15, 0.2), labels=c("0", "0.05", "0.10", "0.15", "0.20"))
segments(1:4,1-est.cdf.orm.cloglog.g1000$lb, 1:4, 1-est.cdf.orm.cloglog.g1000$ub)
points(1:4,1-est.cdf.orm.cloglog.g1000$est, pch="-", cex=1)


segments((1:4)+0.1,1-est.cdf.orm.logit.g1000$lb, (1:4)+0.1, 1-est.cdf.orm.logit.g1000$ub, col=4)
points((1:4)+0.1,1-est.cdf.orm.logit.g1000$est, pch="-", cex=1, col=4)

segments((1:4)+0.2, est.cdf.binary.g1000[,2], (1:4)+0.2, est.cdf.binary.g1000[,3], col=2)
points((1:4)+0.2, est.cdf.binary.g1000[,1], pch="-", cex=1, col=2)


legend("top", c("orm (loglog)","orm (logit)", "logistic regression"), col=c(1,4,2), lty=c(1,1),bty="n")
axis(1, 1:4, labels=c("BPI","NNRTI","UBPI","Other"),
     las=1, tck=-0.01, col=1, cex.axis=0.8)
box()

est.quantile.orm.logit <- quantile.orm(mod=orm.logit,new.data=new.data2, probs=0.95,se=TRUE )

est.quantile.orm.cloglog <- quantile.orm(mod=orm.cloglog,new.data=new.data2, probs=0.95,se=TRUE )
plot(1:4,est.quantile.orm.cloglog$est, type="n",
     axes=FALSE, ylab="viral load", xlab="", xlim=c(0.5, 4.5), ylim=c(300, 10000), main="", log="y")

axis(2, at=c(400, 800,2500, 5000, 10000), labels=c("DL", "800","2500", "5000", "10000"))
segments(1:4,est.quantile.orm.cloglog$lb, 1:4, est.quantile.orm.cloglog$ub)
points(1:4,est.quantile.orm.cloglog$quantile, pch="-", cex=1)


segments((1:4)+0.2, est.quantile.orm.logit$lb, (1:4)+0.2, est.quantile.orm.logit$ub, col=4)
points((1:4)+0.2, est.quantile.orm.logit$quantile, pch="-", cex=1, col=4)

legend("top", c("orm (loglog)","orm (logit)"), col=c(1,4), lty=c(1,1),bty="n")
axis(1, 1:4, labels=c("BPI","NNRTI","UBPI","Other"),
     las=1, tck=-0.01, col=1, cex.axis=0.8)
abline(h=400, col=1, lty=2)
box()

mtext("P[detectable VL |Z]", side=3, at=1/6, outer=TRUE, cex=1)
mtext("P[VL>1000 |Z]", side=3, at=3/6, outer=TRUE, cex=1)
mtext("95th percentile", side=3, at=5/6, outer=TRUE, cex=1)
dev.off()



