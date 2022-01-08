rm(list=ls())
libs <- c("dplyr", "stringr", "readr", "tidyr", "purrr", "ggplot2", "evd", "magrittr")
invisible(lapply(libs, library, character.only = TRUE))

# evd package is used for for Gumbel dist.

#paper directory (sims for other functions)
pdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper")
psdir <- file.path(pdir,"sims","out")

figdir <- file.path(pdir,"fig")

# for reference
if (0){
  
# generate data from standard gumbel, aka standard gumbel or gumbel maximum
# cdf is G(x) = exp{-exp[-(z-a)/b]}
generate.data.3 <- function(seed=1, n=50, p=0.5, alpha=0, beta=c(1.0, -0.5), scale=1){
    set.seed(seed)
    z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
    z2 <- rnorm(n, 0, 1)
    y0 <- rgumbel(n, alpha+beta[1]*z1 + beta[2]*z2, scale)
    # y <- exp(y0) #no transform
    data <- data.frame(y=ordered(y0), z1=z1, z2=z2)
    return(data)
  }
  
cutpoint_est <- function(post,fitdata,y){
    truey<-c(-Inf,fitdata$truey0,Inf)
    
    if (y<=min(fitdata$truey0)){
      # return(-Inf)
      return(c(mean=-Inf,se_mean=NA,sd=NA,`2.5%`=NA,`25%`=NA,`50%`=NA,
               `75%`=NA,`97.5%`=NA, n_eff=NA,  Rhat=NA))
    }
    else if (y > max(fitdata$truey0)) {
      #  return(Inf)
      return(c(mean=Inf,se_mean=NA,sd=NA,`2.5%`=NA,`25%`=NA,`50%`=NA,
               `75%`=NA,`97.5%`=NA, n_eff=NA,  Rhat=NA))
    }
    else {
      post[which(fitdata$truey0>=y)[1]-1,]
    }
  }
  
# get CDF at given y value
# fit_cdf is conditional cdf
# fitdata is data used for BayesCPM fit
cdf_val_est <- function(fit_cdf,fitdata,y){
    truey<-c(-Inf,fitdata$truey0,Inf)
    
    if (y<=min(fitdata$truey0)){
      fit_cdf %>% filter(Var2==-Inf) %>% mutate(yin=y)
    }
    else if (y > max(fitdata$truey0)) {
      fit_cdf %>% filter(Var2==Inf) %>% mutate(yin=y)
    }
    else {
      idx0<-which(fitdata$truey0>=y)[1]
      V2<-fit_cdf %>% pull(Var2)
      fit_cdf %>% filter(Var2==V2[idx0]) %>% mutate(yin=y)
    }
  }
  
# get contrast between two conditional means
getMeanContrast<-function(fit, fitdata, cont.df, sm=TRUE){
    fit_raw <- getMean(fit, fitdata, newdata=cont.df, summ=FALSE)
    
    a<-fit_raw %>% filter(ndrow==1) %>% select(mn)
    b<-fit_raw %>% filter(ndrow==2) %>% select(mn)
    cont.vals<-b-a
    
    if (sm){
      mn_cont_summ <- cont.vals %>%
        dplyr::summarize(mn_cont_mean=mean(mn),
                         mn_cont_med=median(mn),
                         mn_cont_sd=sd(mn),
                         mn_cont_q2.5=quantile(mn,probs=0.025),
                         mn_cont_q5=quantile(mn,probs=0.05),
                         mn_cont_q10=quantile(mn,probs=0.10),
                         mn_cont_q25=quantile(mn,probs=0.25),
                         mn_cont_q75=quantile(mn,probs=0.75),
                         mn_cont_q90=quantile(mn,probs=0.90),
                         mn_cont_q95=quantile(mn,probs=0.95),
                         mn_cont_q97.5=quantile(mn,probs=0.975))
      return(mn_cont_summ)
    } else {
      return(cont.vals)
    }
    
  }
  
# get contrast between two conditional quantiles
getQuantileContrast<-function(fit, fitdata, cont.df, q, sm=TRUE){
    fit_raw <- getQuantile(fit, fitdata, newdata=cont.df, q=q, summ=FALSE)
    
    a<-fit_raw %>% filter(ndrow==1) %>% select(qtile)
    b<-fit_raw %>% filter(ndrow==2) %>% select(qtile)
    cont.vals<-b-a
    
    if (sm){
      qtile_cont_summ <- cont.vals %>%
        dplyr::summarize(qtile_cont_mean=mean(qtile),
                         qtile_cont_med=median(qtile),
                         qtile_cont_sd=sd(qtile),
                         qtile_cont_q2.5=quantile(qtile,probs=0.025),
                         qtile_cont_q5=quantile(qtile,probs=0.05),
                         qtile_cont_q10=quantile(qtile,probs=0.10),
                         qtile_cont_q25=quantile(qtile,probs=0.25),
                         qtile_cont_q75=quantile(qtile,probs=0.75),
                         qtile_cont_q90=quantile(qtile,probs=0.90),
                         qtile_cont_q95=quantile(qtile,probs=0.95),
                         qtile_cont_q97.5=quantile(qtile,probs=0.975))
      return(qtile_cont_summ)
    } else {
      return(cont.vals)
    }
    
  }
  

sim3_coeff.fun <- function(sim=5, seed=1, n=50, p=0.5, alpha=0, beta=c(1,-0.5),
                             scale = 1, yvals=c(-0.3, 0, 0.5, 1.5, 2.5),
                             zvals = data.frame(z1=c(1,1),z2=c(1,0))){
    set.seed(seed)
    seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]
    
    ## preallocate matrices & lists
    
    # betas
    beta.est <- beta.est.cn <- vector("list",sim)
    
    # cutpoints
    gamma.y <- gamma.y.cn <- vector("list",sim)
    
    # conditional cdf at y values
    cond.cdf <- cond.cdf.cn <- vector("list",sim)
    
    # conditional mean & contrast
    cond.mn <- cond.mn.cn <- vector("list",sim)
    cond.mn.cnt <- cond.mn.cnt.cn <- vector("list",sim)
    
    # conditional median & contrast
    cond.med <- cond.med.cn <- vector("list",sim)
    cond.med.cnt <- cond.med.cnt.cn <- vector("list",sim)
    
    # conditional q20 & contrast
    cond.q20 <- cond.q20.cn <- vector("list",sim)
    cond.q20.cnt <- cond.q20.cnt.cn <- vector("list",sim)
    
    # conditional q80 & contrast
    cond.q80 <- cond.q80.cn <- vector("list",sim)
    cond.q80.cnt <- cond.q80.cnt.cn <- vector("list",sim)
    
    for(i in 1:sim){
      print(i)
      try({
        data <- generate.data.3(seed=seeds[i], n=n, p=p,
                                alpha=alpha, beta=beta, scale=scale)
        
        ycens <- as.numeric(levels(data$y)[as.numeric(data$y)])
        ycens[ycens < 0] <- 0
        data$y.cens<-ordered(ycens)
        
        ## fit full outcome model
        mod_data <- mkStanDat(data, outcome="y",
                              preds = c("z1", "z2"),
                              link=3, #loglog link
                              conc=function(n) 1/n)
        
        cpm_fit <- sampling(ord_mod1, data=mod_data, seed=seeds[i],
                            iter=4000, warmup=2000, chains=2, refresh=2000,
                            control = list(adapt_delta = 0.9))
        
        # beta
        beta.est[[i]] <- summary(cpm_fit, pars=c("b[1]","b[2]"), use_cache=FALSE)$summary %>%
          as_tibble() %>% mutate(par=c("b[1]","b[2]"))
        
        # gamma cutpoint estimates at given y values
        cpout <- summary(cpm_fit, pars="cutpoints", use_cache=FALSE)$summary
        gamma.y[[i]] <- lapply(yvals, function(x) cutpoint_est(cpout,mod_data,x)) %>% do.call(rbind,.) %>%
          as_tibble() %>% mutate(par=c("gamma[y1]","gamma[y2]","gamma[y3]","gamma[y4]","gamma[y5]"))
        
        # conditional cdf at y values
        fit_cdf <- try(getCDF(cpm_fit, mod_data, newdata=zvals))
        cond.cdf[[i]] <- try(lapply(yvals,function(x) cdf_val_est(fit_cdf, mod_data,x)) %>% do.call(rbind,.) %>% as_tibble())
        
        # conditional mean (& contrast)
        cond.mn[[i]] <- try(getMean(cpm_fit, mod_data, newdata=zvals))
        # cond.mn.cnt[[i]] <- try(getMeanContrast(cpm_fit, mod_data, cont.df=zvals))
        
        # conditional median (& contrast)
        cond.med[[i]] <- try(getQuantile(cpm_fit, mod_data, newdata=zvals, q=0.5))
        #  cond.med.cnt[[i]] <- try(getQuantileContrast(cpm_fit, mod_data, cont.df=zvals, q=0.5))
        
        # conditional q20 (& contrast)
        cond.q20[[i]] <- try(getQuantile(cpm_fit, mod_data, newdata=zvals, q=0.2))
        # cond.q20.cnt[[i]] <- try(getQuantileContrast(cpm_fit, mod_data, cont.df=zvals, q=0.2))
        
        # conditional q80 (& contrast)
        cond.q80[[i]] <- try(getQuantile(cpm_fit, mod_data, newdata=zvals, q=0.8))
        # cond.q80.cnt[[i]] <- try(getQuantileContrast(cpm_fit, mod_data, cont.df=zvals, q=0.8))
        
        ## fit censored outcome model
        mod_data_cn <- mkStanDat(data, outcome="y.cens",
                                 preds = c("z1", "z2"),
                                 link=3, #loglog link
                                 conc=function(n) 1/n)
        
        cpm_fit_cn <- sampling(ord_mod1, data=mod_data_cn, seed=seeds[i],
                               iter=4000, warmup=2000, chains=2, refresh=2000,
                               control = list(adapt_delta = 0.9))
        
        # beta
        beta.est.cn[[i]] <- summary(cpm_fit_cn, pars=c("b[1]","b[2]"),use_cache=FALSE)$summary %>%
          as_tibble() %>% mutate(par=c("b[1]","b[2]"))
        
        # gamma cutpoint estimates at given y values
        cpout.cn <- summary(cpm_fit_cn, pars="cutpoints",use_cache=FALSE)$summary
        gamma.y.cn[[i]] <- lapply(yvals, function(x) cutpoint_est(cpout.cn,mod_data_cn,x)) %>% do.call(rbind,.) %>%
          as_tibble() %>% mutate(par=c("gamma[y1]","gamma[y2]","gamma[y3]","gamma[y4]","gamma[y5]"))
        
        # conditional cdf at yvals values
        fit_cdf_cn <- try(getCDF(cpm_fit_cn, mod_data_cn, newdata=zvals))
        cond.cdf.cn[[i]] <- try(lapply(yvals,function(x) cdf_val_est(fit_cdf_cn, mod_data_cn, x)) %>% do.call(rbind,.) %>% as_tibble())
        
        
        # conditional mean (& contrast)
        cond.mn.cn[[i]] <- try(getMean(cpm_fit_cn, mod_data_cn, newdata=zvals))
        #  cond.mn.cnt.cn[[i]] <- try(getMeanContrast(cpm_fit_cn, mod_data_cn, cont.df=zvals))
        
        # conditional median (& contrast)
        cond.med.cn[[i]] <- try(getQuantile(cpm_fit_cn, mod_data_cn, newdata=zvals, q=0.5))
        # cond.med.cnt.cn[[i]] <- try(getQuantileContrast(cpm_fit_cn, mod_data_cn, cont.df=zvals, q=0.5))
        
        # conditional q20 (& contrast)
        cond.q20.cn[[i]] <- try(getQuantile(cpm_fit_cn, mod_data_cn, newdata=zvals, q=0.2))
        #cond.q20.cnt.cn[[i]] <- try(getQuantileContrast(cpm_fit_cn, mod_data_cn, cont.df=zvals, q=0.2))
        
        # conditional q80 (& contrast)
        cond.q80.cn[[i]] <- try(getQuantile(cpm_fit_cn, mod_data_cn, newdata=zvals, q=0.8))
        #cond.q80.cnt.cn[[i]] <- try(getQuantileContrast(cpm_fit_cn, mod_data_cn, cont.df=zvals, q=0.8))
      })
    }
    
    
    return(list(full=list(beta.est=beta.est,
                          gamma.y=gamma.y,
                          cond.cdf=cond.cdf,
                          cond.mn=cond.mn,
                          cond.mn.cnt=cond.mn.cnt,
                          cond.med=cond.med,
                          cond.med.cnt=cond.med.cnt,
                          cond.q20=cond.q20,
                          cond.q20.cnt=cond.q20.cnt,
                          cond.q80=cond.q80,
                          cond.q80.cnt=cond.q80.cnt),
                cens=list(beta.est.cn=beta.est.cn,
                          gamma.y.cn=gamma.y.cn,
                          cond.cdf.cn=cond.cdf.cn,
                          cond.mn.cn=cond.mn.cn,
                          cond.mn.cnt.cn=cond.mn.cnt.cn,
                          cond.med.cn=cond.med.cn,
                          cond.med.cnt.cn=cond.med.cnt.cn,
                          cond.q20.cn=cond.q20.cn,
                          cond.q20.cnt.cn=cond.q20.cnt.cn,
                          cond.q80.cn=cond.q80.cn,
                          cond.q80.cnt.cn=cond.q80.cnt.cn)) )
  }
  
}


# Load sim array
simarray <- readRDS(file.path(pdir,"sims","bayes_cpm_simarray2.rds"))

# Load all sim d data

# sim_d0_n - loglog link, no trans, conc=1/ncats
# sim_d1_n - loglog link, log trans, conc=1/(0.8 + 0.35*max(ncats, 3))
# sim_d2_n - loglog link, log trans, conc=1/(2+(ncats/3))
# sim_d3_n - loglog link, log trans, conc=1/2

for (j in 1:50){
  for (f in c("sim_d0_n","sim_d1_n","sim_d2_n","sim_d3_n")){
    nam <- paste0(f, simarray[j,"nsamps"], "_",simarray[j,"rep"])
    fp <- file.path(psdir,paste0(nam,".rds"))
    try(assign(nam, readRDS(fp)))
  }
}

rm(nam,fp,j,f)

# combine 10 reps of 100 sims for each setting
for (n in c(25,50,100,200,400)){
  for (pre in c("sim_d0_n","sim_d1_n","sim_d2_n","sim_d3_n")){
    nm0 <- paste0(pre, n)
    nm <- paste0(nm0, "_", 1:10)
    nmlist <- map(nm,get)
    fl <- map(nmlist, function(x) x[[1]]) %>% pmap(rbind) %>% map(bind_rows)
    cn <- map(nmlist, function(x) x[[2]]) %>% pmap(rbind) %>% map(bind_rows)
    assign(nm0,list(full=fl,cens=cn))
    rm(list=nm)
  }
}

rm(cn,fl,nmlist,nm,nm0,n,pre)

## function to summarize betas

sim_beta_summ <- function(result){
  
  if (names(result[1])=='beta.est'){
    sims_sub <- result[['beta.est']] %>% select(mean, `50%`, par)
  } else {
    sims_sub <- result[['beta.est.cn']] %>% select(mean, `50%`, par)
  }

truebeta <- data.frame(par=c("b[1]","b[2]"),trueval=c(1,-0.5))

mrg_df <- merge(sims_sub,truebeta,by="par",all.y=FALSE) %>% 
  mutate(bias=`50%`-trueval, pct.bias = 100*(bias/abs(trueval))) %>% 
  group_by(par) %>% 
  summarize(avg.pct.bias=mean(pct.bias, na.rm=TRUE),
            avg.bias=mean(bias, na.rm=TRUE),
            mse=mean((bias^2), na.rm=TRUE))

return(mrg_df)
}

# sim_beta_summ(sim_d0_n400[['full']]); sim_beta_summ(sim_d0_n25[['cens']])


## function to summarize gammas

#! pct bias doesn't make sense when denom is 0, show bias or remove

sim_gamma_summ <- function(result){
  
  if (names(result[2])=='gamma.y'){
    sims_sub <- result[['gamma.y']] %>% select(mean, `50%`, par)
  } else {
    sims_sub <- result[['gamma.y.cn']] %>% select(mean, `50%`, par)
  }
  
  truegamma <- data.frame(par=c("gamma[y1]","gamma[y2]","gamma[y3]","gamma[y4]","gamma[y5]"),
                          trueval=c(-0.3, 0, 0.5, 1.5, 2.5))  
  
  mrg_df <- merge(sims_sub,truegamma,by="par",all.y=FALSE) %>% 
    mutate(bias=`50%`-trueval, pct.bias = 100*(bias/abs(trueval))) %>% 
    group_by(par) %>% 
    summarize(avg.pct.bias=mean(pct.bias, na.rm=TRUE),
              avg.bias=mean(bias, na.rm=TRUE),
              mse=mean((bias^2), na.rm=TRUE))
  
return(mrg_df)
}

#sim_gamma_summ(sim_d0_n400[['full']]); sim_gamma_summ(sim_d0_n400[['cens']])
#sim_gamma_summ(sim_d0_n25[['full']]); sim_gamma_summ(sim_d0_n25[['cens']])

## function to summarize conditional CDF

# true cdf values
if(0){
  yvals <- c(-0.3, 0, 0.5, 1.5, 2.5)
  
  # at z1=1, z2=1
  pgumbel(yvals,1*1-0.5*1,1)
  curve(pgumbel(x,1*1-0.5*1,1),-2,5, xlab='y')
  segments(x0=yvals,
           y0=rep(0,5),
           y1=pgumbel(yvals,1*1-0.5*1,1),
           lty=2)
  segments(x0=rep(-3,5),
           x1=yvals,
           y0=pgumbel(yvals,1*1-0.5*1,1),
           lty=2)
  
  # at z1=1, z2=0
  pgumbel(yvals,1*1-0.5*0,1)
  curve(pgumbel(x,1*1-0.5*0,1),-2,5, xlab='y')
  segments(x0=yvals,
           y0=rep(0,5),
           y1=pgumbel(yvals,1*1-0.5*0,1),
           lty=2)
  segments(x0=rep(-3,5),
           x1=yvals,
           y0=pgumbel(yvals,1*1-0.5*0,1),
           lty=2)
  
  qgumbel(c(0.1, 0.4, 0.5, 0.6, 0.9),1*1-0.5*1,1)  
  curve(dgumbel(x,1*1-0.5*1,1),-2,5)

}

sim_cdf_summ <- function(result, beta=c(1, -0.5), scale=1,
                          gamma.y=c(-0.3, 0, 0.5, 1.5, 2.5),
                          zvals = c(1,1)){
  
  if (names(result[3])=='cond.cdf'){
    sims_sub <- result[['cond.cdf']] %>% filter(z1==zvals[1] & z2==zvals[2])
  } else {
    sims_sub <- result[['cond.cdf.cn']] %>% filter(z1==zvals[1] & z2==zvals[2])
    gamma.y<-gamma.y[3:5]
  }
  
  truevals <- pgumbel(gamma.y,beta[1]*zvals[1]+ beta[2]*zvals[2],1)
  true_df <- data.frame(yin=gamma.y, true_cdf=truevals)
  
  mrg_df <- merge(sims_sub, true_df, by='yin',all.y=FALSE) %>% 
    mutate(bias=med_cdf-true_cdf, pct.bias = 100*(bias/abs(true_cdf))) %>% 
    group_by(yin) %>% 
    summarize(avg.pct.bias=mean(pct.bias, na.rm=TRUE),
              avg.bias=mean(bias, na.rm=TRUE),
              mse=mean((bias^2), na.rm=TRUE))
  
  return(mrg_df)
}

# sim_cdf_summ(sim_d0_n400[['full']]); sim_cdf_summ(sim_d0_n400[['cens']])
# sim_cdf_summ(sim_d0_n25[['full']]); sim_cdf_summ(sim_d0_n25[['cens']])

## function to summarize conditional mean

# true mean values
if(0){
#https://stackoverflow.com/questions/18106797/how-to-get-euler-mascheronis-constant-in-r
euler_masch <- -digamma(1)
  
#beta1*z1+beta2*z2
#1*1-0.5*1
mean(rgumbel(1e7, 1*1-0.5*1, 1)) # sim
1*1-0.5*1 + 1*euler_masch # analytical loc + scale*Euler–Mascheroni constant

#1*1-0.5*0
mean(rgumbel(1e7, 1*1-0.5*0, 1)) # sim
1*1-0.5*0 + 1*euler_masch # analytical

}

euler_masch <- -digamma(1)
sim_mn_summ <- function(result, beta=c(1, -0.5),
                         zdf = data.frame(z1=c(1,1),z2=c(1,0))){
  
  true_df <- data.frame(zdf,true_mn=c(beta[1]*zdf[1,1]+beta[2]*zdf[1,2] + 1*euler_masch,
                             beta[1]*zdf[2,1]+beta[2]*zdf[2,2] + 1*euler_masch))
  
  if (names(result[4])=='cond.mn'){
    sims_sub <- result[['cond.mn']]
  } else {
    sims_sub <- result[['cond.mn.cn']]
  }

  mrg_df <- merge(sims_sub, true_df, by=c('z1','z2'), all.y=FALSE) %>%
    mutate(bias=med_mn-true_mn, pct.bias = 100*(bias/abs(true_mn))) %>%
     group_by(ndrow,z1,z2) %>%
     summarize(avg.pct.bias=mean(pct.bias, na.rm=TRUE),
               avg.bias=mean(bias, na.rm=TRUE),
               mse=mean((bias^2), na.rm=TRUE))

  return(mrg_df)
}

#sim_mn_summ(sim_d0_n25[['full']]); sim_mn_summ(sim_d0_n25[['cens']])


# function to summarize conditional quantile

# true quantiles
if(0){
# true median
qgumbel(0.5, 1*1-0.5*1); 1*1-0.5*1 - 1*log(log(2))
qgumbel(0.5, 1*1-0.5*0); 1*1-0.5*0 - 1*log(log(2))

# true 20th qtile
qgumbel(0.2, 1*1-0.5*1)
qgumbel(0.2, 1*1-0.5*0)

# true 80th qtile
qgumbel(0.8, 1*1-0.5*1)
qgumbel(0.8, 1*1-0.5*0)
}

sim_qtile_summ <- function(result, beta=c(1, -0.5),
                        zdf = data.frame(z1=c(1,1),z2=c(1,0)),
                        q = 0.5, statnm='cond.med'){
  
  true_df <- data.frame(zdf, true_qtile = 
                          c(qgumbel(q, beta[1]*zdf[1,1]+beta[2]*zdf[1,2]),
                            qgumbel(q,beta[1]*zdf[2,1]+beta[2]*zdf[2,2])))
  
  if(q==0.5 & statnm != 'cond.med'|
     q==0.2 & statnm != 'cond.q20'|
     q==0.8 & statnm != 'cond.q80') {stop("q must match statnm!")}    
                    
  if (names(result[6])==statnm|names(result[8])==statnm|names(result[10])==statnm){
    sims_sub <- result[[statnm]] 
  } else {
    sims_sub <- result[[paste0(statnm,".cn")]]
  }
  
  mrg_df <- merge(sims_sub, true_df, by=c('z1','z2'), all.y=FALSE) %>% 
    mutate(bias=med_qtile-true_qtile, pct.bias = 100*(bias/abs(true_qtile))) %>% 
    group_by(ndrow,z1,z2) %>% 
    summarize(avg.pct.bias=mean(pct.bias, na.rm=TRUE),
              avg.bias=mean(bias, na.rm=TRUE),
              mse=mean((bias^2), na.rm=TRUE))
  
  return(mrg_df)
}

# sim_qtile_summ(sim_d0_n25[['full']], q = 0.5, statnm='cond.med')
# sim_qtile_summ(sim_d0_n25[['full']], q = 0.2, statnm='cond.q20')
# sim_qtile_summ(sim_d0_n25[['full']], q = 0.8, statnm='cond.q80')

# sim_qtile_summ(sim_d0_n25[['cens']], q = 0.5, statnm='cond.med')
# sim_qtile_summ(sim_d0_n25[['cens']], q = 0.2, statnm='cond.q20')
# sim_qtile_summ(sim_d0_n25[['cens']], q = 0.8, statnm='cond.q80')


# get summaries and convert data for all simulations
nsamps <- c(25,50,100,200,400)
concs <- c('1/J','1/(0.8 + 0.35*J)','1/(2+(J/3))')
pres <- c('sim_d0_n','sim_d1_n','sim_d2_n')

# for alpha=1/2 conc
#concs <- c('1/J','1/(0.8 + 0.35*J)','1/(2+(J/3))','1/2')
#pres <- c('sim_d0_n','sim_d1_n','sim_d2_n','sim_d3_n')
outcome <- c('full.','cens.')

for (k in seq_along(outcome)){
  for (j in seq_along(pres)){
    for (i in nsamps){
      nm0 <- paste0(pres[j], i)
      nm_beta <- paste0(outcome[k],'beta.',nm0)
      nm_gamma <- paste0(outcome[k],'gamma.',nm0)
      nm_cdf <- paste0(outcome[k],'cdf.',nm0)
      nm_mn <- paste0(outcome[k],'mn.',nm0)
      nm_med <- paste0(outcome[k],'med.',nm0)
      nm_q20 <- paste0(outcome[k],'q20.',nm0)
      nm_q80 <- paste0(outcome[k],'q80.',nm0)
      
      beta <- sim_beta_summ(get(nm0)[[k]]) %>% mutate(n=i, conc=concs[j])
      gamma <- sim_gamma_summ(get(nm0)[[k]]) %>% mutate(n=i, conc=concs[j])
      cdf <- sim_cdf_summ(get(nm0)[[k]]) %>% mutate(n=i, conc=concs[j])
      mn <- sim_mn_summ(get(nm0)[[k]]) %>% mutate(n=i, conc=concs[j])
      med <- sim_qtile_summ(get(nm0)[[k]], q=0.5, statnm='cond.med') %>% mutate(n=i, conc=concs[j]) 
      q20 <- sim_qtile_summ(get(nm0)[[k]], q=0.2, statnm='cond.q20') %>% mutate(n=i, conc=concs[j]) 
      q80 <- sim_qtile_summ(get(nm0)[[k]], q=0.8, statnm='cond.q80') %>% mutate(n=i, conc=concs[j]) 
      
      assign(nm_beta,beta)
      assign(nm_gamma,gamma)
      assign(nm_cdf,cdf)
      assign(nm_mn,mn)
      assign(nm_med,med)
      assign(nm_q20,q20)
      assign(nm_q80,q80)
      
      rm(beta, gamma, cdf, mn, med, q20, q80, nm0,
         nm_beta, nm_gamma, nm_cdf, nm_mn, nm_med, nm_q20, nm_q80)
    }
  }
}



# complete datasets for each param/statistic of interest
for (mm in c("full.","cens.")){
  for (nn in c("beta.","gamma.","cdf.","mn.","med.","q20.","q80.")){
    nam <- paste0(mm,nn,"sim")
    dnam <- gsub("\\.","_",nam)
    rmnm <- gsub("\\.","\\\\.",nam)
    out<-lapply(ls()[grep(nam,ls())],get) %>% bind_rows() %>%
      mutate(n=factor(n,levels=c(400,200,100,50,25)))
    assign(paste0(dnam,"_dat"),out)
    rm(list=ls()[grep(rmnm,ls())])
  }
}

# labels and plot params
# https://stackoverflow.com/questions/62662144/conditional-probability-is-not-displaying-properly-in-facet-label-in-ggplot2
nlabs <- paste0('F(y==',c(-0.3, 0, 0.5, 1.5, 2.5),"*'|'*",'~X[1]==1,X[2]==1)')
pltw<-10; plth<-5; atxtsz<-9; fctsiz<-13

### Betas & Gammas ###

#combined plot
full_par<-bind_rows(full_beta_sim_dat,full_gamma_sim_dat) %>% 
  mutate(outcome="uncensored")
cens_par<-bind_rows(cens_beta_sim_dat,cens_gamma_sim_dat) %>% 
  mutate(outcome="censored")

pltw<-10; plth<-7; atxtsz<-9; fctsiz<-13

rbind(full_par,cens_par) %>% 
  mutate(outcome=factor(outcome,levels=c("uncensored","censored"),
                        labels=c("uncensored Y","censored Y")),
         par=case_when(
           par == 'b[1]' ~ "beta[1]",
           par == 'b[2]' ~ "beta[2]",
           par == 'gamma[y1]' ~ "gamma[y[1]]",
           par == 'gamma[y2]' ~ "gamma[y[2]]",
           par == 'gamma[y3]' ~ "gamma[y[3]]",
           par == 'gamma[y4]' ~ "gamma[y[4]]",
           par == 'gamma[y5]' ~ "gamma[y[5]]"
         )) %>% 
  # filter(!is.na(avg.bias)) %>% 
  # ggplot(aes(x=avg.bias,y=n,col=conc,shape=conc)) +
  filter(!is.na(avg.pct.bias)) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75) +
  facet_grid(outcome ~ par,drop=TRUE,labeller = labeller(par=label_parsed),
             switch="y") +
  xlab("average percent bias of posterior parameters") + ylab("sample size") +
  scale_shape_discrete(name=bquote(alpha)) +
  scale_color_discrete(name=bquote(alpha)) +
  theme(axis.title.x = element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))

ggsave(file.path(figdir,"sim_d_pars.png"),width=pltw,height=plth)
ggsave(file.path(figdir,"sim_d_pars.tiff"),width=pltw,height=plth,dpi=600)

### CDF ###

#combined plot
pltw<-10; plth<-7; atxtsz<-10; fctsiz<-9

full_cdf_sim_dat %<>% mutate(outcome="uncensored")
cens_cdf_sim_dat %<>% mutate(outcome="censored")

rbind(full_cdf_sim_dat,cens_cdf_sim_dat) %>% 
  mutate(yin=factor(yin,labels=nlabs),
         outcome=factor(outcome,levels=c("uncensored","censored"),
                        labels=c("uncensored Y","censored Y"))) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(outcome ~ yin, labeller=labeller(yin=label_parsed),
             switch="y")+
  xlab("average percent bias of posterior conditional CDF") + 
  ylab("sample size") +
  scale_shape_discrete(name=bquote(alpha)) +
  scale_color_discrete(name=bquote(alpha)) +
  theme(axis.title.x = element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))

ggsave(file.path(figdir,"sim_d_cdf.png"),width=pltw,height=plth)
ggsave(file.path(figdir,"sim_d_cdf.tiff"),width=pltw,height=plth,dpi=600)

### Mean ###

#combined plot
pltw<-10; plth<-7; atxtsz<-9; fctsiz<-13
full_mn_sim_dat %<>% mutate(outcome="uncensored")
cens_mn_sim_dat %<>% mutate(outcome="censored")

rbind(full_mn_sim_dat, cens_mn_sim_dat) %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "E(Y*'|'*~X[1]==1,X[2]==1)",
                       "E(Y*'|'*~X[1]==1,X[2]==0)"),
         outcome=factor(outcome,levels=c("uncensored","censored"),
                        labels=c("uncensored Y","censored Y"))) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(outcome ~ ndrow, labeller=labeller(ndrow=label_parsed),
             switch="y")+
  xlab("average percent bias of posterior conditional mean") + 
  ylab("sample size")+
  scale_shape_discrete(name=bquote(alpha)) +
  scale_color_discrete(name=bquote(alpha)) +
  theme(axis.title.x = element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))

ggsave(file.path(figdir,"sim_d_mn.png"),width=pltw,height=plth)
ggsave(file.path(figdir,"sim_d_mn.tiff"),width=pltw,height=plth,dpi=600)

### Median ###

# combined plot
full_med_sim_dat %<>% mutate(outcome="uncensored")
cens_med_sim_dat %<>% mutate(outcome="censored")

rbind(full_med_sim_dat,cens_med_sim_dat) %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==0)"),
         outcome=factor(outcome,levels=c("uncensored","censored"),
                        labels=c("uncensored Y","censored Y"))) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(outcome ~ ndrow, labeller=labeller(ndrow=label_parsed),
             switch="y")+
  xlab("average percent bias of posterior conditional median") + 
  ylab("sample size") +
  scale_shape_discrete(name=bquote(alpha)) +
  scale_color_discrete(name=bquote(alpha)) +
  theme(axis.title.x = element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))

ggsave(file.path(figdir,"sim_d_med.png"),width=pltw,height=plth)
ggsave(file.path(figdir,"sim_d_med.tiff"),width=pltw,height=plth,dpi=600)

### 20% quantile ###

# use avg.bias instead of avg.pct.bias?
#combined plot
full_q20_sim_dat %<>% mutate(outcome="uncensored")
cens_q20_sim_dat %<>% mutate(outcome="censored")

cens_q20_sim_dat_mod <- cens_q20_sim_dat %>% filter(ndrow==2)

rbind(full_q20_sim_dat,cens_q20_sim_dat_mod) %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.2}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.2}*'|'*list(X[1]==1,X[2]==0)"),
         outcome=factor(outcome,levels=c("uncensored","censored"),
                        labels=c("uncensored Y","censored Y"))) %>%
#  ggplot(aes(x=avg.bias,y=n,col=conc,shape=conc)) +
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(outcome ~ ndrow, labeller=labeller(ndrow=label_parsed),
             switch="y")+
  xlab("average percent bias of posterior conditional 20th percentile") + 
  ylab("sample size") +
  scale_shape_discrete(name=bquote(alpha)) +
  scale_color_discrete(name=bquote(alpha)) +
  theme(axis.title.x = element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))+
  geom_text(aes(x=150,y=3,label="Q^{0.2}*'|'*list(X[1]==1,X[2]==1)~censored"),parse=TRUE,
          data=data.frame(outcome=factor("censored Y"),ndrow="Q^{0.2}*'|'*list(X[1]==1,X[2]==1)"),
          inherit.aes=FALSE)

ggsave(file.path(figdir,"sim_d_q20.png"),width=pltw,height=plth)
ggsave(file.path(figdir,"sim_d_q20.tiff"),width=pltw,height=plth,dpi=600)


### 80% quantile ###

# combined plot
full_q80_sim_dat %<>% mutate(outcome="uncensored")
cens_q80_sim_dat %<>% mutate(outcome="censored")

rbind(full_q80_sim_dat,cens_q80_sim_dat) %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.8}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.8}*'|'*list(X[1]==1,X[2]==0)"),
         outcome=factor(outcome,levels=c("uncensored","censored"),
                        labels=c("uncensored Y","censored Y"))) %>%
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(outcome ~ ndrow, labeller=labeller(ndrow=label_parsed),
             switch="y")+
  xlab("average percent bias of posterior conditional 80th percentile") + 
  ylab("sample size") + 
  scale_shape_discrete(name=bquote(alpha)) +
  scale_color_discrete(name=bquote(alpha)) +
  theme(axis.title.x = element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))

ggsave(file.path(figdir,"sim_d_q80.png"),width=pltw,height=plth)

