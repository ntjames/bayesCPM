rm(list=ls())
libs <- c("dplyr", "stringr", "readr", "tidyr", "purrr", "ggplot2", "actuar")
invisible(lapply(libs, library, character.only = TRUE))

# actuar package is used for for log logistic dist.
# https://www.rdocumentation.org/packages/actuar/versions/3.0-0/topics/Loglogistic

#paper directory (sims for other functions)
pdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper")
psdir <- file.path(pdir,"sims","out")

figdir <- file.path(pdir,"fig")

#source(file.path(dir,"cpm_functions.r"))

# load sim array
simarray <- readRDS(file.path(pdir,"sims","bayes_cpm_simarray.rds"))

# load all sim c data

# sim_c0_n - logistic link, conc=1/ncats
# sim_c1_n - logistic link, conc=1/(0.8 + 0.35*max(ncats, 3))
# sim_c2_n - logistic link, conc=1/(2+(ncats/3))
# sim_c3_n - logistic link, conc=1/2


if(FALSE){
# reduced loop below
for (j in 1:25){
  nam_c0 <- paste0("sim_c0_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_c0 <- file.path(psdir,paste0(nam_c0,".rds"))
  nam_c1 <- paste0("sim_c1_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_c1 <- file.path(psdir,paste0(nam_c1,".rds"))
  nam_c2 <- paste0("sim_c2_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_c2 <- file.path(psdir,paste0(nam_c2,".rds"))
  nam_c3 <- paste0("sim_c3_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_c3 <- file.path(psdir,paste0(nam_c3,".rds"))
  try(assign(nam_c0, readRDS(fp_c0)))
  try(assign(nam_c1, readRDS(fp_c1)))
  try(assign(nam_c2, readRDS(fp_c2)))
  try(assign(nam_c3, readRDS(fp_c3)))
}

rm(nam_c0,nam_c1,nam_c2,nam_c3,fp_c0,fp_c1,fp_c2,fp_c3,j)
}

for (j in 1:25){
  for (f in c("sim_c0_n","sim_c1_n","sim_c2_n","sim_c3_n")){
    nam <- paste0(f, simarray[j,"nsamps"], "_",simarray[j,"rep"])
    fp <- file.path(psdir,paste0(nam,".rds"))
    try(assign(nam, readRDS(fp)))
  }
}

rm(nam,fp,j,f)

# for reference
if (0){

generate.data.2 <- function(seed=1, n=50, p=0.5, alpha=0, beta=c(1.0, -0.5), scale=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
  y0 <- rlogis(n, alpha+beta[1]*z1 + beta[2]*z2, scale)
  log.y <- exp(y0) #log logistic
  data <- data.frame(y=ordered(log.y), z1=z1, z2=z2)
  return(data)
}


# get_cutpoints <- function(fit,cps=cp_seq,summ_stat='50%',pr=c(0.5)){
#   summary(fit, pars=paste0("cutpoints[",cps,"]"), probs=pr)$summary[,summ_stat]
# }

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
                           scale = 1/3, yvals=exp(c(-0.5, 0, 0.5, 1, 1.5)),
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
      data <- generate.data.2(seed=seeds[i], n=n, p=p,
                              alpha=alpha, beta=beta, scale=scale)
      
      ycens <- as.numeric(levels(data$y)[as.numeric(data$y)])
      ycens[ycens < 1] <- 1
      data$y.cens<-ordered(ycens)
      
      ## fit full outcome model
      mod_data <- mkStanDat(data, outcome="y",
                            preds = c("z1", "z2"),
                            link=1, #logistic link
                            conc=function(n) 1/n)
      
      cpm_fit <- sampling(ord_mod1, data=mod_data, seed=seeds[i],
                          iter=4000, warmup=2000, chains=2, refresh=2000,
                          control = list(adapt_delta = 0.9))
      
      # beta
      beta.est[[i]] <- summary(cpm_fit, pars=c("b[1]","b[2]"), use_cache=FALSE)$summary %>%
        as_tibble() %>% mutate(par=c("b[1]","b[2]"))
      
      # gamma cutpoint estimates at given y values
      cpout <- summary(cpm_fit, pars="cutpoints", use_cache=FALSE)$summary
      #! gamma.y[[i]] <- lapply(yvals, function(x) cutpoint_est(cpout,mod_data,x)) %>% bind_rows() %>%
      #!   mutate(par=c("gamma[y1]","gamma[y2]","gamma[y3]","gamma[y4]","gamma[y5]"))
      gamma.y[[i]] <- lapply(yvals, function(x) cutpoint_est(cpout,mod_data,x)) %>% do.call(rbind,.) %>%
        as_tibble() %>% mutate(par=c("gamma[y1]","gamma[y2]","gamma[y3]","gamma[y4]","gamma[y5]"))
      
      # conditional cdf at y values
      fit_cdf <- try(getCDF(cpm_fit, mod_data, newdata=zvals))
      #! cond.cdf[[i]] <- try(lapply(yvals,function(x) cdf_val_est(fit_cdf, mod_data,x)) %>% bind_rows())
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
                               link=1, #logistic link
                               conc=function(n) 1/n)
      
      cpm_fit_cn <- sampling(ord_mod1, data=mod_data_cn, seed=seeds[i],
                             iter=4000, warmup=2000, chains=2, refresh=2000,
                             control = list(adapt_delta = 0.9))
      
      # beta
      beta.est.cn[[i]] <- summary(cpm_fit_cn, pars=c("b[1]","b[2]"),use_cache=FALSE)$summary %>%
        as_tibble() %>% mutate(par=c("b[1]","b[2]"))
      
      # gamma cutpoint estimates at given y values
      cpout.cn <- summary(cpm_fit_cn, pars="cutpoints",use_cache=FALSE)$summary
      #! gamma.y.cn[[i]] <- lapply(yvals, function(x) cutpoint_est(cpout.cn,mod_data_cn,x)) %>% bind_rows()
      gamma.y.cn[[i]] <- lapply(yvals, function(x) cutpoint_est(cpout.cn,mod_data_cn,x)) %>% do.call(rbind,.) %>%
        as_tibble() %>% mutate(par=c("gamma[y1]","gamma[y2]","gamma[y3]","gamma[y4]","gamma[y5]"))
      
      # conditional cdf at yvals values
      fit_cdf_cn <- try(getCDF(cpm_fit_cn, mod_data_cn, newdata=zvals))
      #! cond.cdf.cn[[i]] <- try(lapply(yvals,function(x) cdf_val_est(fit_cdf_cn, mod_data_cn, x)) %>% bind_rows())
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



# combine 5 reps of 200 sims for each setting
for (n in c(25,50,100,200,400)){
  for (pre in c("sim_c0_n","sim_c1_n","sim_c2_n","sim_c3_n")){
    nm0 <- paste0(pre, n)
    nm <- paste0(nm0, "_", 1:5)
    nmlist <- map(nm,get)
    
    # fl0 <- map(nmlist, function(x) x[[1]]) 
    # fl1 <- pmap(fl0,rbind)
    # fl2 <- lapply(fl0,bind_rows)
    fl <- map(nmlist, function(x) x[[1]]) %>% pmap(rbind) %>% map(bind_rows)
    
    # cn0 <- map(nmlist, function(x) x[[2]])
    # cn1 <- pmap(cn0,rbind)
    # cn2 <- lapply(cn1,bind_rows)
    cn <- map(nmlist, function(x) x[[2]]) %>% pmap(rbind) %>% map(bind_rows)
    
    assign(nm0,list(full=fl,cens=cn))
    rm(list=nm)
  }
}

#rm(cn0,cn1,cn2,fl0,fl1,fl2,nmlist,nm,nm0,n,pre)
rm(cn,fl,nmlist,nm,nm0,n,pre)

## function to summarize betas

sim_beta_summ <- function(result){
  
  if (names(result[1])=='beta.est'){
    sims_sub <- result[['beta.est']] %>% select(mean, `50%`, par)
  } else {
    sims_sub <- result[['beta.est.cn']] %>% select(mean, `50%`, par)
  }

truebeta <- data.frame(par=c("b[1]","b[2]"),trueval=c(1,-0.5))  # scaled truebeta <- c(1/(1/3),-0.5/(1/3))

mrg_df <- merge(sims_sub,truebeta,by="par",all.y=FALSE) %>% 
  mutate(bias=`50%`-trueval, pct.bias = 100*(bias/abs(trueval))) %>% 
  group_by(par) %>% 
  summarize(avg.pct.bias=mean(pct.bias, na.rm=TRUE),
            avg.bias=mean(bias, na.rm=TRUE),
            mse=mean((bias^2), na.rm=TRUE))

return(mrg_df)
}

# sim_beta_summ(sim_c0_n400[['full']])
# sim_beta_summ(sim_c0_n400[['cens']])


## function to summarize gamma

# foo<-sim_c0_n400[['full']]
# 
# goo<- foo[['gamma.y']] %>% select(mean, `50%`, par) %>% group_by(par) %>% 
#   summarize(mn=mean(`50%`))
# truegammas<-c(-0.5,  0.0,  0.5,  1.0,  1.5)
# estimates are scaled by 1/3 i.e estimate=truegammas*3?
# 
# tb <- data.frame(par=c("b[1]","b[2]"),trueval=c(1,-0.5))
# 
# merge(goo,tb,by="par",all.y=FALSE) %>% 
#   mutate(bias=`50%`-trueval, pct.bias = 100*(bias/abs(trueval))) %>% 
#   group_by(par) %>% 
#   summarize(avg.pct.bias=mean(pct.bias))

# pct bias doesn't really make sense when denom is 0, maybe change to show bias??

sim_gamma_summ <- function(result){
  
  if (names(result[2])=='gamma.y'){
    sims_sub <- result[['gamma.y']] %>% select(mean, `50%`, par)
  } else {
    sims_sub <- result[['gamma.y.cn']] %>% select(mean, `50%`, par)
  }
  
  truegamma <- data.frame(par=c("gamma[y1]","gamma[y2]","gamma[y3]","gamma[y4]","gamma[y5]"),
                          trueval=c(-0.5,  0.0,  0.5,  1.0,  1.5))  
  # scaled trueval = c(-0.5,  0.0,  0.5,  1.0,  1.5)*3
  
  mrg_df <- merge(sims_sub,truegamma,by="par",all.y=FALSE) %>% 
    mutate(bias=`50%`-trueval, pct.bias = 100*(bias/abs(trueval))) %>% 
    group_by(par) %>% 
    summarize(avg.pct.bias=mean(pct.bias, na.rm=TRUE),
              avg.bias=mean(bias, na.rm=TRUE),
              mse=mean((bias^2), na.rm=TRUE))
  
return(mrg_df)
}

# sim_gamma_summ(sim_c0_n400[['full']])
# sim_gamma_summ(sim_c0_n400[['cens']])
# 
# sim_gamma_summ(sim_c0_n25[['full']])
# sim_gamma_summ(sim_c0_n25[['cens']])

## function to summarize conditional CDF

# true cdf values
if(0){
  # relationship between logistic and log-logistic
  # https://en.wikipedia.org/wiki/Log-logistic_distribution#Related_distributions
  y0 <- exp(rlogis(1e7, log(1*1-0.5*1), 1/2)); mean(y0)
  y1 <- rllogis(1e7, shape = 2, scale = 1*1-0.5*1); mean(y1)
  ((0.5*pi)/2)/sin(pi/2)
  
  ## setting used for sims
  y0 <- exp(rlogis(1e7, 1*1-0.5*1, 1/3)); mean(y0)
  y1 <- rllogis(1e7, shape = 3, scale = exp(1*1-0.5*1)); mean(y1)
  ((exp(0.5)*pi)/3)/sin(pi/3) # analytic mean
  
  yvals <- c(-0.5, 0, 0.5, 1, 1.5)
  pllogis(exp(yvals),shape = 3, scale = exp(1*1-0.5*1))
  curve(pllogis(x,shape = 3, scale = exp(1*1-0.5*1)),0,6,xlab='y')
  segments(x0=exp(yvals),
           y0=rep(0,5),
           y1=pllogis(exp(yvals),shape = 3, scale = exp(1*1-0.5*1)),
           lty=2)
  segments(x0=rep(-3,5),
           x1=exp(yvals),
           y0=pllogis(exp(yvals),shape = 3, scale = exp(1*1-0.5*1)),
           lty=2)
  
  # at z1=0, z2=0
  pllogis(exp(yvals),shape = 3, scale = exp(0))
  curve(pllogis(x,shape = 3, scale = exp(0)),0,6,xlab='y')
  segments(x0=exp(yvals),
           y0=rep(0,5),
           y1=pllogis(exp(yvals),shape = 3, scale = exp(0)),
           lty=2)
  segments(x0=rep(-3,5),
           x1=exp(yvals),
           y0=pllogis(exp(yvals),shape = 3, scale = exp(0)),
           lty=2)
  
  # log scale
  plogis(yvals,1*1-0.5*1,1/3)
  curve(plogis(x,1*1-0.5*1,1/3),-3,3, xlab='log(y)')
  segments(x0=yvals,
           y0=rep(0,5),
           y1=plogis(yvals,1*1-0.5*1,1/3),
           lty=2)
  segments(x0=rep(-3,5),
           x1=yvals,
           y0=plogis(yvals,1*1-0.5*1,1/3),
           lty=2)
  
  # at z1=0, z2=0
  plogis(yvals,0,1/3)
  curve(plogis(x,0,1/3),-3,3, xlab='log(y)')
  segments(x0=yvals,
           y0=rep(0,5),
           y1=plogis(yvals,0,1/3),
           lty=2)
  segments(x0=rep(-3,5),
           x1=yvals,
           y0=plogis(yvals,0,1/3),
           lty=2)
  
  qlogis(c(0.1, 0.4, 0.5, 0.6, 0.9),1*1-0.5*1,1/3)  
  
  curve(dllogis(x,shape = 3, scale = exp(1*1-0.5*1)),-0.5,6)
}

sim_cdf_summ <- function(result, beta=c(1, -0.5), scale=1/3,
                          gamma.y=c(-0.5, 0, 0.5, 1, 1.5),
                          zvals = c(1,1)){
  
  if (names(result[3])=='cond.cdf'){
    sims_sub <- result[['cond.cdf']] %>% filter(z1==zvals[1] & z2==zvals[2])
  } else {
    sims_sub <- result[['cond.cdf.cn']] %>% filter(z1==zvals[1] & z2==zvals[2])
    gamma.y<-gamma.y[3:5]
  }
  
  truevals <- pllogis(exp(gamma.y),shape = 3, scale = exp(beta[1]*zvals[1]+ beta[2]*zvals[2]))
  true_df <- data.frame(yin=exp(gamma.y), true_cdf=truevals)
  
  mrg_df <- merge(sims_sub, true_df, by='yin',all.y=FALSE) %>% 
    mutate(bias=med_cdf-true_cdf, pct.bias = 100*(bias/abs(true_cdf))) %>% 
    group_by(yin) %>% 
    summarize(avg.pct.bias=mean(pct.bias, na.rm=TRUE),
              avg.bias=mean(bias, na.rm=TRUE),
              mse=mean((bias^2), na.rm=TRUE))
  
  return(mrg_df)
}


# sim_cdf_summ(sim_c0_n400[['full']])
# sim_cdf_summ(sim_c0_n400[['cens']])


## function to summarize conditional mean

# true mean values
if(0){
#beta1*z1+beta2*z2
#1*1-0.5*1
mean(exp(rlogis(1e7, 1*1-0.5*1, 1/3)))
mean(rllogis(1e7, shape = 3, scale = exp(1*1-0.5*1)))
((exp(1*1-0.5*1)*pi)/3)/sin(pi/3)

#1*1-0.5*0
mean(exp(rlogis(1e7, 1*1-0.5*0, 1/3)))
mean(rllogis(1e7, shape = 3, scale = exp(1*1-0.5*0)))
((exp(1*1-0.5*0)*pi)/3)/sin(pi/3)
}

sim_mn_summ <- function(result, beta=c(1, -0.5),
                         zdf = data.frame(z1=c(1,1),z2=c(1,0))){
  
  true_df <- data.frame(zdf,true_mn=c(((exp(beta[1]*zdf[1,1]+beta[2]*zdf[1,2])*pi)/3)/sin(pi/3),
                            ((exp(beta[1]*zdf[2,1]+beta[2]*zdf[2,2])*pi)/3)/sin(pi/3)))
  
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

# sim_mn_summ(sim_c0_n25[['full']])
# sim_mn_summ(sim_c0_n25[['cens']])


# function to summarize conditional quantile

# true quantiles
if(0){
# true median
qllogis(0.5, shape = 3, scale = exp(1*1-0.5*1))
qllogis(0.5, shape = 3, scale = exp(1*1-0.5*0))

# true 20th qtile
qllogis(0.2, shape = 3, scale = exp(1*1-0.5*1))
qllogis(0.2, shape = 3, scale = exp(1*1-0.5*0))

# true 80th qtile
qllogis(0.8, shape = 3, scale = exp(1*1-0.5*1))
qllogis(0.8, shape = 3, scale = exp(1*1-0.5*0))
}

sim_qtile_summ <- function(result, beta=c(1, -0.5),
                        zdf = data.frame(z1=c(1,1),z2=c(1,0)),
                        q = 0.5, statnm='cond.med'){
  
  true_df <- data.frame(zdf, true_qtile = 
                          c(qllogis(q, shape = 3, scale = exp(beta[1]*zdf[1,1]+beta[2]*zdf[1,2])),
                            qllogis(q, shape = 3, scale = exp(beta[1]*zdf[2,1]+beta[2]*zdf[2,2]))))
  
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

# sim_qtile_summ(sim_c0_n25[['full']], q = 0.5, statnm='cond.med')
# sim_qtile_summ(sim_c0_n25[['full']], q = 0.2, statnm='cond.q20')
# sim_qtile_summ(sim_c0_n25[['full']], q = 0.8, statnm='cond.q80')

# sim_qtile_summ(sim_c0_n25[['cens']], q = 0.5, statnm='cond.med')
# sim_qtile_summ(sim_c0_n25[['cens']], q = 0.2, statnm='cond.q20')
# sim_qtile_summ(sim_c0_n25[['cens']], q = 0.8, statnm='cond.q80')

# b<-sim_qtile_summ(sim_b0_n25[[1]], q = 0.5, statnm='cond.med')
# bb<-sim_qtile_summ(sim_b0_n400[[1]], q = 0.2, statnm='cond.q20')
# c<-sim_qtile_summ(sim_b0_n25[[2]], q = 0.5, statnm='cond.med')
# cc<-sim_qtile_summ(sim_b0_n400[[2]], q = 0.2, statnm='cond.q20')


# get summaries and convert data for all simulations
nsamps <- c(25,50,100,200,400)
concs <- c('1/J','1/(0.8 + 0.35*J)','1/(2+(J/3))','1/2')
pres <- c('sim_c0_n','sim_c1_n','sim_c2_n','sim_c3_n')
outcome <- c('full.','cens.')


# nsamps <- c(25,50)
# concs <- c('1/J','1/(0.8 + 0.35*J)')
# pres <- c('sim_c0_n','sim_c1_n')
# outcome <- c('full.','cens.')


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
    nam<-paste0(mm,nn,"sim")
    dnam<-gsub("\\.","_",nam)
    rmnm<-gsub("\\.","\\\\.",nam)
    print(nam); print(dnam); print(rmnm)
    out<-lapply(ls()[grep(nam,ls())],get) %>% bind_rows() %>%
      mutate(n=factor(n,levels=c(400,200,100,50,25)))
    assign(paste0(dnam,"_dat"),out)
    rm(list=ls()[grep(rmnm,ls())])
  }
}

if(0){
## full outcome datasets
full_beta_dat <- bind_rows(lapply(ls()[grep('full.beta.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('full.beta.sim',ls())])

full_gamma_dat <- bind_rows(lapply(ls()[grep('full.gamma.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('full.gamma.sim',ls())])

full_cdf_dat <- bind_rows( lapply(ls()[grep('full.cdf.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('full.cdf.sim',ls())])

full_mn_dat <- bind_rows( lapply(ls()[grep('full.mn.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('full.mn.sim',ls())])

full_med_dat <- bind_rows( lapply(ls()[grep('full.med.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('full.med.sim',ls())])

full_q20_dat <- bind_rows( lapply(ls()[grep('full.q20.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('full.q20.sim',ls())])

full_q80_dat <- bind_rows( lapply(ls()[grep('full.q80.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('full.q80.sim',ls())])

## censored outcome dataset

cens_beta_dat <- bind_rows( lapply(ls()[grep('cens.beta.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('cens.beta.sim',ls())])

cens_gamma_dat <- bind_rows( lapply(ls()[grep('cens.gamma.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('cens.gamma.sim',ls())])

cens_cdf_dat <- bind_rows( lapply(ls()[grep('cens.cdf.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('cens.cdf.sim',ls())])

cens_mn_dat <- bind_rows( lapply(ls()[grep('cens.mn.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('cens.mn.sim',ls())])

cens_med_dat <- bind_rows( lapply(ls()[grep('cens.med.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('cens.med.sim',ls())])

cens_q20_dat <- bind_rows( lapply(ls()[grep('cens.q20.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('cens.q20.sim',ls())])

cens_q80_dat <- bind_rows( lapply(ls()[grep('cens.q80.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))
rm(list=ls()[grep('cens.q80.sim',ls())])
}


## !!! Start HERE !!! ##

# labels and plot params
# https://stackoverflow.com/questions/62662144/conditional-probability-is-not-displaying-properly-in-facet-label-in-ggplot2
nlabs <- paste0('F(y==e^',c(-0.5, 0, 0.5, 1, 1.5),"*'|'*",'~X[1]==1,X[2]==1)')
pltw<-10; plth<-5; atxtsz<-9; fctsiz<-13

### betas & gamma ###

# full outcome plot
bind_rows(full_beta_sim_dat,full_gamma_sim_dat) %>% 
  ggplot(aes(x=avg.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75) +
  facet_grid(. ~ par) +
  xlab("average bias of posterior parameters") + ylab("sample size")

# censored outcome plot
bind_rows(cens_beta_sim_dat,cens_gamma_sim_dat) %>% 
  filter(!is.na(avg.bias)) %>% 
  ggplot(aes(x=avg.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75) +
  facet_grid(. ~ par,drop=TRUE) +
  xlab("average bias of posterior parameters") + ylab("sample size")

### CDF ###

# full outcome plot
full_cdf_sim_dat %>% mutate(yin=factor(yin,labels=nlabs)) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ yin, labeller=labeller(yin=label_parsed))+
  xlab("average percent bias of posterior conditional CDF") + ylab("sample size")

#ggsave(file.path(figdir,"sim_cdf_full.png"),width=pltw,height=plth)

# censored outcome plot
cens_cdf_sim_dat %>% mutate(yin=factor(yin,labels=nlabs[3:5])) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ yin, labeller=labeller(yin=label_parsed))+
  xlab("average percent bias of posterior conditional CDF") + ylab("sample size")

#ggsave(file.path(figdir,"sim_cdf_cens.png"),width=pltw,height=plth)


### Mean ###

# full outcome plot
full_mn_sim_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "E(Y*'|'*~X[1]==1,X[2]==1)",
                       "E(Y*'|'*~X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional mean") + ylab("sample size")

# +coord_cartesian(xlim=c(-4,4))

#ggsave(file.path(figdir,"sim_mn_full.png"),width=pltw,height=plth)

# censored outcome plot
# expected to be biased because of censored y vals, can't really get unbiased est.
cens_mn_sim_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "E(Y*'|'*~X[1]==1,X[2]==1)",
                       "E(Y*'|'*~X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional mean") + ylab("sample size")

#ggsave(file.path(figdir,"sim_mn_cens.png"),width=pltw,height=plth)



### Median ###

# full outcome plot
full_med_sim_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional median") + ylab("sample size") +
  coord_cartesian(xlim=c(-15,15))

#ggsave(file.path(figdir,"sim_med_full.png"),width=pltw,height=plth)

# censored outcome plot
cens_med_sim_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional median") + ylab("sample size") +
  coord_cartesian(xlim=c(-15,15))

#ggsave(file.path(figdir,"sim_med_cens.png"),width=pltw,height=plth)



### 20% quantile ###

# full outcome plot
full_q20_sim_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.2}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.2}*'|'*list(X[1]==1,X[2]==0)")) %>%
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional 20th percentile") + ylab("sample size")

#ggsave(file.path(figdir,"sim_q20_full.png"),width=pltw,height=plth)

# censored outcome plot
# expected to be biased because of censored y vals, can't really get unbiased est.
cens_q20_sim_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.2}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.2}*'|'*list(X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional 20th percentile") + ylab("sample size")

#ggsave(file.path(figdir,"sim_q20_cens.png"),width=pltw,height=plth)



### 80% quantile ###

# full outcome plot
full_q80_sim_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.8}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.8}*'|'*list(X[1]==1,X[2]==0)")) %>%
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional 80th percentile") + ylab("sample size")

#ggsave(file.path(figdir,"sim_q80_full.png"),width=pltw,height=plth)

# censored outcome plot
# expected to be biased because of censored y vals, can't really get unbiased est.
cens_q80_sim_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.8}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.8}*'|'*list(X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional 80th percentile") + ylab("sample size")

#ggsave(file.path(figdir,"sim_q80_cens.png"),width=pltw,height=plth)

