rm(list=ls())
libs <- c("dplyr", "stringr", "readr", "tidyr", "purrr", "ggplot2")
invisible(lapply(libs, library, character.only = TRUE))

#paper directory (sims for other functions)
pdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper")
psdir <- file.path(pdir,"sims","out")

figdir <- file.path(pdir,"fig")

#source(file.path(dir,"cpm_functions.r"))

# load sim array
simarray <- readRDS(file.path(pdir,"sims","bayes_cpm_simarray.rds"))

# load sim data
for (j in 1:25){
  nam_b0 <- paste0("sim_b0_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_b0 <- file.path(psdir,paste0(nam_b0,".rds"))
  nam_b1 <- paste0("sim_b1_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_b1 <- file.path(psdir,paste0(nam_b1,".rds"))
  nam_b2 <- paste0("sim_b2_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_b2 <- file.path(psdir,paste0(nam_b2,".rds"))
  try(assign(nam_b0, readRDS(fp_b0)))
  try(assign(nam_b1, readRDS(fp_b1)))
  try(assign(nam_b2, readRDS(fp_b2)))
}

rm(nam_b0,nam_b1,nam_b2,fp_b0,fp_b1,fp_b2,j)

# for reference
if (0){
# from Liu et al. sim code
generate.data.1 <- function(seed=1, n=50, p=0.5, alpha=0, beta=c(1.0, -0.5), sigma=1){
  set.seed(seed)
  z1 <- sample(c(0,1), size=n, replace=TRUE, prob=c(1-p, p))
  z2 <- rnorm(n, 0, 1)
 # y0 <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
 # y <- exp(y0) #log normal
  log.y <- rnorm(n, alpha+beta[1]*z1 + beta[2]*z2, sigma)
  data <- data.frame(log.y=ordered(log.y), z1=z1, z2=z2)
  return(data)
}

get_cutpoints <- function(fit,cps=cp_seq,summ_stat='50%',pr=c(0.5)){
    summary(fit, pars=paste0("cutpoints[",cps,"]"), probs=pr)$summary[,summ_stat]
  }

cutpoint_est <- function(post,fitdata,y){
    truey<-c(-Inf,fitdata$truey0,Inf)

    if (y<=min(fitdata$truey0)){
      return(-Inf)
    }
    else if (y > max(fitdata$truey0)) {
      return(Inf)
    }
    else {
      post[which(fitdata$truey0>=y)[1]-1]
    }
}

sim1_coeff.fun <- function(sim=5, seed=1, n=50, p=0.5, alpha=0, beta=c(1,-0.5),
                          sigma=1, log.y=c(-1, -0.33, 0.5, 1.33, 2)){
  set.seed(seed)
  seeds <- unique(round(runif(sim*10,0,10)*sim,0))[seq(sim)]

  # cutpoints
  alpha.y <- matrix(NA, ncol=length(log.y), nrow=sim)
  alpha.y.se <- matrix(NA, ncol=length(log.y), nrow=sim)
  alpha.y_cn <- matrix(NA, ncol=length(log.y), nrow=sim)
  alpha.y.se_cn <- matrix(NA, ncol=length(log.y), nrow=sim)

  # betas
  md.beta.est <- matrix(NA, ncol=length(beta), nrow=sim)
  mn.beta.est <- matrix(NA, ncol=length(beta), nrow=sim)
  beta.se <- matrix(NA, ncol=length(beta), nrow=sim)
  md.beta.est_cn <- matrix(NA, ncol=length(beta), nrow=sim)
  mn.beta.est_cn <- matrix(NA, ncol=length(beta), nrow=sim)
  beta.se_cn <- matrix(NA, ncol=length(beta), nrow=sim)

  for(i in 1:sim){
    print(i)
    try({
      data <- generate.data.1(seed=seeds[i], n=n, p=p,
                              alpha=alpha, beta=beta, sigma=sigma)

      ycens <- as.numeric(levels(data$log.y)[as.numeric(data$log.y)])
      ycens[ycens < 0]<- 0
      data$log.y.cens<-ordered(ycens)

      mod_data <- mkStanDat(data, outcome="log.y",
                            preds = c("z1", "z2"),
                            link=2) #probit link

      cpm_fit <- sampling(ord_mod1, data=mod_data, seed=seeds[i],
                 iter=4000, warmup=2000, chains=2, refresh=2000,
                 control = list(adapt_delta = 0.9))

      fit_summ_beta<-summary(cpm_fit, pars=c("b[1]","b[2]"),
                             use_cache=FALSE)$summary

      md.beta.est[i,] <- fit_summ_beta[,'50%']
      mn.beta.est[i,] <- fit_summ_beta[,'mean']
      beta.se[i, ] <- fit_summ_beta[,'sd']

      # median and sd of posterior cutpoints
      cp_seq <- 1:(mod_data$ncat - 1)

      cp_md <- get_cutpoints(cpm_fit, cps=cp_seq)
      cp_sd <- get_cutpoints(cpm_fit, cps=cp_seq, summ_stat='sd')

      # cutpoint estimates at given log(y) values
      alpha.y[i,] <- sapply(log.y, function(x) cutpoint_est(cp_md,mod_data,x))
      alpha.y.se[i,] <- sapply(log.y, function(x) cutpoint_est(cp_sd,mod_data,x))

      mod_data_cn <- mkStanDat(data, outcome="log.y.cens",
                            preds = c("z1", "z2"),
                            link=2) #probit link

      cpm_fit_cn <- sampling(ord_mod1, data=mod_data_cn, seed=seeds[i],
                 iter=4000, warmup=2000, chains=2, refresh=2000,
                 control = list(adapt_delta = 0.9))

      fit_summ_beta_cn<-summary(cpm_fit_cn, pars=c("b[1]","b[2]"),
                             use_cache=FALSE)$summary

      md.beta.est_cn[i,] <- fit_summ_beta_cn[,'50%']
      mn.beta.est_cn[i,] <- fit_summ_beta_cn[,'mean']
      beta.se_cn[i, ] <- fit_summ_beta_cn[,'sd']

      # median and sd of posterior cutpoints
      cp_seq_cn <- 1:(mod_data_cn$ncat - 1)
      cp_md_cn <- get_cutpoints(cpm_fit_cn, cps=cp_seq_cn)
      cp_sd_cn <- get_cutpoints(cpm_fit_cn, cps=cp_seq_cn, summ_stat='sd')

      # cutpoint estimates at given log(y) values
      alpha.y_cn[i,] <- sapply(log.y, function(x) cutpoint_est(cp_md_cn,mod_data_cn,x))
      alpha.y.se_cn[i,] <- sapply(log.y, function(x) cutpoint_est(cp_sd_cn,mod_data_cn,x))

    })
  }

  return(list(full=list(md.beta.est=md.beta.est,
                    mn.beta.est=mn.beta.est,
                    beta.se=beta.se,
                    alpha.y=alpha.y,
                    alpha.y.se=alpha.y.se),
              cens=list(md.beta.est_cn=md.beta.est_cn,
                    mn.beta.est_cn=mn.beta.est_cn,
                    beta.se_cn=beta.se_cn,
                    alpha.y_cn=alpha.y_cn,
                    alpha.y.se_cn=alpha.y.se_cn)))
}

}

# sim_n - probit link, conc=1/ncats
# sim_a1_n - probit link, conc=1/(0.8 + 0.35*max(ncats, 3))
# sim_a2_n - probit link, conc=1/(2+(ncats/3))

# combine 5 reps of 200 sims for each setting
for (n in c(25,50,100,200,400)){
  for (pre in c("sim_b0_n","sim_b1_n","sim_b2_n")){
    nm0 <- paste0(pre, n)
    nm <- paste0(nm0, "_", 1:5)
    nmlist <- map(nm,get)
    
    fl0 <- map(nmlist, function(x) x[[1]])
    fl1 <- pmap(fl0,rbind)
    fl2 <- lapply(fl1,bind_rows)
    
    cn0 <- map(nmlist, function(x) x[[2]])
    cn1 <- pmap(cn0,rbind)
    cn2 <- lapply(cn1,bind_rows)
    
    assign(nm0,list(full=fl2,cens=cn2))
    rm(list=nm)
  }
}

rm(cn0,cn1,cn2,fl0,fl1,fl2,nmlist,nm,nm0,n,pre)

# true cdf values
if(0){
# beta[1]*z1+beta[2]*z2
plnorm(exp(c(-1, -0.33, 0.5, 1.33, 2)),1*1-0.5*1,1)
curve(plnorm(x,1*1-0.5*1,1),0,10,xlab='y')
segments(x0=exp(c(-1, -0.33, 0.5, 1.33, 2)),
         y0=rep(0,5),
         y1=plnorm(exp(c(-1, -0.33, 0.5, 1.33, 2)),1*1-0.5*1,1),
         lty=2)
segments(x0=rep(-3,5),
         x1=exp(c(-1, -0.33, 0.5, 1.33, 2)),
         y0=plnorm(exp(c(-1, -0.33, 0.5, 1.33, 2)),1*1-0.5*1,1),
         lty=2)
  
pnorm(c(-1, -0.33, 0.5, 1.33, 2),1*1-0.5*1,1)
curve(pnorm(x,1*1-0.5*1,1),-3,3, xlab='log(y)')
segments(x0=c(-1, -0.33, 0.5, 1.33, 2),
         y0=rep(0,5),
         y1=pnorm(c(-1, -0.33, 0.5, 1.33, 2),1*1-0.5*1,1),
         lty=2)
segments(x0=rep(-3,5),
         x1=c(-1, -0.33, 0.5, 1.33, 2),
         y0=pnorm(c(-1, -0.33, 0.5, 1.33, 2),1*1-0.5*1,1),
         lty=2)

curve(dnorm(x,1*1-0.5*1,1),-2,2.5)
}

## summarize conditional CDF
sim_cdf_summ <- function(result, beta=c(1, -0.5),
                          alpha.y=c(-1, -0.33, 0.5, 1.33, 2),
                          zvals = c(1,1)){
  
  if (names(result[1])=='cond.cdf'){
    sims_sub <- result[['cond.cdf']] %>% filter(z1==zvals[1] & z2==zvals[2])
  } else {
    sims_sub <- result[['cond.cdf.cn']] %>% filter(z1==zvals[1] & z2==zvals[2])
    alpha.y<-alpha.y[3:5]
  }
  
  truevals <- pnorm(alpha.y,beta[1]*zvals[1]+ beta[2]*zvals[2])
  true_df <- data.frame(yin=alpha.y, true_cdf=truevals)
  
  mrg_df <- merge(sims_sub, true_df, by='yin',all.y=FALSE) %>% 
    mutate(bias=med_cdf-true_cdf, pct.bias = 100*(bias/abs(true_cdf))) %>% 
    group_by(yin) %>% 
    summarize(avg.pct.bias=mean(pct.bias))
  
  return(mrg_df)
}



## summarize conditional mean

# true mean
#beta1*z1+beta2*z2
1*1-0.5*1
1*1-0.5*0

sim_mn_summ <- function(result, beta=c(1, -0.5),
                         zdf = data.frame(z1=c(1,1),z2=c(1,0))){
  
  true_df <- data.frame(zdf,true_mn=c(sum(zdf[1]*beta), sum(zdf[2]*beta)))

  if (names(result[2])=='cond.mn'){
    sims_sub <- result[['cond.mn']] 
  } else {
    sims_sub <- result[['cond.mn.cn']]
  }
  
  mrg_df <- merge(sims_sub, true_df, by=c('z1','z2'), all.y=FALSE) %>% 
    mutate(bias=med_mn-true_mn, pct.bias = 100*(bias/abs(true_mn))) %>% 
     group_by(ndrow,z1,z2) %>% 
     summarize(avg.pct.bias=mean(pct.bias))
  
  return(mrg_df)
}

# summarize conditional quantile

# true median
qnorm(0.5,1*1-0.5*1,1)
qnorm(0.5,1*1-0.5*0,1)

# true 20th qtile
qnorm(0.2,1*1-0.5*1,1)
qnorm(0.2,1*1-0.5*0,1)

sim_qtile_summ <- function(result, beta=c(1, -0.5),
                        zdf = data.frame(z1=c(1,1),z2=c(1,0)),
                        q = 0.5, statnm='cond.med'){
  
  true_df <- data.frame(zdf,true_qtile=c(qnorm(q,sum(zdf[1]*beta),1),
                                       qnorm(q,sum(zdf[2]*beta),1)))
  
  if (names(result[4])==statnm|names(result[6])==statnm){
    sims_sub <- result[[statnm]] 
  } else {
    sims_sub <- result[[paste0(statnm,".cn")]]
  }
  

  mrg_df <- merge(sims_sub, true_df, by=c('z1','z2'), all.y=FALSE) %>% 
    mutate(bias=med_qtile-true_qtile, pct.bias = 100*(bias/abs(true_qtile))) %>% 
    group_by(ndrow,z1,z2) %>% 
    summarize(avg.pct.bias=mean(pct.bias))
  
  return(mrg_df)
}


# b<-sim_qtile_summ(sim_b0_n25[[1]], q = 0.5, statnm='cond.med')
# bb<-sim_qtile_summ(sim_b0_n400[[1]], q = 0.2, statnm='cond.q20')
# c<-sim_qtile_summ(sim_b0_n25[[2]], q = 0.5, statnm='cond.med')
# cc<-sim_qtile_summ(sim_b0_n400[[2]], q = 0.2, statnm='cond.q20')


# get summaries and convert data for all simulations
nsamps <- c(25,50,100,200,400)
concs <- c('1/J','1/(0.8 + 0.35*J)','1/(2+(J/3))')
pres <- c('sim_b0_n','sim_b1_n','sim_b2_n')
outcome <- c('full.','cens.')

for (k in seq_along(outcome)){
  for (j in seq_along(pres)){
    for (i in nsamps){
      nm0 <- paste0(pres[j], i)
      nm_cdf <- paste0(outcome[k],'cdf.',nm0)
      nm_mn <- paste0(outcome[k],'mn.',nm0)
      nm_med <- paste0(outcome[k],'med.',nm0)
      nm_q20 <- paste0(outcome[k],'q20.',nm0)
      
      cdf <- sim_cdf_summ(get(nm0)[[k]]) %>% mutate(n=i, conc=concs[j])
      mn <- sim_mn_summ(get(nm0)[[k]]) %>% mutate(n=i, conc=concs[j])
      med <- sim_qtile_summ(get(nm0)[[k]], q=0.5, statnm='cond.med') %>% mutate(n=i, conc=concs[j]) 
      q20 <- sim_qtile_summ(get(nm0)[[k]], q=0.2, statnm='cond.q20') %>% mutate(n=i, conc=concs[j]) 
      
      assign(nm_cdf,cdf)
      assign(nm_mn,mn)
      assign(nm_med,med)
      assign(nm_q20,q20)
      
      rm(cdf,mn,med,q20,nm0,nm_cdf,nm_mn,nm_med,nm_q20)
    }
  }
}


## full outcome datasets
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

## censored outcome dataset
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

# labels and plot params
# https://stackoverflow.com/questions/62662144/conditional-probability-is-not-displaying-properly-in-facet-label-in-ggplot2
nlabs <- paste0('F(y==e^',c(-1, -0.33, 0.5, 1.33, 2),"*'|'*",'~X[1]==1,X[2]==1)')
pltw<-10; plth<-5; atxtsz<-9; fctsiz<-13

### CDF ###

# full outcome plot
full_cdf_dat %>% mutate(yin=factor(yin,labels=nlabs)) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ yin, labeller=labeller(yin=label_parsed))+
  xlab("average percent bias of posterior conditional CDF") + ylab("sample size")

ggsave(file.path(figdir,"sim_cdf_full.png"),width=pltw,height=plth)

# censored outcome plot
cens_cdf_dat %>% mutate(yin=factor(yin,labels=nlabs[3:5])) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ yin, labeller=labeller(yin=label_parsed))+
  xlab("average percent bias of posterior conditional CDF") + ylab("sample size")

ggsave(file.path(figdir,"sim_cdf_cens.png"),width=pltw,height=plth)


### Mean ###

# full outcome plot
full_mn_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "E(Y*'|'*~X[1]==1,X[2]==1)",
                       "E(Y*'|'*~X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional mean") + ylab("sample size")+
  coord_cartesian(xlim=c(-4,4))

ggsave(file.path(figdir,"sim_mn_full.png"),width=pltw,height=plth)

# censored outcome plot
# expected to be biased because of censored y vals, can't really get unbiased est.
cens_mn_dat %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow)+
  xlab("average percent bias of posterior conditional mean") + ylab("sample size")

#ggsave(file.path(figdir,"sim_mn_cens.png"),width=pltw,height=plth)



### Median ###

# full outcome plot
full_med_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional median") + ylab("sample size") +
  coord_cartesian(xlim=c(-15,15))

ggsave(file.path(figdir,"sim_med_full.png"),width=pltw,height=plth)

# censored outcome plot
cens_med_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.5}*'|'*list(X[1]==1,X[2]==0)")) %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional median") + ylab("sample size") +
  coord_cartesian(xlim=c(-15,15))

ggsave(file.path(figdir,"sim_med_cens.png"),width=pltw,height=plth)



### 20% quantile ###

# full outcome plot
full_q20_dat %>% 
  mutate(ndrow=if_else(ndrow==1,
                       "Q^{0.2}*'|'*list(X[1]==1,X[2]==1)",
                       "Q^{0.2}*'|'*list(X[1]==1,X[2]==0)")) %>%
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow, labeller=labeller(ndrow=label_parsed))+
  xlab("average percent bias of posterior conditional 20th percentile") + ylab("sample size")

ggsave(file.path(figdir,"sim_q20_full.png"),width=pltw,height=plth)

# censored outcome plot
# expected to be biased because of censored y vals, can't really get unbiased est.
cens_q20_dat %>% 
  ggplot(aes(x=avg.pct.bias,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ ndrow)+
  xlab("average percent bias of posterior conditional 20th percentile") + ylab("sample size")

#ggsave(file.path(figdir,"sim_q20_cens.png"),width=pltw,height=plth)

