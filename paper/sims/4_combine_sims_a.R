rm(list=ls())
libs <- c("rstan", "dplyr", "stringr", "readr", "tidyr", "purrr", "ggplot2")
invisible(lapply(libs, library, character.only = TRUE))

#seminar directory (sims for 1/J)
sdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/biostat_sem")
ssdir <- file.path(sdir,"sims","out")

#paper directory (sims for other functions)
pdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper")
psdir <- file.path(pdir,"sims","out")

figdir <- file.path(pdir,"fig")

# load sim array
simarray <- readRDS(file.path(sdir,"sims","bayes_cpm_simarray.rds"))

# load sim data
for (j in 1:25){
  # load sims from seminar 
  nam <- paste0("sim_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp <- file.path(ssdir,paste0(nam,".rds"))
  assign(nam, readRDS(fp))
  
  # load additional sims for paper
  nam_a1 <- paste0("sim_a1_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_a1 <- file.path(psdir,paste0(nam_a1,".rds"))
  nam_a2 <- paste0("sim_a2_n", simarray[j,"nsamps"], "_",simarray[j,"rep"])
  fp_a2 <- file.path(psdir,paste0(nam_a2,".rds"))
  try(assign(nam_a1, readRDS(fp_a1)))
  try(assign(nam_a2, readRDS(fp_a2)))
}

rm(nam,nam_a1,nam_a2,fp,fp_a1,fp_a2,j)

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
  for (pre in c("sim_n","sim_a1_n","sim_a2_n")){
  nm0 <- paste0(pre, n)
  nm <- paste0(nm0, "_", 1:5)
  nmlist <- map(nm,get)
  
  fl <- map(nmlist, function(x) x[[1]]) %>% pmap(rbind)
  cn <- map(nmlist, function(x) x[[2]]) %>% pmap(rbind)
  
  assign(nm0,list(full=fl,cens=cn))
  rm(list=nm)
  }
}

rm(cn,fl,nmlist,nm,nm0,n,pre)


sim_summary <- function(result, md=TRUE, beta=c(1, -0.5),
                              alpha.y=c(-1, -0.33, 0.5, 1.33, 2)){
  
  # mn.beta.est or md.beta.est
  bt <- result[[2]] # result$mn.beta.est
  if(md) bt <- result[[1]] #result$md.beta.est
  
  #alph.cln <- ifelse(is.infinite(result$alpha.y), NA, result$alpha.y)
  alph.cln <- ifelse(is.infinite(result[[4]]), NA, result[[4]])
  
  truebt <- matrix(beta,byrow=TRUE,nrow=nrow(bt),ncol=length(beta))
  truealph <- matrix(alpha.y,byrow=TRUE,nrow=nrow(alph.cln),ncol=length(alpha.y))
  
  # bias in mean or median
  bias.bt <- bt - truebt
  pct.bias.bt <- 100*(bias.bt/abs(truebt))
  
  # bias in cutpoints
  bias.alpha <- alph.cln - truealph
  pct.bias.alpha <- 100*(bias.alpha/abs(truealph))
  
  # average percent bias
  result <- colMeans(cbind(pct.bias.bt,pct.bias.alpha), na.rm=TRUE) %>% 
    t() %>% data.frame()
  
  names(result)<-c("beta[1]","beta[2]",
                   "gamma[y[1]]","gamma[y[2]]",
                   "gamma[y[3]]","gamma[y[4]]","gamma[y[5]]")
  return(result)
}

# sim_summary(sim_a2_n25[[1]]) 



# get summaries and convert data for all simulations
nsamps <- c(25,50,100,200,400)
concs <- c('1/J','1/(0.8 + 0.35*J)','1/(2+(J/3))')
pres <- c('sim_n','sim_a1_n','sim_a2_n')
outcome <- c('full.','cens.')

for (k in seq_along(outcome)){
  for (j in seq_along(pres)){
    for (i in nsamps){
      nm0 <- paste0(pres[j], i)
      dd <- sim_summary(get(nm0)[[k]]) %>% mutate(n=i, conc=concs[j])
      nm <- paste0(outcome[k],nm0)
      assign(nm,dd)
    }
  }
}

#full.sim_a1_n100

## full outcome dataset
full_dat <- bind_rows( lapply(ls()[grep('full.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25)))

# for plotting
full_plt_dat <- full_dat %>% 
  pivot_longer(cols=c("beta[1]","beta[2]",
               "gamma[y[1]]","gamma[y[2]]",
               "gamma[y[3]]","gamma[y[4]]","gamma[y[5]]"))

## censored outcome dataset
cens_dat <- bind_rows( lapply(ls()[grep('cens.sim',ls())],get) ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25))) %>%
  select(-c("gamma[y[1]]","gamma[y[2]]"))

# for plotting
cens_plt_dat <- cens_dat %>%  
  pivot_longer(cols=c("beta[1]","beta[2]",
                      "gamma[y[3]]","gamma[y[4]]","gamma[y[5]]"))


# labels and plot params
#! nlabs <- c("average bias of \nposterior median (%)","average bias of \nposterior se (%)")
#! names(nlabs)<-c("bias.est","bias.se")

pltw<-10; plth<-5; atxtsz<-9; fctsiz<-13

# full outcome plot
full_plt_dat %>% 
  ggplot(aes(x=value,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ name,labeller = labeller(name=label_parsed))+
  xlab("average percent bias of posterior median")+ylab("sample size")+
    theme(axis.title.x = element_text(size=fctsiz),
          axis.title.y = element_text(size=fctsiz),
          axis.text =  element_text(size=atxtsz),
          strip.text = element_text(size=fctsiz),
          strip.text.y = element_text(angle=0))

ggsave(file.path(figdir,"sim_a_pars_full.png"),width=pltw,height=plth)

# censored outcome plot
cens_plt_dat %>% 
  ggplot(aes(x=value,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(. ~ name,labeller = labeller(name=label_parsed))+
  xlab("average percent bias of posterior median")+ylab("sample size")+
  theme(axis.title.x = element_text(size=fctsiz),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))

ggsave(file.path(figdir,"sim_a_pars_cens.png"),width=pltw,height=plth)
