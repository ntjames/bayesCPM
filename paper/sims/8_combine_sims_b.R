rm(list=ls())
#libs <- c("rstan", "dplyr", "stringr", "readr", "tidyr","pammtools","evd","sfsmisc")
libs <- c("rstan", "dplyr", "stringr", "readr", "tidyr","purrr", "ggplot2")
invisible(lapply(libs, library, character.only = TRUE))

# set directories

#paper directory (sims for other functions)
pdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper")
psdir <- file.path(pdir,"sims","out")

figdir <- file.path(pdir,"fig")

#source(file.path(dir,"cpm_functions.r"))

# load sim array
simarray <- readRDS(file.path(sdir,"sims","bayes_cpm_simarray.rds"))

# load addtl sims for paper
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

##!!!!! modify below for sims b (conditional CDF, mean, median, q20 and contrasts)

# sims, probit link, conc=1/ncats
for (i in c(25,50,100,200,400)){
nm0 <- paste0("sim_b0_n", i)
nm <- paste0(nm0, "_", 1:5)
nmlist <- map(nm,get)

fl0 <- map(nmlist, function(x) x[[1]]) #full outcome
cn0 <- map(nmlist, function(x) x[[2]]) #censored outcome

#!!! modify to combine lists
fl <- pmap(fl0,rbind)
cn <- pmap(cn0,rbind)

assign(nm0,list(full=fl,cens=cn))
rm(list=nm)

}

rm(cn,cn0,fl,fl0,nmlist,nm,nm0)

# sims probit link, conc=1/(0.8 + 0.35*max(ncats, 3))
for (i in c(25,50,100,200,400)){
  nm0 <- paste0("sim_b1_n", i)
  nm <- paste0(nm0, "_", 1:5)
  nmlist <- map(nm,get)
  
  fl0 <- map(nmlist, function(x) x[[1]])
  cn0 <- map(nmlist, function(x) x[[2]])
  
  fl <- pmap(fl0,rbind)
  cn <- pmap(cn0,rbind)
  
  assign(nm0,list(full=fl,cens=cn))
  rm(list=nm)
}

rm(cn,cn0,fl,fl0,nmlist,nm,nm0)

# sims probit link, conc=1/(2+(ncats/3))
for (i in c(25,50,100,200,400)){
  nm0 <- paste0("sim_b2_n", i)
  nm <- paste0(nm0, "_", 1:5)
  nmlist <- map(nm,get)
  
  fl0 <- map(nmlist, function(x) x[[1]])
  cn0 <- map(nmlist, function(x) x[[2]])
  
  fl <- pmap(fl0,rbind)
  cn <- pmap(cn0,rbind)
  
  assign(nm0,list(full=fl,cens=cn))
  rm(list=nm)
}

rm(cn,cn0,fl,fl0,nmlist,nm,nm0)


sim_coeff.summary <- function(result, md=FALSE, beta=c(1, -0.5),
                              alpha.y=c(-1, -0.33, 0.5, 1.33, 2)){

  # mn.beta.est or md.beta.est
  bt <- result[[2]] # result$mn.beta.est
  if(md) bt <- result[[1]] #result$md.beta.est

  bt.se <- result[[3]] #result$beta.se
  #alph.cln <- ifelse(is.infinite(result$alpha.y), NA, result$alpha.y)
  alph.cln <- ifelse(is.infinite(result[[4]]), NA, result[[4]])
  #alph.se.cln <- ifelse(is.infinite(result$alpha.y.se), NA, result$alpha.y.se)
  alph.se.cln <- ifelse(is.infinite(result[[5]]), NA, result[[5]])

  mean.est <- c(colMeans(bt, na.rm=TRUE),
                colMeans(alph.cln, na.rm=TRUE))

  mean.se <- c(colMeans(bt.se, na.rm=TRUE),
               colMeans(alph.se.cln, na.rm=TRUE))

  emp.se <- c(apply(bt, 2, function(x) sd(x, na.rm=TRUE)),
              apply(alph.cln, 2, function(x) sd(x, na.rm=TRUE)))

  mse <- c(colMeans((bt-matrix(beta, nrow=dim(bt)[1],
                               ncol=length(beta),byrow=TRUE))^2, na.rm=TRUE),
           colMeans((alph.cln-matrix(alpha.y, nrow=dim(alph.cln)[1],
                                     ncol=length(alpha.y), byrow=TRUE))^2,
                    na.rm=TRUE))

  # coverage - need to get 95% cred interval for Bayes version of this
  if(0){
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
  }

  result <- cbind(format(c(beta, alpha.y), nsmall=1),
                  format(round(mean.est, 2), nsmall=2),
                  format(round(mean.se, 3) , nsmall=3),
                  format(round(emp.se, 3) , nsmall=3),
                  format(round(mse, 4) , nsmall=3))
  return(result)

}

# first 2 rows are betas, next 5 rows are log.y values

conv <- function(summ){
  out0<-sim_coeff.summary(summ) %>%
    as.data.frame(stringsAsFactors=FALSE) %>%
    mutate_all(as.numeric)

  names(out0)<-c("true","mean.est","mean.se","emp.se","mse")

  out<-out0 %>% mutate(bias.est=100*(mean.est-true)/true,
                       bias.se=100*(mean.se-emp.se)/emp.se,
                       param=c("beta[1]","beta[2]",
                             "gamma[y[1]]","gamma[y[2]]",
                             "gamma[y[3]]","gamma[y[4]]","gamma[y[5]]"))

return(out)
}


## full dataset
sim_coeff.a0.n25.full <- conv(sim_n25[[1]]) %>% mutate(n=25,conc='1/J')
sim_coeff.a0.n50.full <- conv(sim_n50[[1]]) %>% mutate(n=50,conc='1/J')
sim_coeff.a0.n100.full <- conv(sim_n100[[1]]) %>% mutate(n=100,conc='1/J')
sim_coeff.a0.n200.full <- conv(sim_n200[[1]]) %>% mutate(n=200,conc='1/J')
sim_coeff.a0.n400.full <- conv(sim_n400[[1]]) %>% mutate(n=400,conc='1/J')

sim_coeff.a1.n25.full <- conv(sim_a1_n25[[1]]) %>% mutate(n=25,conc='1/(0.8 + 0.35*J)')
sim_coeff.a1.n50.full <- conv(sim_a1_n50[[1]]) %>% mutate(n=50,conc='1/(0.8 + 0.35*J)')
sim_coeff.a1.n100.full <- conv(sim_a1_n100[[1]]) %>% mutate(n=100,conc='1/(0.8 + 0.35*J)')
sim_coeff.a1.n200.full <- conv(sim_a1_n200[[1]]) %>% mutate(n=200,conc='1/(0.8 + 0.35*J)')
#sim_coeff.a1.n400.full <- conv(sim_a1_n400[[1]]) %>% mutate(n=400,conc='1/(0.8 + 0.35*J)')

sim_coeff.a2.n25.full <- conv(sim_a2_n25[[1]]) %>% mutate(n=25,conc='1/(2+(J/3))')
sim_coeff.a2.n50.full <- conv(sim_a2_n50[[1]]) %>% mutate(n=50,conc='1/(2+(J/3))')
sim_coeff.a2.n100.full <- conv(sim_a2_n100[[1]]) %>% mutate(n=100,conc='1/(2+(J/3))')
sim_coeff.a2.n200.full <- conv(sim_a2_n200[[1]]) %>% mutate(n=200,conc='1/(2+(J/3))')
sim_coeff.a2.n400.full <- conv(sim_a2_n400[[1]]) %>% mutate(n=400,conc='1/(2+(J/3))')



full_dat <- bind_rows(sim_coeff.a0.n25.full,
                      sim_coeff.a0.n50.full,
                      sim_coeff.a0.n100.full,
                      sim_coeff.a0.n200.full,
                      sim_coeff.a0.n400.full,
                      sim_coeff.a1.n25.full,
                      sim_coeff.a1.n50.full,
                      sim_coeff.a1.n100.full,
                      sim_coeff.a1.n200.full,
                      #sim_coeff.a1.n400.full,
                      sim_coeff.a2.n25.full,
                      sim_coeff.a2.n50.full,
                      sim_coeff.a2.n100.full,
                      sim_coeff.a2.n200.full,
                      sim_coeff.a2.n400.full
                      ) %>%
                      mutate(n=factor(n,levels=c(400,200,100,50,25)))

full_plt_dat <- full_dat %>% select(bias.est,bias.se,param,n,conc) %>%
  pivot_longer(cols=c(bias.est,bias.se))

nlabs <- c("average bias of \nposterior median (%)","average bias of \nposterior se (%)")
names(nlabs)<-c("bias.est","bias.se")

pltw<-10
plth<-5
atxtsz<-9
fctsiz<-13

full_plt_dat %>%
  ggplot(aes(x=value,y=n,col=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(name ~ param,
             labeller = labeller(param=label_parsed,
                                 name=nlabs)) +
  ylab("sample size")+
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))

#ggsave(file.path(figdir,"sim_param_full2.png"),width=pltw,height=plth)



#sim_n100[[1]]["md.beta.est"]

## outcome censored at y=0
# sim_coeff.summary(sim_n25[[2]])
# sim_coeff.summary(sim_n50[[2]])
# sim_coeff.summary(sim_n100[[2]])
# sim_coeff.summary(sim_n200[[2]])
# sim_coeff.summary(sim_n400[[2]])

sim_coeff.a0.n25.cens <- conv(sim_n25[[2]]) %>% mutate(n=25,conc='1/J')
sim_coeff.a0.n50.cens <- conv(sim_n50[[2]]) %>% mutate(n=50,conc='1/J')
sim_coeff.a0.n100.cens <- conv(sim_n100[[2]]) %>% mutate(n=100,conc='1/J')
sim_coeff.a0.n200.cens <- conv(sim_n200[[2]]) %>% mutate(n=200,conc='1/J')
sim_coeff.a0.n400.cens <- conv(sim_n400[[2]]) %>% mutate(n=400,conc='1/J')

sim_coeff.a1.n25.cens <- conv(sim_a1_n25[[2]]) %>% mutate(n=25,conc='1/(0.8 + 0.35*J)')
sim_coeff.a1.n50.cens <- conv(sim_a1_n50[[2]]) %>% mutate(n=50,conc='1/(0.8 + 0.35*J)')
sim_coeff.a1.n100.cens <- conv(sim_a1_n100[[2]]) %>% mutate(n=100,conc='1/(0.8 + 0.35*J)')
sim_coeff.a1.n200.cens <- conv(sim_a1_n200[[2]]) %>% mutate(n=200,conc='1/(0.8 + 0.35*J)')
#sim_coeff.a1.n400.cens <- conv(sim_a1_n400[[2]]) %>% mutate(n=400,conc='1/(0.8 + 0.35*J)')

sim_coeff.a2.n25.cens <- conv(sim_a2_n25[[2]]) %>% mutate(n=25,conc='1/(2+(J/3))')
sim_coeff.a2.n50.cens <- conv(sim_a2_n50[[2]]) %>% mutate(n=50,conc='1/(2+(J/3))')
sim_coeff.a2.n100.cens <- conv(sim_a2_n100[[2]]) %>% mutate(n=100,conc='1/(2+(J/3))')
sim_coeff.a2.n200.cens <- conv(sim_a2_n200[[2]]) %>% mutate(n=200,conc='1/(2+(J/3))')
sim_coeff.a2.n400.cens <- conv(sim_a2_n400[[2]]) %>% mutate(n=400,conc='1/(2+(J/3))')



cens_dat <- bind_rows(sim_coeff.a0.n25.cens,
                      sim_coeff.a0.n50.cens,
                      sim_coeff.a0.n100.cens,
                      sim_coeff.a0.n200.cens,
                      sim_coeff.a0.n400.cens,
                      sim_coeff.a1.n25.cens,
                      sim_coeff.a1.n50.cens,
                      sim_coeff.a1.n100.cens,
                      sim_coeff.a1.n200.cens,
                      #sim_coeff.a1.n400.cens,
                      sim_coeff.a2.n25.cens,
                      sim_coeff.a2.n50.cens,
                      sim_coeff.a2.n100.cens,
                      sim_coeff.a2.n200.cens,
                      sim_coeff.a2.n400.cens
                      ) %>%
  mutate(n=factor(n,levels=c(400,200,100,50,25))) %>% filter(!true %in% c(-1.00,-0.33))


cens_plt_dat <- cens_dat %>% select(bias.est,bias.se,param,n,conc) %>%
  pivot_longer(cols=c(bias.est,bias.se))

cens_plt_dat %>%
  ggplot(aes(x=value,y=n,color=conc,shape=conc)) +
  geom_point(size=3,alpha=0.75)  +
  facet_grid(name ~ param,
             labeller = labeller(param=label_parsed,
                                 name=nlabs)) +
  ylab("sample size")+
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(size=fctsiz),
        axis.text =  element_text(size=atxtsz),
        strip.text = element_text(size=fctsiz),
        strip.text.y = element_text(angle=0))

#ggsave(file.path(figdir,"sim_param_cens2.png"),width=pltw,height=plth)
