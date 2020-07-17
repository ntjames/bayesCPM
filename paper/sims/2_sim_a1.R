rm(list=ls())
#libs <- c("rstan", "dplyr", "stringr", "readr", "tidyr","pammtools","evd","sfsmisc")
libs <- c("rstan", "dplyr", "stringr", "readr", "tidyr","sfsmisc")
invisible(lapply(libs, library, character.only = TRUE)) 

sessionInfo()
head(Sys.cpuinfo(),19)

# set directory
dir <- file.path("~/bayes_cpm_paper1")

source(file.path(dir,"cpm_functions.r"))

# call this once to distribute MCMC chains across cpu cores:
options(mc.cores=parallel::detectCores())

# read in compiled model
ord_mod1 <- readRDS(file.path(dir,"ordinal_model_1.rds"))

# load sim array
simarray <- readRDS(file.path(dir,"bayes_cpm_simarray.rds"))

# select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

sim_pars <- simarray[s_id,]

nsamps <- as.numeric(sim_pars["nsamps"])
rep <- as.numeric(sim_pars["rep"])
seed <- as.numeric(sim_pars["seed"])

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
                            link=2,
                            conc=function(n) 1/(0.8 + 0.35*max(n, 3))) #probit link

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
                            link=2,
                            conc=function(n) 1/(0.8 + 0.35*max(n, 3))) #probit link

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

#p <- 0.5
#alpha <- 0
#beta <- c(1, -0.5)
#sigma <- 1
#log.y <- c(-1, -0.33, 0.5, 1.33, 2)

simout <- paste0("sim_a1_n",nsamps,"_",rep)
assign(simout, sim1_coeff.fun(sim=200, seed=seed, n=nsamps))

saveRDS(get(simout), file = file.path(dir,"out",paste0(simout,".rds")))
