rm(list=ls())
#libs <- c("rstan", "dplyr", "stringr", "readr", "tidyr","pammtools","evd","sfsmisc")
libs <- c("rstan", "dplyr", "stringr", "readr", "tidyr","sfsmisc","evd")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

# set directory
dir <- file.path("~/bayes_cpm_paper1")

#! TEMP
#dir<-file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper/sims")

source(file.path(dir,"cpm_functions.r"))

# call this once to distribute MCMC chains across cpu cores:
options(mc.cores=parallel::detectCores())

# read in compiled model
ord_mod1 <- readRDS(file.path(dir,"ordinal_model_1.rds"))

#! TEMP
#ord_mod1 <- readRDS(file.path(dir,"ordinal_model_1_local.rds"))

# load sim array
simarray <- readRDS(file.path(dir,"bayes_cpm_simarray2.rds"))

# select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

sim_pars <- simarray[s_id,]

nsamps <- as.numeric(sim_pars["nsamps"])
rep <- as.numeric(sim_pars["rep"])
seed <- as.numeric(sim_pars["seed"])

# from Liu et al. sim code


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

if(0){#TEMP test
data <- generate.data.3(seed=12345, n=50, p=0.5)

ycens <- as.numeric(levels(data$y)[as.numeric(data$y)])
ycens[ycens < 0] <- 0
data$y.cens<-ordered(ycens)

## fit full outcome model
mod_data <- mkStanDat(data, outcome="y",
                      preds = c("z1", "z2"),
                      link=3, #loglog link
                      conc=function(n) 1/n)

cpm_fit <- sampling(ord_mod1, data=mod_data, seed=12345,
                    iter=4000, warmup=2000, chains=2, refresh=2000,
                    control = list(adapt_delta = 0.9))

# beta
beta.est <- summary(cpm_fit, pars=c("b[1]","b[2]"), use_cache=FALSE)$summary %>%
  as_tibble() %>% mutate(par=c("b[1]","b[2]"))

# gamma cutpoint estimates at given y values
cpout <- summary(cpm_fit, pars="cutpoints", use_cache=FALSE)$summary
gamma.y <- lapply(yvals, function(x) cutpoint_est(cpout,mod_data,x)) %>% do.call(rbind,.) %>%
  as_tibble() %>% mutate(par=c("gamma[y1]","gamma[y2]","gamma[y3]","gamma[y4]","gamma[y5]"))

yvals <- c(-0.3, 0, 0.5, 1.5, 2.5)

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

qgumbel(c(0.1, 0.4, 0.5, 0.6, 0.9),1*1-0.5*1,1)  
curve(dgumbel(x,1*1-0.5*1,1),-2,5)




# conditional cdf at y values
fit_cdf <- try(getCDF(cpm_fit, mod_data, newdata=zvals))
cond.cdf <- try(lapply(yvals,function(x) cdf_val_est(fit_cdf, mod_data,x)) %>% do.call(rbind,.) %>% as_tibble())


fit_cdf %>% filter(ndrow %in% c(1,2)) %>% ggplot(aes(group=ndrow)) +
  geom_stepribbon(aes(x=yval, ymin=cdf_q5, ymax=cdf_q95, 
                      fill=factor(ndrow)) , alpha=0.4) +
  geom_step(aes(x=yval, y=med_cdf,color=factor(ndrow))) +
  xlab("y") + ylab("Conditional CDF") + 
  scale_fill_discrete(name = "covariate value", 
                      labels=c("z1=1,z2=1", "z1=1,z2=0")) + 
  scale_color_discrete(name = "covariate value", 
                       labels=c("z1=1,z2=1", "z1=1,z2=0")) +
  stat_function(fun=function(x) pgumbel(x,0.5),color="red",
                linetype=2, alpha=0.4) + 
  stat_function(fun=function(x) pgumbel(x,1), color="blue",
                linetype=2, alpha=0.4)


# conditional mean (& contrast)
cond.mn <- try(getMean(cpm_fit, mod_data, newdata=zvals))
# cond.mn.cnt[[i]] <- try(getMeanContrast(cpm_fit, mod_data, cont.df=zvals))

#true mean
#loc + scale*em where em is Euler-Mascheroni constant
em=-digamma(1)
1*1-0.5*1+1*em
1*1-0.5*0+1*em


# conditional median (& contrast)
cond.med <- try(getQuantile(cpm_fit, mod_data, newdata=zvals, q=0.5))
#  cond.med.cnt[[i]] <- try(getQuantileContrast(cpm_fit, mod_data, cont.df=zvals, q=0.5))

#true median
#loc - scale*log(log(2))
1*1-0.5*1 - 1*log(log(2))
1*1-0.5*0 - 1*log(log(2))

# conditional q20 (& contrast)
cond.q20 <- try(getQuantile(cpm_fit, mod_data, newdata=zvals, q=0.2))
# cond.q20.cnt[[i]] <- try(getQuantileContrast(cpm_fit, mod_data, cont.df=zvals, q=0.2))

# true q20
qgumbel(0.2,1*1-0.5*1)
qgumbel(0.2,1*1-0.5*0)

# conditional q80 (& contrast)
cond.q80 <- try(getQuantile(cpm_fit, mod_data, newdata=zvals, q=0.8))
# cond.q80.cnt[[i]] <- try(getQuantileContrast(cpm_fit, mod_data, cont.df=zvals, q=0.8))




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

#p <- 0.5
#alpha <- 0
#beta <- c(1, -0.5)
#scale <- 1/3
#yvals <- c(-0.3, 0, 0.5, 1.5, 2.5)

#! TEMP
# foo<-sim3_coeff.fun(sim=4, seed=12345, n=25)
# foo2<-rbind(foo[[1]])
# foo3<-lapply(foo2, bind_rows)

simout <- paste0("sim_d0_n",nsamps,"_",rep)
assign(simout, sim3_coeff.fun(sim=100, seed=seed, n=nsamps))

saveRDS(get(simout), file = file.path(dir,"out",paste0(simout,".rds")))
