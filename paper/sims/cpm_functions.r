## Functions to make data and process output of Bayes CPM model fit

#! need to make these faster
# https://dtplyr.tidyverse.org/
# https://www.r-bloggers.com/big-data-wrangling-4-6m-rows-with-dtplyr-the-new-data-table-backend-for-dplyr/

#! add roxygen annotation

## -- format data for CPM models using Stan -- ##

mkStanDat<-function(ds, outcome, preds, link, conc=function(n) 1/n){
  require(dplyr)

  N <- nrow(ds)

  y <- ds %>% dplyr::select(outcome) %>% pull()

  if (!is.ordered(y))
    stop("outcome must be ordered factor")

  ncat <- length(unique(y))

  Ylev <- as.numeric(y)

  Q <- ds %>% dplyr::select(preds)

  K <- ncol(Q)

 # alpha <- 1/ncat
  alpha <- conc(ncat)

  ydat <- as.character(y) %>% as.numeric()

  truey0 <- levels(y) %>% as.numeric() %>% sort()

return( list(N=N, ncat=ncat, Ylev=Ylev, link=link, K=K, Q=Q, alpha=alpha,
             ydat=ydat, truey0=truey0) )
}

## -- estimate conditional CDF -- ##

getCDF <- function(fit, fitdata, newdata, summ=TRUE, ...){
  require(dplyr)
  require(stringr)
  require(Matrix)

  #check that cumulative model used?

  #check for brmfit and newdata
  if( !identical(class(fit)[1],"stanfit")) stop("fit must be `stanfit` object")

  #check newdata is a data.frame?

  # check that names in newdata match coefs from model
  if( !identical(sort(names(fitdata$Q)),
                 sort(names(newdata))) ) stop("newdata vars must match model vars")

  # other checks?

  # get data.frame with posterior samples
  posts <- as.data.frame(fit)

  nsamps <- nrow(posts)

  # get ordered values of outcome
  #!truey0 <- fitdata$truey0
  #!truey0levs <- fitdata$ncat

  # prepend value less than min(y) for alpha_0=-Inf intercept
  truey <- c(-Inf,fitdata$truey0)

  # format newdata, betas, and intercepts
  #!npreds <- fit@par_dims$b
  cv_nms <- paste0("b[",1:fit@par_dims$b,"]")
  cut_nms <- paste0("cutpoints[",1:(fitdata$ncat-1),"]")

  ndr <- newdata %>% dplyr::mutate(ndrow=1:n())
  nd <- ndr %>% dplyr::select(-ndrow) %>% as.matrix()

  beta <- posts %>% dplyr::select(cv_nms) %>% as.matrix()
  int <- posts %>% dplyr::select(cut_nms) %>% as.matrix()

  # get matrix of linear predictions Xb
  # (rxp) x (pxs) = rxs
  # r is rows in newdata, p is parameters (cols) in newdata,
  # s is number of MCMC samples
  # below is equivalent to Xb <- nd %*% t(beta), but faster
  Xb <- Matrix::tcrossprod(nd,beta)

  # add Xb to each intercept (4000xints)
  #dim(int) => s x (ints-1)
  #dim(Xb) => r x s

  # use inverse function based on family
  # 1 = logistic; 2 = probit; 3 = loglog; 4 = cloglog; 5 = cauchit !check cauchit is right
  inv_func <- switch(fitdata$link,
                     plogis,
                     pnorm,
                     function(y) exp(-exp(-y)),
                     function(y) 1-exp(-exp(y)),
                     pcauchy)

  #will have 1 for each row of nd
  #! check model/doc to make sure values are being calculated correctly
  #! are cutpoints y<= or y< ??

  nci <- ncol(int)
  nrnd <- nrow(nd)

  nd_lst <- vector("list", length=nrnd)
  for (i in 1:nrnd){
    #  F(y|X)=G^-1(alpha_i-betaX)
    tmpcdf0 <- int - t(Xb[rep(i,nci),, drop=FALSE])

    # add alpha_0=-Inf and alpha_n = Inf, convert to data.frame.table, get cdf
    nd_lst[[i]] <- cbind(`-Inf`=-Inf, tmpcdf0, `Inf`=Inf) %>%
      as.data.frame.table() %>%
      dplyr::mutate(cdf=inv_func(Freq), ndrow=i) %>%
      cbind(nd[i,,drop=FALSE])
  }


  # nd_ds<-vector("character",length=nrnd) # names of conditional cdf datasets
  # for (i in 1:nrnd){
  #   #  F(y_1|X)=G^-1(alpha_i-betaX)
  #   tmpcdf0 <- int - t(Xb[rep(i,nci),, drop=FALSE])
  #   # add alpha_0=-Inf and alpha_n = Inf, convert to data.frame.table, get cdf
  #   tmpcdf <- cbind(`-Inf`=-Inf, tmpcdf0, `Inf`=Inf) %>%
  #     as.data.frame.table() %>%
  #     dplyr::mutate(cdf=inv_func(Freq), ndrow=i) %>%
  #     cbind(nd[i,,drop=FALSE])
  #
  #   ccnm<-paste0("cc",i)
  #   nd_ds[i]<-ccnm
  #   assign(ccnm, tmpcdf)
  # }

  # combine conditional cdfs
  #! cdf_vals <- do.call(rbind, lapply(nd_ds, function(x) get(as.character(x))))
  #! cdf_vals <- bind_rows( lapply(nd_ds, function(x) get(as.character(x))) )
  cdf_vals <- bind_rows( nd_lst )

  if (summ){
    cdf_summ <- cdf_vals %>%
      ungroup() %>%
      group_by(ndrow, Var2) %>%
      dplyr::summarize(mn_cdf=mean(cdf),
                       med_cdf=median(cdf),
                       cdf_q2.5=quantile(cdf,probs=0.025),
                       cdf_q5=quantile(cdf,probs=0.05),
                       cdf_q10=quantile(cdf,probs=0.10),
                       cdf_q25=quantile(cdf,probs=0.25),
                       cdf_q75=quantile(cdf,probs=0.75),
                       cdf_q90=quantile(cdf,probs=0.90),
                       cdf_q95=quantile(cdf,probs=0.95),
                       cdf_q97.5=quantile(cdf,probs=0.975)) %>%
      ungroup() %>% dplyr::mutate(yval=rep(truey,nrow(nd))) %>%
      full_join(., ndr, by="ndrow")
    return(cdf_summ)
  } else {
    cdf_out <- cdf_vals %>%
      ungroup() %>%
      dplyr::arrange(ndrow, Var1) %>%
      dplyr::mutate(yval=rep(truey,nrnd*nsamps))
    return(cdf_out)
  }

}

# version using dtplyr for attempted speedup
getCDF_dt <- function(fit, fitdata, newdata, summ=TRUE, ...){
  require(dtplyr)
  require(dplyr)
#  require(stringr)
  require(Matrix)

  #check that cumulative model used?

  #check for brmfit and newdata
  if( !identical(class(fit)[1],"stanfit")) stop("fit must be `stanfit` object")

  #check newdata is a data.frame?

  # check that names in newdata match coefs from model
  if( !identical(sort(names(fitdata$Q)),
                 sort(names(newdata))) ) stop("newdata vars must match model vars")

  # other checks?

  # get data.frame with posterior samples
  posts <- as.data.frame(fit)
  nsamps <- nrow(posts)

  # prepend value less than min(y) for alpha_0=-Inf intercept to ordered outcome values
  truey <- c(-Inf, fitdata$truey0)

  # format newdata, betas, and intercepts
  cv_nms <- paste0("b[",1:fit@par_dims$b,"]")
  cut_nms <- paste0("cutpoints[",1:(fitdata$ncat-1),"]")

  ndr <- newdata %>% dplyr::mutate(ndrow=1:n())
  nd <- ndr %>% dplyr::select(-ndrow) %>% as.matrix()
  nrnd <- nrow(nd)

  beta <- posts %>% dplyr::select(cv_nms) %>% as.matrix()

  # get matrix of linear predictions Xb
  # (rxp) x (pxs) = rxs
  # r is rows in newdata, p is parameters (cols) in newdata,
  # s is number of MCMC samples
  # below is equivalent to Xb <- nd %*% t(beta), but faster
  Xb <- Matrix::tcrossprod(nd, beta)

  # add Xb to each intercept (4000xints)
  #dim(int) => s x (ints-1)
  #dim(Xb) => r x s

  # use inverse function based on family
  # 1 = logistic; 2 = probit; 3 = loglog; 4 = cloglog; 5 = cauchit !check this is right
  inv_func <- switch(fitdata$link,
                     plogis,
                     pnorm,
                     function(y) exp(-exp(-y)),
                     function(y) 1-exp(-exp(y)),
                     pcauchy)

  #will have 1 for each row of nd
  #! check model/doc to make sure values are being calculated correctly
  #! are cutpoints y<= or y< ??
  int <- posts %>% dplyr::select(cut_nms) %>% as.matrix()
  nci <- ncol(int)

  nd_lst <- vector("list", length=nrnd)
  for (i in 1:nrnd){
    #  F(y|X)=G^-1(alpha_i-betaX)
    tmpcdf0 <- int - t(Xb[rep(i,nci),, drop=FALSE])

    # add alpha_0=-Inf and alpha_n = Inf, convert to data.frame.table, get cdf
    nd_lst[[i]] <- cbind(`-Inf`=-Inf, tmpcdf0, `Inf`=Inf) %>%
      as.data.frame.table() %>%
      dplyr::mutate(cdf=inv_func(Freq), ndrow=i) %>%
      cbind(nd[i,,drop=FALSE])
  }

  # combine conditional cdfs
  cdf_vals <- bind_rows( nd_lst )

  cdf_vals2 <- lazy_dt(cdf_vals)

  if (summ){
    cdf_summ<-cdf_vals2 %>%
      ungroup() %>%
      group_by(ndrow, Var2) %>%
      summarize(mn_cdf=mean(cdf),
               med_cdf=median(cdf),
               cdf_q2.5=quantile(cdf,probs=0.025),
               cdf_q5=quantile(cdf,probs=0.05),
               cdf_q10=quantile(cdf,probs=0.10),
               cdf_q25=quantile(cdf,probs=0.25),
               cdf_q75=quantile(cdf,probs=0.75),
               cdf_q90=quantile(cdf,probs=0.90),
               cdf_q95=quantile(cdf,probs=0.95),
               cdf_q97.5=quantile(cdf,probs=0.975)) %>%
      ungroup() %>% mutate(yval=rep(truey,nrow(nd))) %>%
      full_join(., ndr, by="ndrow") %>% as_tibble()
    return(cdf_summ)
  } else {
    cdf_out <- cdf_vals2 %>%
      ungroup() %>%
      arrange(ndrow, Var1) %>%
      mutate(yval=rep(truey,nrnd*nsamps)) %>% as_tibble()
    return(cdf_out)
  }

}


## -- estimate conditional mean -- ##

getMean<-function(fit, fitdata, newdata, summ=TRUE,...){
  require(dplyr)

  cdf_samps <- getCDF(fit, fitdata, newdata, summ=FALSE)

  mn_vals <- cdf_samps %>% dplyr::filter(cdf!=0) %>% group_by(ndrow, Var1) %>%
    mutate(pdf0=lag(cdf), pdf=ifelse(is.na(pdf0),cdf,cdf-pdf0), fy_Py=pdf*yval) %>%
    dplyr::summarize(n=n(),mn=sum(fy_Py)) %>% ungroup()

  #! pass this from getCDF rather than recreating?
  ndr <- newdata %>% dplyr::mutate(ndrow=1:n())

  if (summ){
    mn_summ<-mn_vals %>%
      ungroup() %>%
      group_by(ndrow) %>%
      dplyr::summarize(mean_mn=mean(mn),
                       med_mn=median(mn),
                       sd_mn=sd(mn),
                       mn_q2.5=quantile(mn,probs=0.025),
                       mn_q5=quantile(mn,probs=0.05),
                       mn_q10=quantile(mn,probs=0.10),
                       mn_q25=quantile(mn,probs=0.25),
                       mn_q75=quantile(mn,probs=0.75),
                       mn_q90=quantile(mn,probs=0.90),
                       mn_q95=quantile(mn,probs=0.95),
                       mn_q97.5=quantile(mn,probs=0.975)) %>%
      full_join(., ndr, by="ndrow")
    return(mn_summ)
  } else {
    #append variable info
    mn_vals_app<-mn_vals %>% full_join(., ndr, by="ndrow")
    return(mn_vals_app)
  }

}


## -- estimate conditional quantiles -- ##

getQuantile<-function(fit, fitdata, newdata, q, summ=TRUE,...){
  require(dplyr)

  cdf_samps <- getCDF(fit, fitdata, newdata, summ=FALSE)

  #! pass these rather than recreating

  #! check if truey is correct
  # need to deal with quantiles below lowest obs and greater than highest obs
  truey0 <- cdf_samps %>% filter(ndrow==1 & Var1=="A") %>% dplyr::select(yval) %>% pull()
  truey <- c(truey0[-1],Inf)

  ndr <- newdata %>% mutate(ndrow=1:n())
  cv_nms <- names(newdata)

  qtile_vals <- cdf_samps %>% group_by(ndrow, Var1) %>%
    mutate( idx.1 = max(which(cdf<=q)), idx.2 = min(which(cdf>=q)),
            cdf.1 = cdf[idx.1], cdf.2 = cdf[idx.2], n=1:n()) %>%
    filter(n==1) %>%
    mutate(idx.y1.cdf=ifelse(idx.1==-Inf,0,cdf.1),
           idx.y2.cdf=ifelse(idx.2==-Inf,0,cdf.2),
           idx.y1=ifelse(idx.1==-Inf,-Inf,truey[idx.1]),
           idx.y2=ifelse(idx.2==-Inf,-Inf,truey[idx.2]),
           qtile=ifelse(idx.1==idx.2,idx.y1,
                        (idx.y2-idx.y1)/(idx.y2.cdf - idx.y1.cdf)*(q-idx.y1.cdf) + idx.y1)) %>%
    ungroup()


  if (summ){
    qtile_summ <- qtile_vals %>%
      group_by(ndrow) %>%
      dplyr::summarize(mean_qtile=mean(qtile),
                       med_qtile=median(qtile),
                       sd_qtile=sd(qtile),
                       qtile_q2.5=quantile(qtile,probs=0.025),
                       qtile_q5=quantile(qtile,probs=0.05),
                       qtile_q10=quantile(qtile,probs=0.10),
                       qtile_q25=quantile(qtile,probs=0.25),
                       qtile_q75=quantile(qtile,probs=0.75),
                       qtile_q90=quantile(qtile,probs=0.90),
                       qtile_q95=quantile(qtile,probs=0.95),
                       qtile_q97.5=quantile(qtile,probs=0.975)) %>%
      full_join(., ndr, by="ndrow")
    return(qtile_summ)
  } else {
    #! append nd vars, keep other intermediary vars??

    qtile_vals_out<-qtile_vals %>% select(ndrow, cv_nms,
                                          idx.y1.cdf, idx.y2.cdf,
                                          idx.1, idx.2, qtile)
    return(qtile_vals_out)
  }

}


## -- Exceedance Probability -- ##

if(0){
?rms::ExProb

set.seed(1)
x1 <- runif(100)
yvar <- x1 + runif(100)
f <- orm(yvar ~ x1)
d <- ExProb(f)
lp <- predict(f, newdata=data.frame(x1=c(.2,.8)))
w <- d(lp)

w

s1 <- abs(x1 - .2) < .1
s2 <- abs(x1 - .8) < .1
plot(w, data=data.frame(x1=c(rep(.2, sum(s1)), rep(.8, sum(s2))),
                        yvar=c(yvar[s1], yvar[s2])))
}
