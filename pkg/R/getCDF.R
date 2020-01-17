#' @title Get conditional CDF from fit Bayes CPM model
#'
#' @description Calculate the posterior conditional CDF from a Bayes CPM fit
#' @param fit 'stanfit' object
#' @param fitdata list containing data used to fit model
#' @param newdata a data frame with columns for each predictor used in the model
#' @param summ logical. Should the function return a summary of the posterior
#'    conditional CDF? (default=TRUE)
#' @return posterior conditional CDF summary (summ=TRUE) or values (summ=FALSE)
#' @examples
#' fit <- bayes_cpm(dat1_stan)
#' getCDF(fit, dat1_stan, data.frame(pred1=c(0,1),pred2=c(1,1)))

#' @export
getCDF <- function(fit, fitdata, newdata, summ=TRUE, ...){
  #! require(stringr) #is stringr used??

  #check that cumulative model used?

  #check for stanfit and newdata
  if( !identical(class(fit)[1],"stanfit")) stop("fit must be `stanfit` object")

  #! check newdata is a data.frame

  # check that names in newdata match coefs from model
  if( !identical(sort(names(fitdata$Q)),
                 sort(names(newdata))) ) stop("newdata vars must match model vars")

  #! other checks

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
  # dim(int) => s x (ints-1)
  # dim(Xb) => r x s

  # use inverse function based on family
  # 1 = logistic; 2 = probit; 3 = loglog; 4 = cloglog; 5 = cauchit !check cauchit
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
