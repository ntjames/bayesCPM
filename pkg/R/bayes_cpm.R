#' @title Bayes CPM model using Stan
#'
#' @description get posterior samples from Bayes CPM model using Stan
#' @param formula an object of class "formula": a symbolic description of the model to be fitted. outcome should be an ordered factor
#' @param data 	a data frame containing the variables in the model
#' @param link the link function to be used. must match one of c("logistic", "probit", "loglog", "cloglog", "cauchit")
#' @param dir_prior_conc concentration parameter for Dirichlet distribution as a function of the number of categories, n. The default is 1/n
#' @param prior NOT USED
#' @param ... other parameters passed to rstan::sampling()
#' @return a 'stanfit' object
#' @examples
#' fit <- bayes_cpm(y~x1+x2, data=dat1, link="logistic", dir_prior_conc=function(n) 1/n)

#' @export
bayes_cpm <- function(formula, data, link, dir_prior_conc = function(n) 1/n,
                          prior=NULL,...){
  mf <- model.frame(formula, data=data)
  mf.nm <- names(mf)
  lnk <- match.arg(link,c("logistic","probit","loglog","cloglog","cauchit"))
  lnkn <- switch(lnk, "logistic"=1, "probit"=2, "loglog"=3,
                 "cloglog"=4, "cauchit"=5)
  standata <- mkStanDat(data, outcome=mf.nm[1], preds = mf.nm[-1], link=lnkn,
                        conc=dir_prior_conc)
  out <- rstan::sampling(stanmodels$bayes_cpm_mod, data = standata, ...)
  return(out)
}
