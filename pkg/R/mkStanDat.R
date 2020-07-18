#' @title Format data for CPM models using Stan
#'
#' @description This function formats data for use with rstan::stan() or rstan::sampling()
#' @param ds a data frame
#' @param outcome a character string containing the model outcome
#' @param preds vector of linear predictors for model
#' @param link the link function to be used (1 = logistic; 2 = probit;
#'    3 = loglog; 4 = cloglog; 5 = cauchit)
#' @return a list containing the data to be used for the model
#' @examples
#' dat1_stan <- mkStanDat(dat1, "outcome_var", c("pred1","pred2"), 2)

#' @export
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

  alpha <- conc(ncat)

  ydat <- as.character(y) %>% as.numeric()

  truey0 <- levels(y) %>% as.numeric() %>% sort()

return( list(N=N, ncat=ncat, Ylev=Ylev, link=link, K=K, Q=Q, alpha=alpha,
             ydat=ydat, truey0=truey0) )
}
