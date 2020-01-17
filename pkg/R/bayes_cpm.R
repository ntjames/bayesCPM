#' @title Bayes CPM model using Stan
#'
#' @description get posterior samples from Bayes CPM model using Stan
#' @param standata a list output from bayesCPM::mkStanDat() containing the data
#' @param ... other parameters passed to rstan::sampling()
#' @return a 'stanfit' object
#' @examples
#' dat1_stan <- mkStanDat(dat1, "outcome_var", c("pred1","pred2"), 2)
#' fit <- bayes_cpm(dat1_stan)

#' @export
bayes_cpm <- function(standata,...){
  out <- rstan::sampling(stanmodels$bayes_cpm_mod, data = standata, ...)
  return(out)
}
