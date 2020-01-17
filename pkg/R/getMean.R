#' @title Get conditional mean from fit Bayes CPM model
#'
#' @description Calculate the posterior conditional mean from a Bayes CPM fit
#' @param fit 'stanfit' object
#' @param fitdata list containing data used to fit model
#' @param newdata a data frame with columns for each predictor used in the model
#' @param summ logical. Should the function return a summary of the posterior
#'    conditional mean? (default=TRUE)
#' @return posterior conditional mean summary (summ=TRUE) or values (summ=FALSE)
#' @examples
#' fit <- bayes_cpm(dat1_stan)
#' getMean(fit, dat1_stan, data.frame(pred1=c(0,1),pred2=c(1,1)))

#' @export
getMean <- function(fit, fitdata, newdata, summ=TRUE,...){

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
