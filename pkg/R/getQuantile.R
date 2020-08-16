#' @title Get conditional quantile from fit Bayes CPM model
#'
#' @description Calculates the posterior conditional quantile
#' @param fit output list from bayes_cpm() containing 'stanfit' object and 'standata' data used to fit model
#' @param newdata a data frame with columns for each predictor used in the model
#' @param q quantile to be estimated (default is median, i.e. q=0.5)
#' @param summ logical. Should the function return a summary of the posterior
#'    conditional mean? (default=TRUE)
#' @return posterior conditional quantile summary (summ=TRUE) or values (summ=FALSE)
#' @examples
#' fit <- bayes_cpm(y~x1+x2, data=dat)
#' getQuantile(fit, data.frame(pred1=c(0,1),pred2=c(1,1)), q=0.5)

#' @export
getQuantile <- function(fit, newdata, q=0.5, summ=TRUE,...){
  
  cdf_samps <- getCDF(fit, newdata, summ=FALSE)
  
  #! pass these rather than recreating
  
  #! check truey is correct
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
                       qtile_q97.5=quantile(qtile,probs=0.975),.groups="keep") %>%
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


# getQuantile <- function(fit, fitdata, newdata, q=0.5, summ=TRUE,...){
# 
#   cdf_samps <- getCDF(fit, fitdata, newdata, summ=FALSE)
# 
#   #! pass these rather than recreating
# 
#   #! check truey is correct
#   # need to deal with quantiles below lowest obs and greater than highest obs
#   truey0 <- cdf_samps %>% filter(ndrow==1 & Var1=="A") %>% dplyr::select(yval) %>% pull()
#   truey <- c(truey0[-1],Inf)
# 
#   ndr <- newdata %>% mutate(ndrow=1:n())
#   cv_nms <- names(newdata)
# 
#   qtile_vals <- cdf_samps %>% group_by(ndrow, Var1) %>%
#     mutate( idx.1 = max(which(cdf<=q)), idx.2 = min(which(cdf>=q)),
#             cdf.1 = cdf[idx.1], cdf.2 = cdf[idx.2], n=1:n()) %>%
#     filter(n==1) %>%
#     mutate(idx.y1.cdf=ifelse(idx.1==-Inf,0,cdf.1),
#            idx.y2.cdf=ifelse(idx.2==-Inf,0,cdf.2),
#            idx.y1=ifelse(idx.1==-Inf,-Inf,truey[idx.1]),
#            idx.y2=ifelse(idx.2==-Inf,-Inf,truey[idx.2]),
#            qtile=ifelse(idx.1==idx.2,idx.y1,
#                         (idx.y2-idx.y1)/(idx.y2.cdf - idx.y1.cdf)*(q-idx.y1.cdf) + idx.y1)) %>%
#     ungroup()
# 
# 
#   if (summ){
#     qtile_summ <- qtile_vals %>%
#       group_by(ndrow) %>%
#       dplyr::summarize(mean_qtile=mean(qtile),
#                        med_qtile=median(qtile),
#                        sd_qtile=sd(qtile),
#                        qtile_q2.5=quantile(qtile,probs=0.025),
#                        qtile_q5=quantile(qtile,probs=0.05),
#                        qtile_q10=quantile(qtile,probs=0.10),
#                        qtile_q25=quantile(qtile,probs=0.25),
#                        qtile_q75=quantile(qtile,probs=0.75),
#                        qtile_q90=quantile(qtile,probs=0.90),
#                        qtile_q95=quantile(qtile,probs=0.95),
#                        qtile_q97.5=quantile(qtile,probs=0.975)) %>%
#       full_join(., ndr, by="ndrow")
#     return(qtile_summ)
#   } else {
#     #! append nd vars, keep other intermediary vars??
# 
#     qtile_vals_out<-qtile_vals %>% select(ndrow, cv_nms,
#                                           idx.y1.cdf, idx.y2.cdf,
#                                           idx.1, idx.2, qtile)
#     return(qtile_vals_out)
#   }
# 
# }
