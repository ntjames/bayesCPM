#' @title Get posterior predictive values from fit Bayes CPM model
#'
#' @description Draw posterior predictive values from a Bayes CPM fit
#' @param fit output list from bayes_cpm() containing 'stanfit' object and 'standata' data used to fit model
#' @param newdata a data frame with columns for each predictor used in the model
#' @param draws number of draws to return
#' @return matrix with posterior predictive values
#' @examples
#' fit <- bayes_cpm(y~x1+x2, data=dat)
#' getPostPred(fit, dat[,c("x1","x2")])

#' @export
getPostPred<-function(fit,newdata,draws=5){

postdf <- as.data.frame(fit$stanfit)
drawsdf <- postdf[sample(nrow(postdf),draws),]
draws_cp <- drawsdf[,grep("cutpoints",names(drawsdf))]

draws_beta <- as.matrix(drawsdf[,grep("b",names(drawsdf)),drop=FALSE])

xx <- as.matrix(newdata)
lp <- draws_beta %*% t(xx)

J <- fit$standata$ncat
N <- fit$standata$N
y_tilde <- y_tilde_cat <- matrix(NA,nrow=N,ncol=draws)

# use inverse function based on family
# 1 = logistic; 2 = probit; 3 = loglog; 4 = cloglog; 5 = cauchit !check cauchit
inv_func <- switch(fit$standata$link,
                   plogis,
                   pnorm,
                   function(y) exp(-exp(-y)),
                   function(y) 1-exp(-exp(y)),
                   pcauchy)

for (n in 1:N){
 theta <- matrix(NA,nrow=J,ncol=draws)
# first <- retCDF(draws_cp[,1]-lp[,n],fit$standata$link) # 2nd arg is link
 first <- inv_func(draws_cp[,1]-lp[,n])
 theta[1,] <- previous <- first
  for (j in 2:(J-1)) {
  #  current <- retCDF(draws_cp[,j]-lp[,n],fit$standata$link) # 2nd arg is link
    current <- inv_func(draws_cp[,j]-lp[,n])
    theta[j,] <- current - previous
    previous <- current
  }
 theta[J,] = 1 - previous
  if (previous <= 0 || previous >= 1) {
    # do nothing
  } else {
    for(d in 1:draws){
      y_tilde_cat[n,d]<-sample.int(N,size=1,replace=TRUE,prob=theta[,d])
    }
  }
}

# associate y_tilde_cat with observed y value
 for(d in 1:draws){
    y_tilde[,d]<-fit$standata$truey0[y_tilde_cat[,d]]
 }

return(y_tilde)
}
