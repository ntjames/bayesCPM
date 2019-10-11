//modified from code sent by Philip G Jones <pgjones@saint-lukes.org>
// itself based on stan_polr code
// concentration (alpha) is estimated with gamma(2,2) prior

functions {
#include /bayes_cpm_funs.stan
}

data {
  int N;
  int ncat;
  int Ylev[N];
  int link;
  int K;
  matrix[N, K] Q;
}

parameters {
  simplex[ncat] pi;
  vector[K] b;
  //estimate alpha parameter
  real<lower=0> alpha; 
}

transformed parameters {
  vector[ncat - 1] cutpoints;
  vector[N] log_lik;
  //vector[ncat] prior_counts;
  
  // repeat alpha for all prior counts (i.e. symmetric Dirichlet)
 // for(i in 1:ncat)
 //   prior_counts[i]=alpha;
   
  cutpoints = make_cutpoints(pi, ncat, link);
  log_lik = loglik(Ylev, N, cutpoints, ncat, Q * b, link);
}

model {
  // priors
  // gamma w/ mean 1, var 2/4=1/2
  // check sensitivity to this choice of hyperprior
  alpha ~ gamma(2,2);
  
  //prior for counts (is it on right scale? transform to scale of cutpoints/intercepts??)
  target += dirichlet_lpdf(pi | rep_vector(alpha, ncat));
  // equivalently
  // pi ~ dirichlet(rep_vector(alpha, ncat));
  
  //prior for betas
  // default weakly inform prior is student t w/ 3 df, centered at 0, and scale 10, check
  // comment out for noninformative, again check this
  //target += student_t_lpdf(b | 3, 0, 10);
  // equivalently
  // b ~ student_t(3, 0, 10);
  
  target += log_lik;
}

generated quantities {
  //vector[K] beta = R_inv * b;
}

