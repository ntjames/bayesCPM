//modified from code sent by Philip G Jones <pgjones@saint-lukes.org>
// itself based on stan_polr code
// model with noninformative prior for alpha

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
   
  // make alpha constant parameter
  //alpha=1;
  
  cutpoints = make_cutpoints(pi, ncat, link);
  log_lik = loglik(Ylev, N, cutpoints, ncat, Q * b, link);
}

model {
  // priors
  
  //prior for counts (is it on right scale? transform to scale of cutpoints/intercepts??)
  // don't specify distribution for alpha so non-informative
  // repeat alpha for all params (i.e. symmetric Dirichlet)
  target += dirichlet_lpdf(pi | rep_vector(alpha, ncat));
  // equivalently
  // pi ~ dirichlet(rep_vector(alpha, ncat));
  
  //prior for betas
  // don't specify for noninformative
  //target += student_t_lpdf(b | 3, 0, 10);
  // equivalently
  // b ~ student_t(3, 0, 10);
  
  target += log_lik;
}

generated quantities {
  //vector[K] beta = R_inv * b;
}

