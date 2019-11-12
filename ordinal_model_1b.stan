// Bayesian CPM model with dirichlet prior
// concentration (alpha) is given as a scalar parameter along with data
// modified in part from code by Philip G Jones <pgjones@saint-lukes.org>
// itself based on stan_polr() code
// includes mixture link

functions {
#include /bayes_cpm_funs_update.stan
}
data {
  int N; // number of observations
  int ncat; // number of unique outcome values
  int Ylev[N]; // ranks of unique outcome values
//  int link; // link function (1=logistic, 2=probit, ...)
  int K; // number of predictors
  matrix[N, K] Q; // N x K design matrix
  real<lower=0> alpha; // concentration parameter
}

transformed data{
// initial guess and empty array for mixture link inverse (algebra_solver)
vector[1] q_guess = [0]'; // check dims on this
real x_r[0];
int x_i[0];
}

parameters {
  simplex[ncat] pi;
  vector[K] b;
  real lambda;
}

transformed parameters {
  vector[ncat - 1] cutpoints;
  vector[N] log_lik;
  
  cutpoints = make_cutpoints_mix(pi, ncat, q_guess, lambda, x_r, x_i);
  log_lik = loglik_mix(Ylev, N, cutpoints, ncat, Q * b, lambda);
}

model {
  // prior for mixture dist. favoring logistic 
  lambda ~ normal(0,3);
  
  //prior for counts (is it on right scale? transform to scale of cutpoints/intercepts??)
  // repeat alpha for all params (i.e. symmetric Dirichlet)
  target += dirichlet_lpdf(pi | rep_vector(alpha, ncat));
  // equivalently
  // pi ~ dirichlet(rep_vector(alpha, ncat));
  
  //prior for betas
  //target += student_t_lpdf(b | 3, 0, 10);
  // equivalently
  // b ~ student_t(3, 0, 10);
  
  target += log_lik;
}

generated quantities {
  //vector[K] beta = R_inv * b;
}

