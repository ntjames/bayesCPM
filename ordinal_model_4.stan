// Bayesian CPM Stan code
// add random intercept

functions {
#include /bayes_cpm_funs.stan
}
data {
  int N; // num individuals
  int ncat; // num categories
  int Ylev[N]; // ordered individual outcome
  int link; // link function
  int K; // num ind predictors
  matrix[N, K] Q; // individual predictors
  int J; // num groups
  int L; //num group predictors
  int<lower=1, upper=J> jj[N]; //group for individual
  matrix[J, L] u; // group predictors
  real<lower=0> alpha;
}
parameters {
  simplex[ncat] pi;
  vector[K] bet; // indiv coeffs by group
  vector[L] b; // group coeffs
}
transformed parameters {
  vector[ncat - 1] cutpoints;
  vector[N] log_lik;
  //row_vector[K] u_b[J];
  vector[J] u_b;
  
  cutpoints = make_cutpoints(pi, ncat, link);
  
  for (j in 1:J)
    u_b[j]=u[j,]*b;

 //   log_lik = loglik(Ylev, N, cutpoints, ncat, Q * bet[jj[n]], u_b[jj[n]], link);
  log_lik = loglik_hier(Ylev, N, cutpoints, ncat, Q * bet, u_b, jj, link);
}

model {
  // priors
  //prior for counts (is it on right scale? transform to scale of cutpoints/intercepts??)
  // repeat alpha for all params (i.e. symmetric Dirichlet)
  target += dirichlet_lpdf(pi | rep_vector(alpha, ncat));
  
  // prior for ind coeffs
  target += student_t_lpdf(bet | 3, 0, 10); 
  
  // prior for group coeffs
  target += student_t_lpdf(b | 3, 0, 10);
  
  target += log_lik;
}

generated quantities {
  //vector[K] beta = R_inv * b;
}

