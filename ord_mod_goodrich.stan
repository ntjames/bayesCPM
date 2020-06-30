// code from benjamin.goodrich@columbia.edu
// posted at https://discourse.datamethods.org/t/fully-sequential-bayesian-clinical-trial-design-for-in-hospital-treatment-of-covid-19/3129/6
// Based on section 13.3.1 of http://hbiostat.org/doc/rms.pdf
functions { // can check these match rms::orm via rstan::expose_stan_functions
  // pointwise log-likelihood contributions
  vector pw_log_lik(vector alpha, vector beta, row_vector[] X, int[] y) {
    int N = size(X);
    vector[N] out;
    int k = max(y); // assumes all possible categories are observed
    for (n in 1:N) {
      real eta = X[n] * beta;
      int y_n = y[n];
      if (y_n == 1) out[n] = log1m_exp(-log1p_exp(-alpha[1] - eta));
      else if (y_n == k) out[n] = -log1p_exp(-alpha[k - 1]  - eta);
      else out[n] = log_diff_exp(-log1p_exp(-alpha[y_n - 1] - eta),
                                 -log1p_exp(-alpha[y_n]     - eta));
    }
    return out;
  }

  // Pr(y == j)
  matrix Pr(vector alpha, vector beta, row_vector[] X, int[] y) {
    int N = size(X);
    int k = max(y); // assumes all possible categories are observed
    matrix[N, k] out;
    for (n in 1:N) {
      real eta = X[n] * beta;
      out[n, 1] = log1m_exp(-log1p_exp(-alpha[1] - eta));
      out[n, k] = -log1p_exp(-alpha[k - 1]  - eta);
      for (y_n in 2:(k - 1))
        out[n, y_n] = log_diff_exp(-log1p_exp(-alpha[y_n - 1] - eta),
                                   -log1p_exp(-alpha[y_n]     - eta));
    }
    return exp(out);
  }
}
data {
  int<lower = 1> N;   // number of observations
  int<lower = 1> p;   // number of predictors
  row_vector[p] X[N]; // columnwise CENTERED predictors
  int<lower = 3> k;   // number of outcome categories (7)
  int<lower = 1, upper = k> y[N]; // outcome on 1 ... k

  // prior standard deviations
  vector<lower = 0>[p] sds;
}
parameters {
  vector[p] beta; // coefficients
  simplex[k] pi;  // category probabilities for a person w/ average predictors
}
transformed parameters {
  vector[k - 1] alpha;                               // intercepts
  vector[N] log_lik;                                 // log-likelihood pieces
  for (j in 2:k) alpha[j - 1] = logit(sum(pi[j:k])); // predictors are CENTERED
  log_lik = pw_log_lik(alpha, beta, X, y);
}
model {
  target += log_lik;
  target += normal_lpdf(beta | 0, sds);
  // implicit: pi ~ dirichlet(ones)
}
generated quantities {
  vector[p] OR = exp(beta);
}
