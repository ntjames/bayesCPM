functions {
  real CDF_polr(real x, int link) {
    if (link == 1) return(inv_logit(x));
    else if (link == 2) return(Phi(x));
    else if (link == 3) return(gumbel_cdf(x, 0, 1));
    else if (link == 4) return(inv_cloglog(x));
    else if (link == 5) return(cauchy_cdf(x, 0, 1));
    else reject("Invalid link");
    return x;
  }
  vector make_cutpoints(vector probabilities, int ncat, int link) {
    vector[ncat - 1] cutpoints;
    real running_sum = 0;
    if (link == 1) for(i in 1:ncat - 1) {
      running_sum += probabilities[i];
      cutpoints[i] = logit(running_sum);
    }
    else if (link == 2) for(i in 1:ncat - 1) {
      running_sum += probabilities[i];
      cutpoints[i] = inv_Phi(running_sum);
    }
    else if (link == 3) for(i in 1:ncat - 1) {
      running_sum += probabilities[i];
      cutpoints[i] = -log(-log(running_sum));
    }
    else if (link == 4) for(i in 1:ncat - 1) {
      running_sum += probabilities[i];
      cutpoints[i] = log(-log1m(running_sum));
    }
    else if (link == 5) for(i in 1:ncat - 1) {
      running_sum += probabilities[i];
      cutpoints[i] = tan(pi() * (running_sum - 0.5));
    }
    else reject("invalid link");
    return cutpoints;
  }
  vector loglik(int[] Ylev, int N, vector cutpoints, int ncat, vector eta, int link) {
    vector[N] ll;
    for (n in 1:N) {
      if (Ylev[n] == 1)
        ll[n] = log(CDF_polr(cutpoints[1] - eta[n], link));
      else if (Ylev[n] == ncat)
        ll[n] = log(1 - CDF_polr(cutpoints[ncat - 1] - eta[n], link));
      else
        ll[n] = log(CDF_polr(cutpoints[Ylev[n]] - eta[n], link) - CDF_polr(cutpoints[Ylev[n] - 1] - eta[n], link));
    }
    return ll;
  }
}

data {
  int N;
  int ncat;
  int Ylev[N];
  vector[ncat] Ysupp;
  int link;
  int K;
  matrix[N, K] Q;
  matrix[K, K] R_inv;
}

parameters {
  simplex[ncat] pi;
  vector[K] b;
}

transformed parameters {
  vector[ncat - 1] cutpoints;
  vector[N] log_lik;
  cutpoints = make_cutpoints(pi, ncat, link);
  log_lik = loglik(Ylev, N, cutpoints, ncat, Q * b, link);
}

model {
  target += student_t_lpdf(b | 3, 0, 10);
  target += log_lik;
}

generated quantities {
  vector[K] beta = R_inv * b;
}

