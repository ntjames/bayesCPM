// functions to fit Bayesian CPM 

// define link function
real CDF_polr(real x, int link) {
  if (link == 1) return(inv_logit(x));
  else if (link == 2) return(Phi(x));
  else if (link == 3) return(gumbel_cdf(x, 0, 1));
  else if (link == 4) return(inv_cloglog(x));
  else if (link == 5) return(cauchy_cdf(x, 0, 1));
  else reject("Invalid link");
  return x;
}

//mixture link
//real CDF_mix(real x, ...)
//pis<-c(0.1,0.2,0.13,0.4,0.02,0.12,0.03)
//qlogis(cumsum(pis))
//qgumbel(cumsum(pis))

// make cutpoints 
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
// log-likelihood 
// same as eq (4) from Liu et al
vector loglik(int[] Ylev, int N, vector cutpoints, int ncat, vector eta, int link) {
  vector[N] ll;
  for (n in 1:N) {
    if (Ylev[n] == 1)
      ll[n] = log(CDF_polr(cutpoints[1] - eta[n], link));
    else if (Ylev[n] == ncat)
    // change to use log1m() ??
      ll[n] = log(1 - CDF_polr(cutpoints[ncat - 1] - eta[n], link));
    else
    // change to use log_diff_exp(log(CDF_polr(.)), log(CDF_polr(.))) ??
      ll[n] = log(CDF_polr(cutpoints[Ylev[n]] - eta[n], link) - CDF_polr(cutpoints[Ylev[n] - 1] - eta[n], link));
  }
  return ll;
}

// log-likelihood for hierarchical model
vector loglik_hier(int[] Ylev, int N, vector cutpoints, int ncat, vector eta, vector zeta, int[] jj, int link) {
  vector[N] ll;
  for (n in 1:N) {
    if (Ylev[n] == 1)
      ll[n] = log(CDF_polr(cutpoints[1] - eta[n] - zeta[jj[n]], link));
    else if (Ylev[n] == ncat)
    // change to use log1m()??
      ll[n] = log(1 - CDF_polr(cutpoints[ncat - 1] - eta[n] - zeta[jj[n]], link));
    else
    // change to use log_diff_exp(log(CDF_polr(.)), log(CDF_polr(.))) ??
      ll[n] = log(CDF_polr(cutpoints[Ylev[n]] - eta[n] - zeta[jj[n]], link) - CDF_polr(cutpoints[Ylev[n] - 1] - eta[n] - zeta[jj[n]], link));
  }
  return ll;
}
