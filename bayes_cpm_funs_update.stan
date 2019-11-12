// functions to fit Bayesian CPM 

/// functions for individual links

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

// log-likelihood; same as eq (4) from Liu et al
// possible to vectorize ??
vector loglik(int[] Ylev, int N, vector cutpoints, int ncat, vector eta, int link) {
  vector[N] ll;
  for (n in 1:N) {
    if (Ylev[n] == 1)
      ll[n] = log(CDF_polr(cutpoints[1] - eta[n], link));
    else if (Ylev[n] == ncat)
    // changed to use log1m() was ll[n] = log(1 - CDF_polr(cutpoints[ncat - 1] - eta[n], link));
    ll[n] = log1m(CDF_polr(cutpoints[ncat - 1] - eta[n], link));
    else
    // change to use log_diff_exp(log(CDF_polr(.)), log(CDF_polr(.)))??
      ll[n] = log(CDF_polr(cutpoints[Ylev[n]] - eta[n], link) - CDF_polr(cutpoints[Ylev[n] - 1] - eta[n], link));
  }
  return ll;
}

//// functions for mixture link from Lang 1999

//mixture link CDF 
real mix_cdf(real x, real lambda){
  real m1 = exp(-exp(3.5*lambda+2));
  real m3 = exp(-exp(-3.5*lambda+2));
  real m2 = 1 - m1 - m3;
return( m1*inv_cloglog(x) + m2*inv_logit(x) + m3*gumbel_cdf(x, 0, 1) );
}

// set-up algebra_system for solver to numerically find inverse CDF
vector mix_cdf_minp(vector q,     //unknown 
                    vector theta, // parameters
                    real[] x_r,     // data (real)
                    int[] x_i) {    // data (int) 
  vector[1] z; 
  z[1] = mix_cdf(q[1],theta[1]) - theta[2];
return( z ); 
}

vector make_cutpoints_mix(vector probabilities, int ncat, 
                      vector q_guess, real lambda, real[] x_r, int[] x_i) {
    vector[ncat - 1] cutpoints;
    vector[2] theta;
    vector[1] cp0;
    real running_sum = 0;
    for(i in 1:ncat - 1) {
      running_sum += probabilities[i];
      // these are defined in transformed data
      //vector[1] q_guess = [0]';
      //real x_r[0];
      //int x_i[0];
      //mixture link inverse CDF 
       theta = [lambda, running_sum]';
       // move this outside loop? evaluate on fine grid & reference within loop?
       cp0 = algebra_solver(mix_cdf_minp, q_guess, theta, x_r, x_i);
       cutpoints[i] = cp0[1];
    }
    return cutpoints;
  }

vector loglik_mix(int[] Ylev, int N, vector cutpoints, int ncat, vector eta,
                  real lambda) {
  vector[N] ll;
  for (n in 1:N) {
    if (Ylev[n] == 1)
      ll[n] = log(mix_cdf(cutpoints[1] - eta[n], lambda));
    else if (Ylev[n] == ncat)
      ll[n] = log1m(mix_cdf(cutpoints[ncat - 1] - eta[n], lambda));
    else
    // change to use log_diff_exp(log(mix_cdf(.)), log(mix_cdf(.)))??
      ll[n] = log(mix_cdf(cutpoints[Ylev[n]] - eta[n], lambda) - mix_cdf(cutpoints[Ylev[n] - 1] - eta[n], lambda));
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
