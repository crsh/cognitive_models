# Generalized Context Model
data {
  //
  real rho;
  int p;
  int ntests;
  int nmemory;
  int ndim;
  int ntrials[ntests];
  int y[ntests];
  matrix[ntests, ndim] tests;
  matrix[nmemory, ndim] memory;
}

transformed data {
  vector[ndim] alpha;

  for(i in 1:ndim) {
    alpha[i] = 1;
  }
}


parameters {
  real<lower=0> c;
  real w_phi[ndim];
  real<lower=0> k;
}

transformed parameters {
  vector<lower=0,upper=1>[ntests] r;
  vector<lower=0>[ndim] exp_w_phi;
  vector<lower=0,upper=1>[ndim] w;

  for(i in 1:ndim) { # see http://andrewgelman.com/2009/04/29/conjugate_prior/
    exp_w_phi[i] = exp(w_phi[i]);
  }
  for(i in 1:ndim) {
    w[i] = exp_w_phi[i] / sum(exp_w_phi);
  }

  # Decision Probabilities
  for (i in 1:ntests) {
    vector[nmemory] s;
    real f;

    for (j in 1:nmemory) {
      vector[ndim] d;

      // Similarities
      for(l in 1:ndim) {
        d[l] = w[l] * fabs(tests[i, l] - memory[j, l])^rho;
      }
      s[j] = exp(-c * ((sum(d) + 0.000001)^(1/rho))^p);
    }

    f = sum(s);
    r[i] = f / (f + k);
  }
}


model {
  # Priors
  c ~ uniform(0, 10);
  k ~ uniform(0, 5);

  for(i in 1:ndim) { # see http://andrewgelman.com/2009/04/29/conjugate_prior/
    w_phi[i] ~ normal(0, 1);
  }

  # Decision Data
  y ~ binomial(ntrials, r);
}

generated quantities {
  int pred_y[ntests];

  for (i in 1:ntests) {
    pred_y[i] = binomial_rng(ntrials[i], r[i]);
  }
}
