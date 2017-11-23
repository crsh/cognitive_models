// // Generalized Context Model
data {
  real rho;
  int p;
  int ntests;
  int nmemory;
  int ndim;
  int ntrials[ntests];
  int category[nmemory];
  int y[ntests];
  matrix[ntests,2] tests;
  matrix[nmemory,2] memory;
}

parameters {
  real<lower=0> c;
  real<lower=0,upper=1> w;
  real<lower=0,upper=1> b;
}

transformed parameters {
  vector<lower=0,upper=1>[ntests] r;
  vector<lower=0,upper=1>[ndim] wk;
  real tmp1[ntests,nmemory,2];
  real tmp2[ntests,nmemory,2];
  real numerator[ntests,nmemory];
  real denominator[ntests,nmemory];

  wk[1] = w;
  wk[2] = 1 - w;

  for (i in 1:ntests) {
    for (j in 1:nmemory) {
      real s;
      vector<lower=0>[ndim] d;

      // Similarities
      for(k in 1:ndim) {
        d[k] = wk[k] * fabs(tests[i, k] - memory[j, k])^rho;
      }
      s = exp(-c * ((sum(d) + 0.000001)^(1/rho))^p);

      // Decision Probabilities
      tmp1[i,j,1] = b * s;
      tmp1[i,j,2] = 0;
      tmp2[i,j,1] = 0;
      tmp2[i,j,2] = (1 - b) * s;

      numerator[i, j] = tmp1[i,j,category[j]];
      denominator[i, j] = tmp1[i,j,category[j]] + tmp2[i,j,category[j]];
    }
    r[i] = sum(numerator[i, ]) / sum(denominator[i, ]);
  }
}

model {
  // Prior
  c ~ uniform(0, 5);
  w ~ beta(1, 1);
  b ~ beta(1, 1);

  // Decision Data
  y ~ binomial(ntrials, r);
}

generated quantities {
  vector[ntests] pred_y;

  for (i in 1:ntests) {
    pred_y[i] = binomial_rng(ntrials[i], r[i]);
  }
}

