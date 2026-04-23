data {
  real<lower=0> a;  // Prior parameter for theta
  real<lower=0> b;  
  int<lower=0> K; // number of historical datas
  vector[K] z0; // Observed historical successes
  vector[K] n0; // Observed historical trials
  vector[K+1] alpha; // prior parameter for gamma
  int<lower=0, upper=1> seq;
}

parameters {
  simplex[K+1] gamma; 
}

transformed parameters {
  vector[K] delta = cumulative_sum(gamma)[1:K];
}

model {
  target += dirichlet_lpdf(gamma | alpha);  // Prior
  if (seq == 1) {
    for (i in 1:K) {
      target += -lbeta(a + z0[i] * delta[i], 
                          b + (n0 - z0)[i] * delta[i]);
    }
    target += lbeta(K*(a-1) + 1 + sum(z0.*delta), K*(b-1) + 1 + sum((n0 - z0).*delta));
  }
}
