data {
  real<lower=0> a;  // Prior parameter for theta
  real<lower=0> b;  
  int<lower=0> K; // number of historical datas
  vector[K] z0; // Observed historical successes
  vector[K] n0; // Observed historical trials
  real<lower=0> al;  // Prior parameters for eta
  real<lower=0> bl;
}

parameters {
  array[K] real<lower=0, upper=1> delta;
}

transformed parameters {
  vector[K] deltal;
  for (i in 1:K){
    deltal[i] = delta[i];
  }
}

model {
  for (i in 1:K) {
    target += beta_lpdf(delta[i] | al, bl ) -
                lbeta(a + z0[i] * deltal[i], 
                         b + (n0 - z0)[i] * deltal[i]);
  }
  target += lbeta(K*(a-1) + 1 + sum(z0.*deltal), K*(b-1) + 1 + sum((n0 - z0).*deltal));
}
