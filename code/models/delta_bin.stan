data {
  real<lower=0, upper=1> a;  // Prior parameter for theta
  real<lower=0, upper=1> b;  
  int<lower=0> K; // number of historical datas
  int<lower=0> z;  // Observed current successes
  int<lower=0> n;  // Observed current trials
  vector[K] z0; // Observed historical successes
  vector[K] n0; // Observed historical trials
  real<lower=0, upper=1> al;  // Prior parameters for beta
  real<lower=0, upper=1> bl;
  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
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
    target += beta_lpdf(delta[i] | al, bl );
  }
  
  // Likelihood term 
  if (post == 1){
    if (seq == 1){
      target += lbeta(z + K*(a-1) + 1 + sum(z0 .* deltal), 
                      n - z + K*(b-1) + 1 + sum( (n0 - z0) .* deltal));
      for (i in 1:K) {
        target += -lbeta(a + z0[i] * deltal[i], 
                         b + (n0 - z0)[i] * deltal[i]);
      }
    } else{
      target += lbeta(z + a + sum(z0 .* deltal), n - z + b + sum( (n0 - z0) .* deltal)) - 
        lbeta(a + sum(z0 .* deltal), b + sum( (n0 - z0) .* deltal));
    }
  }
}