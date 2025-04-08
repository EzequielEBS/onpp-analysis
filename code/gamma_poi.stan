data {
  real<lower=0, upper=1> a;  // Prior parameter for theta
  real<lower=0, upper=1> b;  
  int<lower=0> K; // number of historical datas
  int<lower=0> z;  // Observed current successes
  real<lower=0> t;
  vector[K] z0; // Observed historical successes
  vector[K] t0;
  vector[K+1] alpha; // prior parameter for gamma
  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

parameters {
  simplex[K+1] gamma; 
}

transformed parameters {
  array[K] real<lower=0, upper=1> gamma_;
  
  vector[K] delta; // Cumulative sum
  
  for (i in 1:K) {
    gamma_[i] = gamma[i];
    delta[i] = sum(gamma_[1:i]);
  }
}

model {
  target += dirichlet_lpdf(gamma | alpha);  // Prior
  
  // Likelihood term 
  if (post == 1){
    if (seq == 1){
      target += lgamma(z + K*(a-1) + 1 + sum(z0 .* delta)) - 
                (z + K*(a-1) + 1 + sum(z0 .* delta)) * 
                log(K*b + t + sum(t0 .* delta));
      for (i in 1:K) {
        target += -lgamma(a + z0[i] * delta[i]);
        target += (a + z0[i] * delta[i]) * log(b + t0[i] * delta[i]);
      }
    } else{
      target += lgamma(z + a + sum(z0 .* delta)) - 
                lgamma(a + sum(z0 .* delta)) +
                (a + sum(z0 .* delta))*log(b + sum(t0 .* delta)) -
                (a + sum(z0 .* delta) + z)*log(b + sum(t0 .* delta) + t);
    }
  }
}