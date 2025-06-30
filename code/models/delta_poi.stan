data {
  real<lower=0, upper=1> a;  // Prior parameter for theta
  real<lower=0, upper=1> b;  
  int<lower=0> K; // number of historical datas
  int<lower=0> z;  // Observed current successes
  real<lower=0> t;
  vector[K] z0; // Observed historical successes
  vector[K] t0;
  real<lower=0, upper=1> al; // prior parameter for delta
  real<lower=0, upper=1> bl;
  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

parameters {
  array[K] real<lower=0, upper=1> delta;
}

transformed parameters {
  vector[K] deltal; // Cumulative sum
  
  for (i in 1:K) {
    deltal[i] = delta[i];
  }
}

model {
  for (i in 1:K){
    target += beta_lpdf(delta[i] | al, bl);  // Prior 
  }
  
  // Likelihood term 
  if (post == 1){
    if (seq == 1){
      target += lgamma(z + K*(a-1) + 1 + sum(z0 .* deltal)) - 
                (z + K*(a-1) + 1 + sum(z0 .* deltal)) * 
                log(K*b + t + sum(t0 .* deltal));
      for (i in 1:K) {
        target += -lgamma(a + z0[i] * deltal[i]);
        target += (a + z0[i] * deltal[i]) * log(b + t0[i] * deltal[i]);
      }
    } else{
      target += lgamma(z + a + sum(z0 .* deltal)) - 
                lgamma(a + sum(z0 .* deltal)) +
                (a + sum(z0 .* deltal))*log(b + sum(t0 .* deltal)) -
                (a + sum(z0 .* deltal) + z)*log(b + sum(t0 .* deltal) + t);
    }
  }
}