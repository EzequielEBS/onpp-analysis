data {
  real<lower=0, upper=1> a;  // Prior parameter for theta
  real<lower=0, upper=1> b;  
  int<lower=0> K; // number of historical datas
}

parameters {
  real<lower=0, upper=1> theta;
}

model {
  target += K * beta_lpdf(theta | a, b );
}
