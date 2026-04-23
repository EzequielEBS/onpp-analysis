data {
  real mu0;  // Prior parameter for theta
  real<lower=0> sigma0;  
  int<lower=0> K; // number of historical datas
  array[K] int n0; // Observed historical sample sizes
  array[K] real<lower=0> sigmah; // Observed historical standard deviations
  vector[sum(n0)] y0; // Observed historical data
  array[K] int start_idx; // Starting index for each historical data in y0
  real<lower=0> al;  // Prior parameters for eta
  real<lower=0> bl;
}

transformed data {
  vector[K] n0l;
  for (i in 1:K){
    n0l[i] = n0[i];
  }
  vector[K] sigmahl;
  for (i in 1:K){
    sigmahl[i] = sigmah[i];
  }
  vector[K] sum_y0;
  for (i in 1:K){
    sum_y0[i] = sum(y0[start_idx[i]:(start_idx[i]+n0[i]-1)]);
  }
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
    target += beta_lpdf(delta[i] | al, bl ) +
              // add constant
              0.5 * log(1/sigma0^2 + delta[i] * n0[i] / sigmah[i]) -
              (mu0/sigma0^2 + deltal[i] * sum_y0[i] / sigmahl[i])^2 / 
              (2 * (1/sigma0^2 + deltal[i] * n0[i] / sigmah[i]));
  }
  // add constant
  target += -0.5 * log(K/sigma0^2 + sum(deltal .* n0l ./ sigmahl));
  target += (K*mu0/sigma0^2 +  sum(deltal .* sum_y0 ./ sigmahl))^2 / 
            (2 * (K/sigma0^2 + sum(deltal .* n0l ./ sigmahl)));
}