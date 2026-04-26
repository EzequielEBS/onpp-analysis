data {
  real mu0;  // Prior parameter for theta
  real<lower=0> sigma0;  
  int<lower=0> K; // number of historical datas
  array[K] int n0; // Observed historical sample sizes
  int<lower=0> n; // Current sample size
  array[K] real<lower=0> sigmah; // Observed historical standard deviations
  real<lower=0> sigma; // Observed current standard deviation
  vector[sum(n0)] y0; // Observed historical data
  vector[n] y; // Observed current data
  array[K] int start_idx; // Starting index for each historical data in y0
  vector[K+1] alpha;  // Prior parameters for gamma
  int<lower=0> seq; // 1 for sequential, 0 for non-sequential
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
  array[K+1] real<lower=0> raw_gamma; // raw parameters for gamma
  real theta; // parameter of interest
}

transformed parameters {
  simplex[K+1] gamma;
  for (i in 1:(K+1)){
    gamma[i] = raw_gamma[i] / sum(raw_gamma);
  }
  vector[K] delta = cumulative_sum(gamma[1:K]);
}

model {

  // prior for gamma

  target += gamma_lpdf(raw_gamma | alpha, 1);

  if (seq == 1) {
    // real v0;
    // real m0;
    // v0 = 0;
    // m0 = 0;
    for (i in 1:K) {
      real v0k;
      real m0k;
      v0k = (1/sigma0^2 + delta[i] * n0l[i] / sigmahl[i]^2)^(-1);
      m0k = v0k * (mu0/sigma0^2 + delta[i] * sum_y0[i] / sigmahl[i]^2);
      target += normal_lpdf(theta | m0k, sqrt(v0k));
      // v0 += (1/v0k);
      // m0 += m0k / v0k;
      // constant
      // target += -0.5 * log(2 * pi() * v0k);
      // target += -0.5 * m0k^2 / v0k;

    }
    // v0 = 1 / v0;
    // m0 = v0 * m0;
    // real v;
    // real m;
    // v = (1/v0 + n / sigma^2)^(-1);
    // m = v * (m0 / v0 + sum(y) / sigma^2);
    // target += normal_lpdf(theta | m, sqrt(v));
    // // constant
    // target += 0.5 * m0^2 / v0 + 0.5 * m^2 / v;
    // target += -0.5 * log(2 * pi() * v);
    target += normal_lpdf(y | theta, sigma);
  } else {
    // non-sequential model
    real v0;
    real m0;
    // real v;
    // real m;
    v0 = 1 / (1/sigma0^2 + sum(delta .* n0l ./ sigmahl^2));
    m0 = v0 * (mu0/sigma0^2 + sum(delta .* sum_y0 ./ sigmahl));
    target += normal_lpdf(theta | m0, sqrt(v0));
    target += normal_lpdf(y | theta, sigma);
    // v = (1/v0 + n / sigma^2)^(-1);
    // m = v * (m0 / v0 + sum(y) / sigma^2);
    // target += normal_lpdf(theta | m, sqrt(v));
    // target += -0.5 * log(v0);
    // target += 0.5 * m^2/v + 0.5 * log(v);
  }
}