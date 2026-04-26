data {
  real<lower=0, upper=1> a;  // Prior parameter for theta
  real<lower=0, upper=1> b;  
  int<lower=0> K; // number of historical datas
  int<lower=0> z;  // Observed current successes
  int<lower=0> n;  // Observed current trials
  vector[K] z0; // Observed historical successes
  vector[K] n0; // Observed historical trials
  vector[K+1] alpha; // prior parameter for gamma
  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

parameters {
  simplex[K+1] gamma; 
}

transformed parameters {
  vector[K] delta = cumulative_sum(gamma[1:K]);
  
  // delta[1] = gamma[1];
  // for (i in 2:K) {
  //   delta[i] = delta[i - 1] + gamma[i];
  // }
}

model {
  target += dirichlet_lpdf(gamma | alpha);  // Prior
  // gamma ~ dirichlet(alpha);
  
  // Likelihood term 
  if (post == 1){
    if (seq == 1){
      real alpha_post1 = z + K*(a-1) + 1 + sum(z0 .* delta);
      real beta_post1 = n - z + K*(b-1) + 1 + sum((n0 - z0) .* delta);
      target += lbeta(alpha_post1, beta_post1);
      for (i in 1:K) {
        real alpha_post2 = a + z0[i] * delta[i];
        real beta_post2 = b + (n0 - z0)[i] * delta[i];
        target += -lbeta(alpha_post2, beta_post2);
      }
    } else{
      real alpha_post1 = z + a + sum(z0 .* delta);
      real beta_post1 = n - z + b + sum( (n0 - z0) .* delta);
      target += lbeta(alpha_post1, beta_post1);
      real alpha_post2 = a + sum(z0 .* delta);
      real beta_post2 = b + sum( (n0 - z0) .* delta);
      target += -lbeta(alpha_post2, beta_post2);
    }
  }
}