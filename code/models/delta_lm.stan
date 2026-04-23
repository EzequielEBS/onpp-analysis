functions {
  matrix build_matrix_from_vector(vector X, int start, int nrow, int ncol) {
    int len = nrow * ncol;
    vector[len] x_sub = segment(X, start, len);
    return to_matrix(x_sub, nrow, ncol);
  }
}

data {
  int<lower=1> K; // number of historical data
  array[K] int<lower=0> n0; // number of observations of each hist data
  int<lower=0> n; // number of current observations 
  int<lower=0> p; // number of current covariates
  int<lower=0> len_X0;
  int<lower=0> len_y0;
  vector[len_X0] X0;
  array[K] int<lower=0> startid_X0; 
  vector[len_y0] y0;
  array[K] int<lower=0> startid_y0;
  matrix[n,p] X;
  vector[n] y;
  real<lower=0> a; // prior shape
  real<lower=0> b; // prior scale
  matrix[p,p] V0; // prior precision
  vector[p] mu0; // prior mean
  real<lower=0> tilde_a;
  real<lower=0> tilde_b;
  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

transformed data {
  vector[K] n0_vec = to_vector(n0);
  vector[p] hat_beta = mdivide_left_spd(  X' * X,  X' * y);
  real<lower=0> S = (y - X * hat_beta)' * (y - X * hat_beta);
  matrix[p, K] hat_beta0;
  vector<lower=0>[K] S0;
  matrix[p,K] prod_tX0k_y0k;
  for (k in 1:K) {
    matrix[n0[k], p] X0k = build_matrix_from_vector(X0, startid_X0[k], n0[k], p);
    vector[n0[k]] y0k = segment(y0, startid_y0[k], n0[k]);
    hat_beta0[,k] = mdivide_left_spd( X0k' * X0k, X0k' * y0k);
    S0[k] = (y0k - X0k * hat_beta0[,k])' * (y0k - X0k * hat_beta0[,k]);
    prod_tX0k_y0k[,k] = X0k' * y0k;
  }
}

parameters {
  vector<lower=0, upper=1>[K] delta;
}

transformed parameters {
  // Build default parameters
  real<lower=0> nu0 = n0_vec'* delta*0.5 + a;
  matrix[p,p] Lambda0 = V0;
  vector[p] tilde_beta00 = V0 * mu0;
  for (k in 1:K) {
    matrix[n0[k], p] X0k = build_matrix_from_vector(X0, startid_X0[k], n0[k], p);
    Lambda0 += X0k' * X0k * delta[k];
    tilde_beta00 += prod_tX0k_y0k[,k]*delta[k];
  }
  matrix[p,p] inv_Lambda0 = inverse_spd(Lambda0);
  vector[p] tilde_beta0 = inv_Lambda0 * tilde_beta00;
  real<lower=0> H0 = b;
  for (k in 1:K) {
    H0 += 0.5 * ( S0[k]* delta[k] + 
                  (hat_beta0[,k]' - tilde_beta0' ) * 
                  prod_tX0k_y0k[,k] * delta[k] );
  }
  H0 += 0.5 * mu0' * V0 * (mu0 - tilde_beta0);
  matrix[p,p] Lambda = Lambda0 + X' * X; 
  matrix[p,p] inv_Lambda = inverse_spd(Lambda);
  real<lower=0> H = H0 + 0.5 * (S + (tilde_beta0 - hat_beta)' *
           X' * X * inv_Lambda * Lambda0 * (tilde_beta0 - hat_beta));
           
  // Build sequential parameters
  real<lower=0> H0star = 0;
  vector<lower=0>[K] H0kstar;
  vector<lower=0>[K] nu0kstar;
  vector[p] sum_X0ky0kdeltak = rep_vector(0,p);
  matrix[p, p] sum_Lambda0kstar = rep_matrix(0, p, p);
  vector[K] log_det_Lambda0kstar;
  for (k in 1:K) {
    matrix[n0[k], p] X0k = build_matrix_from_vector(X0, startid_X0[k], n0[k], p);
    matrix[p,p] Lambda0kstar = V0 + X0k'*X0k*delta[k];
    matrix[p,p] inv_Lambda0kstar = inverse_spd(Lambda0kstar);
    sum_X0ky0kdeltak += prod_tX0k_y0k[,k]*delta[k];
    sum_Lambda0kstar += Lambda0kstar;
    H0kstar[k] = b + delta[k]*0.5 * (S0[k] + (mu0 - hat_beta0[,k])' *
                                            (X0k' * X0k) *
                                            inv_Lambda0kstar *
                                            V0 *
                                            (mu0 - hat_beta0[,k])
                                    );
    nu0kstar[k] = n0[k] * delta[k] / 2 + a;
    H0star += H0kstar[k] + 0.5 * (V0*mu0 + prod_tX0k_y0k[,k]*delta[k])' *
                                  inv_Lambda0kstar *
                                  (V0*mu0 + prod_tX0k_y0k[,k]*delta[k]);
    log_det_Lambda0kstar[k] = log_determinant_spd(Lambda0kstar);
  }
  sum_Lambda0kstar = 0.5 * (sum_Lambda0kstar + sum_Lambda0kstar');
  matrix[p,p] inv_sum_Lambda0kstar = inverse_spd(sum_Lambda0kstar);
  H0star += -0.5 * (K*V0*mu0 + sum_X0ky0kdeltak)' *
                      inv_sum_Lambda0kstar *
                      (K*V0*mu0 + sum_X0ky0kdeltak);
  vector[p] tilde_beta0star = inv_sum_Lambda0kstar *
                              (K*V0*mu0 + sum_X0ky0kdeltak);
  matrix[p,p] Lambdastar = sum_Lambda0kstar + X'*X;
  matrix[p,p] inv_Lambdastar = inverse_spd(Lambdastar);
  real<lower=0> Hstar = H0star + 0.5 *(S +
                              (tilde_beta0star - hat_beta)'*
                              (X'*X)*
                              inv_Lambdastar*
                              sum_Lambda0kstar*
                              (tilde_beta0star - hat_beta)
                             );
  real<lower=0> nustar = n0_vec'*delta*0.5 + n*0.5 + K*a;
}

model {
  // prior
  for (k in 1:K) {
    target += beta_lpdf(delta[k] | tilde_a, tilde_b);
  }
  
  
  // Likelihood term 
  if (post == 1){
    if (seq == 1){
      // Sequential method
      for (k in 1:K){
        target += nu0kstar[k] * log(H0kstar[k]) +
              0.5*log_det_Lambda0kstar[k] -
              lgamma(nu0kstar[k]);
      }
      target += -0.5*log_determinant(Lambdastar) -
                nustar * log(Hstar) +
                lgamma(nustar);
    } else{
      // posterior distribution
      target += lgamma(n*0.5 + nu0) + 0.5 * log_determinant_spd(Lambda0) + 
                nu0 * log(H0) - 0.5 * log_determinant_spd(Lambda) - lgamma(nu0) -
                (nu0 + n*0.5) * log(H);
    }
  }
}

generated quantities {
  vector[p] tilde_beta = rep_vector(0,p);
  real<lower=0> df;
  matrix[p,p] Sigma = rep_matrix(0, p, p);
  real<lower=0> shape;
  real<lower=0> scale;
  
  // Likelihood term 
  if (post == 1){
    if (seq == 1){
      tilde_beta = inv_Lambdastar *
                   (K*V0*mu0 + X'*y + sum_X0ky0kdeltak);
      df = 2*nustar;
      Sigma = Hstar / nustar * inv_Lambdastar;
      shape = nustar;
      scale = Hstar;
    } else{
      tilde_beta = inv_Lambda *
                   (X'*y + tilde_beta00);
      df = n + 2*nu0;
      Sigma = 2*H / df * inv_Lambda;
      shape = nu0 + n*0.5;
      scale = H;
    }
  }
  
  vector[p] beta = multi_student_t_rng(df, tilde_beta, Sigma);
  real<lower=0> sigma = inv_gamma_rng(shape, scale);
  vector[p+1] theta = append_row(beta, sigma);
}
