functions {
  real log_Z_beta(real alpha, real beta) {
    return lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
  }
}

data {
  int<lower=1> K;
  array[K] int<lower=0> n0;
  array[K] int<lower=0> s0;   // success counts per historical dataset
  int<lower=0> n;
  int<lower=0> s;             // success count in current data
  real<lower=0> a;        // prior Beta shape 1
  real<lower=0> b;     // prior Beta shape 2
  real<lower=0> al;
  real<lower=0> bl;
  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

parameters {
  vector<lower=0, upper=1>[K] eta;
}

model {
  // ----- prior on eta -----
  for (k in 1:K)
    target += beta_lpdf(eta[k] | al, bl);

  if (post == 1) {

    if (seq == 0) {
      // -------------------------------------------------------
      // NPP
      // -------------------------------------------------------

      real alpha_eta = a;
      real beta_eta  = b;

      for (k in 1:K) {
        alpha_eta += eta[k] * s0[k];
        beta_eta  += eta[k] * (n0[k] - s0[k]);
      }

      // log Z(alpha_eta, beta_eta) - log Z(alpha, beta)  [= log c0(eta)]
      // then log Z(tilde_alpha, tilde_beta) - log Z(alpha_eta, beta_eta)  [= log BF for D]
      // net contribution: log Z(tilde) - log Z(prior)
      real tilde_alpha = alpha_eta + s;
      real tilde_beta  = beta_eta  + (n - s);

      target += log_Z_beta(tilde_alpha, tilde_beta)
              - log_Z_beta(alpha_eta, beta_eta);

    } else {
      // -------------------------------------------------------
      // NPP-SEQ
      // -------------------------------------------------------

      real alpha_seq  = -(K - 1);     // correction: subtract (K-1), not (K-1)*alpha
      real beta_seq   = -(K - 1);
      real logZ_k_sum = 0;

      for (k in 1:K) {
        real alphak = a + eta[k] * s0[k];
        real betak  = b + eta[k] * (n0[k] - s0[k]);

        logZ_k_sum += log_Z_beta(alphak, betak);
        alpha_seq  += alphak;
        beta_seq   += betak;
      }

      // alpha_seq = sum_k alpha_k* - (K-1)
      // beta_seq  = sum_k beta_k*  - (K-1)

      real tilde_alpha = alpha_seq + s;
      real tilde_beta  = beta_seq  + (n - s);

      target += log_Z_beta(tilde_alpha, tilde_beta)
              - logZ_k_sum;
    }
  } else {
    if (seq == 1) {
      // -------------------------------------------------------
      // NPP-SEQ
      // -------------------------------------------------------

      real alpha_seq  = -(K - 1);     // correction: subtract (K-1), not (K-1)*alpha
      real beta_seq   = -(K - 1);
      real logZ_k_sum = 0;

      for (k in 1:K) {
        real alphak = a + eta[k] * s0[k];
        real betak  = b + eta[k] * (n0[k] - s0[k]);

        logZ_k_sum += log_Z_beta(alphak, betak);
        alpha_seq  += alphak;
        beta_seq   += betak;
      }

      // alpha_seq = sum_k alpha_k* - (K-1)
      // beta_seq  = sum_k beta_k*  - (K-1)
      real logZ_seq = log_Z_beta(alpha_seq, beta_seq);

      target += logZ_seq - logZ_k_sum;
    }
  }
}

generated quantities {
  real tilde_alpha;
  real tilde_beta;

  {
    if (post == 1) {
      if (seq == 0) {

        real alpha_eta = a;
        real beta_eta  = b;

        for (k in 1:K) {
          alpha_eta += eta[k] * s0[k];
          beta_eta  += eta[k] * (n0[k] - s0[k]);
        }

        tilde_alpha = alpha_eta + s;
        tilde_beta  = beta_eta  + (n - s);

      } else {

        real alpha_seq = -(K - 1);
        real beta_seq  = -(K - 1);

        for (k in 1:K) {
          alpha_seq += a + eta[k] * s0[k];
          beta_seq  += b + eta[k] * (n0[k] - s0[k]);
        }

        tilde_alpha = alpha_seq + s;
        tilde_beta  = beta_seq  + (n - s);
      }
    } else {
      if (seq == 1) {
        tilde_alpha = K*(a-1) + 1;
        tilde_beta  = K*(b-1) + 1;
      } else {
        tilde_alpha = a;
        tilde_beta  = b;
      }
    }
  }

  // posterior draw: theta | eta, D, D0 ~ Beta(tilde_alpha, tilde_beta)
  real theta = beta_rng(tilde_alpha, tilde_beta);
}
