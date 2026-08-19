functions {
  // log Z for Normal up to constants: (1/2)*log(tau2) + (1/2)*mu^2/tau2
  // (the 2*pi constant cancels in all ratios)
  real log_Z_normal1d(real mu, real tau2) {
    return 0.5 * log(tau2) + 0.5 * mu^2 / tau2;
  }

  // Pooled NPP posterior: [mu_eta, tau_eta_inv2]
  vector npp_mu_tau(vector eta, array[] int n0, array[] real ybar0, real sigma2,
                     real mu0, real tau0_inv2, int K) {
    real tau_eta_inv2 = tau0_inv2;
    real rhs          = mu0 * tau0_inv2;

    for (k in 1:K) {
      tau_eta_inv2 += eta[k] * n0[k] / sigma2;
      rhs          += eta[k] * n0[k] * ybar0[k] / sigma2;
    }

    real tau_eta2 = 1.0 / tau_eta_inv2;
    real mu_eta   = tau_eta2 * rhs;
    return [mu_eta, tau_eta_inv2]';
  }

  // Pooled NPP-SEQ posterior: [mu_seq, tau_seq_inv2, logZ_k_sum]
  vector seq_mu_tau(vector eta, array[] int n0, array[] real ybar0, real sigma2,
                     real mu0, real tau0_inv2, int K) {
    real tau_seq_inv2 = 0;
    real rhs_seq      = 0;
    real logZ_k_sum   = 0;

    for (k in 1:K) {
      real tauk_inv2 = tau0_inv2 + eta[k] * n0[k] / sigma2;
      real tauk2     = 1.0 / tauk_inv2;
      real muk       = tauk2 * (mu0 * tau0_inv2 + eta[k] * n0[k] * ybar0[k] / sigma2);

      logZ_k_sum    += log_Z_normal1d(muk, tauk2);
      tau_seq_inv2  += tauk_inv2;
      rhs_seq       += muk / tauk2;
    }

    real tau_seq2 = 1.0 / tau_seq_inv2;
    real mu_seq   = tau_seq2 * rhs_seq;
    return [mu_seq, tau_seq_inv2, logZ_k_sum]';
  }
}

data {
  int<lower=1> K;
  array[K] int<lower=0> n0;
  array[K] real         ybar0;   // per-dataset sample means
  int<lower=0>  n;
  real          ybar;            // current data sample mean
  real<lower=0> sigma2;          // known variance
  real          mu0;             // prior mean
  real<lower=0> tau02;           // prior variance
  vector[K+1] alpha;
  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

transformed data {
  real tau0_inv2 = 1.0 / tau02;
}

parameters {
  array[K+1] real<lower=0> raw_gamma; // raw parameters for gamma
}

transformed parameters {
  simplex[K+1] gamma;
  for (i in 1:(K+1)){
    gamma[i] = raw_gamma[i] / sum(raw_gamma);
  }
  vector[K] eta = cumulative_sum(gamma[1:K]);
}

model {
  // prior for gamma

  target += gamma_lpdf(raw_gamma | alpha, 1);

  if (post == 1) {

    if (seq == 0) {
      // -------------------------------------------------------
      // NPP
      // -------------------------------------------------------

      vector[2] m = npp_mu_tau(eta, n0, ybar0, sigma2, mu0, tau0_inv2, K);
      real mu_eta       = m[1];
      real tau_eta_inv2 = m[2];
      real tau_eta2     = 1.0 / tau_eta_inv2;

      // update with D
      real tilde_tau_inv2 = tau_eta_inv2 + n / sigma2;
      real tilde_tau2     = 1.0 / tilde_tau_inv2;
      real tilde_mu       = tilde_tau2 * (tau_eta_inv2 * mu_eta + n * ybar / sigma2);

      // log BF for D
      target += log_Z_normal1d(tilde_mu, tilde_tau2)
              - log_Z_normal1d(mu_eta, tau_eta2);

    } else {
      // -------------------------------------------------------
      // NPP-SEQ
      // -------------------------------------------------------

      vector[3] m = seq_mu_tau(eta, n0, ybar0, sigma2, mu0, tau0_inv2, K);
      real mu_seq       = m[1];
      real tau_seq_inv2 = m[2];
      real logZ_k_sum   = m[3];

      // update with D
      real tilde_tau_inv2 = tau_seq_inv2 + n / sigma2;
      real tilde_tau2     = 1.0 / tilde_tau_inv2;
      real tilde_mu       = tilde_tau2 * (tau_seq_inv2 * mu_seq + n * ybar / sigma2);

      // log BF for D
      target += log_Z_normal1d(tilde_mu, tilde_tau2)
              - logZ_k_sum;
    }
  } else {
    if (seq == 1) {
      vector[3] m = seq_mu_tau(eta, n0, ybar0, sigma2, mu0, tau0_inv2, K);
      real mu_seq       = m[1];
      real tau_seq_inv2 = m[2];
      real logZ_k_sum   = m[3];
      real tau_seq2     = 1.0 / tau_seq_inv2;

      // log BF for D
      target += log_Z_normal1d(mu_seq, tau_seq2)
              - logZ_k_sum;
    }
  }
}

generated quantities {
  real tilde_mu;
  real tilde_tau2;

  {
    if (post == 1) {
      if (seq == 0) {

        vector[2] m = npp_mu_tau(eta, n0, ybar0, sigma2, mu0, tau0_inv2, K);
        real mu_eta       = m[1];
        real tau_eta_inv2 = m[2];

        real tilde_tau_inv2 = tau_eta_inv2 + n / sigma2;
        tilde_tau2 = 1.0 / tilde_tau_inv2;
        tilde_mu   = tilde_tau2 * (tau_eta_inv2 * mu_eta + n * ybar / sigma2);

      } else {

        vector[3] m = seq_mu_tau(eta, n0, ybar0, sigma2, mu0, tau0_inv2, K);
        real mu_seq       = m[1];
        real tau_seq_inv2 = m[2];

        real tilde_tau_inv2 = tau_seq_inv2 + n / sigma2;
        tilde_tau2 = 1.0 / tilde_tau_inv2;
        tilde_mu   = tilde_tau2 * (tau_seq_inv2 * mu_seq + n * ybar / sigma2);
      }
    } else {
      if (seq == 1) {
        tilde_mu = mu0;
        tilde_tau2 = tau02/K;
      } else {
        tilde_mu = mu0;
        tilde_tau2 = tau02;
      }
    }
  }

  // posterior draw: theta | eta, D, D0 ~ N(tilde_mu, tilde_tau2)
  real theta = normal_rng(tilde_mu, sqrt(tilde_tau2));
}
