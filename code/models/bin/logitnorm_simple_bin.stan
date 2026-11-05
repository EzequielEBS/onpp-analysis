// =============================================================================
// Bernoulli NPP with a SIMPLE (non-hierarchical) logistic-normal prior on
// eta -- the "just for comparison" baseline against logitnorm_bin.stan's
// compound-symmetry version.
//
//     psi[k] = logit(eta[k])
//     psi[k] ~iid Normal(mu0, sigma0),   k = 1,...,K
//
// mu0 and sigma0 are FIXED hyperparameters passed in as data (not given
// their own hyperpriors, and not shared/estimated across datasets) -- same
// role as eta_bin.stan's al/bl for its iid Beta(al, bl) prior. There is no
// mu/tau/rho/z/e machinery here: each eta[k] is an independent draw, with
// no coupling between the K historical datasets at all. This is to
// logitnorm_bin.stan what eta_bin.stan (iid Beta) is to hier_beta_bin.stan
// (hierarchical Beta) -- the un-pooled logit-normal counterpart, useful as
// a comparison point for how much the compound-symmetry version's rho
// coupling actually buys you.
//
// Everything below the prior block (functions, likelihood, generated
// quantities) is untouched from eta_bin.stan -- only the prior on eta
// changes.
// =============================================================================

functions {
  real log_Z_beta(real alpha, real beta) {
    return lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
  }

  vector npp_ab(vector eta, array[] int s0, array[] int n0, real a, real b, int K) {
    real alpha_eta = a;
    real beta_eta  = b;
    for (k in 1:K) {
      alpha_eta += eta[k] * s0[k];
      beta_eta  += eta[k] * (n0[k] - s0[k]);
    }
    return [alpha_eta, beta_eta]';
  }

  vector seq_ab(vector eta, array[] int s0, array[] int n0, real a, real b, int K) {
    real alpha_seq  = -(K - 1);
    real beta_seq   = -(K - 1);
    real logZ_k_sum = 0;

    for (k in 1:K) {
      real alphak = a + eta[k] * s0[k];
      real betak  = b + eta[k] * (n0[k] - s0[k]);

      logZ_k_sum += log_Z_beta(alphak, betak);
      alpha_seq  += alphak;
      beta_seq   += betak;
    }
    return [alpha_seq, beta_seq, logZ_k_sum]';
  }
}

data {
  int<lower=1> K;
  array[K] int<lower=0> n0;
  array[K] int<lower=0> s0;
  int<lower=0> n;
  int<lower=0> s;
  real<lower=0> a;
  real<lower=0> b;

  // fixed hyperparameters for the simple (non-hierarchical) logistic-normal
  // prior on eta -- same role as al/bl in eta_bin.stan's Beta(al, bl)
  real mu0;             // psi[k] ~ Normal(mu0, sigma0), logit scale
  real<lower=0> sigma0;

  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

parameters {
  vector[K] psi;
}

transformed parameters {
  vector<lower=0, upper=1>[K] eta = inv_logit(psi);
}

model {
  // ----- simple iid logistic-normal prior on eta -----
  psi ~ normal(mu0, sigma0);
  // no extra Jacobian term needed: eta is a deterministic transform of psi,
  // which already carries its own sampling statement above -- eta itself
  // is never declared with a sampling statement.

  if (post == 1) {

    if (seq == 0) {
      // -------------------------------------------------------
      // NPP (logistic-normal prior on eta)
      // -------------------------------------------------------
      vector[2] ab = npp_ab(eta, s0, n0, a, b, K);
      real tilde_alpha = ab[1] + s;
      real tilde_beta  = ab[2] + (n - s);

      target += log_Z_beta(tilde_alpha, tilde_beta)
              - log_Z_beta(ab[1], ab[2]);

    } else {
      // -------------------------------------------------------
      // NPP-SEQ (logistic-normal prior on eta)
      // -------------------------------------------------------
      vector[3] ab = seq_ab(eta, s0, n0, a, b, K);
      real tilde_alpha = ab[1] + s;
      real tilde_beta  = ab[2] + (n - s);

      target += log_Z_beta(tilde_alpha, tilde_beta)
              - ab[3];
    }
  } else {
    if (seq == 1) {
      vector[3] ab = seq_ab(eta, s0, n0, a, b, K);
      target += log_Z_beta(ab[1], ab[2]) - ab[3];
    }
  }
}

generated quantities {
  real tilde_alpha;
  real tilde_beta;

  {
    if (post == 1) {
      if (seq == 0) {
        vector[2] ab = npp_ab(eta, s0, n0, a, b, K);
        tilde_alpha = ab[1] + s;
        tilde_beta  = ab[2] + (n - s);
      } else {
        vector[3] ab = seq_ab(eta, s0, n0, a, b, K);
        tilde_alpha = ab[1] + s;
        tilde_beta  = ab[2] + (n - s);
      }
    } else {
      if (seq == 1) {
        tilde_alpha = K * (a - 1) + 1;
        tilde_beta  = K * (b - 1) + 1;
      } else {
        tilde_alpha = a;
        tilde_beta  = b;
      }
    }
  }

  // posterior draw: theta | eta, D, D0 ~ Beta(tilde_alpha, tilde_beta)
  real theta = beta_rng(tilde_alpha, tilde_beta);
}
