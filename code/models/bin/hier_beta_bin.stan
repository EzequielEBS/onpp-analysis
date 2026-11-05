// =============================================================================
// Bernoulli NPP with a hierarchical Beta prior on eta (no order restriction).
//
// Instead of eta[1],...,eta[K] ~iid Beta(al, bl) (eta_bin.stan), each eta[k]
// is drawn from a Beta with a SHARED mean mu and SHARED concentration tau:
//
//     eta[k] | mu, tau ~ Beta(mu*tau, (1-mu)*tau),   k = 1,...,K
//     mu      ~ Beta(mu_a, mu_b)
//     log_tau ~ Normal(log_tau_mean, log_tau_sd),   tau = exp(log_tau)
//
// This prior is exchangeable across the K historical datasets (no ranking or
// order assumption -- permuting D0_1,...,D0_K leaves it unchanged), but it is
// NOT independent across them: integrating mu out induces positive
// correlation between the eta[k]'s, so evidence that one historical dataset
// should be discounted pulls the (posterior) shared mean mu down, which in
// turn pulls the OTHER eta[k]'s down too. tau controls how strong this
// pooling is (tau -> large: eta[k]'s clustered tightly around mu; tau ->
// small: eta[k]'s only weakly coupled through mu). This is the partial-
// pooling analogue of ONPP's order-restricted shrinkage, without asserting
// which historical dataset is more/less congruent.
//
// PARAMETERIZATION NOTE (why tau lives on the log scale here): putting tau
// directly behind a Gamma(shape, rate) prior -- the "obvious" choice for a
// positive concentration parameter -- creates a Neal's-funnel problem the
// moment shape <= 1: a Gamma(1, rate) (Exponential) prior's density is
// MAXIMAL exactly at tau = 0, so the sampler is pulled straight toward the
// one region where Beta(mu*tau, (1-mu)*tau) degenerates into a bimodal
// spike at {0, 1} (shape parameters -> 0) -- extremely high curvature, and
// exactly where the divergences show up. log_tau ~ Normal(...) instead puts
// ZERO density in the limit tau -> 0 (log_tau -> -Inf is infinitely far out
// in the Normal's tail), so the sampler is never pushed toward that
// pathological corner, while tau can still range as widely as before
// (log_tau_sd controls the spread on the log scale). This is the same
// "sample the log of a scale/concentration parameter" trick recommended
// throughout the Stan User's Guide's reparameterization chapter for
// variance/precision hyperparameters -- no Jacobian correction is needed
// here since we're simply choosing a different hyperprior shape for tau
// directly (via log_tau), not re-deriving an existing Gamma target.
//
// Everything below the prior block (functions, likelihood, generated
// quantities) is untouched from eta_bin.stan -- only the prior on eta
// changes.
// =============================================================================

functions {
  real log_Z_beta(real alpha, real beta) {
    return lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta);
  }

  // Pooled NPP parameters: [alpha_eta, beta_eta]
  vector npp_ab(vector eta, array[] int s0, array[] int n0, real a, real b, int K) {
    real alpha_eta = a;
    real beta_eta  = b;
    for (k in 1:K) {
      alpha_eta += eta[k] * s0[k];
      beta_eta  += eta[k] * (n0[k] - s0[k]);
    }
    return [alpha_eta, beta_eta]';
  }

  // Pooled NPP-SEQ parameters: [alpha_seq, beta_seq, logZ_k_sum]
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
  array[K] int<lower=0> s0;   // success counts per historical dataset
  int<lower=0> n;
  int<lower=0> s;             // success count in current data
  real<lower=0> a;            // prior Beta shape 1 for theta
  real<lower=0> b;            // prior Beta shape 2 for theta

  // hyperpriors for the hierarchical Beta prior on eta
  real<lower=0> mu_a;         // mu ~ Beta(mu_a, mu_b): shared mean of eta[k]'s
  real<lower=0> mu_b;
  real log_tau_mean;          // log_tau ~ Normal(log_tau_mean, log_tau_sd):
  real<lower=0> log_tau_sd;   // shared concentration, sampled on the log scale (see header)

  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

parameters {
  real<lower=0, upper=1> mu;
  real log_tau;
  vector<lower=0, upper=1>[K] eta;
}

transformed parameters {
  real<lower=0> tau = exp(log_tau);
}

model {
  // ----- hierarchical prior on eta -----
  mu      ~ beta(mu_a, mu_b);
  log_tau ~ normal(log_tau_mean, log_tau_sd);
  for (k in 1:K)
    target += beta_lpdf(eta[k] | mu * tau, (1 - mu) * tau);

  if (post == 1) {

    if (seq == 0) {
      // -------------------------------------------------------
      // NPP (hierarchical-Beta prior on eta)
      // -------------------------------------------------------
      vector[2] ab = npp_ab(eta, s0, n0, a, b, K);
      real tilde_alpha = ab[1] + s;
      real tilde_beta  = ab[2] + (n - s);

      target += log_Z_beta(tilde_alpha, tilde_beta)
              - log_Z_beta(ab[1], ab[2]);

    } else {
      // -------------------------------------------------------
      // NPP-SEQ (hierarchical-Beta prior on eta)
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
