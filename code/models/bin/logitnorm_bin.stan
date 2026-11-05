// =============================================================================
// Bernoulli NPP with a compound-symmetry logistic-normal prior on eta
// (no order restriction).
//
//     psi[k]  = logit(eta[k])
//     psi     ~ Normal_K(mu * 1_K, Sigma),   Sigma = tau^2 * [(1-rho) I + rho * 11']
//
// Sigma is exchangeable (invariant to permuting the K historical datasets --
// no dataset is asserted to be more/less congruent a priori) but NOT
// diagonal: rho controls how strongly the eta[k]'s move together on the
// logit scale. rho = 0 recovers K independent logit-normal marginals; rho ->
// 1 pulls all eta[k]'s toward moving as one (near-total pooling). This is the
// direct hypercube analogue of the hierarchical-Beta prior in
// hier_beta_bin.stan; tau plays the same "how tight" role there did.
//
// We restrict to NON-NEGATIVE coupling (rho in (0,1)), which is the regime
// of interest here (shrinking the eta[k]'s toward each other, not away from
// each other). This lets us use a single shared latent factor z to build
// psi, which only reaches Sigma's non-negative-rho region:
//
//     psi[k] = mu + tau * ( sqrt(rho)*z + sqrt(1-rho)*e[k] ),  z, e[k] ~ N(0,1) iid
//
// (The full compound-symmetry PD region also allows small negative rho, but
// representing that requires a Cholesky factor of Sigma instead of this
// one-factor trick; skipped here since negative coupling isn't the goal.)
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

  // hyperpriors for the compound-symmetry logistic-normal prior on eta
  real mu_mean;             // mu ~ Normal(mu_mean, mu_sd), logit scale
  real<lower=0> mu_sd;
  real<lower=0> tau_sd;     // tau ~ Normal+(0, tau_sd) (half-normal)
  real<lower=0> rho_a;      // rho ~ Beta(rho_a, rho_b), coupling strength in (0,1)
  real<lower=0> rho_b;

  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

parameters {
  real mu;
  real<lower=0> tau;
  real<lower=0, upper=1> rho;
  real z;               // shared latent factor across all K datasets
  vector[K] e;           // idiosyncratic (per-dataset) latent factors
}

transformed parameters {
  vector[K] psi = mu + tau * (sqrt(rho) * z + sqrt(1 - rho) * e);
  vector<lower=0, upper=1>[K] eta = inv_logit(psi);
}

model {
  // ----- compound-symmetry logistic-normal prior on eta -----
  mu  ~ normal(mu_mean, mu_sd);
  tau ~ normal(0, tau_sd);   // constrained to (0, Inf) -> half-normal
  rho ~ beta(rho_a, rho_b);
  z   ~ std_normal();
  e   ~ std_normal();
  // no extra Jacobian term needed: eta is a deterministic transform of
  // (mu, tau, rho, z, e), which already carry their own sampling statements
  // above -- eta itself is never declared with a sampling statement.

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
