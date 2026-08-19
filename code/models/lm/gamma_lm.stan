functions {
  matrix build_matrix_from_vector(vector X, int start, int nrow, int ncol) {
    int len = nrow * ncol;
    vector[len] x_sub = segment(X, start, len);
    return to_matrix(x_sub, nrow, ncol);
  }

  // log Z(mu, V, a, b) for NIG normalizing constant
  real log_Z_NIG(vector mu, matrix V, real a, real b, int p) {
    return 0.5 * p * log(2 * pi())
         + 0.5 * log_determinant(V)
         + lgamma(a)
         - a * log(b);
  }

  // Pooled NPP effective-prior NIG parameters, packed as
  // [vec(V_eta_inv) (p*p, column-major), mu_eta (p), a_eta (1), b_eta (1)]
  vector npp_params(vector eta, vector X0, array[] int startid_X0,
                     vector y0, array[] int startid_y0, array[] int n0,
                     matrix V0inv, vector mu0, real a, real b, int K, int p) {
    matrix[p,p] V_eta_inv = V0inv;
    vector[p]   rhs       = V0inv * mu0;
    real        a_eta     = a;
    real        b_eta_sum = 0;  // accumulates eta_k * y0k'y0k

    for (k in 1:K) {
      matrix[n0[k],p] X0k = build_matrix_from_vector(X0, startid_X0[k], n0[k], p);
      vector[n0[k]]   y0k = segment(y0, startid_y0[k], n0[k]);
      V_eta_inv += eta[k] * X0k' * X0k;
      rhs       += eta[k] * X0k' * y0k;
      a_eta     += 0.5 * eta[k] * n0[k];
      b_eta_sum += eta[k] * dot_self(y0k);
    }

    matrix[p,p] V_eta  = inverse_spd(V_eta_inv);
    vector[p]   mu_eta = V_eta * rhs;
    real b_eta = b + 0.5 * (b_eta_sum
                             + quad_form(V0inv, mu0)
                             - quad_form(V_eta_inv, mu_eta));

    return append_row(append_row(to_vector(V_eta_inv), mu_eta), [a_eta, b_eta]');
  }

  // Pooled NPP-SEQ NIG parameters, packed as
  // [vec(V_seq_inv) (p*p, column-major), mu_seq (p), a_seq (1), b_seq (1), logZ_k_sum (1)]
  vector seq_params(vector eta, vector X0, array[] int startid_X0,
                     vector y0, array[] int startid_y0, array[] int n0,
                     matrix V0inv, vector mu0, real a, real b, int K, int p) {
    matrix[p,p] V_seq_inv   = rep_matrix(0, p, p);
    vector[p]   rhs_seq     = rep_vector(0, p);
    real        a_seq       = (K-1)*(0.5*p+1);
    real        b_seq_quad  = 0;
    real        logZ_k_sum  = 0;

    for (k in 1:K) {
      matrix[n0[k],p] X0k = build_matrix_from_vector(X0, startid_X0[k], n0[k], p);
      vector[n0[k]]   y0k = segment(y0, startid_y0[k], n0[k]);

      // per-dataset NIG update
      matrix[p,p] Vk_inv = V0inv + eta[k] * X0k' * X0k;
      matrix[p,p] Vk     = inverse_spd(Vk_inv);
      vector[p]   muk    = Vk * (V0inv * mu0 + eta[k] * X0k' * y0k);
      real        ak     = a + 0.5 * eta[k] * n0[k];
      real        bk     = b + 0.5 * (eta[k] * dot_self(y0k)
                             + quad_form(V0inv, mu0)
                             - quad_form(Vk_inv, muk));

      logZ_k_sum += log_Z_NIG(muk, Vk, ak, bk, p);

      // accumulate V_seq_inv = sum_k Vk_inv - (K-1)*V0inv
      V_seq_inv  += Vk_inv;
      rhs_seq    += Vk_inv * muk;
      a_seq      += ak;
      b_seq_quad += bk + 0.5 * quad_form(Vk_inv, muk);
    }

    matrix[p,p] V_seq  = inverse_spd(V_seq_inv);
    vector[p]   mu_seq = V_seq * rhs_seq;
    real b_seq = b_seq_quad - 0.5 * quad_form(V_seq_inv, mu_seq);

    return append_row(append_row(append_row(to_vector(V_seq_inv), mu_seq), [a_seq, b_seq]'), [logZ_k_sum]');
  }
}

data {
  int<lower=1> K;
  array[K] int<lower=0> n0;
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> len_X0;
  int<lower=0> len_y0;
  vector[len_X0] X0;
  array[K] int<lower=0> startid_X0;
  vector[len_y0] y0;
  array[K] int<lower=0> startid_y0;
  matrix[n,p] X;
  vector[n] y;
  real<lower=0> a;
  real<lower=0> b;
  matrix[p,p] V0;
  vector[p] mu0;
  vector[K+1] alpha;
  int<lower=0, upper=1> post;
  int<lower=0, upper=1> seq;
}

transformed data {
  matrix[p,p] V0inv = inverse_spd(V0);
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
  // ----- prior on eta -----
  target += gamma_lpdf(raw_gamma | alpha, 1);

  if (post == 1) {

    if (seq == 0) {
      // -------------------------------------------------------
      // NPP: log Z(mu_eta, V_eta, a_eta, b_eta)
      //      then if current data, log Z(tilde_mu_eta, ...) - log Z(mu_eta, ...)
      // -------------------------------------------------------

      vector[p*p + p + 2] pk = npp_params(eta, X0, startid_X0, y0, startid_y0, n0,
                                           V0inv, mu0, a, b, K, p);
      matrix[p,p] V_eta_inv = to_matrix(head(pk, p*p), p, p);
      vector[p]   mu_eta    = segment(pk, p*p + 1, p);
      real        a_eta     = pk[p*p + p + 1];
      real        b_eta     = pk[p*p + p + 2];

      matrix[p,p] V_eta = inverse_spd(V_eta_inv);
      real logZ_eta = log_Z_NIG(mu_eta, V_eta, a_eta, b_eta, p);

      // now add log BF for current data D
      // log Z(tilde_mu_eta, tilde_V_eta, tilde_a_eta, tilde_b_eta) - logZ_eta
      matrix[p,p] tilde_V_eta_inv = V_eta_inv + X' * X;
      matrix[p,p] tilde_V_eta     = inverse_spd(tilde_V_eta_inv);
      vector[p]   tilde_mu_eta    = tilde_V_eta * (V_eta_inv * mu_eta + X' * y);
      real        tilde_a_eta     = a_eta + 0.5 * n;
      real        tilde_b_eta     = b_eta + 0.5 * (dot_self(y)
                                      + quad_form(V_eta_inv, mu_eta)
                                      - quad_form(tilde_V_eta_inv, tilde_mu_eta));

      real logZ_tilde = log_Z_NIG(tilde_mu_eta, tilde_V_eta, tilde_a_eta, tilde_b_eta, p);
      target += logZ_tilde - logZ_eta;

    } else {
      // -------------------------------------------------------
      // NPP-SEQ: sum_k log Z(mu_k*, V_k*, a_k*, b_k*) terms
      //          then log Z(mu_seq,...) - sum_k log Z(mu_k*,...)
      //          then log BF for current data
      // -------------------------------------------------------

      vector[p*p + p + 3] pk = seq_params(eta, X0, startid_X0, y0, startid_y0, n0,
                                           V0inv, mu0, a, b, K, p);
      matrix[p,p] V_seq_inv  = to_matrix(head(pk, p*p), p, p);
      vector[p]   mu_seq     = segment(pk, p*p + 1, p);
      real        a_seq      = pk[p*p + p + 1];
      real        b_seq      = pk[p*p + p + 2];
      real        logZ_k_sum = pk[p*p + p + 3];

      // log BF for current data D
      matrix[p,p] tilde_V_seq_inv = V_seq_inv + X' * X;
      matrix[p,p] tilde_V_seq     = inverse_spd(tilde_V_seq_inv);
      vector[p]   tilde_mu_seq    = tilde_V_seq * (V_seq_inv * mu_seq + X' * y);
      real        tilde_a_seq     = a_seq + 0.5 * n;
      real        tilde_b_seq     = b_seq + 0.5 * (dot_self(y)
                                      + quad_form(V_seq_inv, mu_seq)
                                      - quad_form(tilde_V_seq_inv, tilde_mu_seq));

      real logZ_tilde = log_Z_NIG(tilde_mu_seq, tilde_V_seq, tilde_a_seq, tilde_b_seq, p);
      target += logZ_tilde - logZ_k_sum;
    }
  } else {
    if (seq == 1) {
      vector[p*p + p + 3] pk = seq_params(eta, X0, startid_X0, y0, startid_y0, n0,
                                           V0inv, mu0, a, b, K, p);
      matrix[p,p] V_seq_inv  = to_matrix(head(pk, p*p), p, p);
      vector[p]   mu_seq     = segment(pk, p*p + 1, p);
      real        a_seq      = pk[p*p + p + 1];
      real        b_seq      = pk[p*p + p + 2];
      real        logZ_k_sum = pk[p*p + p + 3];

      matrix[p,p] V_seq = inverse_spd(V_seq_inv);
      real logZ_seq = log_Z_NIG(mu_seq, V_seq, a_seq, b_seq, p);
      target += logZ_seq - logZ_k_sum;
    }
  }
}

generated quantities {
  // -------------------------------------------------------
  // Posterior predictive parameters (conditional on eta)
  // returned for downstream use and marginal reconstruction
  // -------------------------------------------------------

  vector[p]   tilde_mu;
  matrix[p,p] tilde_V;
  real        tilde_a;
  real        tilde_b;

  {
    if (post == 1) {
      if (seq == 0) {
        // -------------------------------------------------------
        // NPP: rebuild effective prior (mu_eta, V_eta, a_eta, b_eta)
        //      then update with current data D
        // -------------------------------------------------------

        vector[p*p + p + 2] pk = npp_params(eta, X0, startid_X0, y0, startid_y0, n0,
                                             V0inv, mu0, a, b, K, p);
        matrix[p,p] V_eta_inv = to_matrix(head(pk, p*p), p, p);
        vector[p]   mu_eta    = segment(pk, p*p + 1, p);
        real        a_eta     = pk[p*p + p + 1];
        real        b_eta     = pk[p*p + p + 2];

        // update with D
        matrix[p,p] tilde_V_inv = V_eta_inv + X' * X;
        tilde_V  = inverse_spd(tilde_V_inv);
        tilde_mu = tilde_V * (V_eta_inv * mu_eta + X' * y);
        tilde_a  = a_eta + 0.5 * n;
        tilde_b  = b_eta + 0.5 * (dot_self(y)
                    + quad_form(V_eta_inv, mu_eta)
                    - quad_form(tilde_V_inv, tilde_mu));

      } else {
        // -------------------------------------------------------
        // NPP-SEQ: rebuild (mu_seq, V_seq, a_seq, b_seq)
        //          then update with current data D
        // -------------------------------------------------------

        vector[p*p + p + 3] pk = seq_params(eta, X0, startid_X0, y0, startid_y0, n0,
                                             V0inv, mu0, a, b, K, p);
        matrix[p,p] V_seq_inv = to_matrix(head(pk, p*p), p, p);
        vector[p]   mu_seq    = segment(pk, p*p + 1, p);
        real        a_seq     = pk[p*p + p + 1];
        real        b_seq     = pk[p*p + p + 2];

        // update with D
        matrix[p,p] tilde_V_inv = V_seq_inv + X' * X;
        tilde_V  = inverse_spd(tilde_V_inv);
        tilde_mu = tilde_V * (V_seq_inv * mu_seq + X' * y);
        tilde_a  = a_seq + 0.5 * n;
        tilde_b  = b_seq + 0.5 * (dot_self(y)
                    + quad_form(V_seq_inv, mu_seq)
                    - quad_form(tilde_V_inv, tilde_mu));
      }
    } else {
      if (seq == 1) {
        tilde_mu = mu0;
        tilde_V = V0/K;
        tilde_a = K * (a + 0.5 * p + 1) - 0.5 * p - 1;
        tilde_b = K*b;
      } else {
        tilde_mu = mu0;
        tilde_V = V0;
        tilde_a = a;
        tilde_b = b;
      }
    }
  }

  // -------------------------------------------------------
  // Marginal posterior for beta: multivariate t
  // beta | eta, D, D0 ~ t_{2*tilde_a}(tilde_mu, tilde_b/tilde_a * tilde_V)
  // Draw via: beta = tilde_mu + L * z / sqrt(u / (2*tilde_a))
  // where z ~ N(0,I), u ~ chi2(2*tilde_a), L = cholesky(tilde_b/tilde_a * tilde_V)
  // -------------------------------------------------------

  vector[p] beta_draw;
  real      sigma2_draw;
  vector[p+1] theta;

  {
    // draw sigma2 | eta, D, D0 ~ Inv-Gamma(tilde_a, tilde_b)
    sigma2_draw = inv_gamma_rng(tilde_a, tilde_b);

    // draw beta | sigma2, eta, D, D0 ~ N(tilde_mu, sigma2 * tilde_V)
    matrix[p,p] L = cholesky_decompose(tilde_V);
    beta_draw = multi_normal_cholesky_rng(tilde_mu, sqrt(sigma2_draw) * L);
    theta = append_row(beta_draw, sigma2_draw);
  }

  // -------------------------------------------------------
  // Marginal moments of beta (integrating out sigma2)
  // E[beta | eta, D, D0]   = tilde_mu          (tilde_a > 1/2, always true)
  // Var[beta | eta, D, D0] = tilde_b/(tilde_a - 1) * tilde_V  (tilde_a > 1)
  // -------------------------------------------------------
}
