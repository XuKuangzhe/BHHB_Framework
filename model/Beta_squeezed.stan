// =============================================================================
// Model: Beta Regression (Squeezed)
// Purpose: Competitor model ignoring structural zeros.
// Note: Requires data transformation (squeezing) to avoid 0 and 1, 
//       as standard Beta is defined on (0, 1).
// =============================================================================

data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1,upper=S> SID[N];
  int<lower=1,upper=P> PID[N];
  vector<lower=0,upper=1>[N] Y; // Response variable (squeezed)
  matrix[N,D] X;                // Design matrix
}

parameters {
  vector[D] b_beta;       // Fixed effects
  real<lower=0> phi;      // Precision parameter
  
  // Random effects (Non-centered parameterization)
  vector[S] z_s; 
  vector[P] z_p;
  real<lower=0> sigma_s;
  real<lower=0> sigma_p;
}

transformed parameters {
  vector[S] r_s = z_s * sigma_s;
  vector[P] r_p = z_p * sigma_p;
  vector<lower=0,upper=1>[N] mu;
  vector<lower=0>[N] a;
  vector<lower=0>[N] b;
  
  // Inverse Logit Link for mean
  mu = inv_logit(X * b_beta + r_s[SID] + r_p[PID]);
  
  // Reparameterization
  a = mu * phi;
  b = (1 - mu) * phi;
}

model {
  // --- Priors ---
  b_beta ~ normal(0, 5);
  phi ~ normal(0, 100);
  
  z_s ~ std_normal();
  z_p ~ std_normal();
  sigma_s ~ student_t(3, 0, 2.5);
  sigma_p ~ student_t(3, 0, 2.5);

  // --- Likelihood (Standard Beta) ---
  Y ~ beta(a, b);
}

generated quantities {
  vector[N] Y_rep;
  for (n in 1:N) {
    Y_rep[n] = beta_rng(a[n], b[n]);
  }
}
