// =============================================================================
// Model: Linear Mixed Model (LMM)
// Purpose: Baseline competitor model assuming Gaussian errors.
// Note: This model often fails for zero-inflated bounded data (e.g., proportions).
// =============================================================================

data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1,upper=S> SID[N];
  int<lower=1,upper=P> PID[N];
  vector[N] Y;      // Response variable
  matrix[N,D] X;    // Design matrix
}

parameters {
  vector[D] beta;       // Fixed effects
  real<lower=0> sigma;  // Residual standard deviation
  
  // Random effects (Non-centered parameterization)
  vector[S] z_s; 
  vector[P] z_p;
  real<lower=0> sigma_s;
  real<lower=0> sigma_p;
}

transformed parameters {
  vector[S] r_s = z_s * sigma_s;
  vector[P] r_p = z_p * sigma_p;
  vector[N] mu;
  
  // Linear predictor
  mu = X * beta + r_s[SID] + r_p[PID];
}

model {
  // --- Priors ---
  beta ~ normal(0, 5);
  sigma ~ student_t(3, 0, 2.5);
  
  z_s ~ std_normal();
  z_p ~ std_normal();
  sigma_s ~ student_t(3, 0, 2.5);
  sigma_p ~ student_t(3, 0, 2.5);

  // --- Likelihood (Gaussian) ---
  Y ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] Y_rep;
  
  for (n in 1:N) {
    // Generate posterior predictions for PPC
    Y_rep[n] = normal_rng(mu[n], sigma);
    // Calculate log-likelihood for LOO
    log_lik[n] = normal_lpdf(Y[n] | mu[n], sigma);
  }
}
