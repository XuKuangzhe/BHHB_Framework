// =============================================================================
// Model: Zero-Inflated Beta (ZIB) with Mixed Effects
// Purpose: Main analysis model for eye-tracking data (Gaze vs. Mouse).
// Structure: 
//    1. Bernoulli process (Decision): Probability of non-zero behavior.
//    2. Beta process (Intensity): Magnitude of behavior (given non-zero).
//    - Both processes include fixed effects (X) and random intercepts (Subjects/Items).
// =============================================================================

data {
  int<lower=1> N;                   // Total number of observations
  int<lower=1> D;                   // Number of predictors (Design matrix columns)
  int<lower=1> S;                   // Number of subjects
  int<lower=1> P;                   // Number of items (pictures)
  int<lower=1,upper=S> SID[N];      // Subject IDs
  int<lower=1,upper=P> PID[N];      // Item IDs
  vector<lower=0,upper=1>[N] Y;     // Response variable (e.g., proportion/weight)
  matrix[N,D] X;                    // Design matrix for fixed effects
  
  // Data split for efficiency and zero-inflation handling
  int<lower=0,upper=1> Y_binary[N]; // Binary indicator: 1 if Y > 0, else 0
  int<lower=0> N_pos;               // Number of positive (non-zero) observations
  int<lower=1,upper=N> pos_idx[N_pos]; // Indices of positive observations
}

parameters {
  // Fixed effects coefficients
  // bz[1]: Coefficients for Bernoulli process (Decision)
  // bz[2]: Coefficients for Beta process (Intensity)
  vector[D] bz[2];
  
  // Random effects (Non-centered parameterization)
  // z_*: Standard normal deviates
  // s_*: Standard deviations (scales)
  // Suffix 'B' = Bernoulli process, 'G' = Beta (Gamma-like) process
  vector[P] z_pB;
  vector[P] z_pG;
  vector[S] z_sB;
  vector[S] z_sG;
  
  real<lower=0> s_pB;
  real<lower=0> s_sB;
  real<lower=0> s_pG;
  real<lower=0> s_sG;
  
  // Precision parameter for Beta distribution
  real<lower=0> phi;
}

transformed parameters {
  vector[N] q1;       // Probability of non-zero response (Bernoulli theta)
  vector[N] q2;       // Mean of the Beta distribution (Beta mu)
  vector<lower=0>[N] a; // Beta shape parameter alpha
  vector<lower=0>[N] b; // Beta shape parameter beta
  
  // Recover random effects (Scale * Z)
  vector[P] r_pB = z_pB * s_pB;
  vector[S] r_sB = z_sB * s_sB;
  vector[P] r_pG = z_pG * s_pG;
  vector[S] r_sG = z_sG * s_sG;
  
  // Linear predictors with inverse-logit link
  // Note: Vector operations imply element-wise addition
  q1 = inv_logit(X * bz[1] + r_pB[PID] + r_sB[SID]);
  q2 = inv_logit(X * bz[2] + r_pG[PID] + r_sG[SID]);
  
  // Reparameterization of Beta distribution
  // mean = a / (a + b) = q2
  // precision = a + b = phi
  a = q2 * phi;
  b = (1 - q2) * phi;
}

model {
  // --- Priors ---
  // Weakly informative priors for structure parameters
  z_pB ~ std_normal();
  z_sB ~ std_normal();
  z_pG ~ std_normal();
  z_sG ~ std_normal();
  
  s_pB ~ student_t(3, 0, 2.5);
  s_sB ~ student_t(3, 0, 2.5);
  s_pG ~ student_t(3, 0, 2.5);
  s_sG ~ student_t(3, 0, 2.5);
  
  phi ~ normal(0, 100);
  bz[1] ~ normal(0, 5); 
  bz[2] ~ normal(0, 5);

  // --- Likelihood ---
  // 1. Bernoulli part (Zero vs. Non-zero)
  target += bernoulli_lpmf(Y_binary | q1);
  
  // 2. Beta part (Magnitude given non-zero)
  // Applied only to non-zero observations using pos_idx
  target += beta_lpdf(Y[pos_idx] | a[pos_idx], b[pos_idx]);
}

generated quantities {
  vector[N] log_lik;
  vector[N] Y_rep; 
  int is_zero;
  real q1_safe;
  real beta_dens;

  for (n in 1:N) {
    // 1. Posterior Predictive Check (Generate Y_rep)
    is_zero = bernoulli_rng(1 - q1[n]);
    if (is_zero == 1) {
      Y_rep[n] = 0;
    } else {
      Y_rep[n] = beta_rng(a[n], b[n]);
    }

    // 2. Compute Log-Likelihood (with numerical stability protection)
    // Clamp q1 within [1e-16, 1 - 1e-16] to prevent log(0)
    q1_safe = fmax(1e-16, fmin(1 - 1e-16, q1[n]));

    if (Y_binary[n] == 0) { 
      // Case Y=0: log(P(Y=0)) = log(1 - q1)
      log_lik[n] = log1m(q1_safe);
    } else {
      // Case Y>0: log(P(Y>0)) + log(Beta_density)
      // Calculate Beta log-density
      beta_dens = beta_lpdf(Y[n] | a[n], b[n]);
      
      // Protection against infinite density (model over-certainty)
      if (is_inf(beta_dens)) {
         beta_dens = 100; // Cap at a large finite value
      }
      
      log_lik[n] = log(q1_safe) + beta_dens;
    }
  }
}
