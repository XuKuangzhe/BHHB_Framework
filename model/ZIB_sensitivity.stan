// =============================================================================
// Model: ZIB Sensitivity Analysis
// Purpose: Identical to the main ZIB model, but accepts custom hyperparameters 
//          via the data block to test prior sensitivity (robustness check).
// =============================================================================

data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> S;
  int<lower=1> P;
  int<lower=1,upper=S> SID[N];
  int<lower=1,upper=P> PID[N];
  vector<lower=0,upper=1>[N] Y; 
  matrix[N,D] X;
  int<lower=0,upper=1> Y_binary[N]; 
  int<lower=0> N_pos;                
  int<lower=1,upper=N> pos_idx[N_pos]; 

  // === Hyperparameters for Sensitivity Analysis ===
  real<lower=0> prior_beta_sd;   // Prior SD for fixed effects (e.g., 5 vs 10)
  real<lower=0> prior_phi_sd;    // Prior SD for Phi
  real<lower=0> prior_sd_scale;  // Prior Scale for Random Effects SD (e.g., 2.5 vs 5)
}

parameters {
  vector[D] bz[2];
  
  vector[P] z_pB;
  vector[P] z_pG;
  vector[S] z_sB;
  vector[S] z_sG;
  
  real<lower=0> s_pB;
  real<lower=0> s_sB;
  real<lower=0> s_pG;
  real<lower=0> s_sG;
  
  real<lower=0> phi;
}

transformed parameters {
  vector[N] q1; 
  vector[N] q2;
  vector<lower=0>[N] a;
  vector<lower=0>[N] b;
  
  vector[P] r_pB = z_pB * s_pB;
  vector[S] r_sB = z_sB * s_sB;
  vector[P] r_pG = z_pG * s_pG;
  vector[S] r_sG = z_sG * s_sG;
  
  q1 = inv_logit(X * bz[1] + r_pB[PID] + r_sB[SID]);
  q2 = inv_logit(X * bz[2] + r_pG[PID] + r_sG[SID]);
  
  a = q2 * phi;
  b = (1 - q2) * phi;
}

model {
  // Non-centered parameterization (Structural, fixed to std_normal)
  z_pB ~ std_normal();
  z_sB ~ std_normal();
  z_pG ~ std_normal();
  z_sG ~ std_normal();
  
  // === Modified: Use passed hyperparameters ===
  // Random Effect SDs: use prior_sd_scale
  s_pB ~ student_t(3, 0, prior_sd_scale);
  s_sB ~ student_t(3, 0, prior_sd_scale);
  s_pG ~ student_t(3, 0, prior_sd_scale);
  s_sG ~ student_t(3, 0, prior_sd_scale);
  
  // Fixed Effects: use prior_beta_sd
  bz[1] ~ normal(0, prior_beta_sd); 
  bz[2] ~ normal(0, prior_beta_sd);

  // Phi: use prior_phi_sd
  phi ~ normal(0, prior_phi_sd);

  // Likelihood (Same as main model)
  target += bernoulli_lpmf(Y_binary | q1);
  target += beta_lpdf(Y[pos_idx] | a[pos_idx], b[pos_idx]);
}
