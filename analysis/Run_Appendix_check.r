# ==============================================================================
# Appendix B. Forensic Analysis
# ==============================================================================
# Setup: Define Targets
# ==============================================================================
library(tidyverse)
library(glmmTMB)

# Ensure data is loaded: sumGDTs
if(!exists("sumGDTs")) stop("Please load sumGDTs first.")

# Target definitions
exp_list <- c("S8", "S12")
imp_list <- unique(sumGDTs$infID)
aoi_list <- c("GazeGL", "GazeFH", "GazeEB", "GazeEY", "GazeNO", "GazeMO")

# Initialize result container
forensic_report <- tibble()

print(">>> Starting Forensic Analysis of MLE Instability...")

# ==============================================================================
# Main Loop
# ==============================================================================
for(exp in exp_list) {
  for(imp in imp_list) {
    for(aoi in aoi_list) {
      
      # 1. Data Preparation (With strict hygiene)
      # -------------------------------------------------------
      target_data <- sumGDTs %>% 
        filter(Valp >= 0.6,          # CRITICAL: Consistency with main analysis
               expName == exp, 
               infID == imp) %>%
        select(participant, picname, all_of(aoi)) %>%
        rename(Y = !!sym(aoi))
      
      # Skip empty or insufficient datasets
      if(nrow(target_data) < 10) next
      
      # Calculate Sparsity
      zero_prop <- mean(target_data$Y == 0)
      
      # OPTIONAL: Acceleration Strategy
      # Only run full diagnostics on low-sparsity regimes (< 10% zeros)
      # Comment out this line if you want to check EVERYTHING
      if(zero_prop > 0.10) next 
      
      label <- paste(exp, imp, aoi, sep = " | ")
      print(paste0("Processing: ", label, " (Zeros: ", round(zero_prop*100, 1), "%)"))
      
      # 2. Force Fit glmmTMB
      # -------------------------------------------------------
      df_freq <- target_data %>%
        mutate(sub = as.factor(participant), item = as.factor(picname))
      
      # Use try() to catch hard crashes
      m_attempt <- try(glmmTMB(Y ~ 1 + (1|sub) + (1|item), 
                               ziformula = ~ 1,
                               family = ordbeta(), 
                               data = df_freq), silent = TRUE)
      
      # 3. Forensic Diagnosis
      # -------------------------------------------------------
      diag_row <- tibble(
        Exp = exp,
        Imp = imp,
        AOI = aoi,
        Zero_Prop = zero_prop,
        Status = "OK",
        Hessian_OK = NA,
        Conv_Msg = NA,
        ZI_Intercept = NA,
        ZI_SE = NA
      )
      
      if(inherits(m_attempt, "try-error")) {
        # Case A: Hard Crash
        diag_row$Status <- "CRASH"
        diag_row$Conv_Msg <- as.character(m_attempt)
        
      } else {
        # Case B: Model returned object (Check for soft failure)
        sdr <- m_attempt$sdr
        fit_summ <- summary(m_attempt)
        
        # Check 1: Hessian Positive Definite?
        diag_row$Hessian_OK <- sdr$pdHess
        
        # Check 2: Optimizer Message
        diag_row$Conv_Msg <- m_attempt$fit$message
        
        # Check 3: Parameter Explosion (Zero-Inflation Part)
        # Extract ZI Intercept
        zi_coefs <- fit_summ$coefficients$zi
        if(nrow(zi_coefs) > 0) {
          est <- zi_coefs[1, "Estimate"]
          se  <- zi_coefs[1, "Std. Error"]
          
          diag_row$ZI_Intercept <- est
          diag_row$ZI_SE <- se
          
          # Diagnosis Logic
          is_exploded <- abs(est) > 10
          is_nan_se <- is.nan(se) || is.na(se)
          is_bad_hess <- isFALSE(sdr$pdHess)
          
          if(is_exploded || is_nan_se || is_bad_hess) {
            diag_row$Status <- "UNSTABLE"
          }
        }
      }
      
      # 4. Store Evidence
      # -------------------------------------------------------
      forensic_report <- bind_rows(forensic_report, diag_row)
    }
  }
}

#write_csv(forensic_report,"forensic_report.csv")
# ==============================================================================
# Output & Review
# ==============================================================================
print(">>> Forensic Analysis Complete.")

# Filter for the "Guilty" cases to put in Appendix
evidence_table <- forensic_report %>%
  filter(Status != "OK") %>%
  arrange(Zero_Prop) %>%
  select(Exp, Imp, AOI, Zero_Prop, Status, Hessian_OK, ZI_Intercept, ZI_SE)

# Display Top Evidences
print(head(evidence_table, 10))

# Save for Appendix creation
# write_csv(evidence_table, "Appendix_B_Evidence.csv")


# ==============================================================================
# Appendix C. Robustness_Dual
# ==============================================================================
# Simulation: Robustness Stress Test (Logit-Normal Process)
# Purpose: Generate the "Robustness Horizon" plot for Appendix D.
#          Tests parameter recovery for BOTH Decision and Intensity processes.
# ==============================================================================
library(tidyverse)
library(rstan)
library(patchwork)

# 0. Load Helper Functions
if(file.exists("analysis/HB_Functions.R")) {
  source("analysis/HDI.R")
  source("analysis/HB_Functions.R")
} else {
  try(source("HB_Functions.R"))
  try(source("HDI.R"))
}

# 1. Define Data Generation Function (Logit-Normal Process)
generate_stress_data <- function(sigma, N_subj=40, N_pic=40) {
  
  # Ground Truth Parameters
  # Decision Process (Bernoulli)
  b_bern_true <- c(1.0, -0.5) # Slope Target: -0.5
  
  # Intensity Process (Logit-Normal Misspecification)
  b_int_true  <- c(-0.5, 0.6) # Slope Target: 0.6
  
  # Design Matrix
  grid <- expand.grid(SID = 1:N_subj, PID = 1:N_pic)
  N <- nrow(grid)
  X_vec <- sample(c(0, 1), N, replace = TRUE)
  X_matrix <- cbind(1, X_vec)
  
  # Random Effects
  r_s <- rnorm(N_subj, 0, 0.3); r_p <- rnorm(N_pic, 0, 0.3)
  
  # A. Hurdle Process (Bernoulli)
  logit_p <- X_matrix %*% b_bern_true + r_s[grid$SID] + r_p[grid$PID]
  is_nonzero <- rbinom(N, 1, plogis(logit_p))
  
  # B. Intensity Process (Logit-Normal with varying sigma)
  latent_y <- X_matrix %*% b_int_true + r_s[grid$SID] + r_p[grid$PID] + rnorm(N, 0, sigma)
  y_continuous <- plogis(latent_y)
  
  # Combine Processes & Clamp
  Y <- y_continuous * is_nonzero
  Y[Y > 0.999] <- 0.999; Y[Y > 0 & Y < 0.001] <- 0.001
  
  list(dat = list(N=N, D=2, S=N_subj, P=N_pic, SID=grid$SID, PID=grid$PID, 
                  Y=as.vector(Y), X=X_matrix, Y_binary=as.integer(Y>0), 
                  N_pos=sum(Y>0), pos_idx=which(Y>0)))
}

# 2. Run 10-Level Spectrum Loop
sigma_seq <- seq(0.2, 3.0, length.out = 10) 
res_all <- tibble()

if(!exists("hb_model")) hb_model <- stan_model("models/HB_main.stan")

print(">>> Starting Full-Spectrum Stress Test (Dual Process)...")

pb <- txtProgressBar(min = 0, max = length(sigma_seq), style = 3)

for(i in 1:length(sigma_seq)) {
  sig <- sigma_seq[i]
  sim <- generate_stress_data(sigma = sig)
  
  # Sampling
  fit <- sampling(hb_model, data = sim$dat, chains = 2, iter = 2000, refresh = 0)
  post <- rstan::extract(fit)
  
  # --- Extract Decision Slope (Bernoulli) ---
  # True Value: -0.5
  stats_bern <- summary_MCMC(post$bz[,1,2])
  res_all <- bind_rows(res_all, tibble(
    Sigma = sig,
    Process = "Decision Process (Bernoulli)",
    Parameter = "Slope",
    True_Value = -0.5,
    Estimated = stats_bern["mean"],
    Lower = stats_bern["HDI_low"], 
    Upper = stats_bern["HDI_high"]
  ))
  
  # --- Extract Intensity Slope (Beta) ---
  # True Value: 0.6
  stats_beta <- summary_MCMC(post$bz[,2,2])
  res_all <- bind_rows(res_all, tibble(
    Sigma = sig,
    Process = "Intensity Process (Beta)",
    Parameter = "Slope",
    True_Value = 0.6,
    Estimated = stats_beta["mean"],
    Lower = stats_beta["HDI_low"], 
    Upper = stats_beta["HDI_high"]
  ))
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 3. Plot Dual-Process Robustness Horizon
p_curve <- ggplot(res_all, aes(x = Sigma, y = Estimated)) +
  # Dynamic Reference Line (mapped to True_Value column)
  geom_hline(aes(yintercept = True_Value), linetype = "dashed", color = "#c0392b", linewidth = 0.8) +
  geom_hline(aes(yintercept = 0), linetype = "solid", color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#3498db", alpha = 0.2) +
  geom_line(color = "#2980b9", linewidth = 1) +
  geom_point(size = 3, color = "#2c3e50") +
  # Facet by Process
  facet_wrap(~Process, scales = "free_y") +
  labs(
    title = "Appendix D: Robustness Horizon (Dual Process)",
    #subtitle = "Stability of parameter recovery across a continuous spectrum of distributional violation.",
    y = "Posterior Estimate (Slope)",
    x = "Degree of Misspecification (Logit-Normal Noise Sigma)",
    #caption = "Dashed red lines indicate True Parameter Values. Shaded areas represent 95% HDI."
  ) +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = seq(0.5, 3.0, 0.5))



print(p_curve)

# ggsave("results/Appendix_C_Robustness_Dual.png", p_curve, width = 10, height = 5)
