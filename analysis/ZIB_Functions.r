# ==============================================================================
# Project: Zero-Inflated Beta (ZIB) Model for Eye-Tracking Data
# Purpose: Comprehensive analysis pipeline including parameter recovery simulation, 
#          empirical fitting, posterior predictive checks (PPC), model comparison 
#          (CRPS/RMSE/MAE), and prior sensitivity analysis.
# Author: Kuangzhe Xu
# Date: 2026-01-09
# ==============================================================================

# ==============================================================================
# 00. Dependencies & Helper Functions
# ==============================================================================
required_packages <- c("tidyverse", "rstan", "bayesplot", "patchwork", "ggrepel", "scoringRules")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstan)
  library(bayesplot)
  library(patchwork)
  library(ggrepel)
  library(scoringRules)
})

# ==============================================================================
# 01. Parameter Recovery Simulation
# ==============================================================================

#' Generate Synthetic Data for ZIB Gradient Check
#' @param N_subj Number of subjects
#' @param N_pic Number of items/pictures
#' @param target_zero_prop Target proportion of zeros to simulate sparsity levels
generate_zib_gradient <- function(N_subj=40, N_pic=20, target_zero_prop = 0.5) {
  
  # 1. Inverse Logit to find base intercept for target zero proportion
  base_intercept <- qlogis(1 - target_zero_prop)
  
  # 2. Set Ground Truth Parameters
  b_bern_true = c(base_intercept, -0.5) 
  b_beta_true = c(-0.5, 0.3)
  phi_true = 10
  sd_subj = 0.5
  sd_pic = 0.5
  
  # 3. Generate Design Matrix
  grid <- expand.grid(SID = 1:N_subj, PID = 1:N_pic)
  N <- nrow(grid)
  X_vec <- sample(c(0, 1), N, replace = TRUE)
  X_matrix <- cbind(1, X_vec)
  
  # 4. Generate Random Effects
  r_s_bern <- rnorm(N_subj, 0, sd_subj); r_p_bern <- rnorm(N_pic, 0, sd_pic)
  r_s_beta <- rnorm(N_subj, 0, sd_subj); r_p_beta <- rnorm(N_pic, 0, sd_pic)
  
  # 5. Generate Response Variable Y
  logit_q1 <- (X_matrix %*% b_bern_true) + r_s_bern[grid$SID] + r_p_bern[grid$PID]
  q1 <- plogis(logit_q1) # Probability of Non-Zero
  
  logit_q2 <- (X_matrix %*% b_beta_true) + r_s_beta[grid$SID] + r_p_beta[grid$PID]
  q2 <- plogis(logit_q2) # Mean of Beta distribution
  
  Y <- numeric(N)
  # Vectorized generation is possible, but loop is kept for clarity in logic
  for(i in 1:N) {
    if(rbinom(1, 1, q1[i]) == 0) {
      Y[i] <- 0
    } else {
      shape1 <- q2[i] * phi_true
      shape2 <- (1 - q2[i]) * phi_true
      Y[i] <- rbeta(1, shape1, shape2)
      # Numerical stability truncation for Stan
      if(Y[i] > 0.999) Y[i] <- 0.999; if(Y[i] < 0.001) Y[i] <- 0.001
    }
  }
  
  # 6. Return Data List
  actual_zero <- mean(Y == 0)
  return(list(
    dat = list(N=N, D=2, S=N_subj, P=N_pic, SID=grid$SID, PID=grid$PID, 
               Y=Y, X=X_matrix, Y_binary=as.integer(Y>0), N_pos=sum(Y>0), pos_idx=which(Y>0)),
    truth = list(b_bern=b_bern_true, b_beta=b_beta_true, phi=phi_true),
    scenario = paste0("Sparsity: ", round(actual_zero*100), "%")
  ))
}

#' Run Gradient Simulation Calculation
#' @param stan_model_obj Compiled Stan model object
run_grad_sim_calc <- function(stan_model_obj, 
                              target_props = c(0.10, 0.30, 0.60),
                              iter = 2000, chains = 4, seed = 1234) {
  
  results_list <- list()
  print(">>> Starting Gradient Simulation Check...")
  
  for(tp in target_props) {
    # 1. Generate Data
    sim <- generate_zib_gradient(target_zero_prop = tp)
    scenario_label <- sim$scenario
    print(paste("Running Scenario:", scenario_label))
    
    # 2. Fit Model
    fit <- sampling(stan_model_obj, data = sim$dat, 
                    seed = seed, chains = chains, iter = iter, refresh = 0)
    
    # 3. Extract Results
    post <- rstan::extract(fit)
    
    # Helper to format stats
    get_sum <- function(vec) {
      st <- get_mcmc_stats(vec)
      list(mean = st["mean"], lower = st["HDI_low"], upper = st["HDI_high"])
    }
    
    # 4. Organize Parameters
    params <- list(
      "Bern_Int" = post$bz[,1,1], "Bern_Slope" = post$bz[,1,2],
      "Beta_Int" = post$bz[,2,1], "Beta_Slope" = post$bz[,2,2],
      "Phi" = post$phi
    )
    true_vals <- c(sim$truth$b_bern, sim$truth$b_beta, sim$truth$phi)
    
    for(i in 1:5) {
      p_name <- names(params)[i]
      s <- get_sum(params[[i]])
      results_list[[length(results_list)+1]] <- tibble(
        Scenario = scenario_label,
        Parameter = p_name,
        True_Value = true_vals[i],
        Recovered_Mean = s$mean,
        Lower = s$lower,
        Upper = s$upper
      )
    }
  }
  return(bind_rows(results_list))
}

#' Plot Parameter Recovery Results
plot_grad_sim <- function(recovery_data) {
  ggplot(recovery_data, aes(x = True_Value, y = Recovered_Mean, color = Parameter)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, alpha = 0.6) +
    geom_point(size = 3) +
    facet_wrap(~Scenario, nrow = 1) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Gradient Simulation: Parameter Recovery across Sparsity Levels",
      subtitle = "ZIB accurately recovers parameters regardless of zero-inflation rates (10% - 60%).",
      x = "True Parameter Value",
      y = "Recovered Posterior Mean",
      caption = "Error bars represent 95% HDI."
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
}

# ==============================================================================
# 02. The "Four-Area" Plot Analysis (Empirical)
# ==============================================================================

#' Run Empirical ZIB Analysis Loop
#' @description Iterates through AOIs and Impressions to calculate Decision vs. Intensity effects.
#' @note Requires specific column structure in 'dataset'. Ensure columns align with indices used.
run_empirical_zib_calc <- function(dataset, exp_label, stan_model_obj,
                                   GP_list = c("Glabella","Forehead","Eyebrow","Eye","Nose","Mouth"),
                                   inf_list = NULL, 
                                   seed = 1234) {
  
  perZIB <- tibble()
  rhatsum <- tibble()
  
  CN <- c("Intercept", "Type") 
  disT <- c("Bernoulli(Choice)", "Beta(Weight)")
  
  if(is.null(inf_list)) inf_list <- unique(dataset[[26]]) # Assumes infID is at col 26
  
  print(paste(">>> Starting ZIB Loop for Experiment:", exp_label))
  
  for(g in 1:length(GP_list)) { 
    for(f in 1:length(inf_list)) { 
      
      current_aoi <- GP_list[g]
      current_imp <- inf_list[f]
      print(paste("   -> Processing:", current_aoi, "|", current_imp))
      
      # === Data Preparation ===
      # WARNING: Hardcoded indices based on specific dataset structure
      # 1:participant, 2:expName, 3:picname, 13:Valp, g+13:GazeTarget, g+19:MouseTarget, 26:infID
      target_cols <- c(1, 2, 3, 13, g + 13, g + 19, 26)
      sub_dat <- dataset[, target_cols]
      colnames(sub_dat) <- c("participant", "expName", "picname", "Valp", 
                             "GazeTarget", "MouseTarget", "infID")
      
      TDt <- sub_dat %>% 
        filter(Valp >= 0.6, expName == exp_label, infID == current_imp)
      
      tarDat <- TDt %>% 
        dplyr::select(participant, picname, GazeTarget, MouseTarget) %>%
        pivot_longer(cols = c("GazeTarget", "MouseTarget"), 
                     names_to = "OB", values_to = "weight")
      
      # === Vectorization Prep ===
      Y_vec <- tarDat$weight
      # Beta safety truncation
      Y_vec_safe <- Y_vec
      Y_vec_safe[Y_vec_safe >= 0.9999] <- 0.9999
      
      Y_binary <- as.integer(Y_vec_safe > 0)
      pos_idx <- which(Y_vec_safe > 0)
      
      if(length(pos_idx) < 5) {
        warning(paste("Not enough non-zero data for", current_aoi, current_imp))
        next
      }
      
      datalist <- list(
        N = nrow(tarDat), D = 2,
        S = length(unique(tarDat$participant)),
        P = length(unique(tarDat$picname)),
        SID = as.numeric(as.factor(tarDat$participant)),
        PID = as.numeric(as.factor(tarDat$picname)),
        Y = Y_vec_safe,
        X = tibble(a=1, EV=as.numeric(as.factor(tarDat$OB))), 
        Y_binary = Y_binary,
        N_pos = length(pos_idx),
        pos_idx = pos_idx
      )
      
      # === Sampling ===
      step_result <- tryCatch({
        fitvf <- sampling(stan_model_obj, data=datalist, seed=seed, refresh=0, chains=4, cores=4)
        
        # A. Rhat Check
        rhat_vals <- summary(fitvf)$summary[,"Rhat"]
        is_converged <- all(rhat_vals <= 1.10, na.rm=T)
        
        current_rhat <- tibble(
          Exp = exp_label, Part = current_aoi, Impf = current_imp, 
          Rhat = is_converged, Max_Rhat = max(rhat_vals, na.rm=T)
        )
        
        # B. Parameter Extraction
        dft <- rstan::extract(fitvf)
        current_params <- tibble()
        
        for(i in 1:2) { # Bern vs Beta
          for(n in 1:2) { # Intercept vs Slope
            b_summary <- get_mcmc_stats(dft$bz[,i,n]) 
            
            tb <- tibble(
              Exp = exp_label, GazeP = current_aoi, Impf = current_imp, 
              DisType = disT[i], Para = CN[n], 
              mean = b_summary["mean"], 
              HDI_low = b_summary["HDI_low"], HDI_high = b_summary["HDI_high"]
            )
            current_params <- bind_rows(current_params, tb)
          }
        }
        
        # Clean up large objects
        rm(fitvf, dft); gc()
        
        list(status = "success", rhat = current_rhat, params = current_params)
        
      }, error = function(e) {
        message(paste("Error in", exp_label, ":", current_aoi, current_imp, "-", e$message))
        return(list(status = "error"))
      })
      
      # === Update ===
      if(step_result$status == "success") {
        rhatsum <- bind_rows(rhatsum, step_result$rhat)
        perZIB <- bind_rows(perZIB, step_result$params)
      }
    }
  }
  return(list(params = perZIB, rhat = rhatsum))
}

#' Plot the 4-Area Plot (Single Experiment)
plot_four_area_single <- function(data_file) {
  
  ZIBsum <- data_file
  
  data_wide <- ZIBsum %>% filter(Para != "Intercept") %>%
    dplyr::select(GazeP, Impf, DisType, Para, mean, HDI_low, HDI_high, Exp) %>%
    pivot_wider(names_from = DisType, values_from = c(mean, HDI_low, HDI_high)) %>%
    filter(!is.na(`mean_Bernoulli(Choice)`) & !is.na(`mean_Beta(Weight)`)) %>%
    mutate(
      sig_bern = `HDI_low_Bernoulli(Choice)` * `HDI_high_Bernoulli(Choice)` > 0,
      sig_beta = `HDI_low_Beta(Weight)` * `HDI_high_Beta(Weight)` > 0
    )
  
  ggplot(data_wide, aes(x = `mean_Bernoulli(Choice)`, y = `mean_Beta(Weight)`)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbarh(aes(xmin = `HDI_low_Bernoulli(Choice)`, xmax = `HDI_high_Bernoulli(Choice)`), height = 0, alpha = 0.3, color = "gray60") +
    geom_errorbar(aes(ymin = `HDI_low_Beta(Weight)`, ymax = `HDI_high_Beta(Weight)`), width = 0, alpha = 0.3, color = "gray60") +
    geom_point(aes(color = sig_bern, shape = sig_beta), size = 4, alpha = 0.8) +
    geom_text_repel(aes(label = paste(GazeP, Impf)), size = 3, max.overlaps = 20, box.padding = 0.3) +
    labs(
      title = "Decoupling Decision and Intensity: The 4-Area Plot",
      subtitle = "X-axis: Propensity to use Mouse (Decision). Y-axis: Duration of usage (Intensity).",
      x = "Decision Slope (Bernoulli Log-Odds)",
      y = "Intensity Slope (Beta Log-Odds)",
      color = "Decision Effect\n(Significant?)",
      shape = "Intensity Effect\n(Significant?)"
    ) +
    scale_color_manual(values = c("FALSE" = "#95a5a6", "TRUE" = "#e74c3c")) + 
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 19)) +
    theme_bw(base_size = 14) + theme(legend.position = "right")
}

#' Plot Combined 4-Area Plot for Two Experiments
#' 
#' @description 
#' Visualizes the dissociation between Decision (Bernoulli) and Intensity (Beta) 
#' processes across two experiments.
#'
#' @param data_s8 Data frame containing parameter estimates from Experiment 1 (e.g., S8).
#' @param data_s12 Data frame containing parameter estimates from Experiment 2 (e.g., S12).
#' @return A ggplot object.
plot_four_area_sum <- function(data_s8, data_s12) {
  
  # 1. Label and Merge Data
  # -------------------------------------------------------
  # Ensure 'Exp' column exists; assign default labels if missing
  if(!"Exp" %in% names(data_s8)) data_s8$Exp <- "Experiment 1 (S8)"
  if(!"Exp" %in% names(data_s12)) data_s12$Exp <- "Experiment 2 (S12)"
  
  # Force standardized display labels for the plot facets
  data_s8$Exp_Label <- "Experiment 1 (S8)"
  data_s12$Exp_Label <- "Experiment 2 (S12)"
  
  # Combine datasets
  ZIBsum <- bind_rows(data_s8, data_s12)
  
  # 2. Data Transformation (Long to Wide)
  # -------------------------------------------------------
  data_wide <- ZIBsum %>% 
    filter(Para != "Intercept") %>%
    # Select relevant columns
    dplyr::select(GazeP, Impf, DisType, Para, mean, HDI_low, HDI_high, Exp_Label) %>%
    # Pivot to wide format to have separate columns for Bernoulli and Beta estimates
    pivot_wider(names_from = DisType, values_from = c(mean, HDI_low, HDI_high)) %>%
    # Filter out incomplete entries (NA checks)
    filter(!is.na(`mean_Bernoulli(Choice)`) & !is.na(`mean_Beta(Weight)`)) %>%
    # Calculate Significance based on HDI (Does 0 fall outside the interval?)
    mutate(
      sig_bern = `HDI_low_Bernoulli(Choice)` * `HDI_high_Bernoulli(Choice)` > 0,
      sig_beta = `HDI_low_Beta(Weight)` * `HDI_high_Beta(Weight)` > 0,
      # Define grouping logic for visualization (Color/Shape mapping)
      sig_group = case_when(
        sig_bern & sig_beta ~ "Both Significant",
        sig_bern ~ "Choice Only",
        sig_beta ~ "Weight Only",
        TRUE ~ "Non-Significant"
      )
    )
  
  # 3. Plotting
  # -------------------------------------------------------
  p <- ggplot(data_wide, aes(x = `mean_Bernoulli(Choice)`, y = `mean_Beta(Weight)`)) +
    # A. Quadrant reference lines (x=0, y=0)
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    
    # B. Error bars (Horizontal for Decision/Bernoulli)
    geom_errorbarh(
      aes(xmin = `HDI_low_Bernoulli(Choice)`, xmax = `HDI_high_Bernoulli(Choice)`),
      height = 0, alpha = 0.3, color = "gray60"
    ) +
    # C. Error bars (Vertical for Intensity/Beta)
    geom_errorbar(
      aes(ymin = `HDI_low_Beta(Weight)`, ymax = `HDI_high_Beta(Weight)`),
      width = 0, alpha = 0.3, color = "gray60"
    ) +
    
    # D. Scatter points 
    # Color maps to Decision significance (Primary interest)
    # Shape maps to Intensity significance (Secondary interest)
    geom_point(aes(color = sig_bern, shape = sig_beta), size = 3, alpha = 0.8) +
    
    # E. Text labels (with auto-repulsion to avoid overlap)
    geom_text_repel(aes(label = paste(GazeP, Impf)), 
                    size = 2.5, max.overlaps = 20, box.padding = 0.3) +
    
    # F. Faceting by Experiment
    facet_wrap(~Exp_Label) +
    
    # G. Aesthetics and Labels
    labs(
      title = "Decoupling Decision and Intensity: The 4-Area Plot",
      subtitle = "Comparison of Decision (X) vs. Intensity (Y) processes across two experiments.",
      x = "Decision Slope (Bernoulli Log-Odds)",
      y = "Intensity Slope (Beta Log-Odds)",
      color = "Decision Significant?",
      shape = "Intensity Significant?"
    ) +
    scale_color_manual(values = c("FALSE" = "#95a5a6", "TRUE" = "#e74c3c")) + # Grey vs Red
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 19)) + # Open Circle vs Solid Dot
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  
  return(p)
}

# ==============================================================================
# 03. Posterior Predictive Checks (PPC) Comparison
# ==============================================================================

#' Run PPC Calculation (LMM vs Beta vs ZIB)
run_ppc_calc <- function(dataset, exp_str, imp_str, aoi_str, mouse_str, models_list, seed = 1234, iter = 2000) {
  
  print(paste(">>> Preparing Data for PPC:", exp_str, imp_str, aoi_str))
  
  TDt <- dataset %>% 
    filter(Valp >= 0.6, expName == exp_str, infID == imp_str) %>%
    dplyr::select(participant, picname, all_of(aoi_str), all_of(mouse_str))
  
  tarDat <- TDt %>% 
    pivot_longer(cols = c(all_of(aoi_str), all_of(mouse_str)), names_to = "OB", values_to = "weight")
  
  Y_vec <- tarDat$weight
  N_obs <- length(Y_vec)
  if(N_obs < 10) stop("Not enough data points.")
  
  Y_binary <- as.integer(Y_vec > 0)
  pos_idx <- which(Y_vec > 0)
  
  # --- Stan Data Prep ---
  data_common <- list(
    N = N_obs, D = 2,
    S = length(unique(tarDat$participant)), P = length(unique(tarDat$picname)),
    SID = as.numeric(as.factor(tarDat$participant)), PID = as.numeric(as.factor(tarDat$picname)),
    Y = Y_vec, X = tibble(a=1, EV=as.numeric(as.factor(tarDat$OB)))
  )
  
  data_zib <- c(data_common, list(Y_binary = Y_binary, N_pos = length(pos_idx), pos_idx = pos_idx))
  
  # Beta Transformation (Squeezing)
  Y_trans <- (Y_vec * (N_obs - 1) + 0.5) / N_obs
  data_beta <- data_common; data_beta$Y <- Y_trans
  
  print(">>> Fitting Models (LMM, Beta, ZIB)...")
  fit_lmm <- sampling(models_list$LMM, data = data_common, seed = seed, iter = iter, refresh = 0)
  fit_beta <- sampling(models_list$Beta, data = data_beta, seed = seed, iter = iter, refresh = 0)
  fit_zib <- sampling(models_list$ZIB, data = data_zib, seed = seed, iter = iter, refresh = 0)
  
  extract_rep <- function(fit) rstan::extract(fit)$Y_rep
  
  return(list(
    meta = list(exp = exp_str, imp = imp_str, aoi = aoi_str, zero_prop = mean(Y_vec == 0)),
    obs = Y_vec,
    y_rep = list(lmm = extract_rep(fit_lmm), beta = extract_rep(fit_beta), zib = extract_rep(fit_zib))
  ))
}

#' Plot PPC Comparison
plot_ppc_compare <- function(ppc_data, plot_xlim = c(-0.05, 0.08), n_draws = 200) {
  
  Y_vec <- ppc_data$obs
  y_rep <- ppc_data$y_rep
  meta <- ppc_data$meta
  samp_idx <- sample(nrow(y_rep$zib), min(n_draws, nrow(y_rep$zib)))
  
  # Density Plots
  p_lmm <- ppc_dens_overlay(Y_vec, y_rep$lmm[samp_idx, ]) + coord_cartesian(xlim = plot_xlim) + 
    labs(title = "A. LMM (Gaussian)", subtitle = "Failure: Leakage into negative values") + theme_bw() + theme(legend.position = "none")
  
  p_beta <- ppc_dens_overlay(Y_vec, y_rep$beta[samp_idx, ]) + coord_cartesian(xlim = plot_xlim) + 
    labs(title = "B. Beta (Transformed)", subtitle = "Failure: Cannot produce true zeros") + theme_bw() + theme(legend.position = "none")
  
  p_zib <- ppc_dens_overlay(Y_vec, y_rep$zib[samp_idx, ]) + coord_cartesian(xlim = plot_xlim) + 
    labs(title = "C. ZIB Model", subtitle = "Success: Captures bimodal structure") + theme_bw() + theme(legend.position = c(0.8, 0.7))
  
  # Zero Stats
  p_stat_lmm <- ppc_stat(Y_vec, y_rep$lmm, stat = function(y) mean(y <= 0)) + labs(title = "D. LMM Zero/Neg Pred", x = "Proportion <= 0") + theme_bw()
  p_stat_beta <- ppc_stat(Y_vec, y_rep$beta, stat = function(y) mean(y == 0)) + labs(title = "E. Beta Zero Pred", x = "Proportion == 0") + theme_bw()
  p_stat_zib <- ppc_stat(Y_vec, y_rep$zib, stat = function(y) mean(y == 0)) + labs(title = "F. ZIB Zero Pred", x = "Proportion == 0") + theme_bw()
  
  (p_lmm + p_beta + p_zib) / (p_stat_lmm + p_stat_beta + p_stat_zib) +
    plot_annotation(
      title = paste0("Methodological Breakdown: '", meta$aoi, "' (", round(meta$zero_prop*100, 1), "% Zeros)"),
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
}

# ==============================================================================
# 04. Model Performance Scoring (CRPS, RMSE, MAE)
# ==============================================================================

#' Calculate Performance Metrics for All Models
run_performance_calc <- function(dataset, models_list, exp_list = c("S8", "S12"),
                                 aoi_list = c("GazeGL", "GazeFH", "GazeEB", "GazeEY", "GazeNO", "GazeMO"),
                                 seed = 1234, iter = 2000) {
  
  comparison_results <- tibble()
  imp_list <- unique(dataset$infID)
  
  print(">>> Starting Universal Model Comparison (CRPS & RMSE & MAE)...")
  
  for(exp in exp_list) {
    for(aoi_idx in 1:length(aoi_list)) {
      target_aoi <- aoi_list[aoi_idx]
      short_name <- substr(target_aoi, 5, 6) 
      target_mouse <- paste0("Mouse", short_name)
      
      for(imp in imp_list) {
        label <- paste(exp, target_aoi, imp, sep="_")
        print(paste("Processing:", label))
        
        tryCatch({
          TDt <- dataset %>% filter(Valp >= 0.6, expName == exp, infID == imp) %>%
            dplyr::select(participant, picname, all_of(target_aoi), all_of(target_mouse))
          
          tarDat <- TDt %>% pivot_longer(cols = 3:4, names_to = "OB", values_to = "weight")
          Y_vec <- tarDat$weight; N_obs <- length(Y_vec)
          if(N_obs < 10) next
          
          # Stan Data Lists
          data_common <- list(
            N = N_obs, D = 2, S = length(unique(tarDat$participant)), P = length(unique(tarDat$picname)),
            SID = as.numeric(as.factor(tarDat$participant)), PID = as.numeric(as.factor(tarDat$picname)),
            Y = Y_vec, X = tibble(a=1, EV=as.numeric(as.factor(tarDat$OB)))
          )
          data_zib <- c(data_common, list(Y_binary = as.integer(Y_vec > 0), N_pos = sum(Y_vec > 0), pos_idx = which(Y_vec > 0)))
          data_beta <- data_common; data_beta$Y <- (Y_vec * (N_obs - 1) + 0.5) / N_obs
          
          # Fitting
          fit_lmm <- sampling(models_list$LMM, data = data_common, seed = seed, chains=2, cores=2, iter=iter, refresh=0)
          fit_beta <- sampling(models_list$Beta, data = data_beta, seed = seed, chains=2, cores=2, iter=iter, refresh=0)
          fit_zib <- sampling(models_list$ZIB, data = data_zib, seed = seed, chains=2, cores=2, iter=iter, refresh=0)
          
          # Extraction & Calculation
          y_rep_lmm <- t(rstan::extract(fit_lmm)$Y_rep)
          y_rep_beta <- t(rstan::extract(fit_beta)$Y_rep)
          y_rep_zib <- t(rstan::extract(fit_zib)$Y_rep)
          
          get_metrics <- function(truth, pred_mat) {
            pred_mean <- colMeans(pred_mat)
            resid <- truth - pred_mean
            c(rmse = sqrt(mean(resid^2)), mae = mean(abs(resid)))
          }
          m_lmm <- get_metrics(Y_vec, t(y_rep_lmm))
          m_beta <- get_metrics(Y_vec, t(y_rep_beta))
          m_zib <- get_metrics(Y_vec, t(y_rep_zib))
          
          this_res <- tibble(
            Exp = exp, AOI = target_aoi, Imp = imp, Zero_Prop = mean(Y_vec == 0),
            Model = c("LMM", "Beta", "ZIB"),
            CRPS = c(mean(crps_sample(y = Y_vec, dat = y_rep_lmm)), 
                     mean(crps_sample(y = Y_vec, dat = y_rep_beta)), 
                     mean(crps_sample(y = Y_vec, dat = y_rep_zib))),
            RMSE = c(m_lmm["rmse"], m_beta["rmse"], m_zib["rmse"]),
            MAE  = c(m_lmm["mae"], m_beta["mae"], m_zib["mae"])
          )
          comparison_results <- bind_rows(comparison_results, this_res)
          rm(fit_lmm, fit_beta, fit_zib, y_rep_lmm, y_rep_beta, y_rep_zib); gc()
          
        }, error = function(e) message(paste("Error in", label, ":", e$message)))
      }
    }
  }
  return(comparison_results)
}

plot_crps_trend <- function(results_df) {
  plot_data <- results_df %>% filter(Model != "LMM") %>% 
    left_join(results_df %>% filter(Model == "LMM") %>% dplyr::select(AOI, Imp, Exp, CRPS_LMM = CRPS), by = c("AOI", "Imp", "Exp")) %>%
    mutate(Improvement_CRPS = (CRPS_LMM - CRPS) / CRPS_LMM, Cluster = ifelse(Zero_Prop < 0.3, "Low_Sparsity", "High_Sparsity"))
  
  ggplot(plot_data, aes(x = Zero_Prop, y = Improvement_CRPS, color = Model, fill = Model)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_smooth(aes(group = interaction(Model, Cluster)), method = "lm", alpha = 0.2) +
    geom_point(alpha = 0.7, size = 3, shape = 21, color = "white", stroke = 0.5) +
    scale_y_continuous(labels = scales::percent, name = "Probabilistic Advantage (CRPS vs. LMM)") +
    scale_x_continuous(labels = scales::percent, name = "Data Sparsity (% Zeros)") +
    scale_color_manual(values = c("Beta" = "#e74c3c", "ZIB" = "#3498db")) +
    scale_fill_manual(values = c("Beta" = "#e74c3c", "ZIB" = "#3498db")) +
    labs(title = "Probabilistic Model Superiority (CRPS)", caption = "Trend lines fitted by sparsity cluster.") +
    theme_bw(base_size = 14) + theme(legend.position = "bottom")
}

# ==============================================================================
# 05. Sensitivity Analysis
# ==============================================================================

#' Run Sensitivity Analysis (Original vs Diffuse Priors)
run_sensitivity_calc <- function(dataset, stan_model_obj, exp_list, aoi_list, seed = 1234, iter = 2000) {
  
  sens_results <- tibble()
  imp_list <- unique(dataset$infID)
  print(">>> Starting Full-Scale Prior Sensitivity Analysis...")
  
  for(exp in exp_list) {
    for(aoi in aoi_list) {
      short_name <- substr(aoi, 5, 6); mouse_var <- paste0("Mouse", short_name)
      for(imp in imp_list) {
        label <- paste(exp, aoi, imp, sep="_")
        print(paste("Processing:", label))
        
        tryCatch({
          TDt <- dataset %>% filter(Valp >= 0.6, expName == exp, infID == imp) %>%
            dplyr::select(participant, picname, all_of(aoi), all_of(mouse_var))
          tarDat <- TDt %>% pivot_longer(cols = 3:4, names_to = "OB", values_to = "weight")
          Y_vec <- tarDat$weight; N_obs <- length(Y_vec)
          if(N_obs < 10) next 
          
          # Stan Data
          data_common <- list(N=N_obs, D=2, S=length(unique(tarDat$participant)), P=length(unique(tarDat$picname)),
                              SID=as.numeric(as.factor(tarDat$participant)), PID=as.numeric(as.factor(tarDat$picname)),
                              Y=Y_vec, X=tibble(a=1, EV=as.numeric(as.factor(tarDat$OB))))
          data_zib <- c(data_common, list(Y_binary=as.integer(Y_vec>0), N_pos=sum(Y_vec>0), pos_idx=which(Y_vec>0)))
          
          # Priors
          data_orig <- c(data_zib, list(prior_beta_sd=5, prior_phi_sd=100, prior_sd_scale=2.5))
          data_diff <- c(data_zib, list(prior_beta_sd=10, prior_phi_sd=200, prior_sd_scale=10))
          
          # Fit
          fit_orig <- sampling(stan_model_obj, data=data_orig, seed=seed, chains=2, cores=2, iter=iter, refresh=0)
          fit_diff <- sampling(stan_model_obj, data=data_diff, seed=seed, chains=2, cores=2, iter=iter, refresh=0)
          
          get_stats <- function(fit) {
            post <- rstan::extract(fit)
            list(mean_bern=mean(post$bz[,1,2]), mean_beta=mean(post$bz[,2,2]), max_rhat=max(rstan::summary(fit)$summary[,"Rhat"], na.rm=T))
          }
          s_orig <- get_stats(fit_orig); s_diff <- get_stats(fit_diff)
          
          sens_results <- bind_rows(sens_results, tibble(
            Exp=exp, AOI=aoi, Imp=imp, Parameter=c("Bernoulli Slope", "Beta Slope"),
            Estimate_Original=c(s_orig$mean_bern, s_orig$mean_beta), Rhat_Original=s_orig$max_rhat,
            Estimate_Diffuse=c(s_diff$mean_bern, s_diff$mean_beta), Rhat_Diffuse=s_diff$max_rhat
          ))
          rm(fit_orig, fit_diff); gc()
        }, error = function(e) message(paste("Error:", label)))
      }
    }
  }
  return(sens_results)
}

plot_sensitivity_check <- function(sens_results) {
  ggplot(sens_results, aes(x = Estimate_Original, y = Estimate_Diffuse)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(aes(color = AOI), alpha = 0.6, size = 2) +
    facet_wrap(~Parameter, scales = "free") +
    labs(
      title = "Prior Sensitivity Analysis",
      subtitle = "Bernoulli estimates show regularization in extreme regimes; Beta estimates remain robust.",
      x = "Weakly Informative Priors", y = "Diffuse Priors"
    ) + theme_bw(base_size = 14) + theme(legend.position = "bottom")
}

# ==============================================================================
# 06. General Diagnostics & Final Plotting
# ==============================================================================

#' Universal ZIB Analysis Pipeline for Single Condition
run_zib_analysis <- function(dataset, stan_model_obj, 
                             exp_str = "S8", condition_col = "infID", condition_val = "Ext", 
                             target_y = "GazeMO", target_y_pair = NULL, predictors = NULL, 
                             valp_thresh = 0.6, seed = 1234, chains = 4, iter = 2000) {
  
  mode <- ifelse(is.null(predictors), "Category_Compare", "Regression")
  print(paste(">>> Analysis Mode:", mode))
  
  base_data <- dataset %>% filter(Valp >= valp_thresh, expName == exp_str) %>% filter(!!sym(condition_col) == condition_val)
  if(nrow(base_data) == 0) stop("No data found.")
  
  if(mode == "Category_Compare") {
    if(is.null(target_y_pair)) stop("Provide 'target_y_pair' for Comparison Mode.")
    tarDat <- base_data %>% dplyr::select(participant, picname, all_of(target_y_pair)) %>%
      pivot_longer(cols = all_of(target_y_pair), names_to = "Predictor_Label", values_to = "weight")
    Y_vec <- tarDat$weight
    X_matrix <- tibble(Intercept=1, Slope=as.numeric(as.factor(tarDat$Predictor_Label)))
    pred_names <- c("Intercept", "Gaze_vs_Mouse")
  } else {
    tarDat <- base_data %>% dplyr::select(participant, picname, all_of(target_y), all_of(predictors))
    Y_vec <- tarDat[[target_y]]
    X_matrix <- tibble(Intercept=1, tarDat %>% dplyr::select(all_of(predictors)))
    pred_names <- c("Intercept", predictors)
  }
  
  N_obs <- length(Y_vec)
  Y_vec_safe <- Y_vec; Y_vec_safe[Y_vec_safe >= 0.9999] <- 0.9999
  Y_binary <- as.integer(Y_vec_safe > 0); pos_idx <- which(Y_vec_safe > 0)
  
  datalist <- list(N=N_obs, D=ncol(X_matrix), S=length(unique(tarDat$participant)), P=length(unique(tarDat$picname)),
                   SID=as.numeric(as.factor(tarDat$participant)), PID=as.numeric(as.factor(tarDat$picname)),
                   Y=Y_vec_safe, X=X_matrix, Y_binary=Y_binary, N_pos=length(pos_idx), pos_idx=pos_idx)
  
  print(">>> Running Stan Sampling...")
  fit <- sampling(stan_model_obj, data = datalist, seed = seed, chains = chains, cores = chains, iter = iter, refresh = 0)
  
  print(">>> Extracting Results...")
  post <- rstan::extract(fit)
  res_list <- list()
  for(i in 1:length(pred_names)) {
    p_name <- pred_names[i]
    s_bern <- get_mcmc_stats(post$bz[, 1, i])
    s_beta <- get_mcmc_stats(post$bz[, 2, i])
    
    res_list[[length(res_list)+1]] <- tibble(Process="Decision", Predictor=p_name, Mean=s_bern["mean"], Lower=s_bern["HDI_low"], Upper=s_bern["HDI_high"])
    res_list[[length(res_list)+1]] <- tibble(Process="Intensity", Predictor=p_name, Mean=s_beta["mean"], Lower=s_beta["HDI_low"], Upper=s_beta["HDI_high"])
  }
  
  res_df <- bind_rows(res_list) %>% mutate(Is_Significant = ifelse(Lower > 0 | Upper < 0, "Significant", "Not Significant"), Exp = exp_str, Condition = condition_val)
  return(list(fit = fit, results = res_df, mode = mode, meta = list(predictors = pred_names)))
}

plot_zib_diagnostics <- function(analysis_obj) {
  fit <- analysis_obj$fit
  fit_sum <- rstan::summary(fit)$summary
  rhats <- fit_sum[, "Rhat"]; neff <- fit_sum[, "n_eff"] / (fit@sim$iter * fit@sim$chains / 2)
  
  # Regex to keep relevant parameters
  pars_regex <- c("bz", "phi", "s_", "z_", "r_"); all_pars <- rownames(fit_sum)
  kept <- unlist(lapply(pars_regex, function(x) grep(paste0("^", x), all_pars, value = TRUE)))
  
  p_rhat <- mcmc_rhat(rhats[kept]) + labs(title = "A. Convergence (R-hat)") + theme_bw()
  p_neff <- mcmc_neff(neff[kept]) + labs(title = "B. Sampling Efficiency") + theme_bw()
  
  # Traceplots for first 2 predictors
  trace_pars <- c("bz[1,2]", "bz[2,2]", "phi")
  p_trace <- mcmc_trace(fit, pars = trace_pars[trace_pars %in% names(fit)], facet_args = list(ncol = 1)) + labs(title = "C. Traceplots") + theme_bw()
  
  (p_rhat / p_neff) | p_trace
}

plot_zib_forest <- function(analysis_obj) {
  res_df <- analysis_obj$results %>% filter(Predictor != "Intercept")
  ggplot(res_df, aes(y = Predictor, x = Mean, xmin = Lower, xmax = Upper, color = Is_Significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(height = 0.2, linewidth = 1) + geom_point(size = 4) +
    facet_wrap(~Process, scales = "free_x") +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "#E74C3C")) +
    labs(title = paste("ZIB Effects:", unique(res_df$Exp), "|", unique(res_df$Condition)), x = "Effect Size (Log-Odds)") +
    theme_bw(base_size = 14) + theme(legend.position = "bottom")
}
