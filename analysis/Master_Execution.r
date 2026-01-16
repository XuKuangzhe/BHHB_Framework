################################################################################
# Project: BHHB Framework
# Script: Master Analysis Execution Pipeline
# Date: 2026-01
################################################################################

# ==============================================================================
# 00. Setup & Dependencies
# ==============================================================================
# Clean workspace
rm(list = ls())

# Load Libraries
library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(ggrepel)
library(gridExtra)
library(scoringRules)

# Parallel Processing Settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# --- Set Paths (Modify these for your local environment) ---
# Tip for BRM: Using relative paths or the 'here' package is recommended
path_code <- "./"          # Location of R scripts
path_data <- "./data/"     # Location of data files
path_model <- "./models/"  # Location of Stan files
path_out  <- "./results/"  # Output directory

# Create output dir if not exists
if(!dir.exists(path_out)) dir.create(path_out)

# --- Source Helper Functions ---
# Note: Ensure HB_Functions.R is the one we defined previously
source(file.path(path_code, 'HDI.R')) 
source(file.path(path_code, 'HB_Functions.R'))

# ==============================================================================
# 01. Data Loading & Preprocessing
# ==============================================================================
print(">>> Loading Data...")
# Adjust file name as needed
sumGDT <- read_csv(file.path(path_data, "sumGDT.csv"), show_col_types = FALSE)

threshold <- 1e-5

# Preprocessing: Extract infID and expName
sumGDTs <- sumGDT %>%
  mutate(
    infID = str_sub(picname, 1, 3),
    expName = str_sub(expName, 13, 15) # Adjust indices based on your string format
  )%>%mutate(across(starts_with("Gaze"), ~ ifelse(. < threshold, 0, .)),
             across(starts_with("Mouse"), ~ ifelse(. < threshold, 0, .)))

# Inspection
print(head(sumGDTs))

# ==============================================================================
# 02. Model Compilation
# ==============================================================================
print(">>> Compiling Stan Models (This may take a minute)...")

# Main HB Model
HBmodvef  <- stan_model(file.path(path_model, "HB_main.stan"))

# Competitor Models (for PPC & Comparison)
LMMmod     <- stan_model(file.path(path_model, "LMM_baseline.stan"))
BetaMod    <- stan_model(file.path(path_model, "Beta_squeezed.stan"))

# Sensitivity Analysis Model (Allows custom priors)
model_sens <- stan_model(file.path(path_model, "HB_sensitivity.stan"))

# Bundle for passing to functions
models_list <- list(HB = HBmodvef, LMM = LMMmod, Beta = BetaMod)

# ==============================================================================
# 03. Parameter Recovery Check (Simulation)
# ==============================================================================
print(">>> Running Section 03: Parameter Recovery...")

# Run Simulation
sim_results_df <- run_grad_sim_calc(
  stan_model_obj = HBmodvef, 
  target_props = c(0.10, 0.30, 0.60, 0.90),
  iter = 2000
)

# Plot Results
p_sim <- plot_grad_sim(sim_results_df)
print(p_sim)

# Save
# ggsave(file.path(path_out, "Appendix_Gradient_Simulation.png"), p_sim, width = 12, height = 6)


# ==============================================================================
# 04. The 4-Area Plot (Empirical Analysis)
# ==============================================================================
print(">>> Running Section 04: Empirical 4-Area Plot...")

# --- Experiment 1 (S8) ---
res_s8 <- run_empirical_hb_calc(
  dataset = sumGDTs, 
  exp_label = "S8", 
  stan_model_obj = HBmodvef
)
# Optional: Save intermediate results
# write_csv(res_s8$params, file.path(path_out, "HB_S8_Params.csv"))

# --- Experiment 2 (S12) ---
res_s12 <- run_empirical_hb_calc(
  dataset = sumGDTs, 
  exp_label = "S12", 
  stan_model_obj = HBmodvef
)
# write_csv(res_s12$params, file.path(path_out, "HB_S12_Params.csv"))

# --- Plotting ---
# Single Plot Example (S8)
p_4area_s8 <- plot_four_area_single(res_s8$params)
print(p_4area_s8)
p_4area_s12 <- plot_four_area_single(res_s12$params)
print(p_4area_s12)

# Combined Plot (If you have the wrapper for both, otherwise plot individually)
p_4area_combined <- plot_four_area_sum(res_s8$params, res_s12$params) 
print(p_4area_combined)

# Save
#ggsave(file.path(path_out, "Fig4_Four_Area_Plot.png"), p_4area_combined, width = 10)


# ==============================================================================
# 05. Posterior Predictive Checks (PPC)
# ==============================================================================
print(">>> Running Section 05: PPC...")

# Target Condition: S12, Ext, Mouth
ppc_data <- run_ppc_calc(
  dataset = sumGDTs, 
  exp_str = "S12", 
  imp_str = "Ext", 
  aoi_str = "GazeMO", 
  mouse_str = "MouseMO",
  models_list = models_list,
  iter = 2000 
)

# Plot
p_ppc <- plot_ppc_compare(ppc_data)
print(p_ppc)

# ggsave(file.path(path_out, "Fig3_PPC_Method_Breakdown.png"), p_ppc, width = 12, height = 8)


# ==============================================================================
# 06. Model Performance Scoring (Universal Comparison)
# ==============================================================================
print(">>> Running Section 06: Model Scoring (CRPS/RMSE/MAE)...")

# Calculate Metrics (Time Consuming!)
perf_results <- run_performance_calc(
  dataset = sumGDTs, 
  models_list = models_list, 
  iter = 2000
)

# write_csv(perf_results, file.path(path_out, "Universal_Model_Comparison_Metrics.csv"))

perf_results%>%group_by(Model)%>%
  summarise(MeanCRPS=mean(CRPS),MeanRMSE=mean(RMSE),MeanMAE=mean(MAE))%>%
  arrange(MeanCRPS)

# Plot Trends
p_rmse <- plot_rmse_trend(perf_results)
p_mae  <- plot_mae_trend(perf_results)
p_crps <- plot_crps_trend(perf_results)

# Arrange and Print
grid.arrange(p_rmse, p_mae, p_crps, ncol=1)

# ggsave(file.path(path_out, "Final_CRPS_Trend.png"), p_crps, width = 10, height = 6)


# ==============================================================================
# 07. Sensitivity Analysis
# ==============================================================================
print(">>> Running Section 07: Prior Sensitivity Analysis...")

# Run Check (Time Consuming!)
sens_data <- run_sensitivity_calc(
  dataset = sumGDTs, 
  stan_model_obj = model_sens, 
  exp_list = c("S8", "S12"),
  aoi_list = c("GazeGL", "GazeFH", "GazeEB", "GazeEY", "GazeNO", "GazeMO"),
  iter = 2000
)

# write_csv(sens_data, file.path(path_out, "Full_Sensitivity_With_Diagnostics.csv"))

# Plot
p_sens <- plot_sensitivity_check(sens_data)
print(p_sens)

# ggsave(file.path(path_out, "Appendix_Full_Sensitivity.png"), p_sens, width = 8, height = 5)


# ==============================================================================
# 08. Detailed Diagnostics & Forest Plots
# ==============================================================================
print(">>> Running Section 08: Diagnostics...")

# --- Mode A: Category Comparison (Gaze vs Mouse) ---
res_cat <- run_hb_analysis(
  dataset = sumGDTs,
  stan_model_obj = HBmodvef,
  exp_str = "S12",
  condition_col = "infID", condition_val = "Ext",
  target_y_pair = c("GazeMO", "MouseMO"), 
  predictors = NULL # NULL triggers comparison mode
)

print(plot_hb_diagnostics(res_cat))
print(plot_hb_forest(res_cat))

# --- Mode B: Personality Regression ---
res_reg <- run_hb_analysis(
  dataset = sumGDTs,
  stan_model_obj = HBmodvef,
  exp_str = "S12",
  condition_col = "infID", condition_val = "Agr",
  target_y = "GazeMO",
  predictors = c("Ext", "Agr", "Con", "Neu", "Ope") # Regression mode
)

# Specific Traceplots
print(plot_hb_diagnostics(res_reg, trace_pars=c("bz[1,2]", "bz[2,2]","bz[2,4]", "phi")))
print(plot_hb_forest(res_reg))

print(">>> Analysis Complete.")
