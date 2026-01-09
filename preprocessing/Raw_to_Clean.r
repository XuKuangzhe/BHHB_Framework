# ==============================================================================
# Script: Raw Data Preprocessing & Feature Extraction
# Purpose: 
#   1. Synchronize timestamps between PsychoPy (Stimulus) and Tobii (Gaze).
#   2. Convert raw gaze/mouse coordinates into heatmaps using Gaussian smoothing.
#   3. Calculate weighted density scores for specific AOIs (Glabella, Eyes, etc.).
# Author: Kuangzhe Xu
# Note: Raw data files are required to run this script. Paths must be adjusted.
# ==============================================================================

library(tidyverse)
library(openxlsx)
library(progress)
library(data.table)
library(mmand) # For Gaussian smoothing

# ==============================================================================
# 01. Documentation: Tobii Validity Codes (Tobii Fusion 120)
# ==============================================================================
# Reference: https://connect.tobii.com/s/article/What-Do-Validity-Codes-Mean
# 0: The system is certain that it has recorded all relevant data for a specific eye,
#    and that the data belongs to that specific eye (No risk of confusion).
# 1: The system has only recorded one eye and has made assumptions/estimations 
#    about whether it is the left or right eye. The estimation is likely correct.
# 2: The system has recorded one eye but cannot determine if it is left or right.
# 3/4: Eye not found or signal lost.

# ==============================================================================
# 02. Load Metadata & Masks
# ==============================================================================
# Load AOI Mask Data
# Note: Ensure the path points to your local mask file
load("/masklist.Rdata")
maskN <- mask_list %>% names()

# Load and Process Personality Data (Facesheet)
perdt <- read_csv("/ParFacesheet.csv", show_col_types = FALSE)

subSta <- perdt %>%
  mutate(
    subjID = paste0("Subj", c(1:84)),
    gender = str_replace_all(gender, c("男"="Male", "女"="Female")),
    # Calculate Big Five Traits (Reverse scoring applied where necessary)
    Ext = (Q1 + (8 - Q6)) / 2,
    Agr = ((8 - Q2) + Q7) / 2,
    Con = (Q3 + (8 - Q8)) / 2,
    Neu = (Q4 + (8 - Q9)) / 2,
    Ope = (Q5 + (8 - Q10)) / 2
  ) %>%
  dplyr::select(-c(Email, Name, Q1:Q10))

# ==============================================================================
# 03. Define File Paths
# ==============================================================================
# Path to raw data folder (Modify as needed)
data_root <- "/datacontr"

tsv_path <- list.files(path = data_root, pattern = "*.tsv", full.names = TRUE)
csv_path <- list.files(path = data_root, pattern = "*.csv", full.names = TRUE)

# ==============================================================================
# 04. Main Processing Loop
# ==============================================================================

# Initialize container for summarized data
sumGDT <- c()

# Setup progress bar
pb <- progress_bar$new(format = "Processing [:bar] :percent eta: :eta", total = 84 * 50)

print(">>> Starting Data Preprocessing...")

for(n in 1:84){
  
  # Load raw files for current subject
  dfcsv <- read_csv(csv_path[n], show_col_types = FALSE)
  dftsv <- fread(tsv_path[n])
  dftsv$TimeStamp <- as.numeric(dftsv$TimeStamp)
  
  # ----------------------------------------------------------------------------
  # A. Timestamp Synchronization (PsychoPy <-> Tobii)
  # ----------------------------------------------------------------------------
  # Extract PsychoPy stimulus start times
  vx <- dfcsv$imageStimulus.started %>% na.omit() %>% as.vector()
  
  # Extract unique Event IDs from Tobii stream
  pidn <- dftsv$Event %>% unique() %>% na.omit()
  pidm <- dfcsv$picname %>% na.omit()
  
  # Match events to build synchronization model
  tvdf <- c()
  for(p in 1:50){
    tv <- dftsv %>% filter(Event == pidn[p+1])
    if(nrow(tv) > 0){
      tvd <- tibble(vX = vx[p], vY = as.numeric(tv$TimeStamp[1]))
      tvdf <- rbind(tvdf, tvd)
    }
  }
  
  # Build Linear Model for Time Correction
  lmv <- lm(vY ~ vX, tvdf)
  
  # ----------------------------------------------------------------------------
  # B. Coordinate Processing & Feature Extraction
  # ----------------------------------------------------------------------------
  # Select relevant behavioral columns
  td <- dfcsv %>% 
    dplyr::select(participant, expName, picname, 
                  key_resp_2.keys, key_resp_confi.keys, mouse.x, mouse.y) %>% 
    na.omit()
  
  for(v in 1:50){
    pb$tick()
    
    # --- 1. Process Mouse Coordinates ---
    # Parse string "[x,y]", rescale to screen resolution (1920x1080)
    MX <- round(td$mouse.x[v] %>% str_replace_all("\\[","") %>% str_replace_all("]","") %>% 
                  str_split(",") %>% as_vector() %>% as.numeric() * 1200 + 960, 0)
    MY <- round(td$mouse.y[v] %>% str_replace_all("\\[","") %>% str_replace_all("]","") %>% 
                  str_split(",") %>% as_vector() %>% as.numeric() * 1080 + 540, 0)
    
    # Retrieve subject personality
    perDT <- subSta %>% filter(subjID == td$participant[1])
    
    # --- 2. Process Gaze Coordinates ---
    # Check if image was presented before (using pre-calculated checkfile 'chesum' logic if needed)
    # Here we assume standard flow using the synchronization model
    lmdf <- lmv$coefficients %>% as_tibble()
    
    # Predict Tobii Start/End times based on PsychoPy time
    stT <- lmdf$value[1] + lmdf$value[2] * vx[v]
    enT <- lmdf$value[1] + lmdf$value[2] * vx[v] + 3500 # 3500ms duration
    
    # Filter Valid Gaze Data (Validity code < 3 means reliable)
    tardf <- dftsv %>% 
      filter(TimeStamp > stT & TimeStamp < enT,
             ValidityLeft < 3 | ValidityRight < 3, 
             GazePointX != "NaN")
    
    # Calculate Data Quality (Validity Percentage)
    validperc <- nrow(tardf) / nrow(dftsv %>% filter(TimeStamp > stT & TimeStamp < enT))
    
    # --- 3. Generate Heatmaps (Gaussian Smoothing) ---
    # Rescale Gaze to Screen Resolution
    GX <- round(tardf$GazePointX * 1200 + 960, 0)
    GY <- round(tardf$GazePointY * 1080 + 540, 0)
    
    # Initialize empty matrices for heatmaps (1080p)
    hitdata_g <- hitdata_m <- matrix(0, 1080, 1920)
    
    # Map Gaze Points
    n_r <- length(GX)
    if(n_r > 0){
      for(i in 1:n_r){
        if(GY[i] > 0 & GY[i] < 1080 & GX[i] > 0 & GX[i] < 1920){
          hitdata_g[GY[i], GX[i]] = hitdata_g[GY[i], GX[i]] + 1
        }
      }
    }
    
    # Map Mouse Points
    n_r <- length(MX)
    if(n_r > 0){
      for(i in 1:n_r){
        if(MY[i] > 0 & MY[i] < 1080 & MX[i] > 0 & MX[i] < 1920){
          hitdata_m[MY[i], MX[i]] = hitdata_m[MY[i], MX[i]] + 1
        }
      }
    }
    
    # Apply Gaussian Smoothing to create density maps
    heatmap_g <- gaussianSmooth(hitdata_g, c(10, 10))
    heatmap_m <- gaussianSmooth(hitdata_m, c(10, 10))
    
    # --- 4. Calculate AOI Weights ---
    # Identify the mask index corresponding to the current picture
    num <- which(str_detect(maskN, str_sub(pidm[v], 1, 5)))
    
    # Normalize and Calculate Density for each AOI
    # Formula: Sum(Heatmap * Mask) / Sum(Mask Area)
    matchdat <- td[v, 1:5] %>%
      mutate(perDT[, 2:8], 
             Valp = validperc,
             
             # Gaze Densities
             GazeGL = sum(heatmap_g * mask_list[[num[1]]]) / sum(mask_list[[num[1]]]), # Glabella
             GazeFH = sum(heatmap_g * mask_list[[num[2]]]) / sum(mask_list[[num[2]]]), # Forehead
             GazeEB = (sum(heatmap_g * mask_list[[num[4]]]) + sum(heatmap_g * mask_list[[num[8]]])) /
               (sum(mask_list[[num[4]]]) + sum(mask_list[[num[8]]])),           # Eyebrow (L+R)
             GazeEY = (sum(heatmap_g * mask_list[[num[3]]]) + sum(heatmap_g * mask_list[[num[7]]])) /
               (sum(mask_list[[num[3]]]) + sum(mask_list[[num[7]]])),           # Eye (L+R)
             GazeNO = sum(heatmap_g * mask_list[[num[6]]]) / sum(mask_list[[num[6]]]), # Nose
             GazeMO = sum(heatmap_g * mask_list[[num[5]]]) / sum(mask_list[[num[5]]]), # Mouth
             
             # Mouse Densities
             MouseGL = sum(heatmap_m * mask_list[[num[1]]]) / sum(mask_list[[num[1]]]),
             MouseFH = sum(heatmap_m * mask_list[[num[2]]]) / sum(mask_list[[num[2]]]),
             MouseEB = (sum(heatmap_m * mask_list[[num[4]]]) + sum(heatmap_m * mask_list[[num[8]]])) /
               (sum(mask_list[[num[4]]]) + sum(mask_list[[num[8]]])),
             MouseEY = (sum(heatmap_m * mask_list[[num[3]]]) + sum(heatmap_m * mask_list[[num[7]]])) /
               (sum(mask_list[[num[3]]]) + sum(mask_list[[num[7]]])),
             MouseNO = sum(heatmap_m * mask_list[[num[6]]]) / sum(mask_list[[num[6]]]),
             MouseMO = sum(heatmap_m * mask_list[[num[5]]]) / sum(mask_list[[num[5]]])
      )
    
    sumGDT <- rbind(sumGDT, matchdat)
    
    # Slight pause to update progress bar smoothly
    Sys.sleep(1/84*50)
  }
}

# ==============================================================================
# 05. Save Final Dataset
# ==============================================================================
# write_csv(sumGDT, "/sumGDT.csv")
print(">>> Preprocessing Complete.")

