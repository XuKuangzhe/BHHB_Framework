# Replication Materials: The BHHB Framework

**Authors:** Kuangzhe Xu   
**License:** MIT License  

## 📌 Overview

This repository contains the data, R code, and Stan model specifications for the manuscript *"The BHHB Framework: Decoupling Decision and Intensity in Visual Attention"*.

The core contribution is the **Bayesian Hierarchical Hurdle Beta (BHHB)** framework. Unlike previous Zero-Inflated Beta (ZIB) approaches that assume a mixture process, the BHHB framework adopts a two-part Hurdle assumption. It statistically decouples the decision process (whether to orient) from the intensity process (how long to maintain attention). This distinction provides improved ecological validity for high-frequency eye-tracking data where zeros represent structural non-engagement rather than measurement artifacts.


## 📂 Repository Structure

The project is organized into four main directories:

* **`analysis/`**: Contains the executable R scripts.
    * `Master_Execution.R`: The main pipeline. **Run this script to reproduce results.**
    * `HB_Functions.R`: Custom functions for simulation, HMC sampling, and PPC calculations.
    * `HDI.R`: Utility functions for Highest Density Interval calculations.
* **`models/`**: Contains the Bayesian model specifications written in Stan.
    * `HB_main.stan`: The primary Bayesian Hierarchical Hurdle Beta model.
    * `LMM_baseline.stan` & `Beta_squeezed.stan`: Competitor models for comparison.
    * `HB_sensitivity.stan`: Model for prior sensitivity analysis.
* **`data/`**: Contains the pre-processed dataset.
    * `sumGDT.csv`: Analysis-ready data.
    * `variable_codebook.txt`: Detailed description of all variables.
* **`preprocessing/`**: Documentation of the data cleaning process.
    * `00_Raw_to_Clean.R`: Logic for converting raw Tobii/PsychoPy outputs to the analysis dataset.

## 💻 Computational Environment

The analysis was originally performed and tested on **macOS Tahoe 26.2 (ARM64/Apple Silicon)** using **R version 4.5.1 (2025-06-13)**.

To ensure exact reproducibility, please ensure you have the following key software and packages installed.

### System Requirements
* **R**: Version 4.x.x (Tested on 4.5.1)
* **Rstudio** Version 025.09.0+387 (2025.09.0+387)
* **Xcode (macOS)**: Required for compiling Stan models.

### Key R Dependencies
The analysis relies heavily on the `rstan` ecosystem and the `tidyverse`. Below are the specific versions used in the manuscript generation:

| Category | Package | Version |
| :--- | :--- | :--- |
| **Bayesian Modeling** | `rstan` | 2.32.7 |
| | `bayesplot` | 1.14.0 |
| | `loo` | 2.8.0 |
| **Scoring & Metrics** | `scoringRules` | 1.1.3 |
| **Data Manipulation** | `tidyverse` | 2.0.0 |
| | `dplyr` | 1.1.4 |
| | `readr` | 2.1.5 |
| **Visualization** | `ggplot2` | 3.5.2 |
| | `ggrepel` | 0.9.6 |
| | `gridExtra` | 2.3 |
| | `patchwork` | 1.3.2 |

*Note: A full snapshot of the session information, including all dependencies, is provided in `session_info.txt`.*

## 🚀 Usage Instructions

1.  **Clone or Download** this repository to your local machine.
2.  Open `analysis/Master_Execution.R` in RStudio.
3.  **Set Working Directory**: Ensure your working directory is set to the project root.
    ```r
    # Example inside RStudio
    setwd("/path/to/Project_BHHB")
    ```
4.  **Install Missing Packages**: Run the dependency check at the top of `Master_Execution.R`.
5.  **Run the Pipeline**: Execute the script line-by-line or source the entire file. The script is divided into sections corresponding to the paper's results:
    * Simulation (Parameter Recovery): Verifies parameter recovery across sparsity levels (10%–90% zeros).
    * Empirical Analysis: Generates the "Decision-Intensity State Space"(The 4-Area Plot).
    * Model Comparison: Calculates CRPS/RMSE/MAE comparing BHHB, Logit-LMM, and Freq-HB.
    * Diagnostics: Generates PPC plots for density collapse checks.
    * Another: General Diagnostics / Sensitivity Analysis

## 🛡️ Data Availability & Privacy

**Analysis Data:**
The dataset `data/sumGDT.csv` is provided for reproduction purposes. It contains anonymized, feature-extracted behavioral metrics.

**Originality Statement:** 
This dataset is multimodal. While the eye-tracking component utilizes a previously established experimental design, the mouse-tracking data are original to this study and have not been published elsewhere. The raw coordinate logs are excluded due to size constraints but are available upon reasonable request.

**Raw Data:**
The raw eye-tracking data files (high-frequency coordinate logs) are **not included** in this repository due to their large file size and privacy considerations regarding participant metadata. However, the script `preprocessing/00_Raw_to_Clean.R` is included to transparently demonstrate the logic used to derive the analysis dataset from the raw inputs.

## 📄 Citation

If you use this code or model in your research, please cite the following paper:

> Xu, K. (Under Review). The BHHB Framework: Decoupling Decision and Intensity in Visual Attention. Behavior Research Methods.

## 📞 Contact

For questions regarding the code or data, please create a GitHub Issue or contact the corresponding author.
