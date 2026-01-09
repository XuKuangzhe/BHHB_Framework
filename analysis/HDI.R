empHDI <- function(dat, confidence = 0.95) {
  n <- length(dat)
  sorted_dat <- sort(dat)
  ci_idx_inc <- ceiling(confidence * n)
  n_intervals <- n - ci_idx_inc
  # width[i] = sorted[i + gap] - sorted[i]
  lower_idxs <- 1:n_intervals
  upper_idxs <- lower_idxs + ci_idx_inc
  
  widths <- sorted_dat[upper_idxs] - sorted_dat[lower_idxs]
  
  min_idx <- which.min(widths)
  
  # reback HDI
  return(c(HDI_low = sorted_dat[min_idx], 
           HDI_high = sorted_dat[min_idx + ci_idx_inc]))
}

summary_MCMC <- function(dat, confidence = 0.95) {
  s_mean <- mean(dat)
  s_median <- median(dat)
  s_sd <- sd(dat)
  
  hdi <- empHDI(dat, confidence)
  
  return(c(mean = s_mean, 
           median = s_median, 
           sd = s_sd, 
           HDI_low = hdi[["HDI_low"]], 
           HDI_high = hdi[["HDI_high"]]))
}
