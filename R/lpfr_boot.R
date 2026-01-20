#' Helper function to run bootstrapped marginal LPFR models
lpfr_boot <- function(s_l, b, boot_dat, temp_dir, kz, kb, smooth.cov, family,
                      scalar_varnames, ff_varnames, Y_label){
  p.ff <- length(ff_varnames)
  p.scalar <- length(scalar_varnames)

  # Create matrix of covariates
  covariate_matrix <- as.matrix(boot_dat[,scalar_varnames])

  # Get list of functional predictors
  func_list <- lapply(ff_varnames, function(ff_name) as.matrix(boot_dat[[ff_name]]))

  # Fit pointwise LPFR model
  fit.sl <- refund::lpfr(Y = boot_dat[[Y_label]][,s_l], subj = boot_dat[[subj_var]],
                         covariates = covariate_matrix,
                         funcs = func_list,
                         method = "REML", kz = kz,
                         kb = kb, smooth.cov = smooth.cov)

  lpfr_data <- list(s_l, fit.sl$beta.covariates, fit.sl$BetaHat, fit.sl$ranef, b)

  saveRDS(lpfr_data, file = paste(temp_dir, "/lpfr_res_", "sl", s_l, "_b", b, ".RDS", sep = ""))
  return()
}
