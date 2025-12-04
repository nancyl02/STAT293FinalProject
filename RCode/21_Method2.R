# --------------------------------------------------
# Fit Lasso for Individual VAR(p)
# --------------------------------------------------
fit_lasso_var_p <- function(X, p = 1, lambda_ = NULL) {
  
  d <- ncol(X)
  
  Z_mat <- gen_design_mat(X, p)
  
  # The first d columns are Y_t (Current)
  Y_mat <- Z_mat[, 1:d]
  
  # The remaining columns are the stacked lags (Past)
  X_lags <- Z_mat[, -c(1:d)]
  
  # Fit LASSO for each variable j
  fit_list <- vector("list", d)
  
  for (j in 1:d) {
    fit_list[[j]] <- cv.glmnet(
      x = X_lags,
      y = Y_mat[, j],
      alpha = 1,       # alpha=1 specifies the LASSO penalty (L1 norm)
      lambda = lambda_, # If NULL, glmnet selects its own sequence
      nfolds = 5
    )
  }
  
  est_B_list <- extract_B_hat(fit_list, d, p)
  est_B_array <- reshape_to_array(est_B_list, d, p)
  
  return(est_B_array)
  
}