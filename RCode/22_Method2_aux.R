# --------------------------------------------------
# Construct the VAR(p) Design Matrix
# Inputs:
#   X: The raw time series matrix [T x d]
#   p: Number of lags
# --------------------------------------------------
gen_design_mat <- function(X, p = 1) {
  
  t <- nrow(X) # Total time points (T)
  
  # Create Target Matrix Y (Current time t) 
  Y_mat <- X[(p + 1):t, ]
  
  # Create Predictor Matrix X_lags (Past times)
  # Stack [ X_{t-1} | X_{t-2} | ... | X_{t-p} ] horizontally.
  lag_list <- vector("list", p)
  
  for (lag in 1:p) {
    # If Y ranges from [p+1] to [t],
    # then X_{t-lag} ranges from [p+1-lag] to [t-lag].
    lag_list[[lag]] <- X[(p + 1 - lag):(t - lag), ]
  }
  
  X_lags <- do.call(cbind, lag_list)
  
  Z_mat <- cbind(Y_mat, X_lags)
  
  return(Z_mat)
  
}


# --------------------------------------------------
# Extract Estimated Coefficients from LASSO Results
# Output Dimension: [d] x [d * p]
# Structure: [ Phi_1 | Phi_2 | ... | Phi_p ] (Horizontally stacked)
# --------------------------------------------------
extract_B_hat <- function(fit_list, d, p) {
  
  # Initialize an empty matrix to hold coefficients
  # Rows = Equations (Variables)
  # Cols = All predictors across all lags
  B_hat <- matrix(0, nrow = d, ncol = d * p)
  
  # Iterate through each variable j (each equation)
  for (j in 1:d) {
    
    model <- fit_list[[j]]
    
    # Extract coefficients at the optimal lambda "lambda.min"
    raw_coefs <- coef(model, s = "lambda.min")
    
    # Remove Intercept
    beta_vec <- as.vector(raw_coefs)[-1] 
    
    # Safety Check: Ensure the length matches our expectation
    if (length(beta_vec) == d * p) {
      B_hat[j, ] <- beta_vec
    } else {
      warning(paste("Coefficient length mismatch for variable", j))
    }
  }
  
  return(B_hat)
  
}



# --------------------------------------------------
# Reshape Flat Matrix to 3D Array
# True B Array [d, d, p]
# --------------------------------------------------
reshape_to_array <- function(B_flat, d, p) {
  
  B_array <- array(0, dim = c(d, d, p))
  
  for (lag in 1:p) {
    # Calculate column indices for this specific lag in the flat matrix
    col_idx <- ((lag - 1) * d + 1) : (lag * d)
    
    # Slice the columns and store them in the 3rd dimension of the array
    B_array[,,lag] <- B_flat[, col_idx]
  }
  
  return(B_array)
}
