# --------------------------------------------------
# Estimation Errors
# --------------------------------------------------
diff_func <- function(est_B, B_true){
  
  sum((est_B - B_true)^2)

}


# --------------------------------------------------
# Sensitivity and Specificity
# --------------------------------------------------
sens_spec_func <- function(B_true, B_est, thr = 1e-4) {
  true_nonzero <- abs(B_true) != 0
  est_nonzero <- abs(B_est) > thr
  
  TP <- sum(est_nonzero & true_nonzero)
  FN <- sum(!est_nonzero & true_nonzero)
  TN <- sum(!est_nonzero & !true_nonzero)
  FP <- sum(est_nonzero & !true_nonzero)
  
  if ((TP + FN) > 0) {
    sensitivity <- TP / (TP + FN)
  } else {
    sensitivity <- NA 
  }
  
  if ((TN + FP) > 0) {
    specificity <- TN / (TN + FP)
  } else {
    specificity <- NA
  }
  
  return(list(
    sensitivity = sensitivity, 
    specificity = specificity,
    
    TP = TP, FP = FP, TN = TN, FN = FN 
  ))
  
}


# --------------------------------------------------
# Calculate Mean Sensitivity and Specificity Metrics Across All Subjects
# Inputs:
#   B_true : 3/4 Array [d, d, (p), K] containing true values
#   B_est  : 3/4D Array [d, d, (p), K] containing estimates
#   thr       : Threshold (default 1e-4)
#   output_type
# --------------------------------------------------
calc_metrics <- function(B_true, B_est, thr = 1e-4,
                         output_type = c("mean", "subject", "both")) {
  
  # output_type = "mean" (default)
  output_type <- match.arg(output_type)
  
  # Check dimensions
  if (!all(dim(B_true) == dim(B_est))) {
    stop("Error: Dimensions of True and Estimated arrays do not match!")
  }
  
  dims <- dim(B_true)
  ndim <- length(dims)
  
  if (ndim == 3) {
    # CASE A: Input is 3D [d, d, K] (Implies VAR(1), p=1)
    d <- dims[1]
    K <- dims[3]
    
    # "Reshape" to 4D by adding a singleton dimension for p=1
    dim(B_true) <- c(d, d, 1, K)
    dim(B_est)  <- c(d, d, 1, K)
    
  } else if (ndim == 4) {
    K <- dims[4]
    
  } else {
    stop("Error: Input arrays must be 3D [d,d,K] or 4D [d,d,p,K]")
  }
  
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Loop through each subject
  for (k in 1:K) {
    # Slice the 4D array to get the 3D block [d, d, p] for subject k
    # This block contains all lags for this person
    B_true_k <- B_true[, , , k]
    B_est_k  <- B_est[, , , k]
    
    # Calculate metrics for this person
    res <- sens_spec_func(B_true_k, B_est_k, thr = thr)
    
    sens_vec[k] <- res$sensitivity
    spec_vec[k] <- res$specificity
  }
  
  # Return the Mean across all subjects
  if (output_type == "mean") {
    
    return(list(
      mean_sens = mean(sens_vec),
      mean_spec = mean(spec_vec)
    ))
    
  } else if (output_type == "subject") {
    
    return(list(
      sensitivity = sens_vec,
      specificity = spec_vec
    ))
    
  } else { # "both"
    
    return(list(
      summary = list(mean_sens = mean(sens_vec), mean_spec = mean(spec_vec)),
      raw_data = list(sens_vec = sens_vec, spec_vec = spec_vec)
    ))
  }
}


# --------------------------------------------------
# Calculate RMSFE (Root Mean Square Forecast Error)
# --------------------------------------------------
compute_rmsfe_p <- function(B_est, test_data_matrix) {
  
  # B_est: Coefficient array for a single subject.
  #        Expected dimensions: [d, d, p] or [d, d] (if p=1)
  # test_data_matrix: The test data including the initial lags.
  #                   Structure: Rows 1 to p are the "Lags" (from Training end).
  #                              Rows (p+1) to End are the "Targets" (Test data).
  
  # Ensure B_est is a 3D array [d, d, p]
  # If B_est is a 2D matrix (implying p=1), convert it to 3D array [d, d, 1]
  if (is.matrix(B_est) || length(dim(B_est)) == 2) {
    d_dim <- nrow(B_est)
    B_est <- array(B_est, dim = c(d_dim, d_dim, 1))
  }
  
  d <- dim(B_est)[1]      # Number of variables
  p <- dim(B_est)[3]      # Lag order
  n_total <- nrow(test_data_matrix) # Total rows in test set (including lags)
  
  # Define Actuals (Targets)
  # The actual values predicted start from row (p + 1)
  actuals <- test_data_matrix[(p + 1):n_total, ]
  
  # Initialize Prediction Container
  # Size matches the number of actual target observations
  preds <- matrix(0, nrow = nrow(actuals), ncol = ncol(actuals))
  
  # Rolling One-Step Ahead Forecast
  # Predict X_t based on X_{t-1}, ..., X_{t-p}
  # Loop t runs through the indices of the original test_data_matrix
  for (t in (p + 1):n_total) {
    
    # Initialize prediction vector for time t
    pred_t <- numeric(d)
    
    # Accumulate contributions from each lag: Sum(A_lag * X_{t-lag})
    for (lag in 1:p) {
      A_l <- B_est[, , lag]                 # Coefficient matrix for lag 'l'
      x_lag <- test_data_matrix[t - lag, ]  # Data vector at time t-l
      
      pred_t <- pred_t + (A_l %*% x_lag)
    }
    
    # Store prediction
    # Note: preds index starts at 1, corresponding to t-(p)
    preds[t - p, ] <- pred_t
  }
  
  # Calculate Root Mean Squared Forecast Error
  return(sqrt(mean((preds - actuals)^2)))
}



get_rmsfe_p <- function(B_est_obj, test_data_list, K) {
  
  # Initialize vector to store RMSFE for each subject
  rmsfe_vec <- numeric(K)
  
  for (k in 1:K) {
    
    # --- Logic to extract coefficients for Subject k ---
    
    # Case 1: Input is a 4D Array [d, d, p, K] (e.g., Multi-VAR output)
    if (is.array(B_est_obj) && length(dim(B_est_obj)) == 4) {
      B_k <- B_est_obj[, , , k]
      
      # Case 2: Input is a 3D Array [d, d, K] (e.g., Multi-VAR output when p=1)
    } else if (is.array(B_est_obj) && length(dim(B_est_obj)) == 3) {
      # Extract k-th slice and ensure it keeps array structure for compute_rmsfe_p
      B_temp <- B_est_obj[, , k] 
      # Convert to [d, d, 1] so dimension checks pass
      B_k <- array(B_temp, dim = c(nrow(B_temp), ncol(B_temp), 1))
      
      # Case 3: Input is a List (e.g., LASSO output)
    } else if (is.list(B_est_obj)) {
      B_k <- B_est_obj[[k]]
      
    } else {
      # Case 4: Pooled/Shared model (Same matrix for everyone)
      B_k <- B_est_obj
    }
    
    # --- Compute RMSFE ---
    # Pass the extracted coefficients and the test data (which contains lags)
    rmsfe_vec[k] <- compute_rmsfe_p(B_k, test_data_list[[k]])
  }
  
  # Return the average RMSFE across all subjects
  return(rmsfe_vec)
}