# --------------------------------------------------
# Generate Coefficient Matrix (Can be used for common or unique)
# --------------------------------------------------
gen_coef_matirx <- function(d, idx, min_val = phi_min, max_val = phi_max){
  
  M <- matrix(0, d, d)
  n_nonzero <- length(idx)
  M[idx] <- runif(n_nonzero, min_val, max_val)
  # M[idx] <- runif(n_nonzero, min_val, max_val) * sample(c(-1,1), n_nonzero, replace = TRUE)
  
  return(M)
}


# --------------------------------------------------
# Determine Stability
# --------------------------------------------------
check_phi_stability <- function(phi_array) {
  
  # Require input dimensions to be [row, col, lag]
  # Ensure input is a 3D array [d, d, p].
  if (!is.array(phi_array)) {
    phi_array <- array(phi_array, dim = c(d, d, 1))
  }
  
  d <- dim(phi_array)[1]
  p <- dim(phi_array)[3]
  
  # Convert 3D Array slices into a List
  phi_list <- lapply(1:p, function(i) phi_array[,,i])
  
  # --- the Companion Matrix (F) ---
  # The structure is a (p*d) x (p*d) matrix:
  # F = [ Phi1  Phi2 ... Phip ]
  #     [  I     0   ...   0  ]
  #     [  0     I   ...   0  ]
  #     ...
  
  # 1. Construct Companion Matrix (F)
  # Create Top Block: [Phi1, Phi2, ..., Phip]
  top_block <- do.call(cbind, phi_list)
  
  if (p > 1) {
    # Identity matrix on the sub-diagonal
    ident_block <- diag(1, d * (p - 1))
    # Zero column on the right
    zero_col <- matrix(0, nrow = d * (p - 1), ncol = d)
    bottom_block <- cbind(ident_block, zero_col)
    
    # Combine
    companion_matrix <- rbind(top_block, bottom_block)
  } else {
    # VAR(1) case: Companion matrix is just Phi1 itself
    companion_matrix <- top_block
  }
  
  # 2. Calculate Eigenvalues & Determine Stability
  eigen_values <- eigen(companion_matrix)$values
  max_eigen <- max(Mod(eigen_values))
  
  return(list(
    is_stable = max_eigen < 1,
    max_eigen = max_eigen
  ))
}


# --------------------------------------------------
# Generate Var(p) ILD for each individual
# --------------------------------------------------
gen_VAR_p <- function(phi_array, T_len, p, d, sigma_ep, burn_in = 50) {
  
  # Require input dimensions to be [row, col, lag]
  # Ensure input is a 3D array [d, d, p].
  if (!is.array(phi_array)) {
    phi_array <- array(phi_array, dim = c(d, d, 1))
  }
  p <- dim(phi_array)[3]
  
  # Total Simulation Length (Desired Length + Burn-in)
  total_T <- T_len + burn_in
  
  # Initialize
  X <- matrix(0, nrow = total_T, ncol = d)
  
  # Initialize the first 'p' steps to start the process
  X[1:p, ] <- matrix(rnorm(p * d), nrow = p, ncol = d)
  
  # Generate Data Iteratively
  for (t in (p + 1):total_T) {
    
    # Calculate the Autoregressive component: Sum(Phi_p * X_{t-p})
    ar_term <- numeric(d) # Initialize zero vector
    
    for (lag in 1:p) {
      # Retrieve coefficient matrix for this lag
      phi_p <- phi_array[,,lag]
      
      # Retrieve past data at t-lag
      past_X <- X[t - lag, ]
      
      # Matrix multiplication and accumulation
      ar_term <- ar_term + (phi_p %*% past_X)
    }
    
    # Add random error term
    X[t, ] <- ar_term + rnorm(d, 0, sigma_ep)
  }
  
  # Remove Burn-in Period
  # Return only the valid data sequence of length T_len
  return(X[(burn_in + 1):total_T, ])
}


# --------------------------------------------------
# Generate individual_B (4D array [d, d, p, K])
# --------------------------------------------------
gen_individual_B <- function(K, d, p, n_common, n_unique, max_attempts = 100){
  
  # ==== Generate Common Coefficients Matrix and Indices ==== #
  
  # Position Indices
  # Dimensions: [Number of Unique Elements, Subject, Lag]
  # This stores the specific grid locations for each subject's unique effects
  unique_indices_array <- array(0, dim = c(n_unique, K, p))
  
  # Common Coefficients Matrix (mu)
  # Dimensions: [Row, Col, Lag]
  # This stores the shared coefficient matrices for the population
  common_mu_array <- array(0, dim = c(d, d, p))
  
  # Individual Matrix (B)
  # Dimensions: [Row, Col, Lag, Subject]
  # This will store the final stable transition matrices for every subject
  individual_B <- array(0, dim = c(d, d, p, K))
  
  
  for (lag in 1:p) {
    
    # 1. Generate Common Coefficients Matrix
    # Randomly select grid points for the common structure
    common_idx <- sample(1:(d*d), n_common)
    common_mu_array[,,lag] <- gen_coef_matirx(d, idx = common_idx)
    
    # 2. Allocate Unique Locations
    available_pool <- setdiff(1:(d*d), common_idx)
    
    # Sample indices for all K subjects at once from the remaining pool
    unique_idx <- sample(available_pool, n_unique * K, replace = FALSE)
    
    # Reshape into a matrix where each column represents one subject
    unique_indices_array[,,lag] <- matrix(unique_idx, nrow = n_unique, ncol = K)
    
  }
  
  # ==== Generation Individual Matrix and Stability Check ==== #
  
  for (k in 1:K) {
    
    is_stable <- FALSE
    attempt <- 0
    max_attempts <- 100
    
    while(!is_stable) {
      
      attempt <- attempt + 1
      
      # Temporary container for this subject's p matrices
      current_subject_phi <- array(0, dim = c(d, d, p))
      
      for (lag in 1:p) {
        # Retrieve Common Matrix (mu)
        common_mu <- common_mu_array[,,lag]
        
        # Looks up the specific locations pre-allocated for subject k at this lag
        unique_idx <- unique_indices_array[, k, lag]
        
        # Generate Unique Delta Matrix
        unique_Delta <- gen_coef_matirx(d, idx = unique_idx)
        
        # Phi Matrix at this lag
        current_subject_phi[,,lag] <- common_mu + unique_Delta
      }
      
      # Stability Check
      check_result <- check_phi_stability(current_subject_phi)
      
      if (check_result$is_stable) {
        # Success
        is_stable <- TRUE
        individual_B[,,,k] <- current_subject_phi
        
      } else {
        # Failure: Repeat 'while' loop and regenerate unique_Delta
        # Note: 'unique_idx' and 'common_mu' remain fixed.
        
        # Fallback mechanism if stability is hard to achieve
        if (attempt > max_attempts) {
          
          # Calculate Scaling Factor
          # The current max eigenvalue is > 1 (unstable).
          # Formula: scaling_factor = current_eigen / target_radius
          target_radius <- 0.95
          scaling_factor <- check_result$max_eigen / target_radius
          
          # Apply Scaling
          current_subject_phi <- current_subject_phi / scaling_factor
          
          # Force Stability Flag
          is_stable <- TRUE
          individual_B[,,,k] <- current_subject_phi 
          
          warning(sprintf("Subject %d unstable after %d attempts. Forced scaling applied (Factor=%.2f).", k, max_attempts, scaling_factor))
        }
      }
    }
  }
  
  return(individual_B)
  
}

