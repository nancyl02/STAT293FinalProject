# --------------------------------------------------
# Reshape List to 3D Array
# Estimated B Array [d, d, K] (Only works for p=1)
# --------------------------------------------------
multivar_B <- function(est_multivar_B, d, K){
  
  est_multivar_B <- array(unlist(est_multivar_B), dim = c(d, d, K))
  
  return(est_multivar_B)
}