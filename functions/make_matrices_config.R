# Documentation
make_matrices_config_help <- function(){
  cat("r(
## FUNCTION: make_matrices_config
# This function can be used to create the mxAlgebra objects necessary for a configural CFA model.
# It uses make_alg wrapper functions for mxALgebraFromString that are also part of the set of helper functions
# stored in the `parmod` repository. The functions all rely on previously specified baseline matrices and 
# deviation matrices for each of the subsequent primary studies.

# By providing the number of latent factors (n_lv) in your model, this function generates (mostly, mxAlgebra) objects for
# the following objects:

# When n_lv == 1
#   matT (linear algebra for Tau vector of indicator intercepts)
#   matL (linear algebra for Lambda matrix of indicator loadings)
#   matK (linear algebra for Kappa matrix of factor means)
#   matE (log-linear algebra for Theta matrix of indicator residual variances)
#   matP (log-linear algebra for Phi matrix of factor variance)

# When n_lv > 1
#   matT (linear algebra for Tau vector of indicator intercepts)
#   matL (linear algebra for Lambda matrix of indicator loadings)
#   matK (linear algebra for Kappa matrix of factor means)
#   matE (log-linear algebra for Theta matrix of indicator residual variances)
#   matP is now made up from the following series of algebraic formulas
#   matVar (log-linear algebra for factor variance matrix)
#   matR (Fisher' z-transformation of factor covariances) 
#   matIa (n_lv by n_lv identity matrix)
#   matIb (n_lv by n_lv matrix with zero on diagonal and one in off-diagonal cells)
#   matCov ( (matIa * sqrt(matVar)) %*% matR %*% (matIa * sqrt(matVar) )
#   matP ( matIa * matVar + matIb * matCov )

# See publication [reference] for more background.
)")
  }

# Function to create the matrix algebra objects for a configural model with either 1 or more than 1 latent factor 
make_matrices_config <- function(n_lv = get("n_lv", envir = .GlobalEnv), k = get("k", envir = .GlobalEnv)) {
  # source functions for making mxAlgebra objects from char snippets
  source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_alg.R")

  # Helper: assign and return (so you can track what’s created)
  assign_and_return <- function(obj) {
    assign(obj$name, obj, envir = .GlobalEnv)
    obj
  }
  
  # --- Core matrices: created for all models ---
  matT <- assign_and_return(make_alg_linear("T", "dT", k = k))
  matL <- assign_and_return(make_alg_linear("L", "dL", k = k))
  matK <- assign_and_return(make_alg_linear("K", "dK", k = k))
  matE <- assign_and_return(make_alg_log("E", "dE", matrix_name = "matE", k = k))
  
  # --- Conditional parts ---
  if (n_lv == 1) {
    
    matP <- assign_and_return(make_alg_log("P", "dP", matrix_name = "matP", k = k))
    
    invisible(list(
      matT = matT,
      matL = matL,
      matK = matK,
      matE = matE,
      matP = matP
    ))
    
  } else {
    
    matVar <- assign_and_return(make_alg_log("P", "dP", matrix_name = "matVar", k = k))
    matR   <- assign_and_return(make_alg_z("P", "dP", matrix_name = "matR", k = k))
    
    # Identity matrices
    matIa <- assign_and_return(mxMatrix(
      type = "Diag", nrow = n_lv, ncol = n_lv, values = 1, name = "matIa"))
    matIb <- assign_and_return(mxMatrix(
      type = "Full", nrow = n_lv, ncol = n_lv,
      values = as.numeric(!(diag(nrow=n_lv))), name = "matIb"))
    
    # Latent covariance algebra
    matCov <- assign_and_return(mxAlgebra(
      expression = (matIa * sqrt(matVar)) %*% matR %*% (matIa * sqrt(matVar)),
      name = "matCov"))
    
    # Common factor covariance matrix
    matP <- assign_and_return(mxAlgebra(
      expression = matIa * matVar + matIb * matCov,
      name = "matP"))
    
    invisible(list(
      matT = matT,
      matL = matL,
      matK = matK,
      matE = matE,
      matVar = matVar,
      matR = matR,
      matIa = matIa,
      matIb = matIb,
      matCov = matCov,
      matP = matP
    ))
  }
}
