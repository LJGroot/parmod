# This function can be used to add dummy variables to you IPD data set. 
  
# It assumes that your IPD is stored as raw data in long format with a study identifier (preferably labeled as "study").
# The function has the following arguments:
# df = provide data frame
# study_col = provide column that holds study ID, default is "study"
# ref = reference group, first study by default
# prefix = preferred prefix for names of dummy variables, defaults to "d"
# name_by_value = use numbers 1 through k as dummy suffix (TRUE), 
#                 or use original study ID as suffix (FALSE)
  
add_dummies <- function(df, 
                        study_col = "study", 
                        ref = NULL, 
                        prefix = "d", 
                        name_by_value = FALSE) {
  # error for wrong column name
  if (!study_col %in% names(df)) {
    stop("Column '", study_col, "' not found in data frame.")
  }
  # full-length vector of values (not unique!)
  studies_vec <- as.character(df[[study_col]])
  # unique values as vector
  unique_studies <- unique(studies_vec)
    # error for when column has no values
  if (length(unique_studies) == 0) {
    stop("Study column contains no values")
  }
  # if ref argument is null, take first unique entry as ref
  # if ref is provided but not present, give error
  if (is.null(ref)) {
    ref <- unique_studies[1]
  } else if (!ref %in% unique_studies) {
    stop("Reference '", ref, "' not found among study values:",
         paste(head(unique_studies, 10), collapse = ", "))
  }
  # take all levels minus reference group
  levels_to_dummy <- unique_studies[unique_studies != ref]
  # throw error if remaining levels is zero (only 1 unique study in data)
  if (length(levels_to_dummy) == 0) {
    message("only the reference level '", ref, "' present - no dummy columns created.")
    return(df)
  }
  # generate matrix of dummies
  dummy_mat <- vapply(levels_to_dummy,
                      FUN = function(val) as.integer(studies_vec == val),
                      FUN.VALUE = integer(length(studies_vec)),
                      USE.NAMES = FALSE)
  # convert to data.frame
  dummy_df <- as.data.frame(dummy_mat, stringsAsFactors = FALSE)
  # give names to dummies using provided of default prefix, and using
  # numbers or original study identifiers
  if (name_by_value) {
    colnames(dummy_df) <- paste0(prefix, make.names(levels_to_dummy))
  } else {
    colnames(dummy_df) <- paste0(prefix, seq_along(levels_to_dummy))
  }
  
  out <- df
  out[names(dummy_df)] <- dummy_df
  out
}

# Create a boolean mask matrix for factor loadings

# Args:
#   n_indicators: integer, number of observed indicators (rows)
#   n_factors: integer, number of latent factors (columns)
#   assign: "balanced" (default) or numeric vector length n_indicators with factor indices (1..n_factors)

# Returns:
#   logical matrix n_indicators x n_factors (TRUE = freely estimated)

factor_mask <- function(n_indicators, n_factors, assign = "balanced"){
  if(!is.numeric(n_indicators) || n_indicators <= 0) stop("n_indicators must be positive integer")
  if(!is.numeric(n_factors) || n_factors <= 0) stop("n_factors must be positive integer")
  n_indicators <- as.integer(n_indicators); n_factors <- as.integer(n_factors)

  # determine assignment of each indicator to a factor
  if(is.character(assign) && assign == "balanced"){
    # assign indicators consecutively to factors as evenly as possible
    base <- floor(n_indicators / n_factors)
    extra <- n_indicators %% n_factors
    sizes <- rep(base, n_factors)
    if(extra > 0) sizes[1:extra] <- sizes[1:extra] + 1
    inds <- rep(seq_len(n_factors), times = sizes)
    if(length(inds) != n_indicators) stop("assignment error")
  } else if(is.numeric(assign)){
    if(length(assign) != n_indicators) stop("assign vector must have length n_indicators")
    if(any(assign < 1 | assign > n_factors)) stop("assign values must be between 1 and n_factors")
    inds <- as.integer(assign)
  } else stop("assign must be 'balanced' or numeric vector of factor indices")

  mask <- matrix(FALSE, nrow = n_indicators, ncol = n_factors)
  for(i in seq_len(n_indicators)){
    mask[i, inds[i]] <- TRUE
  }
  rownames(mask) <- paste0("y", seq_len(n_indicators))
  colnames(mask) <- paste0("F", seq_len(n_factors))
  return(mask)
}

# # examples:
# factor_mask(4, 1) creates 
#     F1
# y1 TRUE
# y2 TRUE
# y3 TRUE
# y4 TRUE

# factor_mask(6, 2)
#       F1    F2
# y1  TRUE FALSE
# y2  TRUE FALSE
# y3  TRUE FALSE
# y4 FALSE  TRUE
# y5 FALSE  TRUE
# y6 FALSE  TRUE

#factor_mask(10, 3, assign = c(rep(1,3), rep(2,5), rep(3,2)))
#        F1    F2    F3
# y1   TRUE FALSE FALSE
# y2   TRUE FALSE FALSE
# y3   TRUE FALSE FALSE
# y4  FALSE  TRUE FALSE
# y5  FALSE  TRUE FALSE
# y6  FALSE  TRUE FALSE
# y7  FALSE  TRUE FALSE
# y8  FALSE  TRUE FALSE
# y9  FALSE FALSE  TRUE
# y10 FALSE FALSE  TRUE

# Functions to make mx algebras from string with desired base symbols and deviation letters
# for linear algebra
make_alg_linear <- function(base_symbol,   # base symbol (eg. T for tau vector)
                            matrix_letter, # symbol for deviation vector/matrix
                            k              # number of studies
                            ){
  # matrix expression from string elements
  # baseline object
  mat_expr <- paste0("mat", base_symbol, "0")
  # subsequent objects multiplied by dummy
  for (i in 1:(k - 1)) {
    mat_expr[i + 1] <- paste0("mat", matrix_letter, i, "*", dummies[i])
  }
  # collapse into single expression with plus symbols between elements
  mat_expr <- paste(mat_expr, collapse = " + ")
  # store as mx object with name from base_symbol
  mat <- mxAlgebraFromString( name = paste0("mat", base_symbol), mat_expr )
  return(mat) 
  }

# for log transformed algebra
make_alg_log <- function(base_symbol,    # base symbol (eg. T for tau vector)
                         matrix_letter,  # symbol for deviation vector/matrix
                         k,              # number of studies
                         matrix_name     # matrix name (not symbol)
                         ){
  # matrix expression from string elements
  expr <- paste0("mat", base_symbol, "0", "*exp(")
  for (i in 1:(k - 1)) {
    expr[i + 1] <- paste0("mat", matrix_letter, i, "*", dummies[i]) }
  mat_expr <- paste(expr[1], paste(expr[2:length(expr)], collapse = " + "), ")")
  # store as mx object with name from matrix_name arg
  mat <- mxAlgebraFromString( name = matrix_name, mat_expr )
  return(mat) 
  }

# for fisher z transformation
make_alg_z <- function(base_symbol,    # base symbol (eg. T for tau vector)
                       matrix_letter,  # symbol for deviation vector/matrix
                       k,              # number of studies
                       matrix_name     # matrix name (not symbol)
                       ){
  # matrix expression from string elements
  expr <- paste0("(exp( 2 * (matP0 +")
  for (i in 1:(k - 1)) {
    expr[i + 1] <- paste0("mat", matrix_letter, i, "*", dummies[i]) }
  mat_expr <- paste(
    expr[1],
    paste(expr[2:length(expr)], collapse = " + "),
    ") ) - 1) / ",
    expr[1],
    paste(expr[2:length(expr)], collapse = " + "),
    ") ) + 1)" )
  # store as mx object with name from matrix_name arg
  mat <- mxAlgebraFromString( name = matrix_name, mat_expr )
  return(mat) 
  }


# # For example, the following syntax uses mxAlgebraFromString in the background 
# matR <- make_alg_z("P", "dP", matrix_name = "matR", k = 7)

# # which gives the same output as the following manual specification of an mxAlgebra object
# matR <- mxAlgebra(expression =
#                     (exp( 2 * (
#                         matP0 +
#                           matdP1*d1 + matdP2*d2 + matdP3*d3 +
#                           matdP4*d4 + matdP5*d5 + matdP6*d6 ) ) - 1) /
#                     (exp( 2 * (
#                         matP0 +
#                           matdP1*d1 + matdP2*d2 + matdP3*d3 +
#                           matdP4*d4 + matdP5*d5 +matdP6*d6 ) ) + 1),
#                   name = "matR")

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
      values = if (n_lv == 2) c(0, 1, 1, 0) else 0, name = "matIb"))
    
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
