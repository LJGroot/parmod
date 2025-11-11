make_alg_help <- function(){
  help_text <- c("
### Algebra helper functions
# These functions can be used to help create mxAlgebra objects by providing wrappers for mxAlgebraFromString.

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


## FUNCTION: make_alg_linear()
# This function can be used to make mx algebras from string with desired base symbols and deviation letters
# for linear algebra

# arguments: 
  # base_symbol: this is the base symbol for the matrix, e.g., T for tau vectors, L for lambda matrices
  # matrix_letter: this is the name for the deviation matrices, i suggest d for delta followed by the base symbol (e.g., dL)
  # k: number of studies in IPD, used to obtain number of elements necessary in algebraic formulation

## FUNCTION: make_alg_log()
# This function can be used to make mx algebras from string with desired base symbols and deviation letters
# for log-transformed algebra

# arguments: 
  # base_symbol: this is the base symbol for the matrix, e.g., T for tau vectors, L for lambda matrices
  # matrix_letter: this is the name for the deviation matrices, i suggest d for delta followed by the base symbol (e.g., dL)
  # k: number of studies in IPD, used to obtain number of elements necessary in algebraic formulation
  # matrix_name: internal name for mx object, can be specified separately when necessary 

## FUNCTION: make_alg_z()
# This function can be used to make mx algebras from string with desired base symbols and deviation letters
# for fisher z transformation

# arguments: 
  # base_symbol: this is the base symbol for the matrix, e.g., T for tau vectors, L for lambda matrices
  # matrix_letter: this is the name for the deviation matrices, i suggest d for delta followed by the base symbol (e.g., dL)
  # k: number of studies in IPD, used to obtain number of elements necessary in algebraic formulation
  # matrix_name: internal name for mx object, can be specified separately when necessary 
");print(help_text)
  }
  
  
## FUNCTIONS
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
