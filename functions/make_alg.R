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
