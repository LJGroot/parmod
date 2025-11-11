factor_mask_help <- function(){
  cat(r"(
## FUNCTION: factor_mask()
# This function can be used to create a boolean mask matrix for factor loadings
# A boolean mask matrix is a matrix that reflects the factor structure of the model
# that can be used as input for the `free` argument of the mxMatrix objects for the
# Lambda matrices (or delta Lambda matrices).

# Args:
#   n_indicators: integer, number of observed indicators (rows)
#   n_factors: integer, number of latent factors (columns)
#   assign: "balanced" (default) or numeric vector length n_indicators with factor indices (1..n_factors)

# Returns:
#   logical matrix n_indicators x n_factors (TRUE = freely estimated)

# # examples:
# single-factor structure
# factor_mask(4, 1) 
#     F1
# y1 TRUE
# y2 TRUE
# y3 TRUE
# y4 TRUE

# two-factor structure
# factor_mask(6, 2)
#       F1    F2
# y1  TRUE FALSE
# y2  TRUE FALSE
# y3  TRUE FALSE
# y4 FALSE  TRUE
# y5 FALSE  TRUE
# y6 FALSE  TRUE


# For assign != "balanced", provide a vector that represents factor sructure as follows
# c(rep(
#  "factor number 1", "number of indicators", 
# through 
#  "factor number n_xi", "number of indicators"
# )
# where n_xi is the total number of latent factors.
# For example, a three-factor structure with 3, 5, and 2 indicators respectively:

# factor_mask(10, 3, assign = c(rep(1,3), rep(2,5), rep(3,2)))
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
)")
  
factor_mask <- function(n_indicators, n_factors, assign = "balanced"){
  # Sanity check for user input
  # Stop if there is less than one indicator
  if(!is.numeric(n_indicators) || n_indicators <= 0) stop("n_indicators must be positive integer")
  # Stop if there is less than one latent factor
  if(!is.numeric(n_factors) || n_factors <= 0) stop("n_factors must be positive integer")
  
  # Store input from arguments as integers 
  n_indicators <- as.integer(n_indicators); n_factors <- as.integer(n_factors)

  # Parse `assign` argument
  # determine assignment of each indicator to a factor
  if(is.character(assign) && assign == "balanced"){
  # assign indicators consecutively to factors as evenly as possible
  base <- floor(n_indicators / n_factors) # floor = largest even group possible
  extra <- n_indicators %% n_factors # modulus (rest after taking floor)
  sizes <- rep(base, n_factors) # devide base over factors
  if(extra > 0) sizes[1:extra] <- sizes[1:extra] + 1 # add rest to factors until none left
  inds <- rep(seq_len(n_factors), times = sizes) # repeat each factor index according to computed sizes
  if(length(inds) != n_indicators) stop("assignment error") # sanity check: total assigned equals number of indicators
  } else if(is.numeric(assign)){
  if(length(assign) != n_indicators) stop("assign vector must have length n_indicators") # ensure provided vector matches indicator count
  if(any(assign < 1 | assign > n_factors)) stop("assign values must be between 1 and n_factors") # validate factor indices are in range
  inds <- as.integer(assign) # coerce numeric vector to integer factor indices
  } else stop("assign must be 'balanced' or numeric vector of factor indices") # invalid assign argument type or value

  # Create mask matrix
  mask <- matrix(FALSE, nrow = n_indicators, ncol = n_factors)
  for(i in seq_len(n_indicators)){
    mask[i, inds[i]] <- TRUE
  }

  # add row and col namese
  rownames(mask) <- paste0("y", seq_len(n_indicators))
  colnames(mask) <- paste0("F", seq_len(n_factors))
  return(mask)
}
