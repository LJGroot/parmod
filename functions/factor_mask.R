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
y10 FALSE FALSE  TRUE
