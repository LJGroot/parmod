# # Source all helper functions in this folder at once
# adding dummy variables
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/add_dummies.R")

# Factor structure mask
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/factor_mask.R")

# make mxAlgebra objects from string, a) linear, b) log-linear, and c) fisher's Z transformed
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_alg.R")

# make final matrix objects for factor model with configural invariance (fully moderated)
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_matrices_config.R")

# loading message
cat("
Functions:
  add_dummies
  factor_mask
  make_alg_linear
  make_alg_log
  make_alg_z
  make_matrices_config
loaded into workspace.

For instructions on use, run: 
  factor_mask_help(); make_alg_help(); make_dummies_help()
")
