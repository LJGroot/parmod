# # Source all helper functions in this folder at once

# adding dummy variables
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/add_dummies.R")

# Factor structure mask
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/factor_mask.R")

# make mxAlgebra objects from string, a) linear, b) log-linear, and c) fisher's Z transformed
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_alg.R")

# make final matrix objects for factor model with configural invariance (fully moderated)
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_matrices_config.R")
