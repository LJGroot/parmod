# # Source all helper functions in this `functions` folder at once
# adding dummy variables
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/add_dummies.R")
if (!is.function(get0("add_dummies", ifnotfound = NULL))) message(
  "Function 'add_dummies' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/add_dummies.R could not be loaded. Check repo or contact host.")
if (!is.function(get0("add_dummies_help", ifnotfound = NULL))) message(
  "Function 'add_dummies_help' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/add_dummies.R could not be loaded. Check repo or contact host.")

# Factor structure mask
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/factor_mask.R")
if (!is.function(get0("factor_mask", ifnotfound = NULL))) message(
  "Function 'factor_mask' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/factor_mask.R could not be loaded. Check repo or contact host.")
if (!is.function(get0("factor_mask_help", ifnotfound = NULL))) message(
  "Function 'factor_mask_help' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/factor_mask.R could not be loaded. Check repo or contact host.")

# make mxAlgebra objects from string, a) linear, b) log-linear, and c) fisher's Z transformed
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_alg.R")
if (!is.function(get0("make_alg_linear", ifnotfound = NULL))) message(
  "Function 'make_alg_linear' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_alg.R could not be loaded. Check repo or contact host.")
if (!is.function(get0("make_alg_log", ifnotfound = NULL))) message(
  "Function 'make_alg_log' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_alg.R could not be loaded. Check repo or contact host.")
if (!is.function(get0("make_alg_z", ifnotfound = NULL))) message(
  "Function 'make_alg_z' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_alg.R could not be loaded. Check repo or contact host.")
if (!is.function(get0("make_alg_help", ifnotfound = NULL))) message(
  "Function 'make_alg_help' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_alg.R could not be loaded. Check repo or contact host.")

# make final matrix objects for factor model with configural invariance (fully moderated)
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_matrices_config.R")
if (!is.function(get0("make_matrices_config", ifnotfound = NULL))) message(
  "Function 'make_matrices_config' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_matrices_config.R could not be loaded. Check repo or contact host.")
if (!is.function(get0("make_matrices_config_help", ifnotfound = NULL))) message(
  "Function 'make_matrices_config_help' in https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/make_matrices_config.R could not be loaded. Check repo or contact host.")

# loading message
cat("
Functions:
  add_dummies
  factor_mask
  make_alg_linear
  make_alg_log
  make_alg_z
  make_matrices_config
and associated help functions loaded into workspace.

For instructions on use, run: 
  factor_mask_help(); make_alg_help(); add_dummies_help(); make_matrices_config_help()
")
