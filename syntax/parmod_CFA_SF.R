# Header ----
##
## Script name:
## IPD MASEM via Parameter Moderation with a single-factor model
##
## Author: L.J. Groot
## Copyright (c) Lennert J. Groot, 2025
## Email: l.j.groot@uva.nl
##
## Notes:
##

# Single-Factor Model ----
# This script contains syntax to conduct IPD MASEM via Parameter Moderation (PM)
# This consists of the successive specification of models, that are then compared in
# order to evaluate deterioration of fit. For more background, read the publication
# "Making the Most of Minimal Data: IPD MASEM via Parameter Moderation Using OpenMx"
# (Groot, Kan, & Jak, submitted for review).
# This script only applies to the situation where the hypothesized model is a factor
# model with a single latent factor. When you have more than 1 latent factor, or your
# model is a path model or full SEM model, please refer back to the repository and
# select correct script that corresponds to your hypothesized model.

## Load packages
### check if pacman package is installed, install if not
if (!"pacman" %in% installed.packages()) install.packages("pacman"); library(pacman)

### Load rest of packages
pacman::p_load(osfr, haven, psych, readr, readxl, lavaan, metaSEM, OpenMx)

# Load helper functions from repository (please read provided documentation).
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/helper_functions.R")

## This script uses PM as an approach for IPD MASEM when the
## factor model has a single latent factor. In other words, when the hypothesized
## model is a 1-factor model.
##
## We provide example data for a sample of 7 studies. This is in fact data from the OECD
## PISA educational data. We use data that was gathered using Teachers Questionnaire
## TC227, in a subset of 7 high income studies (Germany, Portugal, Hong Kong, Malaysia,
## Macau, South Korea, and the United Arab Emirates).
##
## If you want to use your own data, make sure that the data is in a long-form aggregated
## format (vertically), and there is a study identifier variable present. Then load your
## file instead of the example data, and extract from that data the following objects
## like is done below for the example data:
##  vars = indicator names
##  n_var = number of indicators (can be obtained through length(vars))
##  n_lv = number of latent factors. Set to 1 in the case of a single-factor model
##         see other script in repository for how to evaluate multi-factor models
##  k = number of primary studies
##  n_dum = number of dummy variables, equals k-1
##  dummies = dummy variable names, set via paste0("d", 1:n_dum)
##            this scripts uses prefix "d" followed by numbers 1 through n_dum.
##            can be customized, but check consistency with rest of script.

## Read data ----
## Example data
## From local folder (after downloading, change file path if necessary)
pisadata <- readRDS("pisadata.rds")
## Load data from OSF repository
data_file <- osf_retrieve_file("g4s8e") |>
  osf_download(conflicts = "overwrite"); pisadata <- readRDS(data_file$local_path)
## Or Github
pisadata <- readRDS(url("https://github.com/LJGroot/parmod/raw/refs/heads/main/data/pisa/pisadata.rds"))
# Or use your own data

# Use sourced function `add_dummies` to create data frame `dat` with dummy vars
# Change study_col to match the study identifier from your data if necessary
# For example data, the identifier is `Country`.
dat <- add_dummies(df = pisadata, study_col = "Country")

# Store data characteristics
# Variable names (observed variables)
vars <- names(dat[2:5])
# Number of observed variables/indicators p
n_var <- length(vars)
# Number of latent variables q
n_lv <- 1
# Number of primary studies k
k <- length(unique(dat$Country))
# Number of dummy variables k-1
n_dum <- k - 1
# dummy variable names
dummies <- paste0("d", 1:n_dum)

## When data is loaded, objects are stored, and dummy variables are made we can proceed
## to store the data as an mxData object to be used in OpenMx. We select all the rows
## and the columns which include the observed variables for the model, and those with
## the newly created dummy variables.
factor_data <- mxData(dat[,c(vars, dummies)], type = "raw")


## Naive Model ----
## Factor Structure
free_lambdas <- factor_mask(n_var, n_lv)
## Matrices
matT <- mxMatrix(
  type = "Full",  nrow = 1,  ncol = n_var,
  free = TRUE,  values = 0.5,  name = "matT")
matL <- mxMatrix(
  type = "Full",  nrow = n_var,  ncol = n_lv,
  free = free_lambdas,
  values = as.numeric(free_lambdas), byrow = TRUE,
  name = "matL")
matE <- mxMatrix(
  type = "Diag", nrow = n_var, free = TRUE, values = 1, name = "matE")
matK <- mxMatrix(
  type = "Full",  nrow = n_lv,  ncol = 1, name = "matK")
matP <- mxMatrix(
  type = "Symm", nrow = n_lv,
  free = !as.logical(diag(nrow = n_lv)),
  values = {v <- matrix(0.5, nrow = n_lv, ncol = n_lv);diag(v) <- 1;v},
  byrow = TRUE, name = "matP")
### Model-implied matrices
matS <- mxAlgebra(
  expression =
    matL %*% matP %*% t(matL) + matE,
  name = "matS")
matM <- mxAlgebra(
  expression =
    matT + t(matL %*% matK),
  name = "matM")
expF <- mxExpectationNormal(covariance = "matS",
                            means = "matM",
                            dimnames = vars)
fitF <- mxFitFunctionML()

modNaive <- mxModel(
  model = "modNaive",   # Model Name
  matT,# Intercepts
  matL, # Loadings
  matE, # Resid vars
  matP,  # Common factor cor
  matK,  # Factor means
  matM,  matS,            # Model algebraic expressions
  expF,                   # Expectation function
  fitF,                   # Fit function
  factor_data )           # IPD
fitNaive <- mxTryHard(modNaive)
summary(fitNaive, fitIndices = TRUE)

## Configural Model --------------------------
# T_0 matrix: baseline intercepts for study 1
matT0 <- mxMatrix(
  type = "Full", nrow = 1, ncol = n_var, free = TRUE, values = 1, name = "matT0")
dTmats <-  setNames(lapply(1:(k - 1), function(i) {
  mxMatrix(
    type = "Full", nrow = 1, ncol = n_var, free = TRUE, name = paste0("matdT", i)
  )}), paste0("matdT", 1:(k - 1)))

# lambda pattern
# this function makes a true/false mask for the pattern of factor loadings to esitmate.
free_lambdas <- factor_mask(n_var, n_lv)
# Lambda_0: baseline factor loadings for reference group
matL0 <- mxMatrix(
  type = "Full",  nrow = n_var,  ncol = n_lv,
  free = free_lambdas,
  values = as.numeric(free_lambdas),
  name = "matL0")
# C1–C6: deviations in loadings for studies 2–7
dLmats <- setNames(lapply(1:(k - 1), function(i) {
  mxMatrix(
    type = "Full", nrow = n_var, ncol = n_lv, free = free_lambdas,
    name = paste0("matdL", i)
  )}), paste0("matdL", 1:(k - 1)))

# E0 Baseline residual variances for reference group
matE0 <- mxMatrix(
  type = "Diag", nrow = n_var, free = TRUE, values = 1, name = "matE0")
dEmats <- setNames(lapply(1:(k - 1), function(i) {
  mxMatrix(
    type = "Diag", nrow = n_var, free = TRUE, name = paste0("matdE", i)
  )}), paste0("matdE", 1:(k - 1)))

# P0 latent variable var/cov matrix for reference group
matP0 <- mxMatrix(
  type = "Symm", nrow = n_lv,
  free = !as.logical(diag(nrow = n_lv)),
  values = {v <- matrix(0.5, nrow = n_lv, ncol = n_lv);diag(v) <- 1;v},
  name = "matP0")
# Deviation of lv cor from baseline in study 2-7
dPmats <- setNames(lapply(1:(k - 1), function(i) {
  mxMatrix(
    type = "Symm", nrow = n_lv, free = matP0$free, name = paste0("matdP", i)
  )}), paste0("matdP", 1:(k - 1)))

# A0 latent variable means in reference group
matK0 <- mxMatrix(
  type = "Full",  nrow = n_lv,  ncol = 1, name = "matK0")
dKmats <- setNames(lapply(1:(k - 1), function(i) {
  mxMatrix(
    type = "Full", nrow = n_lv, ncol = 1, name = paste0("matdK", i)
  )}), paste0("matdK", 1:(k - 1)))

# Define definition variables
Dmats <- list()
for (i in 1:n_dum) {
  Dmats[[i]] <- mxMatrix(
    type = "Full", nrow = 1, ncol = 1, labels = paste0("data.d", i),
    name = paste0("d", i) ) }

# Algebra
# 2 options:
#
# Option 1:
# We need to make matrix algebra objects for matrices tau, Lambda, Kappa, Phi, and Theta.
# We can do this using for each of these matrices separately using the helper functions
# that were sourced in line 31-32. Run make_alg_help() to learn more.
# matT <- make_alg_linear(base_symbol = "T", matrix_letter = "dT", k = k)
# matL <- make_alg_linear(base_symbol = "L", matrix_letter = "dL", k = k)
# matK <- make_alg_linear("K", "dK", k = k)
# matP <- make_alg_log("P", "dP", matrix_name = "matP", k = k)
# matE <- make_alg_log("E","dE", matrix_name = "matE", k = k)

# Option 2:
# Alternatively, we can use the another helper function. This function makes the
# definitive matrix algebra objects necessary for the configural model.
# It does so by obtaining the number of latent variables as stored in the env.
# and makes the appropriate matrices, either without covariances when n_lv == 1, or with
# the appropriate transformation of the common factor covariances when n_lv != 1
make_matrices_config(n_lv)

## Specify model-implied matrices
matS <- mxAlgebra(
  expression =
    matL %*% matP %*% t(matL) + matE,
  name = "matS")
matM <- mxAlgebra(
  expression =
    matT + t(matL %*% matK),
  name = "matM")

## Specify expectation and fit function
expF <- mxExpectationNormal(covariance = "matS", means = "matM", dimnames = vars)
fitF <- mxFitFunctionML()

## Make mxModel object and run the model
## 1 factor model (no common factor covariances)
modConfig <- mxModel(
  model = "Configural",   # Model Name
  matT,  matT0,  dTmats,  # Intercepts
  matL,  matL0,  dLmats,  # Loadings
  matE,  matE0,  dEmats,  # Resid vars
  matK,  matK0,  dKmats,  # Factor means
  matP,  matP0,  dPmats,  # Common factor var-cov
  Dmats,                  # Dummy variables
  matM,  matS,            # Model algebraic expressions
  expF,                   # Expectation function
  fitF,                   # Fit function
  factor_data )           # IPD
fitConfig <- mxTryHard(modConfig)

## Metric Model ----

### Loadings ----
matL <- mxMatrix(
  type = "Full",  nrow = n_var,  ncol = n_lv,
  free = matL0$free,  values = matL0$values,
  name = "matL")

### Factor Covariance ----

# Deviation of lv cor from baseline in study 2-7
for (i in 1:(k-1)) dPmats[[i]]$free <- TRUE

## Make mxModel object and run the model
modMetric <- mxModel(
  model = "Metric",   # Model Name
  matT,  matT0,  dTmats,  # Intercepts
  matL,                   # Loadings
  matE,  matE0,  dEmats,  # Resid vars
  matK,  matK0,  dKmats,  # Factor means
  matP,  matP0,  dPmats,  # Common factor var-cov
  Dmats,                  # Dummy variables
  matM,  matS,            # Model algebraic expressions
  expF,                   # Expectation function
  fitF,                   # Fit function
  factor_data )           # IPD
fitMetric <- mxTryHard(modMetric)
summary(fitMetric, fitIndices = TRUE)

## Scalar model ----

### Intercepts ----
matT <- mxMatrix(
  type = "Full", nrow = 1, ncol = n_var, free = TRUE, values = 1, name = "matT")

### Factor Means ----
for (i in 1:n_dum) dKmats[[i]]$free <- TRUE

### Build and run model ----
modScalar <- mxModel(
  model = "Scalar",   # Model Name
  matT,                   # Intercepts
  matL,                   # Loadings
  matE,  matE0,  dEmats,  # Resid vars
  matP,  matP0,  dPmats,  # Common factor cor
  matK,  matK0,  dKmats,  # Factor means
  Dmats,                  # Dummy variables
  matM,  matS,            # Model algebraic expressions
  expF,                   # Expectation function
  fitF,                   # Fit function
  factor_data )           # IPD
fitScalar <- mxRun(modScalar)
summary(fitScalar, fitIndices = TRUE)
mxCompare(fitConfig, c(fitMetric, fitScalar))
