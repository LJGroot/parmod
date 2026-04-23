# Header ----
##
## Script name:
## IPD MASEM via Parameter Moderation with a path model
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
# This script only applies to the situation where the hypothesized model is a path model.
# When you want to fit a factor model or a full SEM model, please refer back to the
# repository and select the correct script that corresponds to your hypothesized model.

## Load packages
### check if pacman package is installed, install if not
if (!"pacman" %in% installed.packages()) install.packages("pacman"); library(pacman)

### Load rest of packages
pacman::p_load(osfr, haven, psych, readr, readxl, lavaan, metaSEM, OpenMx)

# Load helper functions from repository (please read provided documentation).
source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/helper_functions.R")

## This script uses PM as an approach for IPD MASEM when the hypothesized
## model is a path model.
##
## We provide example data for a sample of x studies. This is in fact data from
##
## If you want to use your own data, make sure that the data is in a long-form aggregated
## format (vertically), and there is a study identifier variable present. Then load your
## file instead of the example data, and extract from that data the following objects
## like is done below for the example data:
##  vars = indicator names
##  n_var = number of indicators (can be obtained through length(vars))
##  k = number of primary studies
##  n_dum = number of dummy variables, equals k-1
##  dummies = dummy variable names, set via paste0("d", 1:n_dum)
##            this scripts uses prefix "d" followed by numbers 1 through n_dum.
##            can be customized, but check consistency with rest of script.

## Read data ----
## Example data
## From local folder (after downloading, change file path if necessary)
integrate.df <- readRDS("integrate.RData")
## Load data from OSF repository
data_file <- osf_retrieve_file("9znd5") |>
  osf_download(conflicts = "overwrite"); load(data_file$local_path)
## Or Github
load(url("https://github.com/LJGroot/parmod/raw/refs/heads/main/data/integrate/integrate.RData"))
# Or use your own data

# Data prep for example data (skip if using own data)
# single dummy for control or intervention condition
integrate.df <- within(integrate.df, {
  tx_int <- as.numeric(tx!="Ctrl")
});integrate.df$study <- as.character(integrate.df$study)


# Use sourced function `add_dummies` to create data frame `dat` with dummy vars
# Change study_col to match the study identifier from your data if necessary
# For example data, the identifier is `study`.
dat <- add_dummies(df = integrate.df, study_col = "study")

# Store data characteristics
# Variable names (insert the names of your observed variables)
vars <- c("pbs0", "alcprob0", "tx_int", "pbs1", "alcprob1")
# Number of observed variables/indicators p
n_var <- length(vars)
# Number of primary studies k
k <- length(unique(dat$study))
# Number of dummy variables k-1
n_dum <- k - 1
# dummy variable names
dummies <- paste0("d", 1:n_dum)

## When data is loaded, objects are stored, and dummy variables are made we can proceed
## to store the data as an mxData object to be used in OpenMx. We select all the rows
## and the columns which include the observed variables for the model, and those with
## the newly created dummy variables.
path_data <- mxData(dat[,c(vars, dummies)], type = "raw")

## Unconstrained model (all parameters moderated by study membership) ----
### Intercepts ----
# Intercept vector A for reference study
matA0 <- mxMatrix(
  type = "Full", nrow = n_var, ncol = 1,
  labels = vars, free = TRUE, values = 1,
  name = "matA0")
# Deviations of intercepts for subsequent studies
dAmats <-  setNames(lapply(1:n_dum, function(i) {
  mxMatrix(
    type = "Full", nrow = n_var, ncol = 1, free = TRUE,
    name = paste0("matdA", i)
  )}), paste0("matdA", 1:n_dum))

### Path Coefficients ----
# matrix B reference study
matB0 <- mxMatrix(
  type = "Full",  nrow = n_var,  ncol = n_var,
  dimnames = list(vars, vars), name = "matB0")

# Setting paths to reflect your path model
# by setting cells with paths in free to TRUE

# For example data
# paths to pbs1
matB0$free["pbs1",                         # arrows going TO
           c("pbs0", "alcprob0", "tx_int") # arrows coming FROM
           ] <- TRUE
# Paths into alcprob1
matB0$free["alcprob1",                             # arrows going TO
           c("pbs0", "alcprob0", "pbs1", "tx_int") # arrows coming FROM
           ] <- TRUE

# deviation matrices for the other studies
dBmats <-  setNames(lapply(1:n_dum, function(i) {
  mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, dimnames = list(vars, vars),
    free = matB0$free, name = paste0("matdB", i)
  )}), paste0("matdB", 1:n_dum))


### Variance-Covariance ----
# Variance-covariance matrix for reference group
matP0 <- mxMatrix(
  type = "Symm", nrow = n_var, dimnames = list(vars, vars), name = "matP0")
# Covariances (exogenous vars) in off-diagonal cells free with starting values 0.5
matP0$free[c("pbs0", "alcprob0", "tx_int"),
           c("pbs0", "alcprob0", "tx_int")] <- TRUE
matP0$values[c("pbs0", "alcprob0", "tx_int"),
             c("pbs0", "alcprob0", "tx_int")] <- 0.5
# Variances (all vars) on diagonal free with starting value 1
diag(matP0$free) <- TRUE
diag(matP0$values) <- 1
dPmats <-  setNames(lapply(1:n_dum, function(i) {
  mxMatrix(
    type = "Symm", nrow = n_var, dimnames = list(vars, vars),
    free = matP0$free, name = paste0("matdP", i)
  )}), paste0("matdP", 1:n_dum))

### Dummy vars ----
# Define definition variables
Dmats <- setNames(lapply(1:n_dum, function(i) {
  mxMatrix(
    type = "Full",    nrow = 1,    ncol = 1,
    labels = paste0("data.d", i),
    name = paste0("d", i)
  )}), paste0("d", 1:n_dum))

### Algebra ----
# alpha column vector
matA <- make_alg_linear(base_symbol = "A", matrix_letter = "dA", k = k)
# algebra for beta matrix
matB <- make_alg_linear(base_symbol = "B", matrix_letter = "dB", k = k)
# psi matrix
matP <- make_alg_linear(base_symbol = "P", matrix_letter = "dP", k = k)

matI <- mxMatrix(
  type = "Iden", nrow = n_var, name = "matI")

### model implied matrices ----
matS <- mxAlgebra(
  expression =
    solve(matI - matB) %*% matP %*% t(solve(matI - matB)),
  name = "matS")
matM <- mxAlgebra(
  expression =
    t(solve(matI - matB) %*% matA),
  name = "matM")

### Expectation and fit function----
exp <- mxExpectationNormal(covariance = "matS",
                           means = "matM",
                           dimnames = vars)
funML <- mxFitFunctionML()

### Build and run model ----
modUnconstr <- mxModel(
  "Unconstrained Model",
  matA, matA0, dAmats,
  matB, matB0, dBmats,
  matP, matP0, dPmats,
  Dmats,
  matS, matM, matI,
  exp, funML,
  path_data)
fitUnconstr <- mxTryHard(modUnconstr)
summary(fitUnconstr, fitIndices = T)

## Constrained path coefficients ----
matB <- mxMatrix(
  type = "Full", nrow = n_var, ncol = n_var, dimnames = list(vars, vars),
  free = matB0$free, name = "matB")

### Build and Run model ----
modConstrB <- mxModel(
  "Constrained B model", # Model name
  matA, matA0, dAmats,   # Intercepts, moderated
  matB,                  # Path coefficients, constrained
  matP, matP0, dPmats,   # (Co)variances, moderated
  Dmats,                 # Dummy variables
  matS, matM, matI,      # Model implied matrices
  exp, funML,            # Expectation and fit function
  path_data)             # IPD
fitConstrB <- mxTryHard(modConstrB)
summary(fitConstrB)

mxCompare(fitUnconstr, fitConstrB)

## Moderation of intercepts ---------
# Both path coefficients and (co)variances are constrained in this model,
# only intercepts are moderated

# Variance-covariance matrix for reference group
matP <- mxMatrix(
  type = "Symm", nrow = n_var, dimnames = list(vars, vars),
  free = matP0$free, values = matP0$values, name = "matP")

### Build and Run model ----
modConstrBP <- mxModel(
  "Constrained B and P model",
  matA, matA0, dAmats,
  matB,
  matP,
  Dmats,
  matS, matM, matI,
  exp, funML,
  path_data)
fitConstrBP <- mxTryHard(modConstrBP, extraTries = 25)
summary(fitConstrBP)

mxCompare(fitUnconstr, c(fitConstrB, fitConstrBP))


