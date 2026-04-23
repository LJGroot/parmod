## Header ----
##
## Script name:
## IPD MASEM via Parameter Moderation Tutorial - Path Model
##
## Author: L.J. Groot
## Copyright (c) Lennert J. Groot, 2025
## Email: l.j.groot@uva.nl
##
## Notes:
##

# IPD MASEM via Parameter Moderation - Path Model Tutorial ----
# This script contains syntax to conduct IPD MASEM via Parameter Moderation (PM). In this
# script, we fit a path model to IPD from 10 primary studies. The data originates from
# [insert publication].
# The steps of this procedure are explained in detail in the publication
# "Making the Most of Minimal Data: IPD MASEM via Parameter Moderation Using OpenMx"
# (Groot, Kan, & Jak, submitted for review).
# Readers can run the syntax to follow along the tutorial. In this tutorial, we chose
# to write code that requires as little prior knowledge of programming in R as possible.
# If, after reading the tutorial, you want to run your own analyses on your own data,
# you can choose one of the additional syntax file in the repository [insert file names].
# These syntax files use code that is more flexible, and easily allows for running the
# analysis steps on different data and a different hypothesized model. Please make sure
# to have read the tutorial carefully and choose the appropriate syntax file when you
# want to do your own analyses.

# Path Model ----
## Libraries and Data ----------
## Load packages
### check if pacman package is installed, install if not
if (!"pacman" %in% installed.packages()) install.packages("pacman"); library(pacman)

### Load rest of packages
pacman::p_load(osfr, haven, psych, readr, readxl, lavaan, metaSEM, OpenMx)

# Load data from online repository
data_file <- osf_retrieve_file("9znd5") |>
  osf_download(conflicts = "overwrite"); load(data_file$local_path)

# single dummy for control or intervention condition
integrate.df <- within(integrate.df, {
  tx_int <- as.numeric(tx!="Ctrl")
})
# Coerce study identifier to character string
integrate.df$study <- as.character(integrate.df$study)
# Dummy vars for study
dat <- cbind.data.frame(
  integrate.df,
  data.frame(
    # st_2 as reference group
    "d1" = as.numeric(integrate.df$study == "8.1"),
    "d2" = as.numeric(integrate.df$study == "8.2"),
    "d3" = as.numeric(integrate.df$study == "8.3"),
    "d4" = as.numeric(integrate.df$study == "9"),
    "d5" = as.numeric(integrate.df$study == "12"),
    "d6" = as.numeric(integrate.df$study == "16"),
    "d7" = as.numeric(integrate.df$study == "18"),
    "d8" = as.numeric(integrate.df$study == "21"),
    "d9" = as.numeric(integrate.df$study == "22")
  ))
vars <- c("pbs0", "alcprob0", "tx_int", "pbs1", "alcprob1")
n_var <- length(vars)
studies <- unique(integrate.df$study)
k <- length(studies)
n_dum <- k - 1
dummies <- paste0("d", 1:n_dum)

# save data as mxData object
path_data <- mxData(dat[,c(vars, dummies)], type = "raw")

## Unconstrained model (all parameters moderated by study membership) ----
### Intercepts ----
# Intercept vector A for reference study
matA0 <- mxMatrix(
  type = "Full", nrow = n_var, ncol = 1,
  labels = vars, free = TRUE, values = 1,
  name = "matA0")
# Deviations of intercepts for subsequent studies
dAmats <- list(
  matdA1 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA1"),
  matdA2 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA2"),
  matdA3 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA3"),
  matdA4 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA4"),
  matdA5 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA5"),
  matdA6 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA6"),
  matdA7 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA7"),
  matdA8 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA8"),
  matdA9 <- mxMatrix(
    type = "Full",  nrow = n_var,  ncol = 1,  free = TRUE, name = "matdA9"))

### Path Coefficients ----
# matrix B reference study
matB0 <- mxMatrix(
  type = "Full",  nrow = n_var,  ncol = n_var,
  dimnames = list(vars, vars), name = "matB0")
# set paths in free to TRUE
# paths to pbs1
matB0$free["pbs1", c("pbs0", "alcprob0", "tx_int")] <- TRUE
# Paths into alcprob1
matB0$free["alcprob1", c("pbs0", "alcprob0", "pbs1", "tx_int")] <- TRUE


# deviation effects for the other studies
dBmats <- list(
  matdB1 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB1"),
  matdB2 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB2"),
  matdB3 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB3"),
  matdB4 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB4"),
  matdB5 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB5"),
  matdB6 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB6"),
  matdB7 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB7"),
  matdB8 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB8"),
  matdB9 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free, name = "matdB9"))

### Variance-Covariance ----
# Variance-covariance matrix for reference group
matP0 <- mxMatrix(
  type = "Symm", nrow = n_var, dimnames = list(vars, vars), name = "matP0")
# Covariances in off-diagonal cells free with starting values 0.5
matP0$free[c("pbs0", "alcprob0", "tx_int"),
           c("pbs0", "alcprob0", "tx_int")] <- TRUE
matP0$values[c("pbs0", "alcprob0", "tx_int"),
             c("pbs0", "alcprob0", "tx_int")] <- 0.5
# Variances on diagonal free with starting value 1
diag(matP0$free) <- TRUE
diag(matP0$values) <- 1

dPmats <- list(
  matdP1 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP1"),
  matdP2 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP2"),
  matdP3 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP3"),
  matdP4 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP4"),
  matdP5 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP5"),
  matdP6 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP6"),
  matdP7 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP7"),
  matdP8 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP8"),
  matdP9 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP9"))

### Dummy vars ----
Dmats <- list(
  matD1 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d1", name = "d1"),
  matD2 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d2", name = "d2"),
  matD3 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d3", name = "d3"),
  matD4 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d4", name = "d4"),
  matD5 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d5", name = "d5"),
  matD6 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d6", name = "d6"),
  matD7 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d7", name = "d7"),
  matD8 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d8", name = "d8"),
  matD9 <- mxMatrix(
    type = "Full",  nrow = 1,  ncol = 1, labels = "data.d9", name = "d9"))

### Algebra ----
# alpha column vector
matA <- mxAlgebra(
  expression =
    matA0 +
    matdA1*d1 + matdA2*d2 + matdA3*d3 +
    matdA4*d4 + matdA5*d5 + matdA6*d6 +
    matdA7*d7 + matdA8*d8 + matdA9*d9,
  name = "matA")
# algebra for beta matrix
matB <- mxAlgebra(
  expression =
    matB0 +
    matdB1*d1 + matdB2*d2 + matdB3*d3 +
    matdB4*d4 + matdB5*d5 + matdB6*d6 +
    matdB7*d7 + matdB8*d8 + matdB9*d9,
  name = "matB")
# psi matrix
matP <- mxAlgebra(
  expression =
    matP0 +
    matdP1*d1 + matdP2*d2 + matdP3*d3 +
    matdP4*d4 + matdP5*d5 + matdP6*d6 +
    matdP7*d7 + matdP8*d8 + matdP9*d9,
  name = "matP")

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
fitUnconstr <- mxTryHard(modUnconstr, extraTries = 50)
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


## Partial synthesis ----
## In this model we want to estimate common effects for the mediation effect that is of
## interest in the meta-analysis. The paths between the baseline measures and the outcome
## variables will be stratified.
##
## Instead of 7 * 10 beta parameters, we now estimate 7 parameters in the baseline study
## and 4 parameters in the nine subsequent studies. That gives 27 degrees of freedom in
## the model comparison.

## Paths
# matrix B reference study
matB0 <- mxMatrix(
  type = "Full",  nrow = n_var,  ncol = n_var,
  dimnames = list(vars, vars), name = "matB0")
# set paths in free to TRUE
# paths to pbs1
matB0$free["pbs1", c("pbs0", "alcprob0", "tx_int")] <- TRUE
# Paths into alcprob1
matB0$free["alcprob1", c("pbs0", "alcprob0", "pbs1", "tx_int")] <- TRUE
matB0$values[matB0$free == TRUE] <- .1

# deviation effects for the other studies
dBmats <- list(
  matdB1 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = list(vars, vars), name = "matdB1"),
  matdB2 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = list(vars, vars), name = "matdB2"),
  matdB3 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = list(vars, vars), name = "matdB3"),
  matdB4 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = list(vars, vars), name = "matdB4"),
  matdB5 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = list(vars, vars), name = "matdB5"),
  matdB6 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = list(vars, vars), name = "matdB6"),
  matdB7 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = dimnames(matB0), name = "matdB7"),
  matdB8 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = dimnames(matB0), name = "matdB8"),
  matdB9 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_var, free = matB0$free,
    dimnames = dimnames(matB0), name = "matdB9"))

# Set paths from baseline variables to outcomes to FREE so they are stratified.
# No deviations from common effect (in matB0) are estimated for the other
# path coefficients.
for (i in 1:n_dum) {
  dBmats[[i]]$free <- FALSE
  dBmats[[i]]$free["pbs1", c("pbs0", "alcprob0")] <- TRUE
  dBmats[[i]]$free["alcprob1", c("pbs0", "alcprob0")] <- TRUE
}

matB <- mxAlgebra(
  expression =
    matB0 +
    matdB1*d1 + matdB2*d2 + matdB3*d3 +
    matdB4*d4 + matdB5*d5 + matdB6*d6 +
    matdB7*d7 + matdB8*d8 + matdB9*d9,
  name = "matB")

## Variances and Covariances
# Variance-covariance matrix for reference group
matP0 <- mxMatrix(
  type = "Symm", nrow = n_var, dimnames = list(vars, vars), name = "matP0")
# Covariances in off-diagonal cells free with starting values 0.5
matP0$free[c("pbs0", "alcprob0", "tx_int"),
           c("pbs0", "alcprob0", "tx_int")] <- TRUE
matP0$values[c("pbs0", "alcprob0", "tx_int"),
             c("pbs0", "alcprob0", "tx_int")] <- 0.5
# Variances on diagonal free with starting value 1
diag(matP0$free) <- TRUE
diag(matP0$values) <- 1

dPmats <- list(
  matdP1 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP1"),
  matdP2 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP2"),
  matdP3 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP3"),
  matdP4 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP4"),
  matdP5 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP5"),
  matdP6 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP6"),
  matdP7 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP7"),
  matdP8 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP8"),
  matdP9 <- mxMatrix(
    type = "Symm",  nrow = n_var, free = matP0$free, name = "matdP9"))

matP <- mxAlgebra(
  expression =
    matP0 +
    matdP1*d1 + matdP2*d2 + matdP3*d3 +
    matdP4*d4 + matdP5*d5 + matdP6*d6 +
    matdP7*d7 + matdP8*d8 + matdP9*d9,
  name = "matP")

### Build and run model ----
modBaselineModerated <- mxModel(
  "Moderated Baseline Measures Model",
  matA, matA0, dAmats,
  matB, matB0, dBmats,
  matP, matP0, dPmats,
  Dmats,
  matS, matM, matI,
  exp, funML,
  path_data)
fitBaselineModerated <- mxTryHard(modBaselineModerated)
summary(fitBaselineModerated, fitIndices = T)

mxCompare(fitUnconstr, fitBaselineModerated)
