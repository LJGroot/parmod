## Header ----
##
## Script name:
## IPD MASEM via Parameter Moderation Tutorial - Confirmatory Factor Model
##
## Author: L.J. Groot
## Copyright (c) Lennert J. Groot, 2025
## Email: l.j.groot@uva.nl
##
## Notes:
##

# IPD MASEM via Parameter Moderation - Factor Model Tutorial ----
# This script contains syntax to conduct IPD MASEM via Parameter Moderation (PM). In this
# script, we fit a factor model to IPD from 7 primary studies. The data originates from
# "Organisation for Economic Co-operation and Development. (2024). PISA 2022 Database
# [Data set] [Retrieved from oecd.org/en/data/datasets/pisa-2022-database.html]"
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

# Factor Model ----
## Libraries and Data ----------
## Load packages
### check if pacman package is installed, install if not
if (!"pacman" %in% installed.packages()) install.packages("pacman"); library(pacman)

### Load rest of packages
pacman::p_load(osfr, haven, psych, readr, readxl, lavaan, metaSEM, OpenMx)

# read data (from online repository)
data_file <- osf_retrieve_file("g4s8e") |>
  osf_download(conflicts = "overwrite"); pisadata <- readRDS(data_file$local_path)

dat <- cbind.data.frame(
  pisadata,
  data.frame(
    # AUS as reference group
    "d1" = as.numeric(pisadata$Country == "DEU"),
    "d2" = as.numeric(pisadata$Country == "HKG"),
    "d3" = as.numeric(pisadata$Country == "KOR"),
    "d4" = as.numeric(pisadata$Country == "MAC"),
    "d5" = as.numeric(pisadata$Country == "PRT"),
    "d6" = as.numeric(pisadata$Country == "ARE")
  )
)
# Store data characteristics
vars <- names(dat[2:5])
n_var <- length(vars)
n_lv <- 1
k <- length(unique(dat$Country))
n_dum <- k - 1
dummies <- paste0("d", 1:n_dum)
# Store data as mxData object
factor_data <- mxData(dat[,c(vars, dummies)], type = "raw")

## Naive Model ----
matT <- mxMatrix(
  type = "Full",  nrow = 1,  ncol = n_var,
  free = TRUE,  values = 0.5,  name = "matT")
free_lambdas <- rep(TRUE, n_var)
matL <- mxMatrix(
  type = "Full",  nrow = n_var,  ncol = n_lv,
  free = free_lambdas,
  values = as.numeric(free_lambdas), byrow = TRUE,
  name = "matL")
matE <- mxMatrix(
  type = "Diag", nrow = n_var, free = TRUE, values = 1, name = "matE")
matK <- mxMatrix(
  type = "Full",  nrow = n_lv,  ncol = 1,
  free = FALSE,  values = 0,
  name = "matK")
matP <- mxMatrix(
  type = "Symm", nrow = n_lv,
  free = !as.logical(diag(nrow = n_lv)),
  values = {v <- matrix(0.5, nrow = n_lv, ncol = n_lv);diag(v) <- 1;v},
  byrow = TRUE,
  name = "matP")
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
fitNaive <- mxTryHard(modNaive, extraTries = 50)
summary(fitNaive, fitIndices = TRUE)


## Configural Model --------------------------
### Intercepts -----
matT0 <- mxMatrix(
  type = "Full",  nrow = 1,  ncol = n_var,
  free = TRUE,  values = 1,  name = "matT0")
# Delta tau matrices for deviation of intercept in studies 2-7
dTmats <- list(
  matdT1 <- mxMatrix(
    type = "Full", nrow = 1, ncol = n_var, free = TRUE, values = 0, name = "matdT1"),
  matdT2 <- mxMatrix(
    type = "Full", nrow = 1, ncol = n_var, free = TRUE, values = 0, name = "matdT2"),
  matdT3 <- mxMatrix(
    type = "Full", nrow = 1, ncol = n_var, free = TRUE, values = 0, name = "matdT3"),
  matdT4 <- mxMatrix(
    type = "Full", nrow = 1, ncol = n_var, free = TRUE, values = 0, name = "matdT4"),
  matdT5 <- mxMatrix(
    type = "Full", nrow = 1, ncol = n_var, free = TRUE, values = 0, name = "matdT5"),
  matdT6 <- mxMatrix(
    type = "Full", nrow = 1, ncol = n_var, free = TRUE, values = 0, name = "matdT6"))

### Factor Loadings ----
# Lambda_0, baseline factor loadings for reference group
free_lambdas <- rep(TRUE, n_var)
matL0 <- mxMatrix(
  type = "Full",  nrow = n_var,  ncol = n_lv,
  free = free_lambdas,
  values = 1, byrow = TRUE,
  name = "matL0")
# delta lambda matrices for study 2-7
dLmats <- list(
  matdL1 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_lv, free = matL0$free, name = "matdL1"),
  matdL2 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_lv, free = matL0$free, name = "matdL2"),
  matdL3 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_lv, free = matL0$free, name = "matdL3"),
  matdL4 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_lv, free = matL0$free, name = "matdL4"),
  matdL5 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_lv, free = matL0$free, name = "matdL5"),
  matdL6 <- mxMatrix(
    type = "Full", nrow = n_var, ncol = n_lv, free = matL0$free, name = "matdL6"))

### Residual Variances ----
# E0 Baseline residual variances for reference group
matE0 <- mxMatrix(
  type = "Diag", nrow = n_var, free = TRUE, values = 1, name = "matE0")
# Delta theta matrices for difference from baseline residual variances for studies 2-7
dEmats <- list(
  matdE1 <- mxMatrix(
    type = "Diag", nrow = n_var, free = TRUE, values = 0, name = "matdE1"),
  matdE2 <- mxMatrix(
    type = "Diag", nrow = n_var, free = TRUE, values = 0, name = "matdE2"),
  matdE3 <- mxMatrix(
    type = "Diag", nrow = n_var, free = TRUE, values = 0, name = "matdE3"),
  matdE4 <- mxMatrix(
    type = "Diag", nrow = n_var, free = TRUE, values = 0, name = "matdE4"),
  matdE5 <- mxMatrix(
    type = "Diag", nrow = n_var, free = TRUE, values = 0, name = "matdE5"),
  matdE6 <- mxMatrix(
    type = "Diag", nrow = n_var, free = TRUE, values = 0, name = "matdE6"))

### Latent Means ----
# K0 latent variable means in reference group
matK0 <- mxMatrix(
  type = "Full",  nrow = n_lv, ncol = 1,
  free = FALSE,  values = 0,
  name = "matK0")
# delta kappa of latent variable means from baseline in study 2-7
dKmats <- list(
  matdK1 <- mxMatrix(
    type = "Full", nrow = n_lv, ncol = 1, free = FALSE, values = 0, name = "matdK1"),
  matdK2 <- mxMatrix(
    type = "Full", nrow = n_lv, ncol = 1, free = FALSE, values = 0, name = "matdK2"),
  matdK3 <- mxMatrix(
    type = "Full", nrow = n_lv, ncol = 1, free = FALSE, values = 0, name = "matdK3"),
  matdK4 <- mxMatrix(
    type = "Full", nrow = n_lv, ncol = 1, free = FALSE, values = 0, name = "matdK4"),
  matdK5 <- mxMatrix(
    type = "Full", nrow = n_lv, ncol = 1, free = FALSE, values = 0, name = "matdK5"),
  matdK6 <- mxMatrix(
    type = "Full", nrow = n_lv, ncol = 1, free = FALSE, values = 0, name = "matdK6"))

### Factor Variance ----
# P latent variable covariance matrix for reference group
# Fixed to 1 for identification purposes, same fixed value for every study
matP <- mxMatrix(
  type = "Symm",  nrow = n_lv, free = FALSE,  values = 1,  name = "matP")

### Dummy Vars ----
# Define background variables (this is our dummy variable created to differentiate
# between study 1 and study 2)
Dmats <- list(
  matD1 <- mxMatrix(
    type = "Full", nrow = 1, ncol = 1, labels = "data.d1", name = "d1"),
  matD2 <- mxMatrix(
    type = "Full", nrow = 1, ncol = 1, labels = "data.d2", name = "d2"),
  matD3 <- mxMatrix(
    type = "Full", nrow = 1, ncol = 1, labels = "data.d3", name = "d3"),
  matD4 <- mxMatrix(
    type = "Full", nrow = 1, ncol = 1, labels = "data.d4", name = "d4"),
  matD5 <- mxMatrix(
    type = "Full", nrow = 1, ncol = 1, labels = "data.d5", name = "d5"),
  matD6 <- mxMatrix(
    type = "Full", nrow = 1, ncol = 1, labels = "data.d6", name = "d6"))

### Algebra ----
#### Intercepts -----
matT <- mxAlgebra(
  expression =
    matT0 + # Baseline intercepts (reference group, study 1)
    matdT1*d1 + matdT2*d2 + matdT3*d3 +
    matdT4*d4 + matdT5*d5 + matdT6*d6, # difference between study 1 and study 2-7
  name = "matT")
#### Loadings ----
matL <- mxAlgebra(
  expression =
    matL0 + # Baseline factor loadings
    matdL1*d1 + matdL2*d2 + matdL3*d3 +
    matdL4*d4 + matdL5*d5 + matdL6*d6, # difference between study 1 and study 2-7
  name = "matL")

#### Residual Variances ----
matE <- mxAlgebra(
  expression =
    matE0*exp(            # Baseline residual variances
      matdE1*d1 + matdE2*d2 + matdE3*d3 +
      matdE4*d4 + matdE5*d5 + matdE6*d6 ), # difference between study 1 and study 2-7
  name = "matE")

#### Factor means ----
matK <- mxAlgebra(
  expression =
    matK0 + # Baseline factor means
    matdK1*d1 + matdK2*d2 + matdK3*d3 +
    matdK4*d4 + matdK5*d5 + matdK6*d6, # difference between study 1 and study 2-7
  name = "matK")

### Model-implied matrices ----
matS <- mxAlgebra(
  expression =
    matL %*% matP %*% t(matL) + matE,
  name = "matS")
matM <- mxAlgebra(
  expression =
    matT + t(matL %*% matK),
  name = "matM")

### Expectation and fit function ----
expF <- mxExpectationNormal(covariance = "matS",
                            means = "matM",
                            dimnames = vars)
fitF <- mxFitFunctionML()

### Build mxModel and run ----
modConfig <- mxModel(
  model = "Configural",   # Model Name
  matT,  matT0,  dTmats,  # Intercepts
  matL,  matL0,  dLmats,  # Loadings
  matE,  matE0,  dEmats,  # Resid vars
  matK,  matK0,  dKmats,  # Factor means
  matP,                   # Common factor var-cov
  Dmats,                  # Dummy variables
  matM,  matS,            # Model algebraic expressions
  expF,                   # Expectation function
  fitF,                   # Fit function
  factor_data )           # IPD
fitConfig <- mxTryHard(modConfig)
summary(fitConfig, fitIndices = TRUE)

## Metric Model ---------------

### Loadings ----
for (i in 1:n_dum) dLmats[[i]]$free <- FALSE

### Factor Covariance ----
# Baseline matrix Phi_0
matP0 <- mxMatrix(
  type = "Symm", nrow = n_lv, values = 1, byrow = TRUE, name = "matP0")
# Deviation of lv cor from baseline in study 2-7
dPmats <- list(
  matdP1 <- mxMatrix(
    type = "Symm", nrow = n_lv, free = TRUE, name = "matdP1"),
  matdP2 <- mxMatrix(
    type = "Symm", nrow = n_lv, free = TRUE, name = "matdP2"),
  matdP3 <- mxMatrix(
    type = "Symm", nrow = n_lv, free = TRUE, name = "matdP3"),
  matdP4 <- mxMatrix(
    type = "Symm", nrow = n_lv, free = TRUE, name = "matdP4"),
  matdP5 <- mxMatrix(
    type = "Symm", nrow = n_lv, free = TRUE, name = "matdP5"),
  matdP6 <- mxMatrix(
    type = "Symm", nrow = n_lv, free = TRUE, name = "matdP6"))

# Algebra
matP <- mxAlgebra(
  expression =
    matP0*exp(
      matdP1*d1 + matdP2*d2 + matdP3*d3 +
      matdP4*d4 + matdP5*d5 + matdP6*d6),
  name = "matP")

## Make mxModel object and run the model
### Build mxModel and run ----
modMetric <- mxModel(
  model = "Metric",   # Model Name
  matT,  matT0,  dTmats,  # Intercepts
  matL,  matL0,  dLmats,  # Loadings
  matE,  matE0,  dEmats,  # Resid vars
  matP,  matP0,  dPmats,  # Common factor cor
  matK,  matK0,  dKmats,  # Factor means
  Dmats,                  # Dummy variables
  matM,  matS,            # Model algebraic expressions
  expF,                   # Expectation function
  fitF,                   # Fit function
  factor_data )           # IPD
fitMetric <- mxRun(modMetric)
summary(fitMetric, fitIndices = TRUE)

# Compare models
mxCompare(fitConfig,fitMetric)

## Scalar model ----
### Intercepts ----
for (i in 1:n_dum) dTmats[[i]]$free <- FALSE

### Factor Means ----
for (i in 1:n_dum) dKmats[[i]]$free <- TRUE

### Build and run model ----
modScalar <- mxModel(
  model = "Scalar",   # Model Name
  matT,  matT0,  dTmats,  # Intercepts
  matL,  matL0,  dLmats,  # Loadings
  matE,  matE0,  dEmats,  # Resid vars
  matP,  matP0,  dPmats,  # Common factor cor
  matK,  matK0,  dKmats,  # Factor means
  Dmats,                  # Dummy variables
  matM,  matS,            # Model algebraic expressions
  expF,                   # Expectation function
  fitF,                   # Fit function
  factor_data )
fitScalar <- mxRun(modScalar)
summary(fitScalar, fitIndices = TRUE)

# Compare models
mxCompare(fitConfig, c(fitMetric, fitScalar))

## Partial Invariance: All-But-One Procedure ----
# This procedure uses a for-loop to reduce constraints on the indicator variables
# one by one. These models are then compared against the scalar model, to test which
# additional constraint has the smallest negative effect on the model fit.
#
# The necessary objects will be stored in the working environment using shorthand
# code that uses some helper functions and alternative syntax, but has the same output
# as the syntax in the section for the configural model above.
#
# This example works for the single factor model in this tutorial. For an example on a
# two-factor model (albeit in a MI context) see paper by Kolbe et al. (2022):
# https://doi.org/10.1037/met0000501

### ABO Procedure ----
# This syntax runs a `for` loop, constraining the intercept and factor loading of
# one of the indicators at a time, saving the output of the model in a list.
# The ABO models have 60 estimated parameters. There were 48 estimated parameters in
# the scalar model. We are estimating an additional k-1 = 6 intercepts and
# k-1 = 6 factor loadings in each of the ABO models.
# In total there are 12 extra estimated parameters, compared to the scalar model.
# That means that the LRT tests that we conduct later using `mxCompare` should
# have df = 12.

# Create empty list for results
fitABO <- list()

# write for loop
for (i in 1:n_var) {
  # column vector for intercepts
  freeparT <- matrix(data = FALSE,
                     nrow = 1,
                     ncol = n_var)
  # set one value to TRUE
  freeparT[i] <- TRUE
  # Loadings
  freeparL <- matrix(data = FALSE,
                     nrow = n_var,
                     ncol = n_lv)
  # set one value to TRUE
  freeparL[i] <- TRUE
  # Tau deviation matrices
  dTmats <- list(
  matdT1 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparT,
                    values = 0, name = "matdT1"),
  matdT2 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparT,
                    values = 0, name = "matdT2"),
  matdT3 <- mxMatrix(type = "Full", nrow = 1, ncol = n_var,free = freeparT,
                    values = 0, name = "matdT3"),
  matdT4 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparT,
                    values = 0, name = "matdT4"),
  matdT5 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparT,
                    values = 0, name = "matdT5"),
  matdT6 <- mxMatrix(type = "Full", nrow = 1, ncol = n_var,free = freeparT,
                    values = 0, name = "matdT6"))
  # Lambda deviation matrices
  dLmats <- list(
  matdL1 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                    values = 0, byrow = TRUE, name = "matdL1"),
  matdL2 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                    values = 0, byrow = TRUE, name = "matdL2"),
  matdL3 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                    values = 0, byrow = TRUE, name = "matdL3"),
  matdL4 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                    values = 0, byrow = TRUE, name = "matdL4"),
  matdL5 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                    values = 0, byrow = TRUE, name = "matdL5"),
  matdL6 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                    values = 0, byrow = TRUE, name = "matdL6"))
  modABO <- mxModel(
    model = paste0("All_but_", i),
    matT,  matT0,  dTmats,  # Intercepts
    matL,  matL0,  dLmats,  # Loadings
    matE,  matE0,  dEmats,  # Resid vars
    matK,  matK0,  dKmats,  # Factor means
    matP,  matP0,  dPmats,  # Common factor var-cov
    Dmats,                  # Dummy variables
    matM,  matS,            # Model algebraic expressions
    expF,                   # Expectation function
    fitF,                   # Fit function
    factor_data             # IPD
  )
  fitABO[[i]] <- mxTryHard(modABO, extraTries = 25)
}

anchorTest <- mxCompare(fitABO, fitScalar)
anchorOut <- data.frame(
  name = paste0("Indicator", 1:n_var),
  X2 = anchorTest$diffLL[seq(2,2*n_var,2)],
  df = anchorTest$diffdf[seq(2,2*n_var,2)],
  p = anchorTest$p[seq(2,2*n_var,2)]
)

anchorOut[order(anchorOut$X2[1:4]), ]

### Anchors-only model------------
# Now that we have found the items that produce the least amount of misfit, we can
# fit a model where only a common effect is estimated for the anchor item(s). The
# rest of the items are stratified.
#
# It is advised to choose ca. 20% of the items per latent factor as anchor items.
# (see paper for reference). In our case, that means that there is 1 anchor item.

anchors <- c(which.min(anchorOut[, "X2"]))


# which parameters to estimate freely
testIn <- c(1:n_var)[-c(anchors)]
freeparT <- matrix(TRUE, nrow = 1, ncol = n_var)
freeparT[1,c(anchors)] <- FALSE
freeparL <- matrix(TRUE,
                   nrow = n_var,
                   ncol = n_lv,
                   byrow = TRUE)
freeparL[c(anchors),1] <- FALSE

# redefine Tau and Lambda deviation matrices
# Tau deviation matrices
dTmats <- list(
  matdT1 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparT,
                     values = 0, name = "matdT1"),
  matdT2 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparT,
                     values = 0, name = "matdT2"),
  matdT3 <- mxMatrix(type = "Full", nrow = 1, ncol = n_var,free = freeparT,
                     values = 0, name = "matdT3"),
  matdT4 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparT,
                     values = 0, name = "matdT4"),
  matdT5 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparT,
                     values = 0, name = "matdT5"),
  matdT6 <- mxMatrix(type = "Full", nrow = 1, ncol = n_var,free = freeparT,
                     values = 0, name = "matdT6"))
# Lambda deviation matrices
dLmats <- list(
  matdL1 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                     values = 0, byrow = TRUE, name = "matdL1"),
  matdL2 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                     values = 0, byrow = TRUE, name = "matdL2"),
  matdL3 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                     values = 0, byrow = TRUE, name = "matdL3"),
  matdL4 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                     values = 0, byrow = TRUE, name = "matdL4"),
  matdL5 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                     values = 0, byrow = TRUE, name = "matdL5"),
  matdL6 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparL,
                     values = 0, byrow = TRUE, name = "matdL6"))
modAnchor <- mxModel(
  model = paste0("Anchor Only"),
  matT,  matT0,  dTmats,  # Intercepts
  matL,  matL0,  dLmats,  # Loadings
  matE,  matE0,  dEmats,  # Resid vars
  matK,  matK0,  dKmats,  # Factor means
  matP,  matP0,  dPmats,  # Common factor var-cov
  Dmats,                  # Dummy variables
  matM,  matS,            # Model algebraic expressions
  expF,                   # Expectation function
  fitF,                   # Fit function
  factor_data             # IPD
)
fitAnchor <- mxTryHard(modAnchor)

mxCompare(fitConfig, fitAnchor)

# Since the hypothesized model in this example has only 4 items on its single latent
# factor, the model with a common effect for the one anchor item and stratified effects
# for the other items is equivalent to the configural model with the only difference
# begin that now the model is identified through the factor loading of the anchor item
# rather than through the common factor variance. This is not the case when the number of
# indicators per latent factor is sufficiently large so that 20% of the indicators per
# latent factor equals 2 or more.
#


### Anchor-plus-one Model----
fitApo <- list()
for (i in testIn) {
  freeparTa <- freeparT
  freeparLa <- freeparL
  freeparTa[i] <- FALSE
  freeparLa[i] <- FALSE
  dTmats <- list(
    matdT1 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparTa,
                       values = 0, name = "matdT1"),
    matdT2 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparTa,
                       values = 0, name = "matdT2"),
    matdT3 <- mxMatrix(type = "Full", nrow = 1, ncol = n_var,free = freeparTa,
                       values = 0, name = "matdT3"),
    matdT4 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparTa,
                       values = 0, name = "matdT4"),
    matdT5 <- mxMatrix(type = "Full",nrow = 1, ncol = n_var, free = freeparTa,
                       values = 0, name = "matdT5"),
    matdT6 <- mxMatrix(type = "Full", nrow = 1, ncol = n_var,free = freeparTa,
                       values = 0, name = "matdT6"))
  # Lambda deviation matrices
  dLmats <- list(
    matdL1 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparLa,
                       values = 0, byrow = TRUE, name = "matdL1"),
    matdL2 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparLa,
                       values = 0, byrow = TRUE, name = "matdL2"),
    matdL3 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparLa,
                       values = 0, byrow = TRUE, name = "matdL3"),
    matdL4 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparLa,
                       values = 0, byrow = TRUE, name = "matdL4"),
    matdL5 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparLa,
                       values = 0, byrow = TRUE, name = "matdL5"),
    matdL6 <- mxMatrix(type = "Full",nrow = n_var, ncol = n_lv, free = freeparLa,
                       values = 0, byrow = TRUE, name = "matdL6"))
  modApo <- mxModel(
    model = paste0("Anchors_plus_", as.character(i)),
    matT,  matT0,  dTmats,  # Intercepts
    matL,  matL0,  dLmats,  # Loadings
    matE,  matE0,  dEmats,  # Resid vars
    matK,  matK0,  dKmats,  # Factor means
    matP,  matP0,  dPmats,  # Common factor var-cov
    Dmats,                  # Dummy variables
    matM,  matS,            # Model algebraic expressions
    expF,                   # Expectation function
    fitF,                   # Fit function
    factor_data             # IPD
  )
  fitApo[[as.character(i)]] <- mxTryHard(modApo, extraTries = 25)
}

# Test anchors-only model against anchors-plus-one models
piTest <- mxCompare(fitAnchor, fitApo)
piOut <- data.frame(name = paste0
                    ("Indicator", testIn),
                    X2 = piTest$diffLL[2:4],
                    df = piTest$diffdf[2:4],
                    p = piTest$p[2:4],
                    p.bon = p.adjust(p = piTest$p[2:4], method = "bonferroni"),
                    p.BH = p.adjust(p = piTest$p[2:4], method = "BH"))
