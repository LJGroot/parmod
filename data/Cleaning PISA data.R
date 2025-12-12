## Header
##
## Script name: PISA Data cleaning
##
## Purpose of script:
##
## Author: L.J. Groot
##
## Date Created: 2025-10-15
##
## Copyright (c) Lennert J. Groot, 2025
## Email: l.j.groot@uva.nl
##
## Notes:
##   package loader:
if (!"pacman" %in% installed.packages()) install.packages("pacman")
library(pacman)
##

# Factor Model ----
## Libraries and Data ----------
## Load packages
### check if pacman package is installed, install if not
if (!"pacman" %in% installed.packages()) install.packages("pacman"); library(pacman)

### Load rest of packages
pacman::p_load(learningtower, dplyr, haven, psych, readr, readxl)

# read data
pisaTQ <- read_sav('CY08MSP_TCH_QQQ.sav')
# Check countries
countrycode[countrycode$country %in% pisaTQ$CNT, ]
# Select High income countries from data set
hinc_countries <- c(
  "AUS", "DEU", "HKG", "KOR", "MAC", "PRT", "ARE")
# filter high income countries
pisaTQ_hinc <- pisaTQ %>%
  filter(CNT %in% hinc_countries)
countrycode[countrycode$country %in% pisaTQ_hinc$CNT, ]
# select and rename math teaching items
pisaTQ_hinc <- pisaTQ_hinc[,c('CNT', 'TC227Q01JA','TC227Q03JA','TC227Q06JA','TC227Q08JA')]
names(pisaTQ_hinc) <- c('Country', paste0('y',1:4))

# remove cases without the items of interest
missingall <- apply(pisaTQ_hinc[, 2:5], 1, function(x) {
  sum(is.na(x[1:4])) == 4
})
pisadata <- pisaTQ_hinc[missingall == FALSE, ]

saveRDS(pisadata, file = "~/R/parmod/pisadata.rds")
