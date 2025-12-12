## Header ---------------------------
##
## Script name: Data Cleaning "itis_dat"
##
## Purpose of script: create long format IPD for IPD MASEM
##
## Author: L.J. Groot
##
## Date Created: 2025-04-15
##
## Copyright (c) Lennert J. Groot, 2025
## Email: l.j.groot@uva.nl
##
## Notes:
##

# Part 1: MNLFA model ---------------------------------
## Part 1a: prepare data ---------------------------------

## Load packages
### check if pacman package is installed, install if not
if (!"pacman" %in% installed.packages()) install.packages("pacman"); library(pacman)

### Load rest of packages
pacman::p_load(haven, psych, readr, readxl)

dataset1 <- read_sav("Blackwell (2007).sav")
# single sample, published


dataset2 <- read_sav("Mindset Assessment Profile Study Data_for OSF.sav")
# Burgoyne and Macnamara (2021) Single sample, published


dataset3 <- read_excel("Ingebritsen-2018-Data.xlsx")
# Ingebritsen school versus university treated as two separate samples
# https://munin.uit.no/handle/10037/12904  dissertation


dataset4 <- read_excel("MTG-StudentData-2018-21-03032021.xlsx")
# Park 2021 (abc) = grey  Three waves. I cannot find the study
# also three measurements. Let's use only T1


head(dataset1)
head(dataset2)
head(dataset3)
head(dataset4)


summary(dataset1)
summary(dataset2)
summary(dataset3)
summary(dataset4)

# Define rescaling function for likert-scale items with different length
resc <- function(newmax, oldmax, value = 1:oldmax) {
  (newmax - 1) * (value - 1) / (oldmax - 1) + 1
}
# Test function, two scenario's
## From 1-6 (length of 6) to 1-7
x1 <- sample(x = 1:6,
             size = 200,
             replace = TRUE)
range(x1)
x1_r <- resc(newmax = 7,
             oldmax = 6,
             value = x1)
range(x1_r)
## from 0-8 (lenght of 9) to 1-7
x2 <- sample(x = 0:8,
             size = 200,
             replace = TRUE)
range(x2)
# plus 1 to shift scale away from zero
x2_r <- resc(newmax = 7,
             oldmax = 9,
             value = x2 + 1)
range(x2_r)

# check StudyLevelInfo tab 'item mapping' to find the relevant items

# Dweck items
#F1	You have a certain amount of intelligence, and you really can't do much to change it.
#F2	Your intelligence is something about you that you can't change very much.
#F3	You can learn new things, but you can't really change your basic intelligence.
#G4	No matter who you are, you can change your intelligence a lot.
#G5	You can always greatly change how intelligent you are.
#G6	No matter how much intelligence you have, you can always change it quite a bit.
#F7	To be honest, you can't really change how intelligent you are.
#M8	You can change even your basic intelligence level considerably.

vars_Blackwell <- c("ent2","ent1","ent3","inc1","inc2","inc3",NA,NA)
vars_Burgoyne <- c("M1R","M3R","M7R","M2","M4","M6","M5R","M8")
vars_Ingebritsen <- c("GM_1_Neg","GM_2_Neg","GM_5_Neg","GM_3_Pos","GM_4_Pos",NA,NA,"GM_6_Pos")
vars_Park <- c("FM3_T1","FM2_T1","FM1_T1","GM5_T1",NA,"GM4_T1",NA,NA)

varnames <- c("study","F1","F2","F3","G4","G5","G6", "F7", "G8")

# Blackwell
data_Blackwell <- data.frame(
  study = 'Blackwell',
  F1 = dataset1$ent2,
  F2 = dataset1$ent1,
  F3 = dataset1$ent3,
  G4 = dataset1$inc1,
  G5 = dataset1$inc2,
  G6 = dataset1$inc3,
  F7 = NA,
  G8 = NA)

colnames(data_Blackwell) <- varnames
summary(data_Blackwell)
# all 1-6

# Burgoyne
data_Burgoyne <- data.frame(
  study = 'Burgoyne',
  F1 = dataset2$Mindset1R,
  F2 = dataset2$Mindset3R,
  F3 = dataset2$Mindset7R,
  G4 = dataset2$Mindset2,
  G5 = dataset2$Mindset4,
  G6 = dataset2$Mindset6,
  F7 = dataset2$Mindset5R,
  G8 = dataset2$Mindset8)
colnames(data_Burgoyne) <- varnames
summary(data_Burgoyne) # all 1-7

# reduce scale lenghts
data_Burgoyne[,2:9] <-  resc(newmax = 6, oldmax = 7, value = data_Burgoyne[,2:9])
summary(data_Burgoyne) # all 1-6

# Ingebritsen
select_school <- dataset3$Education == "VGS-elev"
select_uni <- dataset3$Education == "Student"

data_Ingebritsen_school <- data.frame(
  study = 'Ingebritsen_school',
  F1 = dataset3[select_school,"GM1Neg"],
  F2 = dataset3[select_school,"GM2Neg"],
  F3 = dataset3[select_school,"GM5Neg"],
  G4 = dataset3[select_school,"GM3Pos"],
  G5 = dataset3[select_school,"GM4Pos"],
  G6 = NA,
  F7 = NA,
  G8 = dataset3[select_school,"GM6Pos"])

colnames(data_Ingebritsen_school) <- varnames
summary(data_Ingebritsen_school)
# some 1-5 but not due to scale length i presume, due to small n

data_Ingebritsen_uni <- data.frame(
  study = 'Ingebritsen_uni',
  F1 = dataset3[select_uni,"GM1Neg"],
  F2 = dataset3[select_uni,"GM2Neg"],
  F3 = dataset3[select_uni,"GM5Neg"],
  G4 = dataset3[select_uni,"GM3Pos"],
  G5 = dataset3[select_uni,"GM4Pos"],
  G6 = NA,
  F7 = NA,
  G8 = dataset3[select_uni,"GM6Pos"])

colnames(data_Ingebritsen_uni) <- varnames
summary(data_Ingebritsen_uni) # all 1-6

# Park
select_a <- dataset4$WAVE == "2018"
select_b <- dataset4$WAVE == "2019"
select_c <- dataset4$WAVE == "2020"

data_Park_a <- data.frame(
  study = 'Park_a',
  F1 = dataset4[select_a,"FM3_T1"],
  F2 = dataset4[select_a,"FM2_T1"],
  F3 = dataset4[select_a,"FM1_T1"],
  G4 = dataset4[select_a,"GM5_T1"],
  G5 = NA,
  G6 = dataset4[select_a,"GM4_T1"],
  F7 = NA,
  G8 = NA)

colnames(data_Park_a) <- varnames
summary(data_Park_a)
data_Park_a[,2:9] <-  resc(newmax = 6, oldmax = 7, value = data_Park_a[,2:9])
summary(data_Park_a)


data_Park_b <- data.frame(
  study = 'Park_b',
  F1 = dataset4[select_b,"FM3_T1"],
  F2 = dataset4[select_b,"FM2_T1"],
  F3 = dataset4[select_b,"FM1_T1"],
  G4 = dataset4[select_b,"GM5_T1"],
  G5 = NA,
  G6 = dataset4[select_b,"GM4_T1"],
  F7 = NA,
  G8 = NA)

colnames(data_Park_b) <- varnames
summary(data_Park_b)
data_Park_b[,2:9] <-  resc(newmax = 6, oldmax = 7, value = data_Park_b[,2:9])
summary(data_Park_b)



data_Park_c <- data.frame(
  study = 'Park_c',
  F1 = dataset4[select_c,"FM3_T1"],
  F2 = dataset4[select_c,"FM2_T1"],
  F3 = dataset4[select_c,"FM1_T1"],
  G4 = dataset4[select_c,"GM5_T1"],
  G5 = NA,
  G6 = dataset4[select_c,"GM4_T1"],
  F7 = NA,
  G8 = NA)

colnames(data_Park_c) <- varnames
summary(data_Park_c)
data_Park_c[,2:9] <-  resc(newmax = 6, oldmax = 7, value = data_Park_c[,2:9])
summary(data_Park_c)

# Merged data
itis_dat <- rbind(data_Blackwell,
                  data_Burgoyne,
                  data_Ingebritsen_school,
                  data_Ingebritsen_uni,
                  data_Park_a,
                  data_Park_b,
                  data_Park_c)

saveRDS(object = itis_dat, file = "~/R/itis_dat.Rds")
itis_dat <- readRDS("itis_dat.rds")
