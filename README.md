# parmod
This repostiory accompanies the IPD MASEM tutorial paper "Making the Most of Minimal Data: IPD MASEM via Parameter Moderation Using OpenMx" by Lennert J. Groot, Kees Jan Kan, and Suzanne Jak (Manuscript under revision).

## Contents of this repository
### Syntax
This syntax folder contains the following R scripts:

-   `parmod_tutorial_CFA`
-   `parmod_tutorial_PA`
-   `parmod_CFA_SF`
-   `parmod_CFA_MF`
-   `parmod_PA`

The use and contents of the different scripts are described below.

#### parmod tutorial CFA

`parmod_tutorial_CFA.R` contains the syntax that is used in the tutorial paper, in particular in the section on evaluating a factor model. The script uses a completely manual specification of each of the models.

#### parmod tutorial PA

`parmod_tutorial_CFA.R` contains the syntax that is used in the tutorial paper, in particular in the section on evaluating a path model. The script uses a completely manual specification of each of the models.


#### parmod CFA SF

parmod_CFA_SF.R, contains syntax to use PM for a single-factor model. This syntax file is flexible and allows for the user to conduct IPD MASEM via PM on their own data and hypothesized model. It uses some more helper functions and shorter programming than the manual specification in the tutorial.

#### parmod CFA MF

`parmod_CFA_MF.R` contains syntax to use PM for a single-factor a multi-factor model. This syntax file is flexible and allows for the user to conduct IPD MASEM via PM on their own data and hypothesized model. It uses some more helper functions and shorter programming than the manual specification in the tutorial.

#### parmod PA

`parmod_PA.R` contains syntax to use PM for a path model. This syntax file is flexible and allows for the user to conduct IPD MASEM via PM on their own data and hypothesized model. It uses some more helper functions and shorter programming than the manual specification in the tutorial.

### Data
#### Integrate
The integrate folder contains the integrate.RData file. This file contains IPD from the publication by Huh et al. (2022). The IPD is in long format with study identifier. The data was cleaned and prepared by Huh et al. (2022). The file integrate.RData is used in the tutorial paper, in particual in the section on evaluating a path model. It is also used as example data in the file parmod_PA.R. This latter syntax file can be used to conduct IPD MASEM via PM on your own data when your hyphothesized model is a path model.

This data was originally published in:

Huh, D., Li, X., Zhou, Z., Walters, S. T., Baldwin, S. A., Tan, Z., Larimer, M. E., & Mun, E.-Y. (2022). A structural equation modeling approach to meta-analytic mediation analysis using individual participant data: Testing protective behavioral strategies as a mediator of brief motivational intervention effects on alcohol-related problems. Prevention Science, 23(3), 390–402. https://doi.org/10.1007/s11121-021-01318-4

#### Itis
The itis folder contains four raw data files: (1)Blackwell (2007).sav; (2) Ingebritsen-2018-Data.xlsx; (3) MTG-StudentData-2018-21-03032021.xlsx; (4) Mindset Assessment Profile Study Data_for OSF.sav), a data cleaning script, and a cleaned data file that was generated from the raw data and the cleaning script. The cleaned data file itis_dat.Rds is used as example data in the parmod_CFA_MF.R syntax file. This latter syntax file can be used to conduct IPD MASEM via PM on your own data when your hyphothesized model is a factor model with more than one latent factor.

This data was originally published in:

Scherer, R., & Campos, D. G. (2022). Measuring those who have their minds set: An item-level meta-analysis of the implicit theories of intelligence scale in education. Educational Research Review, 37, 100479. https://doi.org/https://doi.org/10.1016/j.edurev.2022.100479

#### Pisa
the pisa folder contains a very large, raw data file (CY08MSP_TCH_QQQ.SAV), a data cleaning script, and a cleaned data file that was generated from the raw data and the cleaning script. The cleaned data file Pisadata.Rds is used in the tutorial, in particular in the section on evaluating a factor model. It also serves as example data in the parmod_CFA_SF.R syntax file. This latter syntax file can be used to conduct IPD MASEM via PM on your own data when your hyphothesized model is a factor model with one latent factor.

The PISA data are published as:

Organisation for Economic Co-operation and Development. (2024). PISA 2022 Database [Data set] [Retrieved from https://www.oecd.org/en/data/datasets/pisa-2022-database.html].

### Functions
The functions folder contains a set of helper functions. Please note that these are not used in the tutorial or described in the publication. 
The functions are used in the additional syntax files.

There are four separate syntax files with helper functions found in this repository: (1) add_dummies(); (2) factor_mask(); (3) make_alg(); and (4) make_matrices_config(). They can each be loaded into the workspace separately, or all together by sourcing the helper_functions.R script.

Each helper function also comes with its own _help function that prints the documentation for the function in the console.

The best way to use the helper functions is by loading them all together into the workspace using the package osfr and the source() function.

helper_functions <- osf_download(osf_retrieve_file("c2u78"), conflicts = "overwrite"); source(helper_functions$local_path)
Alternatively, one could download the helper functions from the Github repository via the syntax:

source("https://raw.githubusercontent.com/LJGroot/parmod/refs/heads/main/functions/helper_functions.R")
Please check the documentation in the console using the syntax that is shown as the functions are succesfully loaded.
