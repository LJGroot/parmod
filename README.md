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
The data folder contains data files that are used in the tutorial and as example data in the additional applications of PM.
The following data used in the tutorial and included in the repository:

The data used for the factor model tutorial are published in:

-   Organisation for Economic Co-operation and Development. (2024). PISA 2022 Database [Data
set] [Retrieved from https://www.oecd.org/en/data/datasets/pisa-2022-database.html].
  
The data used for the path model tutorial are published in:

-   Huh, D., Li, X., Zhou, Z., Walters, S. T., Baldwin, S. A., Tan, Z., Larimer, M. E., & Mun, E.-Y.
  (2022). A structural equation modeling approach to meta-analytic mediation analysis
  using individual participant data: Testing protective behavioral strategies as a mediator
  of brief motivational intervention effects on alcohol-related problems. Prevention Science,
  23(3), 390–402. https://doi.org/10.1007/s11121-021-01318-4

Additional data is used as example data in the flexible syntax files for multi-factor CFA models

-   Scherer, R., & Campos, D. G. (2022). Measuring those who have their minds set: An item-level
  meta-analysis of the implicit theories of intelligence scale in education. Educational
  Research Review, 37, 100479. https://doi.org/https://doi.org/10.1016/j.edurev.2022.100479

### Functions
The functions folder contains a set of helper functions. Please note that these are not used in the tutorial or described in the publication. 
The functions are used in the additional syntax files.
