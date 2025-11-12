This repostiory accompanies the IPD MASEM tutorial paper "Making the Most of Minimal Data: IPD MASEM via Parameter Moderation Using OpenMx" by Lennert J. Groot, Kees Jan Kan, and Suzanne Jak (Manuscript in preparation).

The syntax folder contains R scripts. 
  parmod_tutorial.R contains the syntax that is used in the publication and uses a completely manual (but somewhat laborious) specification of each of the models.
  parmod_CFA_SF.R, parmod_CFA_MF.R and parmod_PA.R contain syntax to use PM for, respecitively, a single-factor factor model, a multi-factor factor model, and a path model.
    the syntax in these files is flexible and allows for the user to conduct IPD MASEM via PM on their own data and hypothesized model.
    
The data folder contains data files that are used in the tutorial and as example data in the additional applications of PM.

The functions folder contains a set of helper functions that are not used in the tutorial, but are used in the additional syntax files. The 

The data used for the factor model tutorial are published in:
Scherer, R., & Campos, D. G. (2022). Measuring those who have their minds set: An item-level
  meta-analysis of the implicit theories of intelligence scale in education. Educational
  Research Review, 37, 100479. https://doi.org/https://doi.org/10.1016/j.edurev.2022.100479
  
The data used for the path model tutorial are published in:
Huh, D., Li, X., Zhou, Z., Walters, S. T., Baldwin, S. A., Tan, Z., Larimer, M. E., & Mun, E.-Y.
  (2022). A structural equation modeling approach to meta-analytic mediation analysis
  using individual participant data: Testing protective behavioral strategies as a mediator
  of brief motivational intervention effects on alcohol-related problems. Prevention Science,
  23(3), 390–402. https://doi.org/10.1007/s11121-021-01318-4
