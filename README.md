# Van-Der-Heijden-Log-Linear-EM-Algorithm
Van Der Heijden et al Log Linear Expectation-Maximization Algorithm

# Description

This method allows imputing missing cells of a contingency table and for undergoing multiple system estimation. The algorithm assumes that any missingness is missing at random (MAR).

# Academic Papers
Two papers are provided that reference 
this algorithm in detail:

- **Title:** *Multiple System Estimation using Covariates having Missing Values and Measurement Error: Estimating the Size of the Māori Population in New Zealand*
  - **Citation:** Journal of the Royal Statistical Society Series A: Statistics in Society, Volume 185, Issue 1, January 2022, Pages 156–177, https://doi.org/10.1111/rssa.12731
- **Title:** *An Overview of Population Size Estimation where Linking Registers Results in Incomplete Covariates, with an Application to Mode of Transport of Serious Road Casualties*
  - **Citation:** Heijden, Peter G.M. & Smith, Paul & Cruyff, M.J.L.F & Bakker, Bart. (2017). An Overview of Population Size Estimation where Linking Registers Results in Incomplete Covariates, with an Application to Mode of Transport of Serious Road Casualties. Journal of Official Statistics. 34. 10.1515/jos-2018-0011.


# Methods
- Any model-fitting method that allows a poisson-distributed log-link model can be utilized my recursively predicting cells until cell values converge.
- The R Package [Coarsened Variable Modeling (CVAM)](https://github.com/uscensusbureau/cvam) Version 0.9.4.3+ implements this algorithm natively.
- Manually optimizing the likelihood is another method in lieu of using a model-fitting method for poisson-distributed log-link models.
- For cases where there is no missingness in observed cells (ex Heijden 2018) there exists a closed-form solution. A general closed-form solution may exist for other scenarios (ex Heijden 2021), but I have not derived it at this time.
