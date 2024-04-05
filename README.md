# DGP-LVM: Derivative Gaussian process latent variable model
This is a repository for code availability for DGP-LVM. The codes are arranged as follows.
# Models
DGP-LVM model development files for derivative squared exponential (SE), Matern 3/2 (M32) and Matern 5/2 (M52) covariance functions. Model development is perfomed using Stan (https://mc-stan.org/).
# Simulation study
Simulation study designed to validate and compare performance of DGP-LVM and other GP models in estimating latent input variables. Two simulation scenarios have been considered: data generated from a GP and data generated from a periodic function.
# Model convergence
This is to check model convergence using MCMC convergence diagnostics. The plots have been generated using the associated code.
# Simulation results
Code for generating the plots showing recovery of ground truth latent inputs along with model hyperparameter recovery.
# Case study
DGP-LVM is used to analyze a reduced single-cell RNA sequencing data. The data is a pre-processed Cell Cycle data obtained from Cytopath (https://doi.org/10.1016/j.crmeth.2022.100359).