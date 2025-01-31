# DGP-LVM: Derivative Gaussian process latent variable model
This is a repository for the codes used in the paper: DGP-LVM (https://arxiv.org/abs/2404.04074). The codes are arranged as follows.
## Models
DGP-LVM model development files for derivative squared exponential (SE), Matern 3/2 (M32) and Matern 5/2 (M52) covariance functions can be found using the single Stan file 'DerivGPmodels.stan'. Model development is perfomed using Stan (https://mc-stan.org/).
## Simulation study
Simulation study designed to validate and compare performance of DGP-LVM and other GP models in estimating latent input variables. In the manuscipt, five simulation scenarios have been considered: data generated from SE, Matern 3/2, Matern 5/2, periodic data and periodic with trend data. The general simulation study R script 'dgplvmSimStudy.R' can be used alongside 'DgplvmSimFns.R' to call different data generating processes and other functions involved in the Sim Study.
## Model convergence
This is to check model convergence using MCMC convergence diagnostics. The plots have been generated using the associated R script 'mcmcConvergenceDiag.R'. 
## Simulation results
Code for generating the plots showing recovery of ground truth latent inputs (see 'ModelEvaluationLatentX.R') along with model hyperparameter (see 'ModelEvaluationHyperparameters.R') recovery.
## Case study
DGP-LVM is used to analyze a reduced single-cell RNA sequencing data. The data is a pre-processed Cell Cycle data obtained from Cytopath (https://doi.org/10.1016/j.crmeth.2022.100359). The corresponding codes can be found in 'CaseStudy.R'
