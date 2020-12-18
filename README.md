* Organization
  - There are three directories, each correspond to a subsection in the experiments section of the paper: simulation (simulation experiment), movielens (MovieLens data), us_natality (US Natality data).
  - Each directory has four sub directories: code, data, qsub, and result.
  - Directory 'code' has files of all the code (R source code) that was used in the analysis. 
  - Directory 'data' has (if any) simulated data or real data that were used in the analysis. 
  - Directory 'qsub' has SGE files (.q) that were used to submit jobs on a SGE cluster. 
  - Directory 'result' has a sub directory 'img' and stores the result (if any) produced in the analysis. This directory may be empty or absent.

* Files
  - 'simulate_data.R' contains the code to simulate and partition the data. 
  - 'partition_data.R' contains the code to partition the simulated or real data. 
  - 'analyze_result.R' contains the code for analyzing the results of DPMC, WASP, and competing methods and making plots/tables.
  - 'dls_sampler.R' contains the code for the MCMC/Gibbs sampler of a subset posterior distribution. This is a modified version of the code in 'mcmc_sampler.R' using stochastic approximation.
  - 'dls_comb.R' contains the code for combining the subset posterior distribution.
  - accuracy.R' contains the code for measuring approximation errors. 
  - 'submit.R' contains the  code for the R code for submitting a job on the cluster. The files in 'qsub' directory use this file for running simulations. If you want to replicate our results, then you should see this file. 

* Citation:
  If you use the code, then please cite the following paper:
  - Srivastava, S. and Xu, Y.(2020+). Distributed Bayesian Inference for Linear Mixed-Effects Models. Journal of Computational and Graphical Statistics. 

* Contact:
  Please email Sanvesh Srivastava (<sanvesh-srivastava@uiowa.edu>) or Yixiang Xu (<yixiang-xu@berkeley.edu>) if you have any questions related to the code.

* Acknowledgment
  - Some code for linear mixed effects modeling has been borrowed from Patrick O. Perry (<http://ptrckprry.com/code/>).
  - Please email us if you think that we have missed citations to your paper/work. 
  - This research is partially supported by grants from the Office of Naval Research (ONR-BAA N000141812741) and the National Science Foundation (DMS-1854667/1854662). 
