README for the files in "Mixture M-RA" Folder
Written by: Marco H. Benedetti
Last updated: October 31, 2020.

This README file describes the code provided alongside the paper "Identifying regions of inhomogeneities in spatial processes via an M-RA and mixture priors" by Marco H. Benedetti, Veronica J. Berrocal, and Naveen N. Narisetty.

The code herein is intended to provide users with a straightforward simulation study that allows them to fit the Mixture M-RA to randomly generated spatial data.

The simulatied data corresponds to Simulation Study 2 in the manuscript, in which data are generated as a realization of a mixture of Gaussian Processes:
In the left portion of the domain, the data are a realization of a Gaussian process with Matern covariance and range parameter equal to 1.  This is a process in which spatial dependence persists at large distances.
In the right portion of the domain, the data are a realization of a Gaussian process with Matern covariance and range parameter equal to 0.01.  This is a process in which spatial dependence decays rapidly.
See Figure 2 (a) for an illustration of these data.

We expect basis function weights corresponding to knots in the left part of the domain to shrink towards zero at higher levels.
In other words we expect that posterior expectations of the latent binary variables, i.e. E(Z_{m,j}|Y), will be close to zero.

We expect basis function weights corresponding to knots in the right part of the domain to remain active in the model, even at higher levels.
In other words we expect that posterior expectations of the latent binary variables, i.e. E(Z_{m,j}|Y), will be close to one.

Files:

1) materncovcode.cpp: This is a cpp script used to compute matern covariance between two points as a function of their distance and a set of covariance function parameters.
This program takes a distance matrix as as input.  If you use a different covariance function, please note in the event of any errors that the inputs must be the same as maternCovcpp.

2) Run before MCMC code.R: This R script contains all functions needed to prepare to run the Mixture M-RA.
A number of function were written by Matthias Katzfus for the paper "A multi-resolution approximation for massive spatial datasets."
In addition, we provide MCMC functions to sample from the posterior distribution of various model parameters.

3) Example Simulation.R: Provides data generation and calls MCMC functions to fit the Mixture MRA to randomly generated non-stationary spatial data.


General Notes:
1) Both R files are documented, and we encourage users to refer to these notes as you run the code.
2) You can save time by proposing values of phi and nu (the range and smoothness parameters respectively) from a discrete set of proposal values.
This way, you can pre-compute basis functions before running MCMC, rather than computing new basis functions at each iteration.
Another approach is to fix the covariance function parameters apriori.  However, we have not looked into this in the course of working on this paper, and cannot speak to how this would effect inference.
3) Note that calling partition.2d re-orders the data.  Don't use your original data (y) after running partition.2d.  Instead, use y.new, as we have done in Example Simulation.R
