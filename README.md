# Mamba: Markov chain Monte Carlo (MCMC) for Bayesian analysis in julia

julia 0.4: [![Docs](https://readthedocs.org/projects/mambajl/badge/?version=release-0.6)](http://mambajl.readthedocs.org/en/release-0.6/)
[![Mamba](http://pkg.julialang.org/badges/Mamba_0.4.svg)](http://pkg.julialang.org/?pkg=Mamba&ver=0.4)
[![Build Status](https://travis-ci.org/brian-j-smith/Mamba.jl.svg?branch=release-0.6)](https://travis-ci.org/brian-j-smith/Mamba.jl)  
julia 0.3: [![Docs](https://readthedocs.org/projects/mambajl/badge/?version=release-0.4)](http://mambajl.readthedocs.org/en/release-0.4/)
[![Mamba](http://pkg.julialang.org/badges/Mamba_0.3.svg)](http://pkg.julialang.org/?pkg=Mamba&ver=0.3)
[![Build Status](https://travis-ci.org/brian-j-smith/Mamba.jl.svg?branch=release-0.4)](https://travis-ci.org/brian-j-smith/Mamba.jl)

## Purpose

*Mamba* is a [julia](http://julialang.org/) programming environment and toolset for the implementation and inference of Bayesian models using MCMC sampling.  The package provides a framework for (1) specification of hierarchical models through stated relationships between data, parameters, and statistical distributions; (2) block-updating of parameters with samplers provided, defined by the user, or available from other packages; (3) execution of sampling schemes; and (4) posterior inference.  It is designed to give users access to all levels of the design and implementation of MCMC simulators to particularly aid in the development of complex models.

Several software options are available for MCMC sampling of Bayesian models.  Individuals who are primarily interested in data analysis, unconcerned with the details of MCMC, and have models that can be fit in [JAGS](http://mcmc-jags.sourceforge.net/), [Stan](http://mc-stan.org/), or [OpenBUGS](http://www.openbugs.net/) are encouraged to use those programs.  *Mamba* is intended for individuals who wish to have access to lower-level MCMC tools, are knowledgeable of MCMC methodologies, and have experience, or wish to gain experience, with their application.  The package also provides stand-alone convergence diagnostics and posterior inference tools, which are essential for the analysis of MCMC output regardless of the software used to generate it. 

## Features

* An interactive and extensible interface.
* Support for a wide range of model and distributional specifications.
* An environment in which all interactions with the software are made through a single, interpreted programming language.
    * Any **julia** operator, function, type, or package can be used for model specification.
    * Custom distributions and samplers can be written in **julia** to extend the package.
* Directed acyclic graph representations of models.
* Arbitrary blocking of model parameters and designation of block-specific samplers.
* Samplers that can be used with the included simulation engine or apart from it, including Slice, adaptive multivariate Metropolis, adaptive Metropolis within Gibbs, and No-U-Turn (Hamiltonian Monte Carlo) samplers.
* Automatic parallel execution of parallel MCMC chains on multi-processor systems.
* Restarting of chains.
* Command-line access to all package functionality, including its simulation API.
* Convergence diagnostics: Gelman, Rubin, and Brooks; Geweke; Heidelberger and Welch; Raftery and Lewis.
* Posterior summaries: moments, quantiles, HPD, cross-covariance, autocorrelation, MCSE, ESS.
* [Gadfly](https://github.com/dcjones/Gadfly.jl) plotting: trace, density, running mean, autocorrelation.
* Run-time performance on par with compiled MCMC software.

## Getting Started

The following **julia** command will install the package:

```julia
julia> Pkg.add("Mamba")
```
