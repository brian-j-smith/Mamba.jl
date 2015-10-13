Mamba: Markov chain Monte Carlo (MCMC) for Bayesian analysis in julia
=====================================================================

:Version: |release|
:Requires: julia 0.4
:Date: |today|
:Maintainer: Brian J Smith (brian-j-smith@uiowa.edu)
:Contributors: Benjamin Deonovic (benjamin-deonovic@uiowa.edu), Brian J Smith (brian-j-smith@uiowa.edu), and `others <https://github.com/brian-j-smith/Mamba.jl/contributors>`_
:Web site: https://github.com/brian-j-smith/Mamba.jl
:License: `MIT <https://github.com/brian-j-smith/Mamba.jl/blob/master/LICENSE.md>`_

Overview
--------

Purpose
^^^^^^^

`Mamba` is an open platform for the implementation and application of MCMC methods to perform Bayesian analysis in `julia <http://julialang.org/>`_.  The package provides a framework for (1) specification of hierarchical models through stated relationships between data, parameters, and statistical distributions; (2) block-updating of parameters with samplers provided, defined by the user, or available from other packages; (3) execution of sampling schemes; and (4) posterior inference.  It is intended to give users access to all levels of the design and implementation of MCMC simulators to particularly aid in the development of new methods.

Several software options are available for MCMC sampling of Bayesian models.  Individuals who are primarily interested in data analysis, unconcerned with the details of MCMC, and have models that can be fit in `JAGS <http://mcmc-jags.sourceforge.net/>`_, `Stan <http://mc-stan.org/>`_, or `OpenBUGS <http://www.openbugs.net/>`_ are encouraged to use those programs.  `Mamba` is intended for individuals who wish to have access to lower-level MCMC tools, are knowledgeable of MCMC methodologies, and have experience, or wish to gain experience, with their application.  The package also provides stand-alone convergence diagnostics and posterior inference :ref:`tools <section-Chains>`, which are essential for the analysis of MCMC output regardless of the software used to generate it.

Features
^^^^^^^^

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
    * `Gadfly <https://github.com/dcjones/Gadfly.jl>`_ plotting: trace, density, running mean, autocorrelation.
    * Run-time performance on par with compiled MCMC software.

Getting Started
^^^^^^^^^^^^^^^

The following **julia** command will install the package:

.. code-block:: julia

    julia> Pkg.add("Mamba")


Contents
--------

.. toctree::
    :maxdepth: 2

    intro.rst
    tutorial.rst
    variate.rst
    mcmc.rst
    samplers.rst
    examples.rst
    discussion.rst
    supplement.rst
    xrefs.rst


Indices
^^^^^^^

* :ref:`genindex`
