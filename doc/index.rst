Mamba: Markov chain Monte Carlo for Bayesian analysis in julia (version |release|)
==================================================================================

:Version: |release|
:Date: |today|
:Maintainer: Brian J Smith (brian-j-smith@uiowa.edu)
:Contributors: Benjamin Deonovic (benjamin-deonovic@uiowa.edu), Brian J Smith (brian-j-smith@uiowa.edu)
:Web site: https://github.com/brian-j-smith/Mamba.jl
:License: `MIT <https://github.com/brian-j-smith/Mamba.jl/blob/master/LICENSE.md>`_

Purpose
-------

`Mamba` is a `julia <http://julialang.org/>`_ programming environment and toolset for the implementation and inference of Bayesian models using Markov chain Monte Carlo (MCMC) sampling.  The package provides a framework for (1) specification of hierarchical models through stated relationships between data, parameters, and statistical distributions; (2) block-updating of parameters with samplers provided, defined by the user, or available from other packages; (3) execution of MCMC sampling schemes, and (4) posterior inference.  It is designed to give users access to all levels of the design and implementation of MCMC simulators to particularly aid in the development of complex models.

The package aims to provide:

	* An intuitive, interactive, and extensible interface.
	* Support for a wide range of model and distributional specifications.
	* An environment in which all interactions with the software are made through a single, interpreted programming language.
	
		* Any **julia** operator, function, type, or package can be used for model specification.
		* Custom distributions and samplers can be written in **julia** to extend the package.
		
	* Directed acyclic graph representations of models.
	* Arbitrary blocking of model parameters and designation of block-specific samplers.
	* Samplers that can be used with the included simulation engine or apart from it.
	* Command-line access to all package functionality, including its simulation API.
	* Tools for convergence diagnostics and posterior inference.
	* Run-time performance on par with compiled MCMC software.

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
	conclusion.rst
	supplement.rst
	xrefs.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

