Introduction
============

MCMC Software
-------------

Markov chain Monte Carlo (MCMC) methods are a class of algorithms for simulating autocorrelated draws from probability distributions :cite:`brooks:2011:HMC,gamerman:1997:MCM,gilks:1996:MCP,robert:2004:MCS`.  They are widely used to obtain empirical estimates for and make inference on multidimensional distributions that often arise in Bayesian statistical modeling, computational physics, and computational biology.  Because MCMC provides estimates of *distributions* of interest, rather than *point* estimates and asymptotic standard errors, it facilitates wide ranges of inferences and provides for more realist prediction errors.  An MCMC algorithm can be devised for any probability model.  Implementations of algorithms are computational in nature, with the resources needed to execute algorithms directly related to the dimensionality of their associated problems.  Rapid increases in computing power and emergence of MCMC software have enabled models of increasing complexity to be fit.  For all its advantages, MCMC is considered to be one of the most important developments and powerful tools in modern statistical computing.

Several software tools have been developed for MCMC simulation of the posterior distributions of Bayesian model parameters.  Tools range from those designed for general model fitting to those for specific models.  *WinBUGS*, its open-source incarnation *OpenBUGS*, and the 'BUGS' clone Just Another Gibbs Sampler (*JAGS*) are among the most widely used programming tools for general model fitting :cite:`lunn:2009:BUGS,plummer:2003:JAGS`.  These three provide similar programming syntaxes for users to specify a wide array of statistical models by simply stating relationships between data, parameters, and statistical distributions.  Once a model is specified, the programs automatically formulate an MCMC sampling scheme with which to simulate parameter values from their posterior distribution.  *OpenBUGS* additionally employs an expert system for and provides manual control over its scheme.  All aforementioned tasks can be accomplished with minimal programming and without any specific knowledge of MCMC methodology.  Users who are adept at both and so inclined can write software modules to add new distributions and samplers to *OpenBUGS* and *JAGS* :cite:`thomas:2014:ODM,wabersich:2013:EJT`.  *Stan* is another open-source program worth noting for its accessible syntax and automatically tuned Hamiltonian Monte Carlo sampling scheme :cite:`stan-software:2014`.  Additionally, several **R** packages are available for more specific model fitting.  The *arm* package provides Bayesian inference for generalized linear, ordered logistic or probit, and mixed-effects regression models :cite:`gelman:2014:arm`, *MCMCpack* fits a wide range of models commonly encountered in the social and behavioral sciences :cite:`martin:2013:MCP`, and many others that are more focused on specific classes of models can be found in the "Bayesian Inference" task view on the Comprehensive **R** Archive Network :cite:`park:2014:cran`.

MCMCsim Package
---------------

*MCMCsim* :cite:`smith:2014:MCMC` is a **julia** :cite:`bezanson:2012:JFD,julia:2014` package that supports native Gibbs sampling for general Bayesian model fitting.  Like *OpenBUGS* and *JAGS*, it supports Gibbs sampling for a wide range of statistical models.   *MCMCsim* provides a set of **julia** types and method functions that allow users to create their own MCMC schemes and sampling functions, while handling tasks that are common to all implementations.  The package is similar in purpose to the *GRIMS* (General R Interface for Markov Sampling) program :cite:`neal:2012:grims`, but additionally provides a framework for model specifications with a method-based and lower-level toolset.  Its target audience includes individuals who are familiar with the **julia** programming language and wish to develop MCMC samplers in that environment; who have specific models or classes of models to implement; and, in some cases, are able to derive full conditional distributions of model parameters (up to normalizing constants).  One advantage of the package for **julia** users is that no other languages are needed for development; whereas, extensions to *OpenBUGS* and *JAGS* require compiled languages which present different sets of development, testing, and debugging challenges.  Another advantage is the flexibility of using sampling functions defined by the user, supplied by the package, and available in other packages; thus providing easy access to any current and future **julia** operators, functions, and packages.

*MCMCsim* is intended for the user who has formulated a Bayesian model in terms of parameters :math:`(\theta_1, \ldots, \theta_p)`, and wishes to implement an MCMC sampling scheme to simulate draws from the joint posterior distribution.  The package supports the general Gibbs :cite:`gelfand:1990:SBA,geman:1984:SRG` scheme outlined in Algorithm 1.  In its implementation with the package, the user may specify any blocking :math:`\{\Theta_j\}_{j=1}^{B}` of the parameters and corresponding functions :math:`\{f_j\}_{j=1}^{B}` to sample each :math:`\Theta_j` from its full conditional distribution :math:`p(\Theta_j | \Theta \setminus \Theta_{j})`.  Simulation performance (efficiency and runtime) can be affected greatly by the choice of blocking scheme and sampling functions.  For some models, an optimal choice may not be obvious, and different choices may need to be tried to find one that gives a desired level of performance.  This can be a time-consuming process.  The *MCMCsim* package provides a set of **julia** classes and method functions to facilitate the specification of different schemes and functions.  Supported sampling functions include those provided by the package, user-defined functions, and functions from other packages; thus providing great flexibility with respect to sampling methods.  Furthermore, a sampling engine is provided to save the user from having to implement tasks common to all MCMC simulators.  Therefore, time and energy can be focused on implementation aspects that most directly affect performance.

.. _figure-Gibbs:

.. figure:: images/gibbs.png
	:align: center

	*MCMCsim* Gibbs sampling scheme.
	
A summary of the steps involved in using the package to perform MCMC simulation for a Bayesian model is given below.

	#. Decide on names to use for **julia** objects that will represent the model data structures and parameters (:math:`\theta_1, \ldots, \theta_p`).  For instance, the :ref:`section-Line` section describes a linear regression example in which predictor :math:`\bm{x}` and response :math:`\bm{y}` are represented by objects ``x`` and ``y``, and regression parameters :math:`\beta_0`, :math:`\beta_1`, and :math:`\sigma^2` by objects ``b0``, ``b1``, and ``s2``.

	#. Create a dictionary to store all structures considered to be fixed in the simulation; e.g., the ``line`` dictionary in the regression example.

	#. Specify the model using the constructors described in the :ref:`section-MCMC-Types` section, to create the following:
 
		a. An ``MCMCDependent`` object for each model term that has a distributional specification.  This includes parameters and data, such as the regression parameters ``b0``, ``b1``, and ``s2`` that have prior distributions and ``y`` that has a likelihood specification.

		b. A vector of ``MCMCSampler`` objects containing supplied, user-defined, or external functions :math:`\{f_j\}_{j=1}^{B}` for sampling each parameter block :math:`\Theta_j`.

		c. An ``MCMCModel`` object from the resulting nodes and sampler vector.

	#. Simulate parameter values with the :func:`mcmc` function.
	
	#. Use the MCMC output to perform convergence checks and posterior inference.
