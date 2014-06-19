Conclusion
----------

`Mamba` is a platform for the development and application of MCMC simulators for Bayesian modelling.  Such simulators can be difficult to implement in practice.  `Mamba` eases that task by standardizing and automating the generation of initial values, specification of distributions, iteration of Gibbs steps, updating of parameters, and running of MCMC chains.  It automatically evaluates (unnormalized) full conditionals and allows MCMC simulators to be implemented by simply stating relationships between data, parameters, and statistical distributions, similar to the 'BUGS' clones and Stan program.  In general, the package is designed to give users access to all levels of MCMC design and implementation.  To that end, its toolset includes: 1) a model specification syntax, 2) stand-alone and integrated sampling functions, 3) a simulation engine and API, and 4) functions for convergence diagnostics and posterior inference.  Moreover, its tools are designed to be modular so that they can be combined, extended, and used in ways that best meet users' needs.  

`Mamba` can accommodate a wide range of model formulations as well as combinations of user-defined, supplied, and external samplers and distributions. It handles routine implementation tasks, thus allowing users to focus on design issues and making the package well-suited for

	* developing new Bayesian models,
	* implementing simulators for classes of models,
	* testing new sampling algorithms,
	* prototyping MCMC schemes for software development, and
	* teaching MCMC methods.

Furthermore, output is easily generated with the simulation engine in a standardized format that can be analyzed directly with supplied or user-defined tools for convergence diagnostics and inference.  Use of the package is illustrated with several examples.  Future plans include additional sampler implementations and optimizations, development of alternative model-specification interfaces, and automatic specification of sampling schemes.  The software is freely available and licensed under the open-source MIT license.  It is hoped that the package will foster MCMC methods developed by researchers in the field and make their methods more accessible to a broader scientific community.
