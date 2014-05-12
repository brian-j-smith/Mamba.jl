Conclusion
----------

`MCMCsim` provides a range of tools for the implementation of MCMC simulators in **julia**.  Tools include: 1) a modelling language, 2) sampling functions, 3) a sampling engine and API, and 4) convergence diagnostics and statistical summaries.  The supplied sampling functions ``amm!``, ``amwg!``, ``nuts!``, ``slice!``, and ``slicewg!`` can be used alone, as illustrated in their documentation, or with other the package's simulation engine.  The engine applied to models defined in terms of the ``MCMCDependent``, ``MCMCModel``, and ``MCMCSampler`` types standardizes and automates the generation of initial values, specification of distributions, iteration of Gibbs steps, updating of parameters, and running of MCMC chains.  Automatic calculation of unnormalized full conditions allows MCMC simulators to be implemented by simply stating relationships between data, parameters, and statistical distributions, similar to the 'BUGS' clones and Stan program.  Futhermore, the `MCMCsim` tools are designed to be modular so that they can be combined, extended, and used in ways that best meet users' needs.

The package can accommodate a wide range of model formulations as well as combinations of user-defined, supplied, and external samplers and distributions. It handles routine implementation tasks, thus allowing users to focus on design issues, and making the package well-suited for

	* developing new Bayesian models,
	* implementing samplers for classes of models,
	* testing new sampling algorithms,
	* prototyping MCMC schemes for software development, and
	* teaching courses in MCMC methods.

Moreover, output is easily generated with the simulation engine in a standardized format that can be analyzed directly with supplied or user-defined tools for convergence diagnostics and inference.  Use of the package was illustrated several example.  Future plans for the package include additional samplers, development of alternative model-specification interfaces, and automatic specification of sampling schemes.  The software is freely available and licensed under the open-source MIT license.  It is hoped that the package can serve as a platform to help bring together MCMC methods developed by researchers in the field and make their methods more accessible to a broader scientific community.
