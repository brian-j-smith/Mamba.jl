.. index:: Sampling Functions; Adaptive Metropolis within Gibbs

Adaptive Metropolis within Gibbs (AMWG)
---------------------------------------

Implementation of a Metropolis-within-Gibbs sampler :cite:`metropolis:1953:ESC`, :cite:`robert:2009:EAM`, :cite:`tierney:1994:MCE` for iteratively simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: amwg!(v::VariateAMWG, sigma::Vector{Float64}, logf::Function; \
                adapt::Bool=false, batchsize::Integer=50, target::Real=0.44)

	Simulate one draw from a target distribution using an adaptive Metropolis-within-Gibbs sampler.
	
	**Arguments**
	
		* ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function should contain the ``tune`` field returned by the previous call.
		* ``sigma`` : initial standard deviations for the univariate normal proposal distributions.
		* ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.data``.
		* ``adapt`` : whether to adaptively update the proposal distribution.
		* ``batchsize`` : number of samples that must be accumulated before applying an adaptive update to the proposal distributions.
		* ``target`` : a target acceptance rate for the adaptive algorithm.
		
	**Value**
	
		Returns ``v`` updated with the simulated values and associated tuning parameters.
	
	**Example**

		.. literalinclude:: amwg.jl
			:language: julia
			
.. index:: VariateAMWG

VariateAMWG Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``VariateAMWG <: VariateVector``

Fields
``````

* ``data::Vector{VariateType}`` : vector of sampled values.
* ``tune::TuneAMWG`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: VariateAMWG(x::Vector{VariateType}, tune::TuneAMWG)
              VariateAMWG(x::Vector{VariateType}, tune=nothing)

	Construct a ``VariateAMWG`` object that stores values and tuning parameters for adaptive Metropolis-within-Gibbs sampling.
	
	**Arguments**
	
		* ``x`` : vector of sampled values.
		* ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.
		
	**Value**
	
		Returns a ``VariateAMWG`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

		
.. index:: TuneAMWG

TuneAMWG Type
^^^^^^^^^^^^^

Declaration
```````````

``type TuneAMWG``

Fields
``````

* ``adapt::Bool`` : whether the proposal distribution has been adaptively tuned.
* ``accept::Vector{Integer}`` : number of accepted candidate draws generated for each element of the parameter vector during adaptive updating.
* ``batchsize::Integer`` : number of samples that must be accumulated before applying an adaptive update to the proposal distributions.
* ``m::Integer`` : number of adaptive update iterations that have been performed.
* ``sigma::Vector{Float64}`` : updated values of the proposal standard deviations if ``adapt = true``, and the user-defined values otherwise.
* ``target::Real`` : target acceptance rate for the adaptive algorithm.

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SamplerAMWG(params::Vector{T<:String}, sigma::Vector{U<:Real}; \
				adapt::Symbol=:none, batchsize::Integer=50, target::Real=0.44)

	Construct an ``MCMCSampler`` object for adaptive Metropolis-within-Gibbs sampling.
	
	**Arguments**
	
		* ``params`` : named stochastic nodes to be updated with the sampler.
		* ``sigma`` : initial standard deviations for the univariate normal proposal distributions.
		* ``adapt`` : type of adaptation phase.  Options are
			* ``:all`` : adapt proposals during all iterations.
			* ``:burnin`` : adapt proposals during burn-in iterations.
			* ``:none`` : no adaptation (Metropolis-within-Gibbs sampling with fixed proposals).
		* ``batchsize`` : number of samples that must be accumulated before applying an adaptive update to the proposal distributions.
		* ``target`` : a target acceptance rate for the algorithm.

	**Value**
	
		Returns an ``MCMCSampler`` type object.
