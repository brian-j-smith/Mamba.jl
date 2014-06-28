.. index:: Sampling Functions; Adaptive Mixture Metropolis

Adaptive Mixture Metropolis (AMM)
---------------------------------

Implementation of the Roberts and Rosenthal :cite:`robert:2009:EAM` adaptive (multivariate) mixture Metropolis :cite:`haario:2001:AMA,hastings:1970:MCS,metropolis:1953:ESC` sampler for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: amm!(v::AMMVariate, SigmaF::Cholesky{Float64}, logf::Function; \
                adapt::Bool=true)

	Simulate one draw from a target distribution using an adaptive mixture Metropolis sampler.  Parameters are assumed to be continuous and unconstrained.
	
	**Arguments**
	
		* ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function should contain the ``tune`` field returned by the previous call.
		* ``SigmaF`` : Cholesky factorization of the covariance matrix for the non-adaptive multivariate normal proposal distribution.
		* ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.
		* ``adapt`` : whether to adaptively update the proposal distribution.
		
	**Value**
	
		Returns ``v`` updated with simulated values and associated tuning parameters.
		
	.. _example-amm:
	
	**Example**
	
		.. literalinclude:: amm.jl
			:language: julia
					
.. index:: AMMVariate

AMMVariate Type
^^^^^^^^^^^^^^^

Declaration
```````````

``AMMVariate <: VectorVariate``

Fields
``````

* ``value::Vector{VariateType}`` : vector of sampled values.
* ``tune::AMMTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: AMMVariate(x::Vector{VariateType}, tune::AMMTune)
              AMMVariate(x::Vector{VariateType}, tune=nothing)

	Construct a ``AMMVariate`` object that stores sampled values and tuning parameters for adaptive mixture Metropolis sampling.
	
	**Arguments**
	
		* ``x`` : vector of sampled values.
		* ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.
		
	**Value**
	
		Returns a ``AMMVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.
		

.. index:: AMMTune

AMMTune Type
^^^^^^^^^^^^

Declaration
```````````

``type AMMTune``

Fields
``````
* ``adapt::Bool`` : whether the proposal distribution has been adaptively tuned. 
* ``beta::Real`` : proportion of weight given to draws from the non-adaptive proposal with covariance factorization ``SigmaF``, relative to draws from the adaptively tuned proposal with covariance factorization ``SigmaLm``, during adaptive updating.  Fixed at ``beta = 0.05``.
* ``m::Integer`` : number of adaptive update iterations that have been performed.
* ``Mv::Vector{Float64}`` : running mean of draws ``v`` during adaptive updating.  Used in the calculation of ``SigmaLm``.
* ``Mvv::Vector{Float64}`` : running mean of ``v * v'`` during adaptive updating.  Used in the calculation of ``SigmaLm``.
* ``scale::Real`` : fixed value ``2.38^2`` in the factor (``scale / length(v)``) by which the adaptively updated covariance matrix is scaled---adopted from Gelman, Roberts, and Gilks :cite:`gelman:1996:EMJ`.
* ``SigmaF::Cholesky{Float64}`` : factorization of the non-adaptive covariance matrix.
* ``SigmaLm::Matrix{Float64}`` : lower-triangular factorization of the adaptively tuned covariance matrix.

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: AMM(params::Vector{Symbol}, Sigma::Matrix{T<:Real}; \
				adapt::Symbol=:all)

	Construct an ``MCMCSampler`` object for adaptive mixture Metropolis sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.
	
	**Arguments**
	
		* ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``link()`` function.
		* ``Sigma`` : covariance matrix for the non-adaptive multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.
		* ``adapt`` : type of adaptation phase.  Options are
			* ``:all`` : adapt proposal during all iterations.
			* ``:burnin`` : adapt proposal during burn-in iterations.
			* ``:none`` : no adaptation (multivariate Metropolis sampling with fixed proposal).

	**Value**
	
		Returns an ``MCMCSampler`` type object.

	**Example**
	
		See the :ref:`section-Examples` section.
