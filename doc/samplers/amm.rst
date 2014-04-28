.. index:: Sampling Functions; Adaptive Multivariate Metropolis

Adaptive Multivariate Metropolis (AMM)
--------------------------------------

Implementation of the Roberts and Rosenthal :cite:`robert:2009:EAM` adaptive multivariate Metropolis sampler :cite:`haario:2001:AMA`, :cite:`hastings:1970:MCS`, :cite:`metropolis:1953:ESC` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: amm!(v::VariateAMM, SigmaF::Cholesky{Float64}, logf::Function; \
                adapt::Bool=false)

	Simulate one draw from a target distribution using an adaptive multivariate Metropolis sampler.
	
	**Arguments**
	
		* ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function should contain the ``tune`` field returned by the previous call.
		* ``SigmaF`` : Cholesky factorization of the covariance matrix for the non-adaptive multivariate normal proposal distribution.
		* ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.
		* ``adapt`` : whether to adaptively update the proposal distribution.
		
	**Value**
	
		Returns ``v`` updated with the simulated values and associated tuning parameters.
	
	**Example**
	
		.. literalinclude:: amm.jl
			:language: julia
					
.. index:: VariateAMM

VariateAMM Type
^^^^^^^^^^^^^^^

Declaration
```````````

``VariateAMM <: VariateVector``

Fields
``````

* ``value::Vector{VariateType}`` : vector of sampled values.
* ``tune::TuneAMM`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: VariateAMM(x::Vector{VariateType}, tune::TuneAMM)
              VariateAMM(x::Vector{VariateType}, tune=nothing)

	Construct a ``VariateAMM`` object that stores values and tuning parameters for adaptive multivariate Metropolis sampling.
	
	**Arguments**
	
		* ``x`` : vector of sampled values.
		* ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.
		
	**Value**
	
		Returns a ``VariateAMM`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.
		

.. index:: TuneAMM

TuneAMM Type
^^^^^^^^^^^^

Declaration
```````````

``type TuneAMM``

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

.. function:: SamplerAMM(params::Vector{T<:String}, Sigma::Matrix{U:<Real}; \
				adapt::Symbol=:none)

	Construct an ``MCMCSampler`` object for adaptive multivariate Metropolis sampling.
	
	**Arguments**
	
		* ``params`` : named stochastic nodes to be updated with the sampler.
		* ``Sigma`` : covariance matrix for the non-adaptive multivariate normal proposal distribution.
		* ``adapt`` : type of adaptation phase.  Options are
			* ``:all`` : adapt proposal during all iterations.
			* ``:burnin`` : adapt proposal during burn-in iterations.
			* ``:none`` : no adaptation (multivariate Metropolis sampling with fixed proposal).

	**Value**
	
		Returns an ``MCMCSampler`` type object.