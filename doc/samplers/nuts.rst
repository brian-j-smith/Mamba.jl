.. index:: Sampling Functions; No-U-Turn Sampler

No-U-Turn Sampler (NUTS)
------------------------

Implementation of the NUTS extension :cite:`hoffman:2011:nuts` to Hamiltonian Monte Carlo :cite:`neal:2011:hmc` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: nutseps(v::VariateNUTS, fx::Function)
	
	Generate an initial value for the step size parameter of the No-U-Turn sampler.
	
	**Arguments**
	
		* ``v`` : the current state of parameters to be simulated.
		* ``fx`` : function to compute the log-transformed density (up to a normalizing constant) and gradient vector at ``v.data``, and to return the respective results as a tuple.
		
	**Value**
	
		A numeric step size value.

.. function:: nuts!(v::VariateNUTS, eps::Real, fx::Function; adapt::Bool=false, \
                target::Real=0.6)

	Simulate one draw from a target distribution using the No-U-Turn sampler.
	
	**Arguments**
	
		* ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function should contain the ``tune`` field returned by the previous call.
		* ``eps`` : the NUTS algorithm step size parameter.
		* ``fx`` : function to compute the log-transformed density (up to a normalizing constant) and gradient vector at ``v.data``, and to return the respective results as a tuple.
		* ``adapt`` : whether to adaptively update the ``eps`` step size parameter.
		* ``target`` : a target acceptance rate for the algorithm.
		
	**Value**
	
		Returns ``v`` updated with the simulated values and associated tuning parameters.
	
	**Example**

		.. literalinclude:: nuts.jl
			:language: julia

.. index:: VariateNUTS

VariateNUTS Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``VariateNUTS <: VariateVector``

Fields
``````

* ``data::Vector{VariateType}`` : vector of sampled values.
* ``tune::TuneNUTS`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: VariateNUTS(x::Vector{VariateType}, tune::TuneNUTS)
              VariateNUTS(x::Vector{VariateType}, tune=nothing)

	Construct a ``VariateNUTS`` object that stores values and tuning parameters for No-U-Turn sampling.
	
	**Arguments**
	
		* ``x`` : vector of sampled values.
		* ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.
		
	**Value**
	
		Returns a ``VariateNUTS`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.


.. index:: TuneNUTS

TuneNUTS Type
^^^^^^^^^^^^^

Declaration
```````````

``type TuneNUTS``

Fields
``````
* ``adapt::Bool`` : whether the proposal distribution has been adaptively tuned.
* ``alpha::Float64`` : cumulative acceptance probabilities :math:`\alpha` from leapfrog steps.
* ``eps::Float64`` : updated value of the step size parameter :math:`\epsilon_m = \exp\left(\mu - \sqrt{m} \bar{H}_m / \gamma\right)` if ``adapt = true``, and the user-defined value otherwise.
* ``epsbar::Float64`` : dual averaging parameter, defined as :math:`\bar{\epsilon}_m = \exp\left(m^{-\kappa} \log(\epsilon_m) + (1 - m^{-\kappa}) \log(\bar{\epsilon}_{m-1})\right)`.
* ``gamma::Float64`` : dual averaging parameter, fixed at :math:`\gamma = 0.05`.
* ``Hbar::Float64`` : dual averaging parameter, defied as :math:`\bar{H}_m = \left(1 - \frac{1}{m + t_0}\right) \bar{H}_{m-1} + \frac{1}{m + t_0} \left(\text{target} - \frac{\alpha}{n_\alpha}\right)`.
* ``kappa::Float64`` : dual averaging parameter, fixed at :math:`\kappa = 0.05`.
* ``m::Integer`` : number of adaptive update iterations :math:`m` that have been performed.
* ``mu::Float64`` : dual averaging parameter, defined as :math:`\mu = \log(10 \epsilon_0)`.
* ``nalpha::Integer`` : the total number :math:`n_\alpha` of leapfrog steps performed.
* ``t0::Float64`` : dual averaging parameter, fixed at :math:`t_0 = 10`.
* ``target::Float64`` : target acceptance rate for the adaptive algorithm.

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SamplerNUTS(params::Vector{T<:String}; dtype::Symbol=:forward, \
				target::Real=0.6)

	Construct an ``MCMCSampler`` object for No-U-Turn sampling, with the algorithm's step size parameter adaptively tuned during burn-in iterations.
	
	**Arguments**
	
		* ``params`` : named stochastic nodes to be updated with the sampler.
		* ``dtype`` : type of differentiation for gradient calculations.  Options are
			* ``:central`` : central differencing.
			* ``:forward`` : forward differencing.
		* ``target`` : a target acceptance rate for the algorithm.

	**Value**
	
		Returns an ``MCMCSampler`` type object.
