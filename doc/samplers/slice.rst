.. index:: Sampling Functions; Multivariate Slice

Multivariate Slice (Slice)
--------------------------

Implementation of the multivariate shrinkage slice sampler of Neal :cite:`neal:2003:SS` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: slice!(v::VariateSlice, width::Vector{Float64}, logf::Function)

	Simulate one draw from a target distribution using a multivariate (shrinkage) slice sampler.
	
	**Arguments**
	
		* ``v`` : current state of parameters to be simulated.
		* ``width`` : vector of the same length as ``v``, defining initial widths of a hyperrectangle from which to simulate values.
		* ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.
		
	**Value**
	
		Returns ``v`` updated with the simulated values and associated tuning parameters.
	
	**Example**

		.. literalinclude:: slice.jl
			:language: julia


.. index:: VariateSlice

VariateSlice Type
^^^^^^^^^^^^^^^^^

Declaration
```````````

``VariateSlice <: VariateVector``

Fields
``````

* ``value::Vector{VariateType}`` : vector of sampled values.
* ``tune::TuneSlice`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: VariateSlice(x::Vector{VariateType}, tune::TuneSlice)
              VariateSlice(x::Vector{VariateType}, tune=nothing)

  	Construct a ``VariateSlice`` object that stores values and tuning parameters for slice sampling.
	
	**Arguments**
	
		* ``x`` : vector of sampled values.
		* ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.
		
	**Value**
	
		Returns a ``VariateSlice`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: TuneSlice

TuneSlice Type
^^^^^^^^^^^^^^

Declaration
```````````

``type TuneSlice``

Fields
``````
* ``width::Vector{Float64}`` : vector of initial widths defining hyperrectangles from which to simulate values.


MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: Slice(params::Vector{T<:String}, width::Vector{Float64})

	Construct an ``MCMCSampler`` object for multivariate (shrinkage) slice sampling.
	
	**Arguments**
	
		*  ``params`` : named stochastic nodes to be updated with the sampler.
		* ``width`` : vector of the same length as the combined elements of nodes ``params``, defining initial widths of a hyperrectangle from which to simulate values.

	**Value**
	
		Returns an ``MCMCSampler`` type object.

.. index:: Sampling Functions; Slice within Gibbs

Slice within Gibbs (SliceWG)
----------------------------

Implementation of the univariate shrinkage slice sampler of Neal :cite:`neal:2003:SS` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: slicewg!(v::VariateSlice, width::Vector{Float64}, logf::Function)

	Simulate one draw from a target distribution using a univariate (shrinkage) slice-within-Gibbs sampler.
	
	**Arguments**
	
		* ``v`` : current state of parameters to be simulated.
		* ``width`` : vector of the same length as ``v``, defining initial intervals widths from which to simulate values.
		* ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.
		
	**Value**
	
		Returns ``v`` updated with the simulated values and associated tuning parameters.
	
	**Example**

		.. literalinclude:: slice.jl
			:language: julia
	
MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SliceWG(params::Vector{T<:String}, width::Vector{Float64})

	Construct an ``MCMCSampler`` object for univariate (shrinkage) slice-within-Gibbs sampling.
	
	**Arguments**
	
		*  ``params`` : named stochastic nodes to be updated with the sampler.
		* ``width`` : vector of the same length as the combined elements of nodes ``params``, defining initial interval widths from which to simulate values.

	**Value**
	
		Returns an ``MCMCSampler`` type object.
