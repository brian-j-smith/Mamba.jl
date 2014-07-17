.. index:: Sampling Functions; Shrinkage Slice

Shrinkage Slice (Slice)
-----------------------

Implementation of the shrinkage slice sampler of Neal :cite:`neal:2003:SS` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: slice!(v::SliceVariate, width::Vector{Float64}, logf::Function, stype::Symbol=:multivar)

	Simulate one draw from a target distribution using a shrinkage slice sampler.  Parameters are assumed to be continuous, but may be constrained or unconstrained.
	
	**Arguments**
	
		* ``v`` : current state of parameters to be simulated.
		* ``width`` : vector of the same length as ``v``, defining initial widths of a hyperrectangle from which to simulate values.
		* ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.
		* ``stype`` : sampler type. Options are
			* ``:multivar`` : Joint multivariate sampling of parameters.
			* ``:univar`` : Sequential univariate sampling.
		
	**Value**
	
		Returns ``v`` updated with simulated values and associated tuning parameters.
	
	.. _example-slice:
	
	**Example**

		The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.
		
		.. literalinclude:: slice.jl
			:language: julia


.. index:: SliceVariate

SliceVariate Type
^^^^^^^^^^^^^^^^^

Declaration
```````````

``SliceVariate <: VectorVariate``

Fields
``````

* ``value::Vector{VariateType}`` : vector of sampled values.
* ``tune::SliceTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: SliceVariate(x::Vector{VariateType}, tune::SliceTune)
              SliceVariate(x::Vector{VariateType}, tune=nothing)

  	Construct a ``SliceVariate`` object that stores sampled values and tuning parameters for slice sampling.
	
	**Arguments**
	
		* ``x`` : vector of sampled values.
		* ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.
		
	**Value**
	
		Returns a ``SliceVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: SliceTune

SliceTune Type
^^^^^^^^^^^^^^

Declaration
```````````

``type SliceTune``

Fields
``````
* ``width::Vector{Float64}`` : vector of initial widths defining hyperrectangles from which to simulate values.


Sampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: Slice(params::Vector{Symbol}, width::Vector{T<:Real}, \
                stype::Symbol=:multivar; transform::Bool=false)

	Construct an ``Sampler`` object for shrinkage slice sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.
	
	**Arguments**
	
		*  ``params`` : stochastic nodes to be updated with the sampler.
		* ``width`` : vector of the same length as the combined elements of nodes ``params``, defining initial widths of a hyperrectangle from which to simulate values.
		* ``stype`` : sampler type. Options are
			* ``:multivar`` : Joint multivariate sampling of parameters.
			* ``:univar`` : Sequential univariate sampling.
		* ``transform`` : whether to sample parameters on the link-transformed scale (unconstrained parameter space).  If ``true``, then constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``link()`` function, and ``width`` is interpreted as being relative to the unconstrained parameter space.  Otherwise, sampling is relative to the untransformed space.

	**Value**
	
		Returns an ``Sampler`` type object.

	**Example**
	
		See the :ref:`section-Examples` section.
