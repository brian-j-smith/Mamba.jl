.. index:: Sampling Functions; Slice Simplex Sampler

Slice Simplex
-------------

Implementation of the slice simplex sampler as described by Cowles et. al :cite:`cowles:2009:RAMPS` 
for simulating autocorrelated draws on a simplex from a distribution that can be specified up to a 
constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function slicesimplex!(v::SliceSimplexVariate, width::Float64, logf::Function)

    Simulate one draw from a target distribution using a slice simplex sampler.  Parameters are assumed to be continuous
    and constrained to a simplex (i.e. sum to 1). 

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``width`` : The next value is sampled from a shrunken simplex around the current value. The edge lengths of the 
        shrunken simplex are equal to ``width`` times the edge length of the standard (n-1)-Simplex. 0 <= width <= 1. 
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-slicesimplex:

    **Example**

        .. literalinclude:: slicesimplex.jl
            :language: julia


.. index:: Sampler Types; SliceSimplexVariate

SliceSimplexVariate Type
^^^^^^^^^^^^^^^^^^^^^^^^

Declaration
```````````

``SliceSimplexVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::SliceSimplexTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: SliceSimplexVariate(x::Vector{Float64}, tune::SliceSimplexTune)
              SliceSimplexVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``SliceSimplexVariate`` object that stores sampled values and tuning parameters for slice simplex sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``SliceSimplexVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; SliceSimplexTune

SliceSimplexTune Type
^^^^^^^^^^^^^^^^^^^^^

Declaration
```````````

``type SliceSimplexTune``

Fields
``````
* ``width`` : The next value is sampled from a shrunken simplex around the current value. The edge lengths of the 
shrunken simplex are equal to ``width`` times the edge length of the standard (n-1)-Simplex. 0 <= width <= 1. 


Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: SliceSimplex{T<:Real}(params::Vector{Symbol}, width::T; transform::Bool=false)

    Construct a ``Sampler`` object for slice simplex sampling. Parameters are assumed to be continuous
    and constrained to a simplex (i.e. sum to 1). 

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.
        * ``width`` : The next value is sampled from a shrunken simplex around the current value. The edge lengths of the 
        shrunken simplex are equal to ``width`` times the edge length of the standard (n-1)-Simplex. 0 <= width <= 1.
        * ``transform`` : whether to sample parameters on the link-transformed scale (unconstrained parameter space).  If 
        ``true``, then constrained parameters are mapped to unconstrained space according to transformations defined by the 
        :ref:`section-Stochastic` ``link()`` function, and ``width`` is interpreted as being relative to the unconstrained 
        parameter space.  Otherwise, sampling is relative to the untransformed space.


    **Value**

        Returns a ``Sampler`` type object.