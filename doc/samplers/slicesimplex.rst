.. index:: Sampling Functions; Slice Simplex

.. _section-SliceSimplex:

Slice Simplex (SliceSimplex)
----------------------------

Implementation of the slice simplex sampler as described by Cowles et al. :cite:`cowles:2009:RAMPS` for simulating autocorrelated draws of parameters on the simplex :math:`\{\theta_1, \ldots, \theta_d : \theta_i \ge 0, \sum_{i=1}^d \theta_i = 1\}` and from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: slicesimplex!(v::SliceSimplexVariate, logf::Function; scale::Real=1.0)

    Simulate one draw from a target distribution using a slice simplex sampler.  Parameters are assumed to be continuous and constrained to a simplex.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``scale`` : a value ``0 < scale <= 1`` by which to scale the standard simplex to define an initial space from which to simulate values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).

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

.. function:: SliceSimplexVariate(x::AbstractVector{T<:Real}, tune::SliceSimplexTune)
              SliceSimplexVariate(x::AbstractVector{T<:Real}, tune=nothing)

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

* ``scale`` : a value ``0 < scale <= 1`` by which to scale the standard simplex to define an initial space from which to simulate values.


Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: SliceSimplex(params::Vector{Symbol}; scale::Real=1.0)

    Construct a ``Sampler`` object for which slice simplex sampling is to be applied separately to each of the supplied parameters.  Parameters are assumed to be continuous and constrained to a simplex.

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.
        * ``scale`` : a value ``0 < scale <= 1`` by which to scale the standard simplex to define an initial space from which to simulate values.

    **Value**

        Returns a ``Sampler`` type object.

    **Example**

        See the :ref:`Asthma <example-Asthma>`, :ref:`Eyes <example-Eyes>`, and other :ref:`section-Examples`.
