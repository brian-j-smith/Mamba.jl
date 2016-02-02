.. index:: Sampling Functions; Slice Simplex

.. _section-SliceSimplex:

Slice Simplex (SliceSimplex)
----------------------------

Implementation of the slice simplex sampler as described by Cowles et al. :cite:`cowles:2009:RAMPS` for simulating autocorrelated draws of parameters on the simplex :math:`\{\theta_1, \ldots, \theta_d : \theta_i \ge 0, \sum_{i=1}^d \theta_i = 1\}` and from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SliceSimplex(params::ElementOrVector{Symbol}; args...)

    Construct a ``Sampler`` object for which SliceSimplex sampling is to be applied separately to each of the supplied parameters.  Parameters are assumed to be continuous and constrained to a simplex.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.
        * ``args...`` : additional keyword arguments to be passed to the ``SliceSimplexVariate`` constructor.

    **Value**

        Returns a ``Sampler{SliceSimplexTune}`` type object.

    **Example**

        See the :ref:`Asthma <example-Asthma>`, :ref:`Eyes <example-Eyes>`, and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::SliceSimplexVariate)

    Draw one sample from a target distribution using the SliceSimplex sampler.  Parameters are assumed to be continuous and constrained to a simplex.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

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

``typealias SliceSimplexVariate SamplerVariate{SliceSimplexTune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::SliceSimplexTune`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: SliceSimplexVariate(x::AbstractVector{T<:Real}, logf::Function; \
                                  scale::Real=1.0)

    Construct a ``SliceSimplexVariate`` object that stores simulated values and tuning parameters for slice simplex sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``scale`` : value ``0 < scale <= 1`` by which to scale the standard simplex to define an initial space from which to simulate values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).

    **Value**

        Returns a ``SliceSimplexVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; SliceSimplexTune

SliceSimplexTune Type
^^^^^^^^^^^^^^^^^^^^^

Declaration
```````````

``type SliceSimplexTune <: SamplerTune``

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``scale::Float64`` : value ``0 < scale <= 1`` by which to scale the standard simplex to define an initial space from which to simulate values.
