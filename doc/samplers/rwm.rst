.. index:: Sampling Functions; Random Walk Metropolis

.. _section-RWM:

Random Walk Metropolis (RWM)
----------------------------

Random walk Metropolis-Hastings algorithm :cite:`hastings:1970:MCS,metropolis:1953:ESC` in which parameters are sampled from symmetric distributions centered around the current values.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: RWM(params::ElementOrVector{Symbol}, \
                  scale::ElementOrVector{T<:Real}; \
                  proposal::SymDistributionType=Normal)

    Construct a ``Sampler`` object for RWM sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``scale`` : scaling value or vector of the same length as the combined elements of nodes ``params`` for the ``proposal`` distribution.  Values are relative to the unconstrained parameter space, where candidate draws are generated.
        * ``proposal`` : symmetric distribution of type ``Biweight``, ``Cosine``, ``Epanechnikov``, ``Normal``, ``SymTriangularDist``, ``SymUniform``, or ``Triweight`` to be centered around current parameter values and used to generate proposal draws.  Specified ``scale`` determines the standard deviations of Normal proposals and widths of the others.

    **Value**

        Returns a ``Sampler{RWMTune}`` type object.

    **Example**

        See the :ref:`Dyes <example-Dyes>` and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: rwm!(v::RWMVariate, scale::ElementOrVector{T<:Real}, \
                   logf::Function; proposal::SymDistributionType=Normal)

    Simulate one draw from a target distribution using an RWM sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``scale`` : scalar or vector of the same length as ``v`` for the ``proposal`` distribution.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``proposal`` : symmetric distribution of type ``Biweight``, ``Cosine``, ``Epanechnikov``, ``Normal``, ``SymTriangularDist``, ``SymUniform``, or ``Triweight`` to be centered around current parameter values and used to generate proposal draws.  Specified ``scale`` determines the standard deviations of Normal proposals and widths of the others.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-rwm:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.

        .. literalinclude:: rwm.jl
            :language: julia

.. index:: Sampler Types; RWMVariate

RWMVariate Type
^^^^^^^^^^^^^^^

Declaration
```````````

``typealias RWMVariate SamplerVariate{RWMTune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::RWMTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: RWMVariate(x::AbstractVector{T<:Real})
              RWMVariate(x::AbstractVector{T<:Real}, tune::RWMTune)

    Construct a ``RWMVariate`` object that stores simulated values and tuning parameters for RWM sampling.

    **Arguments**

        * ``x`` : simulated values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If not supplied, parameters are set to their defaults.

    **Value**

        Returns a ``RWMVariate`` type object with fields set to the values supplied to arguments ``x`` and ``tune``.


.. index:: Sampler Types; RWMTune

RWMTune Type
^^^^^^^^^^^^

Declaration
```````````

``type RWMTune <: SamplerTune``

Fields
``````

* ``scale::Union{Real, Vector}`` : scaling for the proposal distribution.
* ``proposal::SymDistributionType`` : proposal distribution.
