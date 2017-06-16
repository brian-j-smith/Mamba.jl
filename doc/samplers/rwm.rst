.. index:: Sampling Functions; Random Walk Metropolis

.. _section-RWM:

Random Walk Metropolis (RWM)
----------------------------

Random walk Metropolis-Hastings algorithm :cite:`hastings:1970:MCS,metropolis:1953:ESC` in which parameters are sampled from symmetric distributions centered around the current values.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: RWM(params::ElementOrVector{Symbol}, \
                  scale::ElementOrVector{T<:Real}; args...)

    Construct a ``Sampler`` object for RWM sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``scale`` : scaling value or vector of the same length as the combined elements of nodes ``params`` for the ``proposal`` distribution.  Values are relative to the unconstrained parameter space, where candidate draws are generated.
        * ``args...`` : additional keyword arguments to be passed to the ``RWMVariate`` constructor.

    **Value**

        Returns a ``Sampler{RWMTune}`` type object.

    **Example**

        See the :ref:`Dyes <example-Dyes>` and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::RWMVariate)

    Draw one sample from a target distribution using the RWM sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

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

``const RWMVariate = SamplerVariate{RWMTune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::RWMTune`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: RWMVariate(x::AbstractVector{T<:Real}, \
                         scale::ElementOrVector{U<:Real}, logf::Function; \
                         proposal::SymDistributionType=Normal)

    Construct a ``RWMVariate`` object that stores simulated values and tuning parameters for RWM sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``scale`` : scalar or vector of the same length as ``x`` for the ``proposal`` distribution.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``proposal`` : symmetric distribution of type ``Biweight``, ``Cosine``, ``Epanechnikov``, ``Normal``, ``SymTriangularDist``, ``SymUniform``, or ``Triweight`` to be centered around current parameter values and used to generate proposal draws.  Specified ``scale`` determines the standard deviations of Normal proposals and widths of the others.

    **Value**

        Returns a ``RWMVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.


.. index:: Sampler Types; RWMTune

RWMTune Type
^^^^^^^^^^^^

Declaration
```````````

``type RWMTune <: SamplerTune``

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``scale::Union{Float64, Vector{Float64}}`` : scaling for the proposal distribution.
* ``proposal::SymDistributionType`` : proposal distribution.
