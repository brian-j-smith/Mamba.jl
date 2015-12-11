.. index:: Sampling Functions; Random Walk Metropolis

.. _section-RWM:

Random Walk Metropolis (RWM)
----------------------------

Simple random walk Metropolis-Hastings algorithm :cite:`hastings:1970:MCS,metropolis:1953:ESC` in which parameters are perturbed from their previous iteration by a user specified amount. Simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: RWM(params::Vector{Symbol}, delta::Real)

    Construct a ``Sampler`` object for random walk Metropolis sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``delta`` : Parameters are perturbed from previous iteration by Uniform(-delta, delta).

    **Value**

        Returns a ``Sampler{RWMTune}`` type object.

    **Example**

        See the :ref:`Dyes <example-Dyes>` and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: rwm!(v::RWMVariate, delta::Real, logf::Function)

    Simulate one draw from a target distribution using an random walk Metropolis sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.  
        * ``delta`` : Parameters are perturbed from previous iteration by Uniform(-delta, delta).
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).

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

    Construct a ``RWMVariate`` object that stores simulated values and tuning parameters for random walk Metropolis sampling.

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

* ``delta`` : Parameters are perturbed from previous iteration by Uniform(-delta, delta).
