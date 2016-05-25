.. index:: Sampling Functions; Individual Adaptation Sampler

.. _section-BIA:

Binary MCMC Model Composition (BIA)
------------------------------------

Implementation of the binary-state Individual adaptation sampler of Griffin, et al. :cite:`griffin:2014:BIA` which adjusts a general proposal to the data. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: BIA(params::ElementOrVector{Symbol}; args...)

    Construct a ``Sampler`` object for BIA sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.
        * ``args...`` : additional keyword arguments to be passed to the ``BIAVariate`` constructor.

    **Value**

        Returns a ``Sampler{BIATune}`` type object.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::BIAVariate)

    Draw one sample from a target distribution using the BIA sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bmc3:

    **Example**

        .. literalinclude:: bia.jl
            :language: julia


.. index:: Sampler Types; BIAVariate

BIAVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``typealias BIAVariate SamplerVariate{BIATune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::BIATune`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: BIAVariate(x::AbstractVector{T<:Real}, logf::Function; k::Integer=1)

    Construct a ``BIAVariate`` object that stores simulated values and tuning parameters for BIA sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``A`` : Vector of probabilities to switch from 0 to 1 (i.e. added).
        * ``D`` : Vector of probabilities to switch from 1 to 0 (i.e. deleted).
        * ``epsilon`` : minimum that ``A`` and ``D`` can take on. Must be 0 < epsilon < 0.5.
        * ``lambda`` : Rate of decay of the adaptation. Must be 0.5 < lambda <= 1.0.
        * ``tau`` : Target mutation rate. Must be 0.0 < tau < 1.0.

    **Value**

        Returns a ``BIAVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; BIATune

BIATune Type
^^^^^^^^^^^^^

Declaration
```````````

``type BIATune <: SamplerTune``

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``A`` : Vector of probabilities to switch from 0 to 1 (i.e. added).
* ``D`` : Vector of probabilities to switch from 1 to 0 (i.e. deleted).
* ``epsilon`` : minimum that ``A`` and ``D`` can take on. Must be 0 < epsilon < 0.5.
* ``lambda`` : Rate of decay of the adaptation. Must be 0.5 < lambda <= 1.0.
* ``tau`` : Target mutation rate. Must be 0.0 < tau < 1.0.
