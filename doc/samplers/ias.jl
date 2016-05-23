.. index:: Sampling Functions; Individual Adaptation Sampler

.. _section-IAS:

Binary MCMC Model Composition (IAS)
------------------------------------

Implementation of the binary-state Individual adaptation sampler of Griffin, et al. :cite:`griffin:2014:IAS` which adjusts a general proposal to the data. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: IAS(params::ElementOrVector{Symbol}; args...)

    Construct a ``Sampler`` object for IAS sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.
        * ``args...`` : additional keyword arguments to be passed to the ``IASVariate`` constructor.

    **Value**

        Returns a ``Sampler{IASTune}`` type object.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::IASVariate)

    Draw one sample from a target distribution using the IAS sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bmc3:

    **Example**

        .. literalinclude:: ias.jl
            :language: julia


.. index:: Sampler Types; IASVariate

IASVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``typealias IASVariate SamplerVariate{IASTune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::IASTune`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: IASVariate(x::AbstractVector{T<:Real}, logf::Function; k::Integer=1)

    Construct a ``IASVariate`` object that stores simulated values and tuning parameters for IAS sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``A`` : Vector of probabilities to switch from 0 to 1 (i.e. added).
        * ``D`` : Vector of probabilities to switch from 1 to 0 (i.e. deleted).
        * ``epsilon`` : minimum that ``A`` and ``D`` can take on. Must be 0 < epsilon < 0.5.
        * ``lambda`` : Rate of decay of the adaptation. Must be 0.5 < lambda <= 1.0.
        * ``tau`` : Target mutation rate. Must be 0.0 < tau < 1.0.

    **Value**

        Returns a ``IASVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; IASTune

IASTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type IASTune <: SamplerTune``

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``A`` : Vector of probabilities to switch from 0 to 1 (i.e. added).
* ``D`` : Vector of probabilities to switch from 1 to 0 (i.e. deleted).
* ``epsilon`` : minimum that ``A`` and ``D`` can take on. Must be 0 < epsilon < 0.5.
* ``lambda`` : Rate of decay of the adaptation. Must be 0.5 < lambda <= 1.0.
* ``tau`` : Target mutation rate. Must be 0.0 < tau < 1.0.
