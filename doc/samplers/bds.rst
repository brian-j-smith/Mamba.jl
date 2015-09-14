.. index:: Sampling Functions; Binary Deterministic Sampler

Binary Deterministic Sampler (BDS)
----------------------------------

Implementation of the binary deterministic sampler of Schafer :cite:`schafer:2012:DIS,schafer:2013:SMCB` (sometimes referred to as the *Metropolized Gibbs* or *Modified Metropolized Gibbs*) for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: bds!(v::BDSVariate, Γ::Vector{Vector{Int}}, logf::Function)

    Simulate one draw from a target distribution using a binary deterministic sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``Γ`` : candidate set of indices of the parameters whose states are to be flipped simultaneously.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bds:

    **Example**

        .. literalinclude:: bds.jl
            :language: julia


.. index:: Sampler Types; BDSVariate

BDSVariate Type
^^^^^^^^^^^^^^^

Declaration
```````````

``BDSVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::BDSTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: BDSVariate(x::Vector{Float64}, tune::BDSTune)
              BDSVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``BDSVariate`` object that stores sampled values and tuning parameters for binary deterministic sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``BDSVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; BDSTune

BDSTune Type
^^^^^^^^^^^^

Declaration
```````````

``type BDSTune``

Fields
``````
* ``Γ::Vector{Vector{Int}}`` : candidate set of indices of the parameters whose states are to be flipped simultaneously.


Sampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: BDS(params::Vector{Symbol}, d::Integer, k::Integer=1)
              BDS(params::Vector{Symbol}, Γ::Vector{Vector{Int}})

    Construct a ``Sampler`` object for binary deterministic sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.
        * ``d`` : total length of the parameters in the combined nodes.
        * ``k`` : generate all combinations of ``k <= d`` candidate indices of the parameters to flip.
        * ``Γ`` : candidate set of indices of the parameters to flip.

    **Value**

        Returns a ``Sampler`` type object.
