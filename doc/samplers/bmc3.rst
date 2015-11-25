.. index:: Sampling Functions; Binary MCMC Model Composition

.. _section-BMC3:

Binary MCMC Model Composition (BMC3)
------------------------------------

Implementation of the binary-state MCMC Model Composition of Madigan and York :cite:`madigan:1995:MC3` in which proposed updates are always state changes. Liu :cite:`liu:1996:MMG` shows this sampler is more efficient than Gibbs sampling for a binary vector. Schafer :cite:`schafer:2012:DIS,schafer:2013:SMCB` proposes a method for block updates of binary vectors using this sampler. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.



Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: bmc3!(v::BMMGVariate, logf::Function; k::Integer=1)
              bmc3!(v::BMMGVariate, indexset::Vector{Vector{Int}}, logf::Function)

    Simulate one draw from a target distribution using the BMC3 sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``k`` : number of parameters, such that ``k <= length(v)``, to select at random for simultaneous updating in each call to the sampler.
        * ``indexset`` : candidate set of indices of the parameters whose states are to be changed simultaneously.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bmc3:

    **Example**

        .. literalinclude:: bmc3.jl
            :language: julia


.. index:: Sampler Types; BMC3Variate

BMC3Variate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``BMC3Variate <: SamplerVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::BMC3Tune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: BMC3Variate(x::AbstractVector{T<:Real}, tune::BMC3Tune)
              BMC3Variate(x::AbstractVector{T<:Real}, tune=nothing)

    Construct a ``BMC3Variate`` object that stores sampled values and tuning parameters for BMC3 sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``BMC3Variate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; BMC3Tune

BMC3Tune Type
^^^^^^^^^^^^^

Declaration
```````````

``type BMC3Tune <: SamplerTune``

Fields
``````

* ``indexset::Vector{Vector{Int}}`` : candidate set of indices of the parameters whose states are to be changed simultaneously.


Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: BMC3(params::Vector{Symbol}; k::Integer=1)
              BMC3(params::Vector{Symbol}, indexset::Vector{Vector{Int}})

    Construct a ``Sampler`` object for BMC3 sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.
        * ``k`` : number of parameters to select at random for updating in each call to the sampler.
        * ``indexset`` : candidate set of indices of the parameters to change.

    **Value**

        Returns a ``Sampler`` type object.
