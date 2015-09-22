.. index:: Sampling Functions; Binary Modified Metropolised Gibbs

Binary Modified Metropolised Gibbs (BMMG)
-----------------------------------------

Implementation of the binary-state modified Metropolised Gibbs sampler of Schafer :cite:`schafer:2012:DIS,schafer:2013:SMCB` in which proposed updates are always state changes.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: bmmg!(v::BMMGVariate, indexset::Vector{Vector{Int}}, logf::Function)

    Simulate one draw from a target distribution using the BMMG sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``indexset`` : candidate set of indices of the parameters whose states are to be changed simultaneously.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bmmg:

    **Example**

        .. literalinclude:: bmmg.jl
            :language: julia


.. index:: Sampler Types; BMMGVariate

BMMGVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``BMMGVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::BMMGTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: BMMGVariate(x::Vector{Float64}, tune::BMMGTune)
              BMMGVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``BMMGVariate`` object that stores sampled values and tuning parameters for BMMG sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``BMMGVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; BMMGTune

BMMGTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type BMMGTune``

Fields
``````
* ``indexset::Vector{Vector{Int}}`` : candidate set of indices of the parameters whose states are to be changed simultaneously.


Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: BMMG(params::Vector{Symbol}, d::Integer, k::Integer=1)
              BMMG(params::Vector{Symbol}, indexset::Vector{Vector{Int}})

    Construct a ``Sampler`` object for BMMG sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.
        * ``d`` : total length of the parameters in the combined nodes.
        * ``k`` : generate all combinations of ``k <= d`` candidate indices of the parameters to change.
        * ``indexset`` : candidate set of indices of the parameters to change.

    **Value**

        Returns a ``Sampler`` type object.
