.. index:: Sampling Functions; Binary Metropolized Gibbs

Binary Metropolized Gibbs (BMG)
---------------------------------------------

Implementation of the binary-state Metropolized Gibbs sampler of Schafer :cite:`schafer:2012:DIS,schafer:2013:SMCB` in which proposed updates are always state changes.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: bmg!(v::BMGVariate, logf::Function)

    Simulate one draw from a target distribution using the BMG sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``x``.

    **Value**

        Returns ``x`` updated.

    .. _example-bmg:

    **Example**

        .. literalinclude:: bmg.jl
            :language: julia

.. index:: Sampler Types; BMGVariate

BMGVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``BMGVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.

Constructors
````````````

.. function:: BMGVariate(x::Vector{Float64})

    Construct a ``BMGVariate`` object that stores sampled values for BMG sampling.

    **Arguments**

        * ``x`` : vector of sampled values.

    **Value**

        Returns a ``BMGVariate`` type object with fields pointing to the values supplied to arguments ``x``.

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: BMG(params::Vector{Symbol})

    Construct a ``Sampler`` object for BMG sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.

    **Value**

        Returns a ``Sampler`` type object.
