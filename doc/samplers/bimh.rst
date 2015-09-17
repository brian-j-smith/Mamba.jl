.. index:: Sampling Functions; Binary Independent Metropolis Hastings

Binary Independent Metropolis Hastings (BIMH)
---------------------------------------------

Implementation of the binary-state indepdnent Metropolis Hastings sampler of Schafer :cite:`schafer:2012:DIS,schafer:2013:SMCB` in which proposed updates are always state changes.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: bimh!(x::Vector{Int}, logf::Function)

    Simulate one draw from a target distribution using the BIMH sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``x`` : current state of parameters to be simulated.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``x``.

    **Value**

        Returns ``x`` updated.

    .. _example-bimh:

    **Example**

        .. literalinclude:: bimh.jl
            :language: julia

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: BIMH(params::Vector{Symbol})

    Construct a ``Sampler`` object for BIMH sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.

    **Value**

        Returns a ``Sampler`` type object.
