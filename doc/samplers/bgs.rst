.. index:: Sampling Functions; Binary Gibbs Sampler

.. _section-BGS:

Binary Gibbs Sampler (BGS)
-------------------------------

Implementation of the binary-state Gibbs sampler described by Dellaportas et al. :cite:`dellaportas:2002:GVS` in which components are drawn sequentially from full conditional marginal
distributions.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: bgs!(v::BGSVariate, logf::Function)

    Simulate one draw from a target distribution using the BGS sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bgs:

    **Example**

        .. literalinclude:: bgs.jl
            :language: julia


.. index:: Sampler Types; BGSVariate

BGSVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``BGSVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::BGSTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: BGSVariate(x::Vector{Float64}, tune::BGSTune)
              BGSVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``BGSVariate`` object that stores sampled values and tuning parameters for BGS sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``BGSVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

BGSTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type BGSTune``

Fields
``````

* ``probs::Vector{Float64}`` : probabilities of the full conditional marginal distributions from which components of the candidate vector were sequentially drawn.


Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: BGS(params::Vector{Symbol})

    Construct a ``Sampler`` object for BGS sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.

    **Value**

        Returns a ``Sampler`` type object.

    **Example**

        See the :ref:`Pollution <example-Pollution>` and other :ref:`section-Examples`.
