.. index:: Sampling Functions; Binary Deterministic

Binary Deterministic (BDS)
-----------------------

Implementation of the binary deterministic sampler of Schafer :cite:`schafer:2013:SMCB` (sometimes referred to as the 
metropolised Gibbs or Modiﬁed metropolised Gibbs) for simulating autocorrelated draws from a distribution that can 
be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: bds!(v::BDSVariate, Γ::Vector{Vector{Int}}, logf::Function)

    Simulate one draw from a target distribution using a binary deterministic sampler.  Parameters are assumed to be binary
    integer values (i.e. 0 or 1). 

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``Γ`` : Indices of parameters to flip to other state.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.
       
    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bds:

    **Example**

        .. literalinclude:: bds.jl
            :language: julia


.. index:: Sampler Types; BDSVariate

SliceVariate Type
^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^

Declaration
```````````

``type BDSTune``

Fields
``````
* ``Γ::Vector{Vector{Int}}`` : Indices of parameters to flip to other state.


Sampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: BDS(params::Vector{Symbol}, [d::Integer, k::Integer=1 | Γ::Vector{Vector{Int}}])

    Construct a ``Sampler`` object for binary deterministic sampling.  Parameters are assumed to be binary integers 
    (i.e. 0 or 1).

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.
        * ``d`` : Integer equal to the number of combined elements of nodes ``params``, used to determine the number of
        parameters to flip at the same time.
        * ``k`` : Number of parameters to flip at the same time. If not provided it will be sampled from a truncated 
        Geometric with mean d/2. 
        * ``Γ::Vector{Vector{Int}}`` : Indices of parameters to flip to other state. If not provided it will be all possible
        combinations of indices of size ``k``. 

    **Value**

        Returns a ``Sampler`` type object.

    **Example**

        See the :ref:`section-Examples` section.
