.. index:: Sampling Functions; Binary MCMC Model Composition

.. _section-BMC3:

Binary MCMC Model Composition (BMC3)
------------------------------------

Implementation of the binary-state MCMC Model Composition of Madigan and York :cite:`madigan:1995:MC3` in which proposed updates are always state changes. Liu :cite:`liu:1996:MMG` shows this sampler is more efficient than Gibbs sampling for a binary vector. Schafer :cite:`schafer:2012:DIS,schafer:2013:SMCB` proposes a method for block updates of binary vectors using this sampler. An adaptive version of this sampler is proposed in `lamnisos:2013:AMC3`. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: BMC3(params::ElementOrVector{Symbol}; args...)

    Construct a ``Sampler`` object for BMC3 sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.
        * ``args...`` : additional keyword arguments to be passed to the ``BMC3Variate`` constructor.

    **Value**

        Returns a ``Sampler{BMC3Tune}`` type object.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::BMC3Variate)

    Draw one sample from a target distribution using the BMC3 sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

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

``typealias BMC3Variate SamplerVariate{BMC3Tune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::BMC3Tune`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: BMC3Variate(x::AbstractVector{T<:Real}, logf::Function; k::Integer=1, indexset::Vector{Vector{Int}}=Vector{Vector{Int}}(), epsilon::Float64=1/length(x))

    Construct a ``BMC3Variate`` object that stores simulated values and tuning parameters for BMC3 sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``k`` : number of parameters to select at random for simultaneous updating in each call of the sampler. Used if ``indexset`` is not specified. 
        * ``indexset`` : A set of the sets of parameters to select at random for simultaneous updating in ecah call of the sampler. 
        * ``epsilon`` : If sampling is specified to be adaptive, a parameter is selected to be updated by a mixture of a discrete distribution depending on the current MCMC sample and a uniform distribution, where ``epsilon`` is the mixture weight. 

    **Value**

        Returns a ``BMC3Variate`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; BMC3Tune

BMC3Tune Type
^^^^^^^^^^^^^

Declaration
```````````

``type BMC3Tune <: SamplerTune``

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``k::Int`` : number of parameters to select at random for simultaneous updating in each call of the sampler. Used if ``indexset`` is not specified. 
* ``indexset::Vector{Vector{Int}}`` : A set of the set of parameters to select at random for simultaneous updating in each call of the sampler. 
* ``n::Int`` : current iteration.
* ``m::Vector{Float64}`` : current inclusion frequency, i.e. mean. 
* ``v::Vector{Float64}`` : current sample variance.
* ``epsilon::Float64`` : If sampling is specified to be adaptive, a parameter is selected to be updated by a mixture of a discrete distribution depending on the current MCMC sample and a uniform distribution, where ``epsilon`` is the mixture weight.
