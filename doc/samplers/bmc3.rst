.. index:: Sampling Functions; Binary MCMC Model Composition

.. _section-BMC3:

Binary MCMC Model Composition (BMC3)
------------------------------------

Implementation of the binary-state MCMC Model Composition of Madigan and York :cite:`madigan:1995:MC3` in which proposed updates are always state changes. Liu :cite:`liu:1996:MMG` shows this sampler is more efficient than Gibbs sampling for a binary vector. Schafer :cite:`schafer:2012:DIS,schafer:2013:SMCB` proposes a method for block updates of binary vectors using this sampler. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: BMC3(params::ElementOrVector{Symbol}; k::BMC3Form=1)

    Construct a ``Sampler`` object for BMC3 sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.
        * ``k`` : number of parameters or vector of parameter indices to select at random for simultaneous updating in each call of the sampler.

    **Value**

        Returns a ``Sampler{BMC3Tune{typeof(k)}}`` type object.

    **Example**

        See the :ref:`Pollution <example-Pollution>` and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::SamplerVariate{BMC3Tune{F<:BMC3Form}})

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

BMC3 Variate Type
^^^^^^^^^^^^^^^^^

Declaration
```````````

``SamplerVariate{BMC3Tune{F<:BMC3Form}}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::BMC3Tune{F}`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: BMC3Variate(x::AbstractVector{T<:Real}, logf::Function; k::BMC3Form=1)

    Construct a ``SamplerVariate`` object that stores simulated values and tuning parameters for BMC3 sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``k`` : number of parameters or vector of parameter indices to select at random for simultaneous updating in each call of the sampler.

    **Value**

        Returns a ``SamplerVariate{BMC3Tune{typeof(k)}}`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; BMC3Form
.. index:: Sampler Types; BMC3Tune

BMC3Tune Type
^^^^^^^^^^^^^

Declaration
```````````

.. code-block:: julia

    const BMC3Form = Union{Int, Vector{Vector{Int}}}
    type BMC3Tune{F<:BMC3Form} <: SamplerTune

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``k::F`` : number of parameters or vector of parameter indices to select at random for simultaneous updating in each call of the sampler.
