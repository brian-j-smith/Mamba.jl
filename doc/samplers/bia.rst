.. index:: Sampling Functions; Binary Individual Adaptation

.. _section-BIA:

Binary Individual Adaptation (BIA)
----------------------------------

Implementation of the binary-state Individual Adaptation sampler of Griffin, et al. :cite:`griffin:2014:BIA` which adjusts a general proposal to the data. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: BIA(params::ElementOrVector{Symbol}; args...)

    Construct a ``Sampler`` object for BIA sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.
        * ``args...`` : additional keyword arguments to be passed to the ``BIAVariate`` constructor.

    **Value**

        Returns a ``Sampler{BIATune}`` type object.

    **Example**

        See the :ref:`Pollution <example-Pollution>` and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::BIAVariate)

    Draw one sample from a target distribution using the BIA sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bia:

    **Example**

        .. literalinclude:: bia.jl
            :language: julia


.. index:: Sampler Types; BIAVariate

BIAVariate Type
^^^^^^^^^^^^^^^

Declaration
```````````

``const BIAVariate = SamplerVariate{BIATune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::BIATune`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: BIAVariate(x::AbstractVector{T<:Real}, logf::Function; \\
                         A::Vector{Float64} = ones(x) / length(x), \\
                         D::Vector{Float64} = ones(x) / length(x), \\
                         epsilon::Real = 0.01 / length(x), decay::Real = 0.55, \\
                         target::Real = 0.45)

    Construct a ``BIAVariate`` object that stores simulated values and tuning parameters for BIA sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``A`` : vector of probabilities to switch the elements of ``x`` from 0 to 1 (i.e. added).
        * ``D`` : vector of probabilities to switch elements from 1 to 0 (i.e. deleted).
        * ``epsilon`` : range ``(epsilon, 1 - epsilon)`` for the elements of ``A`` and ``D``, where ``0 < epsilon < 0.5``.
        * ``decay`` : rate of decay of the adaptation, where ``0.5 < decay <= 1.0``.
        * ``target`` : target mutation rate, where ``0.0 < target < 1.0``.

    **Value**

        Returns a ``BIAVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; BIATune

BIATune Type
^^^^^^^^^^^^^

Declaration
```````````

``type BIATune <: SamplerTune``

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``A::Vector{Float64}`` : vector of probabilities to switch from 0 to 1.
* ``D::Vector{Float64}`` : vector of probabilities to switch from 1 to 0.
* ``epsilon::Float64`` : range ``(epsilon, 1 - epsilon)`` for the elements of ``A`` and ``D``.
* ``decay::Float64`` : rate of decay of the adaptation.
* ``target::Float64`` : target mutation rate.
* ``iter::Int`` : iteration number for adaptive updating.
