.. index:: Sampling Functions; Discrete Gibbs Sampler

.. _section-DGS:

Discrete Gibbs Sampler (DGS)
----------------------------

Implementation of a sampler for the simulation of discrete or discretized model parameters with finite support.  Draws are simulated directly from a probability mass function that can be specified up to a constant of proportionality.  Note that versions of this sampler evaluate the probability function over all points in the parameter space; and, as a result, may be very computationally intensive for large spaces.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: DGS(params::ElementOrVector{Symbol})

    Construct a ``Sampler`` object for which DGS sampling is to be applied separately to each of the supplied parameters.  Parameters are assumed to have discrete univariate distributions with finite supports.

    **Arguments**

        *  ``params`` : stochastic node(s) to be updated with the sampler.

    **Value**

        Returns a ``Sampler{DSTune{Function}}`` type object.

    **Example**

        See the :ref:`Eyes <example-Eyes>`, :ref:`Pollution <example-Pollution>`, and other :ref:`section-Examples`.

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::DGSVariate)
              sample!(v::DiscreteVariate)

    Draw one sample directly from a target probability mass function.  Parameters are assumed to have discrete and finite support.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.


.. index:: Sampler Types; DGSVariate
.. index:: Sampler Types; DiscreteVariate

Discrete Variate Types
^^^^^^^^^^^^^^^^^^^^^^

Declaration
```````````

.. code-block:: julia

    typealias DGSVariate SamplerVariate{DSTune{Function}}
    typealias DiscreteVariate SamplerVariate{DSTune{Vector{Float64}}}

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::DSTune{F<:DSForm}`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: DGSVariate(x::AbstractVector{T<:Real}, support::Matrix{U<:Real}, \
                         mass::Function)
              DiscreteVariate(x::AbstractVector{T<:Real}, support::Matrix{U<:Real}, \
                              mass::Vector{Float64})

    Construct an object that stores simulated values and tuning parameters for discrete sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``support`` : matrix whose columns contain the vector coordinates in the parameter space from which to simulate values.
        * ``mass`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the density (up to a normalizing constant), or a vector of sampling probabilities for the parameter space.

    **Value**

        Returns a ``DGSVariate`` or ``DiscreteVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; DSForm
.. index:: Sampler Types; DSTune

DSTune Type
^^^^^^^^^^^

Declaration
```````````

.. code-block:: julia

    typealias DSForm Union{Function, Vector{Float64}}
    type DSTune{F<:DSForm} <: SamplerTune

Fields
``````

* ``mass::Nullable{F}`` : density mass function or vector supplied to the constructor, or null if not supplied.
* ``support::Matrix{Real}`` : matrix whose columns contain the vector coordinates in the parameter space from which to simulate values.
