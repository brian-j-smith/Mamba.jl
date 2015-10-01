.. index:: Sampling Functions; Discrete Gibbs Sampler

.. _section-DGS:

Discrete Gibbs Sampler (DGS)
----------------------------

Implementation of a sampler for the simulation of discrete or discretized model parameters with finite support.  Draws are simulated directly from a probability mass function that can be specified up to a constant of proportionality.  Note that versions of this sampler evaluate the probability function over all points in the parameter space; and, as a result, may be very computationally intensive for large spaces.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: dgs!(v::DGSVariate, support::Vector, logf::Function)
              dgs!(v::DGSVariate, support::Vector, probs::Vector{Float64})

    Simulate one draw directly from a target probability mass function.  Parameters are assumed to have discrete and finite support.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``support`` : vector of coordinates in the parameter space from which to simulate values.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at each element in the ``support`` vector.
        * ``probs`` : sampling probabilities for the ``support`` vector of coordinates.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.


.. index:: Sampler Types; DGSVariate

DGSVariate Type
^^^^^^^^^^^^^^^

Declaration
```````````

``DGSVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::DGSTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: DGSVariate(x::Vector{Float64}, tune::DGSTune)
              DGSVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``DGSVariate`` object that stores sampled values and tuning parameters for DGS sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``DGSVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; DGSTune

DGSTune Type
^^^^^^^^^^^^^^

Declaration
```````````

``type DGSTune``

Fields
``````

* ``support::Vector`` : vector of coordinates in the parameter space from which to simulate values.
* ``probs::Vector{Float64}`` : corresponding sampling probabilities.


Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: DGS(params::Vector{Symbol})

    Construct a ``Sampler`` object for which DGS sampling is to be applied separately to each of the supplied parameters.  Parameters are assumed to have discrete univariate distributions with finite supports.

    **Arguments**

        *  ``params`` : stochastic nodes to be updated with the sampler.

    **Value**

        Returns a ``Sampler`` type object.

    **Example**

        See the :ref:`Eyes <example-Eyes>` example.
