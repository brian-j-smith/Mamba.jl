.. index:: Sampling Functions; Metropolis-adjusted Langevin algorithm 

.. _section-MALA:

Metropolis-adjusted Langevin algorithm (MALA)
---------------------------------------------

Implementation of the Metropolis-adjusted Langevin algorithm of Robers :cite:`roberts:1996:MALA`.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: mala!(v::MALAVariate, scale::Float64, U::Cholesky{Float64}, fx::Function)

    Simulate one draw from a target distribution using the BHMC sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``scale`` : The scale.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-mala:

    **Example**

        .. literalinclude:: mala.jl
            :language: julia


.. index:: Sampler Types; MALAVariate

MALAVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``MALAVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::MALATune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: MALAVariate(x::Vector{Float64}, tune::MALATune)
              MALAVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``MALAVariate`` object that stores sampled values and tuning parameters for MALA sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``MALAVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; MALATune

MALATune Type
^^^^^^^^^^^^^

Declaration
```````````

``type MALATune``

Fields
``````
* ``U::Cholesky{Float64}`` : Cholesky decomposition of covariance matrix. 
* ``scale::Float64`` : Step size. 

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: MALA(params::Vector{Symbol}, scale::Float64, M::Matrix{Float64}, dtype::Symbol=:forward)

    Construct a ``Sampler`` object for MALA sampling.  Parameters are assumed to be unconstrained.

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``link()`` function.
        * ``scale`` : step size.
      	* ``M`` : Covariance matrix.
      	* ``dtype`` : type of differentiation for gradient calculations. Options are
	          * ``:central`` : central differencing.
	          * ``:forward`` : forward differencing.
        
    **Value**

        Returns a ``Sampler`` type object.
