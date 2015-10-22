.. index:: Sampling Functions; Hamiltonian Monte Carlo

.. _section-HMC:

Hamiltonian Monte Carlo (HMC)
-----------------------------
Implementation of the Hybrid Monte Carlo (also known as Hamiltonian Monte Carlo) of Duane :cite:`duane:1987:hmc`. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: hmc!(v::HMCVariate, M::Matrix{Float64}, epsilon::Float64, L::Int, fx::Function)

    Simulate one draw from a target distribution using the BHMC sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function should contain the ``tune`` field returned by the previous call.
        * ``M`` : Covariance matrix for the auxiliary momentum variables. 
        * ``epsilon`` : the HMC algorithm step size parameter.
        * ``L`` : Number of steps to take in the Leapfrog algorithm. 
        * ``fx`` : function that takes a single ``DenseVector`` argument at which to compute the log-transformed density (up to a normalizing constant) and gradient vector, and returns the respective results as a tuple.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.


    **Example**

        .. literalinclude:: hmc.jl
            :language: julia


.. index:: Sampler Types; HMCVariate

HMCVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``HMCVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::HMCTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: HMCVariate(x::Vector{Float64}, tune::HMCTune)
              HMCVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``HMCVariate`` object that stores sampled values and tuning parameters for HMC sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``HMCVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; HMCTune

HMCTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type HMCTune``

Fields
``````
* ``M::Matrix{Float64}`` : Covariance matrix for the auxiliary momentum variables. 
* ``epsilon::Float64`` : the HMC algorithm step size parameter.
* ``L::Int64`` : Number of steps to take in the Leapfrog algorithm. 

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: HMC{T<:Real}(params::Vector{Symbol}, M::Matrix{T}, epsilon::Float64, L::Int; dtype::Symbol=:forward)

    Construct a ``Sampler`` object for HMC sampling.  Parameters are assumed to be unconstrained.

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``link()`` function.
      	* ``M`` : Covariance matrix.
        * ``epsilon`` : the HMC algorithm step size parameter.
        * ``L`` : Number of steps to take in the Leapfrog algorithm. 
      	* ``dtype`` : type of differentiation for gradient calculations. Options are
	          * ``:central`` : central differencing.
	          * ``:forward`` : forward differencing.
        
    **Value**

        Returns a ``Sampler`` type object.
