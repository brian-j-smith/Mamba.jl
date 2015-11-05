.. index:: Sampling Functions; Hamiltonian Monte Carlo

.. _section-HMC:

Hamiltonian Monte Carlo (HMC)
-----------------------------
Implementation of the Hybrid Monte Carlo (also known as Hamiltonian Monte Carlo) of Duane :cite:`duane:1987:hmc`. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.  Code is derived from Neal's implementation :cite:`neal:2011:hmc`.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: hmc!(v::HMCVariate, epsilon::Real, L::Integer, fx::Function)
              hmc!(v::HMCVariate, epsilon::Real, L::Integer, \
                   SigmaF::Cholesky{Float64}, fx::Function)

    Simulate one draw from a target distribution using the HMC sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``epsilon`` : step size.
        * ``L`` : number of steps to take in the Leapfrog algorithm.
        * ``SigmaF`` : Cholesky factorization of the covariance matrix for the multivariate normal proposal distribution.  If omitted, the identity matrix is assumed.
        * ``fx`` : function that takes a single ``DenseVector`` argument at which to compute the log-transformed density (up to a normalizing constant) and gradient vector, and returns the respective results as a tuple.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.


    **Example**

        .. literalinclude:: hmc.jl
            :language: julia


.. index:: Sampler Types; HMCVariate

HMCVariate Type
^^^^^^^^^^^^^^^

Declaration
```````````

``HMCVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::HMCTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: HMCVariate(x::AbstractVector{T<:Real}, tune::HMCTune)
              HMCVariate(x::AbstractVector{T<:Real}, tune=nothing)

    Construct a ``HMCVariate`` object that stores sampled values and tuning parameters for HMC sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``HMCVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; HMCTune

HMCTune Type
^^^^^^^^^^^^

Declaration
```````````

``type HMCTune``

Fields
``````

* ``epsilon::Float64`` : step size.
* ``L::Int`` : number of steps to take in the Leapfrog algorithm.
* ``SigmaF::Cholesky{Float64}`` : Cholesky factorization of the covariance matrix for the multivariate normal proposal distribution.

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: HMC(params::Vector{Symbol}, epsilon::Real, L::Integer; \
                  dtype::Symbol=:forward)
              HMC(params::Vector{Symbol}, epsilon::Real, L::Integer, \
                  Sigma::Matrix{T<:Real}; dtype::Symbol=:forward)

    Construct a ``Sampler`` object for HMC sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``epsilon`` : step size.
        * ``L`` : number of steps to take in the Leapfrog algorithm.
        * ``Sigma`` : covariance matrix for the multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.  If omitted, the identity matrix is assumed.
        * ``dtype`` : type of differentiation for gradient calculations. Options are
            * ``:central`` : central differencing.
            * ``:forward`` : forward differencing.

    **Value**

        Returns a ``Sampler`` type object.

    **Example**

        See the :ref:`Dyes <example-Dyes>` and other :ref:`section-Examples`.
