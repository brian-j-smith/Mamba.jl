.. index:: Sampling Functions; Metropolis-Adjusted Langevin Algorithm

.. _section-MALA:

Metropolis-Adjusted Langevin Algorithm (MALA)
---------------------------------------------

Implementation of the Metropolis-Adjusted Langevin Algorithm (MALA) of Roberts and Tweedie :cite:`roberts:1996:MALA` and Roberts and Stramer :cite:`roberts:2002:LD`.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.  MALA is related to Hamiltonian Monte Carlo as described thoroughly by Girolami and Calderhead :cite:`girolami:2011:RMHMC`.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: mala!(v::MALAVariate, scale::Real, fx::Function)
              mala!(v::MALAVariate, scale::Real, SigmaF::Cholesky{Float64}, \
                    fx::Function)

    Simulate one draw from a target distribution using the MALA sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``scale`` : factor by which the drift and covariance matrix of the proposal distribution are scaled.
        * ``SigmaF`` : Cholesky factorization of the covariance matrix for the multivariate normal proposal distribution.  If omitted, the identity matrix is assumed.
        * ``fx`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant) and gradient vector, and returns the respective results as a tuple.

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

``MALAVariate <: SamplerVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::MALATune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: MALAVariate(x::AbstractVector{T<:Real}, tune::MALATune)
              MALAVariate(x::AbstractVector{T<:Real}, tune=nothing)

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

``type MALATune <: SamplerTune``

Fields
``````

* ``scale::Float64`` : factor by which the drift and covariance matrix of the proposal distribution are scaled.
* ``SigmaF::Cholesky{Float64}`` : Cholesky factorization of the covariance matrix for the multivariate normal proposal distribution.

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: MALA(params::Vector{Symbol}, scale::Real; dtype::Symbol=:forward)
              MALA(params::Vector{Symbol}, scale::Real, Sigma::Matrix{T<:Real}; \
                   dtype::Symbol=:forward)

    Construct a ``Sampler`` object for MALA sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``scale`` : factor by which the drift and covariance matrix of the proposal distribution are scaled.
        * ``Sigma`` : covariance matrix for the multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.  If omitted, the identity matrix is assumed.
        * ``dtype`` : type of differentiation for gradient calculations. Options are
            * ``:central`` : central differencing.
            * ``:forward`` : forward differencing.

    **Value**

        Returns a ``Sampler`` type object.

    **Example**

        See the :ref:`Dyes <example-Dyes>` and other :ref:`section-Examples`.
