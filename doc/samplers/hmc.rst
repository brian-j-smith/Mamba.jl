.. index:: Sampling Functions; Hamiltonian Monte Carlo

.. _section-HMC:

Hamiltonian Monte Carlo (HMC)
-----------------------------
Implementation of the Hybrid Monte Carlo (also known as Hamiltonian Monte Carlo) of Duane :cite:`duane:1987:hmc`. The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.  Code is derived from Neal's implementation :cite:`neal:2011:hmc`.

Model-Based Constructors
^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: HMC(params::ElementOrVector{Symbol}, epsilon::Real, L::Integer; \
                  dtype::Symbol=:forward)
              HMC(params::ElementOrVector{Symbol}, epsilon::Real, \
                  L::Integer, Sigma::Matrix{T<:Real}; dtype::Symbol=:forward)

    Construct a ``Sampler`` object for HMC sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``epsilon`` : step size.
        * ``L`` : number of steps to take in the Leapfrog algorithm.
        * ``Sigma`` : covariance matrix for the multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.  If omitted, the identity matrix is assumed.
        * ``dtype`` : type of differentiation for gradient calculations. Options are
            * ``:central`` : central differencing.
            * ``:forward`` : forward differencing.

    **Value**

        Returns a ``Sampler{HMCTune}`` type object.

    **Example**

        See the :ref:`Dyes <example-Dyes>` and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::HMCVariate)

    Draw one sample from a target distribution using the HMC sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.


    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.

        .. literalinclude:: hmc.jl
            :language: julia


.. index:: Sampler Types; HMCVariate

HMCVariate Type
^^^^^^^^^^^^^^^

Declaration
```````````

``typealias HMCVariate SamplerVariate{HMCTune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::HMCTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: HMCVariate(x::AbstractVector{T<:Real}, epsilon::Real, L::Integer, \
                         logfgrad::Function)
              HMCVariate(x::AbstractVector{T<:Real}, epsilon::Real, L::Integer, \
                         Sigma::Matrix{U<:Real}, logfgrad::Function)

    Construct an ``HMCVariate`` object that stores simulated values and tuning parameters for HMC sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``epsilon`` : step size.
        * ``L`` : number of steps to take in the Leapfrog algorithm.
        * ``Sigma`` : covariance matrix for the multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.  If omitted, the identity matrix is assumed.
        * ``logfgrad`` : function that takes a single ``DenseVector`` argument at which to compute the log-transformed density (up to a normalizing constant) and gradient vector, and returns the respective results as a tuple.

    **Value**

        Returns an ``HMCVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; HMCTune

HMCTune Type
^^^^^^^^^^^^

Declaration
```````````

``type HMCTune <: SamplerTune``

Fields
``````

* ``logfgrad::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density and gradient vector, or null if not supplied.
* ``epsilon::Float64`` : step size.
* ``L::Int`` : number of steps to take in the Leapfrog algorithm.
* ``SigmaL::Union{UniformScaling{Int}, LowerTriangular{Float64}}`` : Cholesky factorization of the covariance matrix for the multivariate normal proposal distribution.
