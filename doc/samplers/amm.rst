.. index:: Sampling Functions; Adaptive Mixture Metropolis

.. _section-AMM:

Adaptive Mixture Metropolis (AMM)
---------------------------------

Implementation of the Roberts and Rosenthal :cite:`robert:2009:EAM` adaptive (multivariate) mixture Metropolis :cite:`haario:2001:AMA,hastings:1970:MCS,metropolis:1953:ESC` sampler for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: AMM(params::ElementOrVector{Symbol}, Sigma::Matrix{T<:Real}; \
                  adapt::Symbol=:all, args...)

    Construct a ``Sampler`` object for AMM sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``Sigma`` : covariance matrix for the non-adaptive multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.
        * ``adapt`` : type of adaptation phase.  Options are
            * ``:all`` : adapt proposal during all iterations.
            * ``:burnin`` : adapt proposal during burn-in iterations.
            * ``:none`` : no adaptation (multivariate Metropolis sampling with fixed proposal).
        * ``args...`` : additional keyword arguments to be passed to the ``AMMVariate`` constructor.

    **Value**

        Returns a ``Sampler{AMMTune}`` type object.

    **Example**

        See the :ref:`Seeds <example-Seeds>` and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::AMMVariate; adapt::Bool=true)

    Draw one sample from a target distribution using the AMM sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function will contain the ``tune`` field returned by the previous call.
        * ``adapt`` : whether to adaptively update the proposal distribution.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-amm:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.

        .. literalinclude:: amm.jl
            :language: julia

.. index:: Sampler Types; AMMVariate

AMMVariate Type
^^^^^^^^^^^^^^^

Declaration
```````````

``typealias AMMVariate SamplerVariate{AMMTune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::AMMTune`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: AMMVariate(x::AbstractVector{T<:Real}, Sigma::Matrix{U<:Real}, logf::Function; \
                         beta::Real=0.05, scale::Real=2.38)

    Construct an ``AMMVariate`` object that stores simulated values and tuning parameters for AMM sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``Sigma`` : covariance matrix for the non-adaptive multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``beta`` : proportion of weight given to draws from the non-adaptive proposal with covariance factorization ``SigmaL``, relative to draws from the adaptively tuned proposal with covariance factorization ``SigmaLm``, during adaptive updating.
        * ``scale`` : factor (``scale^2 / length(x)``) by which the adaptively updated covariance matrix is scaled---default value adopted from Gelman, Roberts, and Gilks :cite:`gelman:1996:EMJ`.

    **Value**

        Returns an ``AMMVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.


.. index:: Sampler Types; AMMTune

AMMTune Type
^^^^^^^^^^^^

Declaration
```````````

``type AMMTune <: SamplerTune``

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``adapt::Bool`` : whether the proposal distribution is being adaptively tuned.
* ``beta::Float64`` : proportion of weight given to draws from the non-adaptive proposal with covariance factorization ``SigmaL``, relative to draws from the adaptively tuned proposal with covariance factorization ``SigmaLm``, during adaptive updating.
* ``m::Int`` : number of adaptive update iterations that have been performed.
* ``Mv::Vector{Float64}`` : running mean of draws ``v`` during adaptive updating.  Used in the calculation of ``SigmaLm``.
* ``Mvv::Matrix{Float64}`` : running mean of ``v * v'`` during adaptive updating.  Used in the calculation of ``SigmaLm``.
* ``scale::Float64`` : factor (``scale^2 / length(v)``) by which the adaptively updated covariance matrix is scaled.
* ``SigmaL::LowerTriangular{Float64}`` : Cholesky factorization of the non-adaptive covariance matrix.
* ``SigmaLm::Matrix{Float64}`` : pivoted factorization of the adaptively tuned covariance matrix.
