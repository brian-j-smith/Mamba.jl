.. index:: Sampling Functions; Simplified Manifold Metropolis-Adjusted Langevin Algorithm

.. _section-SMMALA:

Simplified Manifold Metropolis-Adjusted Langevin Algorithm (SMMALA)
---------------------------------------------

Implementation of the Simplified Manifold Metropolis-Adjusted Langevin Algorithm of Girolami and Calderhead :cite:`girolami:2011:RMHMC`.

Model-Based Constructors
^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SMMALA(params::ElementOrVector{Symbol}, epsilon::Real)
              SMMALA(params::ElementOrVector{Symbol}, epsilon::Real, \
                   Sigma::Matrix{T<:Real})

    Construct a ``Sampler`` object for SMMALA sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``epsilon`` : factor by which the drift and covariance matrix of the proposal distribution are scaled.
        * ``Sigma`` : covariance matrix for the multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.  If omitted, the identity matrix is assumed.

    **Value**

        Returns a ``Sampler{SMMALATune}`` type object.

    **Example**

        See the :ref:`Dyes <example-Dyes>` and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::SMMALAVariate)

    Draw one sample from a target distribution using the SMMALA sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-smmala:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.

        .. literalinclude:: smmala.jl
            :language: julia


.. index:: Sampler Types; SMMALAVariate

SMMALAVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``const SMMALAVariate = SamplerVariate{SMMALATune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::SMMALATune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: SMMALAVariate(x::AbstractVector{T<:Real}, epsilon::Real, logfgradhess::Function)
              SMMALAVariate(x::AbstractVector{T<:Real}, epsilon::Real, \
                          Sigma::Matrix{U<:Real}, logfgradhess::Function)

    Construct a ``SMMALAVariate`` object that stores simulated values and tuning parameters for SMMALA sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``epsilon`` : factor by which the drift and covariance matrix of the proposal distribution are scaled.
        * ``Sigma`` : covariance matrix for the multivariate normal proposal distribution.  The covariance matrix is relative to the unconstrained parameter space, where candidate draws are generated.  If omitted, the identity matrix is assumed.
        * ``logfgradhess`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant), gradient vector, and hessian matrix. Returns the respective results as a tuple.

    **Value**

        Returns a ``SMMALAVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; SMMALATune

SMMALATune Type
^^^^^^^^^^^^^

Declaration
```````````

``type SMMALATune <: SamplerTune``

Fields
``````

* ``logfgradhess::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, gradient vector, and hessian. Null if not supplied.
* ``epsilon::Float64`` : factor by which the drift and covariance matrix of the proposal distribution are scaled.
* ``SigmaL::Union{UniformScaling{Int}, LowerTriangular{Float64}}`` : Cholesky factorization of the covariance matrix for the multivariate normal proposal distribution.
