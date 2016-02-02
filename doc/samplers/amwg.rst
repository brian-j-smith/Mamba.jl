.. index:: Sampling Functions; Adaptive Metropolis within Gibbs

.. _section-AMWG:

Adaptive Metropolis within Gibbs (AMWG)
---------------------------------------

Implementation of a Metropolis-within-Gibbs sampler :cite:`metropolis:1953:ESC,robert:2009:EAM,tierney:1994:MCE` for iteratively simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: AMWG(params::ElementOrVector{Symbol}, \
                   sigma::ElementOrVector{T<:Real}; adapt::Symbol=:all, args...)

    Construct a ``Sampler`` object for AMWG sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``sigma`` : scaling value or vector of the same length as the combined elements of nodes ``params``, defining initial standard deviations for univariate normal proposal distributions.  Standard deviations are relative to the unconstrained parameter space, where candidate draws are generated.
        * ``adapt`` : type of adaptation phase.  Options are
            * ``:all`` : adapt proposals during all iterations.
            * ``:burnin`` : adapt proposals during burn-in iterations.
            * ``:none`` : no adaptation (Metropolis-within-Gibbs sampling with fixed proposals).
        * ``args...`` : additional keyword arguments to be passed to the ``AMWGVariate`` constructor.

    **Value**

        Returns a ``Sampler{AMWGTune}`` type object.

    **Example**

        See the :ref:`Birats <example-Birats>`, :ref:`Blocker <example-Blocker>`, and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::AMWGVariate; adapt::Bool=true)

    Draw one sample from a target distribution using the AMWG sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function will contain the ``tune`` field returned by the previous call.
        * ``adapt`` : whether to adaptively update the proposal distribution.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-amwg:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.  Also, see the :ref:`example-Line_AMWG_Slice` example.

        .. literalinclude:: amwg.jl
            :language: julia

.. index:: Sampler Types; AMWGVariate

AMWGVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``typealias AMWGVariate SamplerVariate{AMWGTune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::AMWGTune`` : tuning parameters for the sampling algorithm.

Constructor
```````````

.. function:: AMWGVariate(x::AbstractVector{T<:Real}, \
                          sigma::ElementOrVector{U<:Real}, logf::Function; \
                          batchsize::Integer=50, target::Real=0.44)

    Construct an ``AMWGVariate`` object that stores simulated values and tuning parameters for AMWG sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``sigma`` : scaling value or vector of the same length as the combined elements of nodes ``params``, defining initial standard deviations for univariate normal proposal distributions.  Standard deviations are relative to the unconstrained parameter space, where candidate draws are generated.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``batchsize`` : number of samples that must be accumulated before applying an adaptive update to the proposal distributions.
        * ``target`` : target acceptance rate for the algorithm.

    **Value**

        Returns an ``AMWGVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.


.. index:: Sampler Types; AMWGTune

AMWGTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type AMWGTune <: SamplerTune``

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``adapt::Bool`` : whether the proposal distribution is being adaptively tuned.
* ``accept::Vector{Int}`` : number of accepted candidate draws generated for each element of the parameter vector during adaptive updating.
* ``batchsize::Int`` : number of samples that must be accumulated before applying an adaptive update to the proposal distributions.
* ``m::Int`` : number of adaptive update iterations that have been performed.
* ``sigma::Vector{Float64}`` : updated values of the proposal standard deviations if ``m > 0``, and user-supplied values otherwise.
* ``target::Float64`` : target acceptance rate for the adaptive algorithm.
