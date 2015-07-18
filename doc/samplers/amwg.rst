.. index:: Sampling Functions; Adaptive Metropolis within Gibbs

Adaptive Metropolis within Gibbs (AMWG)
---------------------------------------

Implementation of a Metropolis-within-Gibbs sampler :cite:`metropolis:1953:ESC,robert:2009:EAM,tierney:1994:MCE` for iteratively simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: amwg!(v::AMWGVariate, sigma::Vector{Float64}, logf::Function; \
                    adapt::Bool=true, batchsize::Integer=50, target::Real=0.44)

    Simulate one draw from a target distribution using an adaptive Metropolis-within-Gibbs sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function should contain the ``tune`` field returned by the previous call.
        * ``sigma`` : initial standard deviations for the univariate normal proposal distributions.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.
        * ``adapt`` : whether to adaptively update the proposal distribution.
        * ``batchsize`` : number of samples that must be newly accumulated before applying an adaptive update to the proposal distributions.
        * ``target`` : a target acceptance rate for the adaptive algorithm.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-amwg:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.

        .. literalinclude:: amwg.jl
            :language: julia

.. index:: AMWGVariate

AMWGVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``AMWGVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::AMWGTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: AMWGVariate(x::Vector{Float64}, tune::AMWGTune)
              AMWGVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``AMWGVariate`` object that stores sampled values and tuning parameters for adaptive Metropolis-within-Gibbs sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``AMWGVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.


.. index:: AMWGTune

AMWGTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type AMWGTune``

Fields
``````

* ``adapt::Bool`` : whether the proposal distribution has been adaptively tuned.
* ``accept::Vector{Integer}`` : number of accepted candidate draws generated for each element of the parameter vector during adaptive updating.
* ``batchsize::Integer`` : number of samples that must be accumulated before applying an adaptive update to the proposal distributions.
* ``m::Integer`` : number of adaptive update iterations that have been performed.
* ``sigma::Vector{Float64}`` : updated values of the proposal standard deviations if ``adapt = true``, and the user-defined values otherwise.
* ``target::Real`` : target acceptance rate for the adaptive algorithm.

Sampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: AMWG(params::Vector{Symbol}, sigma::Vector{T<:Real}; \
                   adapt::Symbol=:all, batchsize::Integer=50, target::Real=0.44)

    Construct a ``Sampler`` object for adaptive Metropolis-within-Gibbs sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``link()`` function.
        * ``sigma`` : initial standard deviations for the univariate normal proposal distributions.  The standard deviations are relative to the unconstrained parameter space, where candidate draws are generated.
        * ``adapt`` : type of adaptation phase.  Options are
            * ``:all`` : adapt proposals during all iterations.
            * ``:burnin`` : adapt proposals during burn-in iterations.
            * ``:none`` : no adaptation (Metropolis-within-Gibbs sampling with fixed proposals).
        * ``batchsize`` : number of samples that must be accumulated before applying an adaptive update to the proposal distributions.
        * ``target`` : a target acceptance rate for the algorithm.

    **Value**

        Returns a ``Sampler`` type object.

    **Example**

        See the :ref:`section-Examples` section.
