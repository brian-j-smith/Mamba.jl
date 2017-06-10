.. index:: Sampling Functions; No-U-Turn Sampler

.. _section-NUTS:

No-U-Turn Sampler (NUTS)
------------------------

Implementation of the No-U-Turn Sampler extension (algorithm 6) :cite:`hoffman:2014:nuts` to Hamiltonian Monte Carlo :cite:`neal:2011:hmc` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: NUTS(params::ElementOrVector{Symbol}; dtype::Symbol=:forward, args...)

    Construct a ``Sampler`` object for NUTS sampling, with the algorithm's step size parameter adaptively tuned during burn-in iterations.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``dtype`` : type of differentiation for gradient calculations.  Options are
            * ``:central`` : central differencing.
            * ``:forward`` : forward differencing.
        * ``args...`` : additional keyword arguments to be passed to the ``NUTSVariate`` constructor.

    **Value**

        Returns a ``Sampler{NUTSTune}`` type object.

    **Example**

        See the :ref:`Dyes <example-Dyes>`, :ref:`Equiv <example-Equiv>`, and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::NUTSVariate; adapt::Bool=false)

    Draw one sample from a target distribution using the NUTS sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function will contain the ``tune`` field returned by the previous call.
        * ``adapt`` : whether to adaptively update the ``epsilon`` step size parameter.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-nuts:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.

        .. literalinclude:: nuts.jl
            :language: julia

.. index:: Sampler Types; NUTSVariate

NUTSVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``const NUTSVariate = SamplerVariate{NUTSTune}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::NUTSTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: NUTSVariate(x::AbstractVector{T<:Real}, epsilon::Real, logfgrad::Function; \
                          target::Real=0.6)
              NUTSVariate(x::AbstractVector{T<:Real}, logfgrad::Function; target::Real=0.6)

    Construct a ``NUTSVariate`` object that stores simulated values and tuning parameters for NUTS sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``epsilon`` : step size parameter.
        * ``logfgrad`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant) and gradient vector, and returns the respective results as a tuple.  If ``epsilon`` is not specified, the function is used by the constructor to generate an initial step size value.
        * ``target`` : target acceptance rate for the algorithm.

    **Value**

        Returns a ``NUTSVariate`` type object with fields set to the supplied ``x`` and tuning parameter values.


.. index:: Sampler Types; NUTSTune

NUTSTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type NUTSTune <: SamplerTune``

Fields
``````

* ``logfgrad::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density and gradient vector, or null if not supplied.
* ``adapt::Bool`` : whether the proposal distribution is being adaptively tuned.
* ``alpha::Float64`` : cumulative acceptance probabilities :math:`\alpha` from leapfrog steps.
* ``epsilon::Float64`` : updated value of the step size parameter :math:`\epsilon_m = \exp\left(\mu - \sqrt{m} \bar{H}_m / \gamma\right)` if ``m > 0``, and the user-supplied value otherwise.
* ``epsbar::Float64`` : dual averaging parameter, defined as :math:`\bar{\epsilon}_m = \exp\left(m^{-\kappa} \log(\epsilon_m) + (1 - m^{-\kappa}) \log(\bar{\epsilon}_{m-1})\right)`.
* ``gamma::Float64`` : dual averaging parameter, fixed at :math:`\gamma = 0.05`.
* ``Hbar::Float64`` : dual averaging parameter, defied as :math:`\bar{H}_m = \left(1 - \frac{1}{m + t_0}\right) \bar{H}_{m-1} + \frac{1}{m + t_0} \left(\text{target} - \frac{\alpha}{n_\alpha}\right)`.
* ``kappa::Float64`` : dual averaging parameter, fixed at :math:`\kappa = 0.05`.
* ``m::Int`` : number of adaptive update iterations :math:`m` that have been performed.
* ``mu::Float64`` : dual averaging parameter, defined as :math:`\mu = \log(10 \epsilon_0)`.
* ``nalpha::Int`` : the total number :math:`n_\alpha` of leapfrog steps performed.
* ``t0::Float64`` : dual averaging parameter, fixed at :math:`t_0 = 10`.
* ``target::Float64`` : target acceptance rate for the adaptive algorithm.
