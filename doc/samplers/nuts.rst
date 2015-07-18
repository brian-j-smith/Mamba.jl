.. index:: Sampling Functions; No-U-Turn Sampler

No-U-Turn Sampler (NUTS)
------------------------

Implementation of the NUTS extension (algorithm 6) :cite:`hoffman:2014:nuts` to Hamiltonian Monte Carlo :cite:`neal:2011:hmc` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: nutsepsilon(v::NUTSVariate, fx::Function)

    Generate an initial value for the step size parameter of the No-U-Turn sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : the current state of parameters to be simulated.
        * ``fx`` : function to compute the log-transformed density (up to a normalizing constant) and gradient vector at ``v.value``, and to return the respective results as a tuple.

    **Value**

        A numeric step size value.

.. function:: nuts!(v::NUTSVariate, epsilon::Real, fx::Function; adapt::Bool=false, \
                    target::Real=0.6)

    Simulate one draw from a target distribution using the No-U-Turn sampler.  Parameters are assumed to be continuous and unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.  When running the sampler in adaptive mode, the ``v`` argument in a successive call to the function should contain the ``tune`` field returned by the previous call.
        * ``epsilon`` : the NUTS algorithm step size parameter.
        * ``fx`` : function to compute the log-transformed density (up to a normalizing constant) and gradient vector at ``v.value``, and to return the respective results as a tuple.
        * ``adapt`` : whether to adaptively update the ``epsilon`` step size parameter.
        * ``target`` : a target acceptance rate for the algorithm.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-nuts:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.

        .. literalinclude:: nuts.jl
            :language: julia

.. index:: NUTSVariate

NUTSVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``NUTSVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::NUTSTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: NUTSVariate(x::Vector{Float64}, tune::NUTSTune)
              NUTSVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``NUTSVariate`` object that stores sampled values and tuning parameters for No-U-Turn sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``NUTSVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.


.. index:: NUTSTune

NUTSTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type NUTSTune``

Fields
``````
* ``adapt::Bool`` : whether the proposal distribution has been adaptively tuned.
* ``alpha::Float64`` : cumulative acceptance probabilities :math:`\alpha` from leapfrog steps.
* ``epsilon::Float64`` : updated value of the step size parameter :math:`\epsilon_m = \exp\left(\mu - \sqrt{m} \bar{H}_m / \gamma\right)` if ``adapt = true``, and the user-defined value otherwise.
* ``epsbar::Float64`` : dual averaging parameter, defined as :math:`\bar{\epsilon}_m = \exp\left(m^{-\kappa} \log(\epsilon_m) + (1 - m^{-\kappa}) \log(\bar{\epsilon}_{m-1})\right)`.
* ``gamma::Float64`` : dual averaging parameter, fixed at :math:`\gamma = 0.05`.
* ``Hbar::Float64`` : dual averaging parameter, defied as :math:`\bar{H}_m = \left(1 - \frac{1}{m + t_0}\right) \bar{H}_{m-1} + \frac{1}{m + t_0} \left(\text{target} - \frac{\alpha}{n_\alpha}\right)`.
* ``kappa::Float64`` : dual averaging parameter, fixed at :math:`\kappa = 0.05`.
* ``m::Integer`` : number of adaptive update iterations :math:`m` that have been performed.
* ``mu::Float64`` : dual averaging parameter, defined as :math:`\mu = \log(10 \epsilon_0)`.
* ``nalpha::Integer`` : the total number :math:`n_\alpha` of leapfrog steps performed.
* ``t0::Float64`` : dual averaging parameter, fixed at :math:`t_0 = 10`.
* ``target::Float64`` : target acceptance rate for the adaptive algorithm.

Sampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: NUTS(params::Vector{Symbol}; dtype::Symbol=:forward, \
                   target::Real=0.6)

    Construct a ``Sampler`` object for No-U-Turn sampling, with the algorithm's step size parameter adaptively tuned during burn-in iterations.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``link()`` function.
        * ``dtype`` : type of differentiation for gradient calculations.  Options are
            * ``:central`` : central differencing.
            * ``:forward`` : forward differencing.
        * ``target`` : a target acceptance rate for the algorithm.

    **Value**

        Returns a ``Sampler`` type object.

    **Example**

        See the :ref:`section-Examples` section.
