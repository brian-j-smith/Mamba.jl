.. index:: Sampling Functions; Approximate Bayesian Computation

.. _section-ABC:

Approximate Bayesian Computation (ABC)
--------------------------------------

Approximate Bayesian Computation in the framework of MCMC (also known as Likelihood-Free MCMC) as proposed by :cite:`marjoram:2003:abc` for simulating autocorrelated draws from a posterior distribution without evaluating its likelihood.  Also see :cite:`sisson:2011:ABC` for a thorough review of Likelihood-Free MCMC.

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: ABC(params::ElementOrVector{Symbol}, \
                  scale::ElementOrVector{T<:Real}, summary::Function, \
                  epsilon::Real; kernel::KernelDensityType=SymUniform, \
                  dist::Function=(Tsim, Tobs) -> sqrt(sumabs2(Tsim - Tobs)), \
                  proposal::SymDistributionType=Normal, maxdraw::Integer=1, \
                  nsim::Integer=1, randeps::Bool=false, ratio::Float64=1.0, \
                  args...)

    Construct a ``Sampler`` object for ABC sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function.
        * ``scale`` : scaling value or vector of the same length as the combined elements of nodes ``params`` for the ``proposal`` distribution.  Values are relative to the unconstrained parameter space, where candidate draws are generated.
        * ``summary`` : function that takes a vector of observed or simulated data and returns a summary statistic or vector of statistics.
        * ``epsilon`` : target tolerance for determining how similar observed and simulated data summary statistics need to be in order to accept a candidate draw.  Internal tolerances are adaptively tuned at each iteration to decrease monotonically to this target.
        * ``kernel`` : weighting kernel density of type ``Biweight``, ``Cosine``, ``Epanechnikov``, ``Normal``, ``SymTriangularDist``, ``SymUniform``, or ``Triweight`` to use in measuring similarity between observed and simulated data summary statistics.  Specified ``epsilon`` determines the standard deviation of Normal kernels and widths of the others.
        * ``dist`` : positive function for the kernel density to compute distance between vectors of observed (``Tobs``) and simulated (``Tsim``) data summary statistics (default: Euclidean distance).
        * ``proposal`` : symmetric distribution of type ``Biweight``, ``Cosine``, ``Epanechnikov``, ``Normal``, ``SymTriangularDist``, ``SymUniform``, or ``Triweight`` to be centered around current parameter values and used to generate proposal draws.  Specified ``scale`` determines the standard deviations of Normal proposals and widths of the others.
        * ``maxdraw`` : maximum number of unaccepted candidates to draw in each call of the sampler.  Draws are generated until one is accepted or the maximum is reached.  Larger values increase acceptance rates at the expense of longer runtimes.
        * ``nsim`` : number of data sets to simulate in deciding whether to accept a candidate draw.  Larger values lead to closer approximations of the target distribution at the expense of longer runtimes.
        * ``randeps`` : Whether to perturb internal tolerance by random exponential variable. 
        * ``ratio`` : Internal tolerance is monotonically decreased by convex combination of previous tolerance and distance between vectors of observed and simulated data summary statics, with weight ``ratio`` given the distance. 
        * ``args...`` : additional keyword arguments to be passed to the ``dist`` function.

    **Value**

        Returns a ``Sampler{ABCTune}`` type object.

    **Example**

        See the :ref:`ABC Line <example-Line_ABC>` and other :ref:`section-Examples`.

.. index:: Sampler Types; ABCTune

ABCTune Type
^^^^^^^^^^^^

Declaration
```````````

``type ABCTune``

Fields
``````

* ``datakeys::Vector{Symbol}`` : stochastic "data" nodes in the full conditional distribution for parameters to be updated and nodes at which summary statistics are computed separately in the sampling algorithm.
* ``Tsim::Vector{Vector{Float64}}`` : simulated data summary statistics for the ``nsim`` data sets.
* ``epsilon::Vector{Float64}`` : adaptively tuned tolerances for the data sets.
* ``mean_eps::Vector{Float64}`` : mean of adaptively tuned tolerances for the data sets (if ``randeps`` is true). 
* ``curr_eps::Vector{Float64}`` : adaptively tuned tolerances for the data sets.
