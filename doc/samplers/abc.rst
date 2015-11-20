.. index:: Sampling Functions; Approximate Bayesian Computing

.. _section-ABC:

Approximate Bayesian Computing
------------------------------

Approximate Bayesian Computing in the framework of MCMC (also known as Likelihood-Free MCMC) as proposed by :cite:`marjoram:2003:abc` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality. Also see :cite:`sisson:2011:ABC` for a thorough review of Likelihood-Free MCMC. 

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: ABC(params::Vector{Symbol}, sigma::Real, summarize::Vector{Function}, \
                  epsilon::Real; rho::Function = (x, y) -> sqrt(sum((x - y) .^ 2)), \
                  s::Int = 1, n::Int = 50, kernel::Symbol = :uniform)

    Construct a ``Sampler`` object for ABC sampling.  Parameters are assumed to be continuous and unconstrained. 

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``link()`` function.
        * ``sigma`` : Proposal is drawn from normal distribution centered at previous draw with ``sigma`` variance. 
        * ``summarize`` : A vector of functions to compute summary statistics for the data.
        * ``epsilon`` : A cut-off parameter to determine how close simulated data needs to be to accept proposal.
        * ``rho`` : A distance function to compute distance between simulated data and real data. Default is Euclidean distance. 
        * ``s`` : Number of data sets to simulate for each proposal. 
        * ``n`` : Number of simulations to accept new proposal before staying with current value.
        * ``kernel`` : Weighting density kernel. One of :uniform, :epanechnikov, :gaussian

    **Value**

        Returns a ``Sampler`` type object.

