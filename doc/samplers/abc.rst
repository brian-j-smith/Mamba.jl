.. index:: Sampling Functions; Approximate Bayesian Computing

.. _section-ABC:

Approximate Bayesian Computing
------------------------------

Approximate Bayesian Computing in the framework of MCMC as proposed by :cite:`marjoram:2003:abc` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality. 

Sampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: ABC(params::Vector{Symbol}, summarize::Function, 
                  rho::Function, epsilon::Vector{Float64}, sigma::Float64)

    Construct a ``Sampler`` object for ABC sampling.  Parameters are assumed to be continuous and unconstrained. 

    **Arguments**

        * ``params`` : stochastic nodes to be updated with the sampler.  Constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``link()`` function.
        * ``summarize`` : A function to compute summary statistics for the data.
        * ``rho`` : A distance function to compute distance between simulated data and real data. 
        * ``epsilon`` : A cut-off parameter to determine how close simulated data needs to be to accept proposal.
        ** ``sigma`` : Proposal is drawn from normal distribution centered at previous draw with ``sigma`` variance. 

    **Value**

        Returns a ``Sampler`` type object.

