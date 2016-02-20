.. index:: Examples; GK: Approximate Bayesian Computation

.. _example-gk:

GK: Approximate Bayesian Computation
------------------------------------

Approximate Bayesian Computation (ABC) is often useful when the likelihood is costly to compute.  The generalized GK distribution :cite:`rayner:2002:GK` is a distribution defined by the inverse of its cumulative distribution function.  It is therefore easy to sample from, but difficult to obtain analytical likelihood functions.  These properties make the GK distribution well suited for ABC.  The following is a simulation study that mirrors :cite:`allingham:2009:ABCGK`.

Model
^^^^^
.. math::

    x_{i} &\sim \text{GK}\left(A, B, g, k\right) \quad\quad i=1, \ldots, 1000 \\
    A, B, g, k &\sim \text{Uniform}(0, 10)

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: gk.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2501:10000
    Thinning interval = 1
    Chains = 1,2,3
    Samples per chain = 7500

    Empirical Posterior Estimates:
         Mean        SD        Naive SE       MCSE         ESS
    k 0.3510874 0.132239097 0.00088159398 0.0082306353  258.13865
    A 3.0037486 0.068834777 0.00045889851 0.0018699334 1355.07569
    B 1.0575823 0.129916443 0.00086610962 0.0074502510  304.07900
    g 2.0259195 0.311756472 0.00207837648 0.0114243024  744.68323

    Quantiles:
         2.5%       25.0%      50.0%      75.0%     97.5%
    k 0.12902213 0.25821839 0.34355889 0.42948590 0.6759640
    A 2.86683504 2.95778075 3.00296541 3.04971216 3.1358843
    B 0.80207071 0.96839859 1.05571418 1.14595782 1.3123378
    g 1.56533594 1.79945626 1.97432281 2.19565833 2.7212103
