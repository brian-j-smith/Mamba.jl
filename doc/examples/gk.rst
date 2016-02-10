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
         Mean         SD        Naive SE       MCSE        ESS
    k 0.35177956 0.114186840 0.00076124560 0.0070749476 260.48676
    A 3.00044501 0.056131349 0.00037420899 0.0021176913 702.56362
    B 1.05879179 0.104464188 0.00069642792 0.0059092212 312.51753
    g 2.04113191 0.317411213 0.00211607476 0.0134387307 557.86361

    Quantiles:
         2.5%       25.0%      50.0%      75.0%     97.5%
    k 0.14224981 0.27641846 0.33980515 0.41834851 0.6208537
    A 2.88887160 2.96200257 3.00069478 3.03852346 3.1088669
    B 0.84952558 0.99167335 1.06172076 1.12769693 1.2681306
    g 1.56703473 1.81345003 1.99234172 2.21926781 2.8253217
