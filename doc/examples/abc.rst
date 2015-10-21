.. index:: Examples; Approximate Bayesian Computing

.. _example-ABC:

Approximate Bayesian Computing
------------------------------

A simple example to demonstrate the Approximate Bayesian Computing sampler within the MCMC framework based on the linear regression model defined in the :ref:`section-Line` section.

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: abc.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia


    Iterations = 1001:100000
    Thinning interval = 1
    Chains = 1,2,3
    Samples per chain = 99000

    Empirical Posterior Estimates:
               Mean        SD       Naive SE       MCSE        ESS
         s2 24.4752508 64.7170427 0.1187518651 0.598603593 11688.4975
    beta[1]  4.1692294  3.9195078 0.0071920602 0.071406491  3012.9206
    beta[2] -0.6555909  1.0995553 0.0020176176 0.019950447  3037.5883

    Quantiles:
                2.5%       25.0%      50.0%      75.0%       97.5%
         s2  0.29037695  2.0887221  8.1153310 24.14863892 145.1549089
    beta[1] -1.56050545  0.6667326  4.1605370  7.43427276  11.2155627
    beta[2] -2.74257047 -1.2601862 -0.8254984  0.28168134   1.2607243