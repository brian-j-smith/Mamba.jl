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

    Iterations = 1001:10000
    Thinning interval = 1
    Chains = 1,2,3
    Samples per chain = 9000
    
    Empirical Posterior Estimates:
               Mean        SD       Naive SE       MCSE       ESS   
    beta[1] 0.63133698 1.30979294 0.0079711460 0.049024555 713.80234
    beta[2] 0.79995255 0.38169153 0.0023229007 0.012115902 992.46124
         s2 1.18427121 2.49666395 0.0151942129 0.131048462 362.95796
    
    Quantiles:
                2.5%        25.0%       50.0%     75.0%      97.5%  
    beta[1] -1.841697099 -0.09055321 0.59443124 1.2937670  3.3970011
    beta[2]  0.039303260  0.58749126 0.80036075 1.0210728  1.5305396
         s2  0.009203774  0.14104046 0.36967615 0.9555609 13.3210173

