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
               Mean       SD      Naive SE       MCSE        ESS  
    beta[1] 0.7534015 1.3485104 0.0082067729 0.0148617246 8233.229
    beta[2] 0.8143747 0.3953645 0.0024061117 0.0042634381 8599.544
         s2 1.1552724 4.0568132 0.0246889790 0.0473048620 7354.590
    
    Quantiles:
                 2.5%         25.0%       50.0%       75.0%     97.5%  
    beta[1] -1.3437121000 -0.008639872 0.639024602 1.30601035 3.9826175
    beta[2] -0.1752602212  0.712807275 0.866900610 1.00431638 1.4024470
         s2  0.0008981952  0.012789860 0.116695894 0.80486787 8.2198644
