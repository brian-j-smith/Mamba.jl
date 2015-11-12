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
               Mean        SD       Naive SE       MCSE        ESS   
    beta[1] 0.62204275 0.95161227 0.0057913278 0.0122490112 6035.5660
    beta[2] 0.84299144 0.21170344 0.0012883861 0.0018919009 9000.0000
         s2 0.08700958 0.25847592 0.0015730344 0.0047900300 2911.8147
    
    Quantiles:
                 2.5%          25.0%       50.0%       75.0%     97.5%  
    beta[1] -1.20663682420 -0.064975089 0.644823296 1.34406084 2.3047911
    beta[2]  0.44800996151  0.685582443 0.847293522 1.00371419 1.2105923
         s2  0.00085608263  0.005821443 0.020176814 0.07216813 0.5734865
