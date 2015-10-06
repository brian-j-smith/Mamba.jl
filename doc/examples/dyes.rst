.. index:: Examples; Dyes: Variance Components Model

.. _example-Dyes:

Dyes: Variance Components Model
-------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`, Davies :cite:`davies:1967:SMR`, and Box and Tiao :cite:`box:1973:BIS` concerning batch-to-batch variation in yields from six batches and five samples of dyestuff.

Model
^^^^^

Yields are modelled as

.. math::

    y_{i,j} &\sim \text{Normal}(\mu_i, \sigma_\text{within}) \quad\quad i=1,\ldots,6; j=1,\ldots,5 \\
    \mu_i &\sim \text{Normal}(\theta, \sigma_\text{between}) \\
    \theta &\sim \text{Normal}(0, 1000) \\
    \sigma^2_\text{within}, \sigma^2_\text{between} &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`y_{i,j}` is the response for batch :math:`i` and sample :math:`j`.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: dyes.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
                  Mean        SD      Naive SE      MCSE         ESS
         theta 1527.0634   22.934017  0.2648192   0.37925906 3656.69166
    s2_between 2504.9567 2821.190905 32.5763066 297.11608227   90.15974
     s2_within 2870.5356  997.858373 11.5222760  57.91302702  296.88297

    Quantiles:
                  2.5%      25.0%     50.0%     75.0%        97.5%
         theta 1481.22789 1512.8876 1527.6533 1540.3759 1.57203607x103
    s2_between  152.11123  850.5886 1676.4907 3000.7898 1.28745381x104
     s2_within 1532.84786 2175.9060 2665.1653 3328.3947 5.54381144x103
