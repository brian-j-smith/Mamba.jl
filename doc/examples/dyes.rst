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
                  Mean        SD       Naive SE      MCSE         ESS
    s2_between 2504.9567 2821.190905 32.57630657 297.11608227   90.15974
         theta 1527.0634   22.934017  0.26481922   0.37925906 3656.69166
     s2_within 2870.5356  997.858373 11.52227601  57.91302702  296.88297
         mu[1] 1511.6652   20.705097  0.23908186   0.51689818 1604.51773
         mu[2] 1527.5936   19.998432  0.23092200   0.38578338 2687.23016
         mu[3] 1553.3854   22.326708  0.25780662   0.70552900 1001.42781
         mu[4] 1506.8091   21.622864  0.24967933   0.63046700 1176.25563
         mu[5] 1578.9083   24.924672  0.28780532   1.14046097  477.63710
         mu[6] 1487.0292   23.837722  0.27525431   0.99950096  568.80457

    Quantiles:
                  2.5%      25.0%     50.0%     75.0%        97.5%
    s2_between  152.11123  850.5886 1676.4907 3000.7898 1.28745381x104
         theta 1481.22789 1512.8876 1527.6533 1540.3759 1.57203607x103
     s2_within 1532.84786 2175.9060 2665.1653 3328.3947 5.54381144x103
         mu[1] 1470.52305 1498.3758 1511.3962 1525.6220 1.55003711x103
         mu[2] 1487.83858 1514.8982 1527.6126 1540.9084 1.56536498x103
         mu[3] 1507.55073 1539.3348 1553.7100 1567.9356 1.59696999x103
         mu[4] 1463.65177 1492.3893 1507.2429 1520.9127 1.54846053x103
         mu[5] 1528.07326 1562.5512 1579.6135 1596.0115 1.62705081x103
         mu[6] 1441.88604 1470.5346 1486.2646 1503.0577 1.53582165x103
