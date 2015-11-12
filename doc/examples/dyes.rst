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
    s2_between 3192.2614 4456.129102 51.45494673 487.44961486   83.57109
         theta 1526.7186   24.549671  0.28347518   0.37724897 3750.00000
     s2_within 2887.5853 1075.174260 12.41504296  76.89117959  195.52607
         mu[1] 1511.4798   20.819711  0.24040531   0.52158448 1593.30921
         mu[2] 1527.9087   20.344151  0.23491402   0.30199960 3750.00000
         mu[3] 1552.6742   21.293738  0.24587891   0.70276515  918.08605
         mu[4] 1506.6440   21.349176  0.24651905   0.61821290 1192.57616
         mu[5] 1578.6636   25.512471  0.29459264   1.29216105  389.82685
         mu[6] 1487.1934   24.693967  0.28514137   1.23710390  398.44592

    Quantiles:
                  2.5%      25.0%     50.0%     75.0%        97.5%
    s2_between  111.92351  815.9012 1651.4605 3269.5663 1.90261752x104
         theta 1475.08243 1513.5687 1527.1763 1540.0550 1.57444106x103
     s2_within 1566.61796 2160.3251 2654.9358 3324.8621 5.65161862x103
         mu[1] 1469.68693 1498.0985 1511.3084 1525.7344 1.55281296x103
         mu[2] 1486.15990 1514.8487 1527.6843 1541.2155 1.56770562x103
         mu[3] 1512.72912 1537.7046 1552.1976 1566.7188 1.59498301x103
         mu[4] 1463.91701 1492.3206 1506.9856 1521.3449 1.54854165x103
         mu[5] 1528.52480 1562.1831 1579.2515 1596.0167 1.62731017x103
         mu[6] 1440.27721 1470.7983 1486.1844 1502.9911 1.54464614x103
