.. index:: Examples; Rats: A Normal Hierarchical Model

.. _example-Rats:

Rats: A Normal Hierarchical Model
---------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex` and section 6 of Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning 30 rats whose weights were measured at each of five consecutive weeks.

Model
^^^^^

Weights are modeled as

.. math::

    y_{i,j} &\sim \text{Normal}\left(\alpha_i + \beta_i (x_j - \bar{x}), \sigma_c\right) \quad\quad i=1,\ldots,30; j=1,\ldots,5 \\
    \alpha_i &\sim \text{Normal}(\mu_\alpha, \sigma_\alpha) \\
    \beta_i &\sim \text{Normal}(\mu_\beta, \sigma_\beta) \\
    \mu_\alpha, \mu_\beta &\sim \text{Normal}(0, 1000) \\
    \sigma^2_\alpha, \sigma^2_\beta, \sigma^2_c &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`y_{i,j}` is repeated weight measurement :math:`j` on rat :math:`i`, and :math:`x_j` is the day on which the measurement was taken.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: rats.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
                Mean         SD       Naive SE       MCSE        ESS
       s2_c  37.2543133 6.026634572 0.0695895819 0.2337982327  664.4576
    mu_beta   6.1830663 0.108042927 0.0012475723 0.0017921615 3634.4474
     alpha0 106.6259925 3.459210115 0.0399435178 0.0526804390 3750.0000

    Quantiles:
               2.5%      25.0%       50.0%       75.0%       97.5%
       s2_c 27.778388  33.0906026  36.4630047  40.5538472  51.5713716
    mu_beta  5.969850   6.1110307   6.1836454   6.2538831   6.3964953
     alpha0 99.815707 104.3369878 106.6105679 108.9124224 113.5045347
