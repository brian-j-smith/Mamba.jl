.. index:: Examples; Stacks: Robust Regression

.. _example-Stacks:

Stacks: Robust Regression
-------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`, Brownlee :cite:`brownlee:1965:STM`, and Birkes and Dodge :cite:`birkes:1993:AMR` concerning 21 daily responses of stack loss, the amount of ammonia escaping, as a function of air flow, temperature, and acid concentration.

Model
^^^^^

Losses are modelled as

.. math::

    y_i &\sim \text{Laplace}(\mu_i, \sigma^2) \quad\quad i=1,\ldots,21 \\
    \mu_i &= \beta_0 + \beta_1 z_{1i} + \beta_2 z_{2i} + \beta_3 z_{3i} \\
    \beta_0, \beta_1, \beta_2, \beta_3 &\sim \text{Normal}(0, 1000) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`y_i` is the stack loss on day :math:`i`; and :math:`z_{1i}, z_{2i}, z_{3i}` are standardized predictors.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: stacks.jl
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
           b[1]   0.836863707 0.13085145 0.0015109423 0.0027601754 2247.4171
           b[2]   0.744454449 0.33480007 0.0038659382 0.0065756939 2592.3158
           b[3]  -0.116648437 0.12214077 0.0014103601 0.0015143922 3750.0000
             b0 -38.776564595 8.81860433 0.1018284717 0.0979006137 3750.0000
          sigma   3.487643717 0.87610847 0.0101164292 0.0279025494  985.8889
     outlier[1]   0.042666667 0.20211796 0.0023338572 0.0029490162 3750.0000
     outlier[3]   0.054800000 0.22760463 0.0026281519 0.0034398827 3750.0000
     outlier[4]   0.298000000 0.45740999 0.0052817156 0.0089200654 2629.5123
    outlier[21]   0.606400000 0.48858046 0.0056416412 0.0113877443 1840.7583

    Quantiles:
                    2.5%         25.0%        50.0%         75.0%         97.5%
           b[1]   0.57218621   0.75741345   0.834874964   0.918345319   1.101502854
           b[2]   0.16177144   0.52291878   0.714951465   0.933171533   1.476258382
           b[3]  -0.36401372  -0.19028697  -0.113463801  -0.036994963   0.118538277
             b0 -56.70056875 -44.11785905 -38.698338454 -33.409149788 -21.453323631
          sigma   2.17947513   2.86899865   3.348631697   3.953033535   5.592773118
     outlier[1]   0.00000000   0.00000000   0.000000000   0.000000000   1.000000000
     outlier[3]   0.00000000   0.00000000   0.000000000   0.000000000   1.000000000
     outlier[4]   0.00000000   0.00000000   0.000000000   1.000000000   1.000000000
    outlier[21]   0.00000000   0.00000000   1.000000000   1.000000000   1.000000000
