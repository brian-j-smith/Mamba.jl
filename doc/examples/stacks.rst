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
          sigma   3.502016751 0.85180493 0.0098357961 0.0260872478 1066.1634
     outlier[1]   0.038133333 0.19153087 0.0022116080 0.0028072835 3750.0000
     outlier[3]   0.052533333 0.22311481 0.0025763079 0.0033696758 3750.0000
     outlier[4]   0.280666667 0.44935488 0.0051887033 0.0085620856 2754.3543
    outlier[21]   0.592266667 0.49144589 0.0056747284 0.0110715500 1970.3100
             b0 -38.839445204 8.62188712 0.0995569770 0.1229160381 3750.0000
           b[1]   0.837923059 0.13227068 0.0015273303 0.0025240531 2746.1876
           b[2]   0.744640682 0.33529404 0.0038716421 0.0060623848 3058.8958
           b[3]  -0.116620758 0.11670888 0.0013476381 0.0018156465 3750.0000

    Quantiles:
                    2.5%         25.0%        50.0%        75.0%         97.5%
          sigma   2.18670097   2.88866448   3.37259353   3.985459541   5.43759582
     outlier[1]   0.00000000   0.00000000   0.00000000   0.000000000   1.00000000
     outlier[3]   0.00000000   0.00000000   0.00000000   0.000000000   1.00000000
     outlier[4]   0.00000000   0.00000000   0.00000000   1.000000000   1.00000000
    outlier[21]   0.00000000   0.00000000   1.00000000   1.000000000   1.00000000
             b0 -56.72768098 -43.91684537 -38.93457916 -33.448406668 -21.75823545
           b[1]   0.57760961   0.75420399   0.83808746   0.918517384   1.10667879
           b[2]   0.14374739   0.52549311   0.71970820   0.945807004   1.48176196
           b[3]  -0.35969464  -0.18763632  -0.11239930  -0.041941611   0.10131938
