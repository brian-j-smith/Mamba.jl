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
                   Mean       SD      Naive SE     MCSE       ESS
           b[1]   0.83445 0.1297656 0.001498404 0.00227082 3265.534
           b[2]   0.75151 0.3329099 0.003844112 0.00562761 3499.489
           b[3]  -0.11714 0.1196021 0.001381046 0.00160727 5537.350
     outlier[1]   0.03893 0.1934490 0.002233757 0.00285110 4603.713
     outlier[3]   0.05587 0.2296794 0.002652109 0.00357717 4122.536
     outlier[4]   0.30000 0.4582881 0.005291855 0.00861175 2832.010
    outlier[21]   0.60987 0.4878125 0.005632774 0.01106491 1943.615
             b0 -38.73747 8.6797300 0.100224890 0.10329322 7061.042
          sigma   3.46651 0.8548797 0.009871301 0.02758864  960.173

    Quantiles:
                   2.5%       25.0%      50.0%      75.0%      97.5%
           b[1]   0.568641   0.754152   0.835908   0.916990   1.089509
           b[2]   0.178021   0.525812   0.723047   0.947780   1.480699
           b[3]  -0.366040  -0.189466  -0.112007  -0.041001   0.111335
     outlier[1]   0.000000   0.000000   0.000000   0.000000   1.000000
     outlier[3]   0.000000   0.000000   0.000000   0.000000   1.000000
     outlier[4]   0.000000   0.000000   0.000000   1.000000   1.000000
    outlier[21]   0.000000   0.000000   1.000000   1.000000   1.000000
             b0 -56.503611 -44.046775 -38.745722 -33.417515 -21.737427
          sigma   2.177428   2.860547   3.337315   3.941882   5.509634
