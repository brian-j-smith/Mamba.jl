.. index:: Examples; Seeds: Random Effect Logistic Regression

.. _example-Seeds:

Seeds: Random Effect Logistic Regression
----------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`, Crowder :cite:`crowder:1978:BBA`, and Breslow and Clayton :cite:`breslow:1993:AIG` concerning the proportion of seeds that germinated on each of 21 plates arranged according to a 2 by 2 factorial layout by seed and type of root extract.

Model
^^^^^

Germinations are modelled as

.. math::

    r_i &\sim \text{Binomial}(n_i, p_i) \quad\quad i=1,\ldots,21 \\
    \operatorname{logit}(p_i) &= \alpha_0 + \alpha_1 x_{1i} + \alpha_2 x_{2i} + \alpha_{12} x_{1i} x_{2i} + b_i \\
    b_i &\sim \text{Normal}(0, \sigma) \\
    \alpha_0, \alpha_1, \alpha_2, \alpha_{12} &\sim \text{Normal}(0, 1000) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`r_i` are the number of seeds, out of :math:`n_i`, that germinate on plate :math:`i`; and :math:`x_{1i}` and :math:`x_{2i}` are the seed type and root extract.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: seeds.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:12500
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 5000

    Empirical Posterior Estimates:
               Mean        SD       Naive SE       MCSE       ESS
     alpha2  1.3107281 0.26053104 0.0026053104 0.016391301 252.63418
         s2  0.0857053 0.09738014 0.0009738014 0.008167153 142.16720
     alpha0 -0.5561543 0.17595432 0.0017595432 0.011366011 239.65346
    alpha12 -0.7464409 0.43006756 0.0043006756 0.026796428 257.58440
     alpha1  0.0887002 0.26872879 0.0026872879 0.013771345 380.78132

    Quantiles:
                2.5%       25.0%      50.0%      75.0%      97.5%
     alpha2  0.80405938  1.1488818  1.3099477  1.4807632  1.8281561
         s2  0.00117398  0.0211176  0.0593763  0.1114008  0.3523464
     alpha0 -0.91491978 -0.6666323 -0.5512929 -0.4426242 -0.2224477
    alpha12 -1.54570414 -1.0275765 -0.7572503 -0.4914919  0.1702970
     alpha1 -0.42501648 -0.0936379  0.0943906  0.2600758  0.6235393
