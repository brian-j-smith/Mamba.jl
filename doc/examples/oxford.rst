.. index:: Examples; Oxford: Smooth Fit to Log-Odds Ratios

.. _example-Oxford:

Oxford: Smooth Fit to Log-Odds Ratios
-------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex` and Breslow and Clayton :cite:`breslow:1993:AIG` concerning the association between death from childhood cancer and maternal exposure to X-rays, for subjects partitioned into 120 age and birth-year strata.

Model
^^^^^

Deaths are modelled as

.. math::

    r^0_i &\sim \text{Binomial}(n^0_i, p^0_i) \quad\quad i=1,\ldots,120 \\
    r^1_i &\sim \text{Binomial}(n^1_i, p^1_i) \\
    \operatorname{logit}(p^0_i) &= \mu_i \\
    \operatorname{logit}(p^1_i) &= \mu_i + \log(\psi_i) \\
    \log(\psi) &= \alpha + \beta_1 \text{year}_i + \beta_2 (\text{year}^2_i - 22) + b_i \\
    \mu_i &\sim \text{Normal}(0, 1000) \\
    b_i &\sim \text{Normal}(0, \sigma) \\
    \alpha, \beta_1, \beta_2 &\sim \text{Normal}(0, 1000) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`r^0_i` is the number of deaths among unexposed subjects in stratum :math:`i`, :math:`r^1_i` is the number among exposed subjects, and :math:`\text{year}_i` is the stratum-specific birth year (relative to 1954).


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: oxford.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:12500
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 5000

    Empirical Posterior Estimates:
             Mean         SD        Naive SE       MCSE        ESS
    alpha  0.5657848 0.063005090 0.00063005090 0.0034007318 343.24680
       s2  0.0262390 0.030798915 0.00030798915 0.0026076857 139.49555
    beta1 -0.0433363 0.016175426 0.00016175426 0.0011077969 213.20196
    beta2  0.0054771 0.003567575 0.00003567575 0.0002358424 228.82453

    Quantiles:
              2.5%       25.0%      50.0%      75.0%       97.5%
    alpha  0.44382579  0.5238801  0.5675039  0.60514271  0.6959681
       s2  0.00071344  0.0033353  0.0146737  0.03971325  0.1182023
    beta1 -0.07451524 -0.0543180 -0.0434426 -0.03212161 -0.0099208
    beta2 -0.00104990  0.0028489  0.0056500  0.00774736  0.0136309
