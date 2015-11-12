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
              Mean          SD         Naive SE        MCSE         ESS
    beta2  0.005477119 0.0035675748 0.00003567575 0.00033192987 115.519023
    beta1 -0.043336269 0.0161754258 0.00016175426 0.00133361554 147.112695
    alpha  0.565784774 0.0630050896 0.00063005090 0.00468384860 180.944576
       s2  0.026238992 0.0307989154 0.00030798915 0.00302056007 103.967091

    Quantiles:
               2.5%         25.0%         50.0%         75.0%         97.5%
    beta2 -0.0010499046  0.0028489198  0.0056500394  0.0077473623  0.013630865
    beta1 -0.0745152363 -0.0543180318 -0.0434425931 -0.0321216097 -0.009920787
    alpha  0.4438257884  0.5238801187  0.5675039159  0.6051427125  0.695968063
       s2  0.0007134423  0.0033352655  0.0146737037  0.0397132522  0.118202258
