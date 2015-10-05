.. index:: Examples; Epilepsy: Repeated Measures on Poisson Counts

.. _example-Epilepsy:

Epilepsy: Repeated Measures on Poisson Counts
---------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`, Thall and Vail :cite:`thall:1990:SCM` Breslow and Clayton :cite:`breslow:1993:AIG` concerning the effects of treatment, baseline seizure counts, and age on follow-up seizure counts at four visits in 59 patients.

Model
^^^^^

Counts are modelled as

.. math::

    y_{i,j} &\sim \text{Poisson}(\mu_{i,j}) \quad\quad i=1,\ldots,59; j=1,\ldots,4 \\
    \log(\mu_{i,j}) &= \alpha_0 + \alpha_\text{Base} \log(\text{Base}_i / 4) +
      \alpha_\text{Trt} \text{Trt}_i + \alpha_\text{BT} \text{Trt}_i \log(\text{Base}_i / 4) + \\
      & \quad\quad \alpha_\text{Age} \log(\text{Age}_i) + \alpha_\text{V4} \text{V}_4 + \text{b1}_i +
        \text{b}_{i,j} \\
    \text{b1}_i &\sim \text{Normal}(0, \sigma_\text{b1}) \\
    \text{b}_{i,j} &\sim \text{Normal}(0, \sigma_\text{b}) \\
    \alpha_* &\sim \text{Normal}(0, 100) \\
    \sigma^2_\text{b1}, \sigma^2_\text{b} &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`y_{ij}` are the counts on patient :math:`i` at visit :math:`j`, :math:`\text{Trt}` is a treatment indicator, :math:`\text{Base}` is baseline seizure counts, :math:`\text{Age}` is age in years, and :math:`\text{V}_4` is an indicator for the fourth visit.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: epil.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:15000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 6250

    Empirical Posterior Estimates:
                  Mean        SD      Naive SE       MCSE       ESS
     alpha_Age  0.4583090 0.3945362 0.0035288392 0.020336364 376.38053
        alpha0 -1.3561708 1.3132402 0.0117459774 0.072100024 331.75503
      alpha_BT  0.2421700 0.1905664 0.0017044781 0.010759318 313.70641
    alpha_Base  0.9110497 0.1353545 0.0012106472 0.007208435 352.58447
          s2_b  0.1352375 0.0318193 0.0002846002 0.001551352 420.68781
     alpha_Trt -0.7593139 0.3977342 0.0035574432 0.023478808 286.96826
         s2_b1  0.2491188 0.0731667 0.0006544231 0.002900632 636.27088
      alpha_V4 -0.0928793 0.0836669 0.0007483393 0.003604194 538.87837

    Quantiles:
                  2.5%      25.0%      50.0%      75.0%      97.5%
     alpha_Age -0.1966699  0.176356  0.4160869  0.6966479  1.3050754
        alpha0 -4.1688878 -2.157933 -1.2634314 -0.4362265  0.8661958
      alpha_BT -0.0902501  0.108103  0.2265617  0.3583543  0.6578050
    alpha_Base  0.6631882  0.817701  0.9026821  0.9974174  1.2006197
          s2_b  0.0715581  0.112591  0.1362650  0.1580326  0.1937159
     alpha_Trt -1.6368221 -1.011391 -0.7565400 -0.4808709 -0.0161134
         s2_b1  0.1381748  0.197135  0.2376114  0.2896550  0.4228051
      alpha_V4 -0.2550453 -0.148157 -0.0931360 -0.0366814  0.0720990
