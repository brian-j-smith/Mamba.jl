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
                   Mean         SD        Naive SE       MCSE        ESS
        alpha0 -1.35617079 1.313240197 0.01174597740 0.1021568814 165.25442
      alpha_BT  0.24217000 0.190566444 0.00170447809 0.0163585849 135.70673
     alpha_Trt -0.75931393 0.397734236 0.00355744316 0.0337796459 138.63592
         s2_b1  0.24911885 0.073166731 0.00065442313 0.0044942066 265.04599
    alpha_Base  0.91104974 0.135354470 0.00121064718 0.0111438503 147.52807
     alpha_Age  0.45830900 0.394536219 0.00352883922 0.0310012419 161.96291
      alpha_V4 -0.09287934 0.083666872 0.00074833925 0.0051087325 268.21356
          s2_b  0.13523750 0.031819272 0.00028460022 0.0025401820 156.91007

    Quantiles:
                   2.5%        25.0%       50.0%        75.0%        97.5%
        alpha0 -4.16888778 -2.157932918 -1.26343143 -0.436226488  0.866195785
      alpha_BT -0.09025014  0.108102968  0.22656174  0.358354280  0.657804969
     alpha_Trt -1.63682212 -1.011390894 -0.75653998 -0.480870874 -0.016113397
         s2_b1  0.13817484  0.197135457  0.23761145  0.289655043  0.422805121
    alpha_Base  0.66318818  0.817700638  0.90268210  0.997417378  1.200619714
     alpha_Age -0.19666987  0.176356196  0.41608686  0.696647899  1.305075377
      alpha_V4 -0.25504531 -0.148157017 -0.09313603 -0.036681431  0.072099011
          s2_b  0.07155807  0.112591385  0.13626498  0.158032604  0.193715882
