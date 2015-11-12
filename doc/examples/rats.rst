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
                Mean        SD       Naive SE       MCSE       ESS
       s2_c  37.0199560 5.51731572 0.0637084743 0.214955115  658.8099
    mu_beta   6.1862265 0.10722363 0.0012381118 0.002028034 2795.3133
     alpha0 106.5745824 3.65326205 0.0421842366 0.055839206 3750.0000

    Quantiles:
               2.5%       25.0%      50.0%       75.0%       97.5%
       s2_c 27.6932918  33.095218  36.5856460  40.2952967  49.372348
    mu_beta  5.9784465   6.115354   6.1866336   6.2565411   6.401243
     alpha0 99.3655988 104.134835 106.5493857 109.0125662 113.782955
