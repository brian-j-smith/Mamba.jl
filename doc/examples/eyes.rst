.. index:: Examples; Eyes: Normal Mixture Model

.. _example-Eyes:

Eyes: Normal Mixture Model
--------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Bowmaker :cite:`bowmaker:1985:TTT`, and Robert :cite:`robert:1994:MDI` concerning 48 peak sensitivity wavelength measurements taken on a set of monkey's eyes.

Model
^^^^^

Measurements are modelled as the mixture distribution

.. math::

    y_i &\sim \text{Normal}(\lambda_{T_i}, \sigma) \quad\quad i=1,\ldots,48 \\
    T_i &\sim \text{Categorical}(p, 1 - p) \\
    \lambda_1 &\sim \text{Normal}(0, 1000) \\
    \lambda_2 &= \lambda_1 + \theta \\
    \theta &\sim \text{Uniform}(0, 1000) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001) \\
    p &= \text{Uniform}(0, 1)

where :math:`y_i` is the measurement on monkey :math:`i`.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: eyes.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
                 Mean       SD      Naive SE      MCSE       ESS
    lambda[1] 536.81290 0.9883483 0.011412463 0.052323509 356.8012
    lambda[2] 548.75963 1.4660892 0.016928940 0.061653235 565.4693
            p   0.59646 0.0909798 0.001050544 0.003005307 916.4579
           s2  15.37891 7.1915331 0.083040671 0.406764080 312.5775

    Quantiles:
                2.5%    25.0%     50.0%     75.0%     97.5%
    lambda[1] 535.086 536.17470 536.75185 537.36318 538.9397
    lambda[2] 545.199 548.07943 548.86953 549.63026 551.1739
            p   0.413   0.54206   0.60080   0.65619   0.7602
           s2   8.637  11.40687  13.68101  16.78394  39.0676
