.. index:: Examples; Eyes: Normal Mixture Model

.. _example-Eyes:

Eyes: Normal Mixture Model
--------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`, Bowmaker :cite:`bowmaker:1985:TTT`, and Robert :cite:`robert:1994:MDI` concerning 48 peak sensitivity wavelength measurements taken on a set of monkey's eyes.

Model
^^^^^

Measurements are modelled as the mixture distribution

.. math::

    y_i &\sim \text{Normal}(\lambda_{T_i}, \sigma) \quad\quad i=1,\ldots,48 \\
    T_i &\sim \text{Categorical}(P) \\
    \lambda_1 &\sim \text{Normal}(0, 1000) \\
    \lambda_2 &= \lambda_1 + \theta \\
    \theta &\sim \text{Uniform}(0, 1000) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001) \\
    P &= \text{Dirichlet}(1, 1)

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
                  Mean         SD       Naive SE       MCSE         ESS
           s2  14.45234459 4.96689604 0.0573527753 0.1853970211  717.73584
         P[1]   0.60357102 0.08379221 0.0009675491 0.0013077559 3750.00000
         P[2]   0.39642898 0.08379221 0.0009675491 0.0013077559 3750.00000
    lambda[1] 536.75250337 0.88484533 0.0102173137 0.0304232594  845.90821
    lambda[2] 548.98693469 1.18938418 0.0137338256 0.0625489045  361.58042

    Quantiles:
                  2.5%         25.0%        50.0%        75.0%        97.5%
           s2   8.55745263  11.37700782  13.39622105  16.26405571  27.08863166
         P[1]   0.43576584   0.54842468   0.60395494   0.66066383   0.76552882
         P[2]   0.23447118   0.33933617   0.39604506   0.45157532   0.56423416
    lambda[1] 535.08951002 536.16891895 536.73114152 537.30710397 538.61179506
    lambda[2] 546.57859195 548.24188698 548.97732892 549.74377421 551.38780148
