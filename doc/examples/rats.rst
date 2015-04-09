.. index:: Examples; Rats: A Normal Hierarchical Model

.. _example-Rats:

Rats: A Normal Hierarchical Model
---------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and section 6 of Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning 30 rats whose weights were measured at each of five consecutive weeks.

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
               Mean       SD      Naive SE     MCSE       ESS
    mu_beta   6.18623 0.1072236 0.001238112 0.00215008 2486.9858
     alpha0 106.57458 3.6532621 0.042184237 0.05738100 4053.4550
       s2_c  37.01996 5.5173157 0.063708474 0.18973116  845.6260

    Quantiles:
              2.5%     25.0%     50.0%     75.0%     97.5%
    mu_beta  5.97845   6.11535   6.18663   6.25654   6.40124
     alpha0 99.36560 104.13484 106.54939 109.01257 113.78295
       s2_c 27.69329  33.09522  36.58565  40.29530  49.37235
