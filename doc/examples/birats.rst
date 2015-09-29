.. index:: Examples; Birats: A Bivariate Normal Hierarchical Model

.. _example-Birats:

Birats: A Bivariate Normal Hierarchical Model
---------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex` and section 6 of Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning 30 rats whose weights were measured at each of five consecutive weeks.

Model
^^^^^

Weights are modeled as

.. math::

    Y_{i,j} &\sim \text{Normal}\left(\mu_{i,j}, \sigma_C\right) \quad\quad i=1,\ldots,30; j=1,\ldots,5 \\
    \mu_{i,j} &= \beta_{1,i} + \beta_{2,i} x_j \\
    \bm{\beta}_i &\sim \text{MvNormal}\left(\bm{\mu}_\beta, \bm{\Sigma}\right) \\
    \bm{\mu}_\beta &\sim \text{MvNormal}\left(\begin{bmatrix}0 \\ 0\end{bmatrix}, \begin{bmatrix}1e6 & 0 \\ 0 & 1e6\end{bmatrix}\right) \\
    \bm{\Sigma} &\sim \text{InverseWishart}\left(2, \begin{bmatrix}200 & 0 \\ 0 & 0.2 \end{bmatrix}\right) \\
    \sigma^2_C &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`y_{i,j}` is repeated weight measurement :math:`j` on rat :math:`i`, and :math:`x_j` is the day on which the measurement was taken.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: birats.jl
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
        sigmaC   6.165580 0.49032603 0.0056617973 0.026607160 339.60404
    mu_beta[1] 106.550196 2.37595237 0.0274351348 0.104141627 520.50720
    mu_beta[2]   6.190639 0.10438935 0.0012053844 0.004162153 629.03710

    Quantiles:
                   2.5%       25.0%      50.0%       75.0%       97.5%
        sigmaC   5.2974579   5.816391   6.1349725   6.4722057   7.2352780
    mu_beta[1] 101.9060895 104.963676 106.5475350 108.1248562 111.2112533
    mu_beta[2]   5.9858728   6.117944   6.1910308   6.2610129   6.3954136

                PSRF 97.5%
        sigmaC 1.003 1.017
    mu_beta[1] 1.003 1.008
    mu_beta[2] 1.003 1.017
