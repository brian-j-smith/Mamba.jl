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
                   Mean         SD       Naive SE       MCSE       ESS
        sigmaC   6.1431758 0.460583341 0.0053183583 0.021005830 480.76933
    mu_beta[1] 106.7046188 2.258246468 0.0260759841 0.081338282 770.81949
    mu_beta[2]   6.1804557 0.104047928 0.0012014420 0.004102793 643.14317

    Quantiles:
                   2.5%       25.0%       50.0%       75.0%      97.5%
        sigmaC   5.3167022   5.8123935   6.1206859   6.445419   7.0971249
    mu_beta[1] 102.3595659 105.2252185 106.6914834 108.164852 111.2001520
    mu_beta[2]   5.9720904   6.1130035   6.1817455   6.248025   6.3848798
