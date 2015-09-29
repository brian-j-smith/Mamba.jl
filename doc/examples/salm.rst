.. index:: Examples; Salm: Extra-Poisson Variation in a Dose-Response Study

.. _example-Salm:

Salm: Extra-Poisson Variation in a Dose-Response Study
------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex` and Breslow :cite:`breslow:1984:EPV` concerning mutagenicity assay data on salmonella in three plates exposed to six doses of quinoline.

Model
^^^^^

Number of revertant colonies of salmonella are modelled as

.. math::

    y_{i,j} &\sim \text{Poisson}(\mu_{i,j}) \quad\quad i=1,\ldots,3; j=1,\ldots,6 \\
    \log(\mu_{i,j}) &= \alpha + \beta \log(x_j + 10) + \gamma x_j + \lambda_{i,j} \\
    \alpha, \beta, \gamma &\sim \text{Normal}(0, 1000) \\
    \lambda_{i,j} &\sim \text{Normal}(0, \sigma) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`y_i` is the number of colonies in plate :math:`i` and dose :math:`j`.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: salm.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
             Mean        SD       Naive SE     MCSE       ESS
    gamma -0.0011251 0.00034537 0.000003988 0.00002158 256.1352
    alpha  2.0100584 0.26156943 0.003020344 0.02106994 154.1158
       s2  0.0690770 0.04304237 0.000497010 0.00192980 497.4697
     beta  0.3543443 0.07160779 0.000826856 0.00564423 160.9577

    Quantiles:
             2.5%       25.0%     50.0%      75.0%      97.5%
    gamma -0.0017930 -0.0013474 -0.001133 -0.0009096 -0.0003927
    alpha  1.5060295  1.8471115  2.006273  2.1655604  2.5054785
       s2  0.0135820  0.0397881  0.059831  0.0874088  0.1774730
     beta  0.2073238  0.3127477  0.358298  0.4000692  0.4878904
