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
               Mean           SD         Naive SE         MCSE          ESS
    alpha  2.0100584321 0.26156942610 0.003020343571 0.027100522461  93.157673
       s2  0.0690769709 0.04304237136 0.000497010494 0.001985478741 469.961085
     beta  0.3543443166 0.07160779229 0.000826855563 0.007122805926 101.069088
    gamma -0.0011250515 0.00034536546 0.000003987937 0.000025419899 184.590851

    Quantiles:
               2.5%         25.0%        50.0%          75.0%          97.5%
    alpha  1.5060295212  1.847111545  2.0062727893  2.16556036713  2.50547846044
       s2  0.0135820292  0.039788136  0.0598307554  0.08740879470  0.17747300708
     beta  0.2073237562  0.312747679  0.3582979574  0.40006921583  0.48789037325
    gamma -0.0017930379 -0.001347435 -0.0011325733 -0.00090960998 -0.00039271486
