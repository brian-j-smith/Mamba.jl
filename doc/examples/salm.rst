.. index:: Examples; Salm: Extra-Poisson Variation in a Dose-Response Study

.. _example-Salm:

Salm: Extra-Poisson Variation in a Dose-Response Study
------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Breslow :cite:`breslow:1984:EPV` concerning mutagenicity assay data on salmonella in three plates exposed to six doses of quinoline.

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
    5x6 Array{Any,2}:
     ""         "Mean"      "SD"         "Naive SE"   "MCSE"        "ESS"
     "gamma"  -0.00112505  0.000345365  3.98794e-6   2.15796e-5  256.135
     "alpha"   2.01006     0.261569     0.00302034   0.0210699   154.116
     "s2"      0.069077    0.0430424    0.00049701   0.0019298   497.47
     "beta"    0.354344    0.0716078    0.000826856  0.00564423  160.958

    Quantiles:
    5x6 Array{Any,2}:
     ""         "2.5%"       "25.0%"      "50.0%"      "75.0%"      "97.5%"
     "gamma"  -0.00179304  -0.00134744  -0.00113257  -0.00090961  -0.000392715
     "alpha"   1.50603      1.84711      2.00627      2.16556      2.50548
     "s2"      0.013582     0.0397881    0.0598308    0.0874088    0.177473
     "beta"    0.207324     0.312748     0.358298     0.400069     0.48789
