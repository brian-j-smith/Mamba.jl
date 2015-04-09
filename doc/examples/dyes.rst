.. index:: Examples; Dyes: Variance Components Model

.. _example-Dyes:

Dyes: Variance Components Model
-------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Davies :cite:`davies:1967:SMR`, and Box and Tiao :cite:`box:1973:BIS` concerning batch-to-batch variation in yields from six batches and five samples of dyestuff.

Model
^^^^^

Yields are modelled as

.. math::

    y_{i,j} &\sim \text{Normal}(\mu_i, \sigma_\text{within}) \quad\quad i=1,\ldots,6; j=1,\ldots,5 \\
    \mu_i &\sim \text{Normal}(\theta, \sigma_\text{between}) \\
    \theta &\sim \text{Normal}(0, 1000) \\
    \sigma^2_\text{within}, \sigma^2_\text{between} &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`y_{i,j}` is the response for batch :math:`i` and sample :math:`j`.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: dyes.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
    4x6 Array{Any,2}:
     ""                "Mean"      "SD"    "Naive SE"     "MCSE"       "ESS"
     "s2_between"  2504.96     2821.19   32.5763       201.768      195.505
     "s2_within"   2870.54      997.858  11.5223        50.364      392.552
     "theta"       1527.06       22.934   0.264819       0.395415  3363.99

    Quantiles:
    4x6 Array{Any,2}:
     ""                "2.5%"      "25.0%"  .      "75.0%"       "97.5%"
     "s2_between"   152.111     850.589        3000.79      12874.5
     "s2_within"   1532.85     2175.91         3328.39       5543.81
     "theta"       1481.23     1512.89         1540.38       1572.04
