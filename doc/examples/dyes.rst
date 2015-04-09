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
                  Mean       SD     Naive SE    MCSE       ESS
    s2_between 2504.9567 2821.1909 32.576307 201.76845  195.50525
     s2_within 2870.5356  997.8584 11.522276  50.36398  392.55247
         theta 1527.0634   22.9340  0.264819   0.39541 3363.98979

    Quantiles:
                   2.5%      25.0%     50.0%    75.0%       97.5%
    s2_between  152.111234  850.5886 1676.491 3000.7898 1.287454x10^4
     s2_within 1532.847859 2175.9060 2665.165 3328.3947  5.54381x10^3
         theta 1481.227885 1512.8876 1527.653 1540.3759  1.57204x10^3
