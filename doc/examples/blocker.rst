.. index:: Examples; Blocker: Random Effects Meta-Analysis of Clinical Trials

.. _example-Blocker:

Blocker: Random Effects Meta-Analysis of Clinical Trials
--------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex` and Carlin :cite:`carlin:1992:ME2` concerning a meta-analysis of 22 clinical trials to prevent mortality after myocardial infarction.

Model
^^^^^

Events are modelled as

.. math::

    r^c_i &\sim \text{Binomial}(n^c_i, p^c_i) \quad\quad i=1,\ldots,22 \\
    r^t_i &\sim \text{Binomial}(n^t_i, p^t_i) \\
    \operatorname{logit}(p^c_i) &= \mu_i \\
    \operatorname{logit}(p^t_i) &= \mu_i + \delta_i \\
    \mu_i &\sim \text{Normal}(0, 1000) \\
    \delta_i &\sim \text{Normal}(d, \sigma) \\
    d &\sim \text{Normal}(0, 1000) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`r^c_i` is the number of control group events, out of :math:`n^c_i`, in study :math:`i`; and :math:`r^t_i` is the number of treatment group events.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: blocker.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
                  Mean        SD       Naive SE       MCSE         ESS
    delta_new -0.25005767 0.15032528 0.0017358068 0.0042970294 1223.84778
           s2  0.01822186 0.02112127 0.0002438874 0.0012497271  285.63372
            d -0.25563567 0.06184194 0.0007140893 0.0030457284  412.27207

    Quantiles:
                  2.5%       25.0%       50.0%       75.0%       97.5%
    delta_new -0.53854055 -0.32799584 -0.25578493 -0.17758841  0.07986060
           s2  0.00068555  0.00416488  0.01076157  0.02444208  0.07735715
            d -0.37341230 -0.29591698 -0.25818488 -0.21834138 -0.12842580
