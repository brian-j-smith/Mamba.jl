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
                  Mean         SD        Naive SE       MCSE        ESS
           s2  0.01822186 0.021121265 0.00024388736 0.0014150714 222.78358
            d -0.25563567 0.061841945 0.00071408927 0.0040205781 236.58613
    delta_new -0.25005767 0.150325282 0.00173580684 0.0050219145 896.03592

    Quantiles:
                   2.5%         25.0%         50.0%         75.0%       97.5%
           s2  0.0006855452  0.0041648765  0.0107615659  0.024442084  0.07735715
            d -0.3734122953 -0.2959169814 -0.2581848849 -0.218341380 -0.12842580
    delta_new -0.5385405488 -0.3279958446 -0.2557849252 -0.177588413  0.07986060
