.. index:: Examples; Dogs: Loglinear Model for Binary Data

.. _example-Dogs:

Dogs: Loglinear Model for Binary Data
-------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Lindley and Smith :cite:`lindley:1972:BEL`, and Kalbfleisch :cite:`kalbfleisch:1985:PSI` concerning the Solomon-Wynne experiment on dogs.  In the experiment, 30 dogs were subjected to 25 trials.  On each trial, a barrier was raised, and an electric shock was administered 10 seconds later if the dog did not jump the barrier.

Model
^^^^^

Failures to jump the barriers in time are modelled as

.. math::

    y_{i,j} &= \text{Bernoulli}(\pi_{i,j}) \quad\quad i=1,\ldots,30; j=2,\ldots,25 \\
    \log(\pi_{i,j}) &= \alpha x_{i,j-1} + \beta (j - 1 - x_{i,j-1}) \\
    \alpha, \beta &\sim \text{Flat}(-\infty, -0.00001) \\

where :math:`y_{i,j} = 1` if dog :math:`i` fails to jump the barrier before the shock on trial :math:`j`, and 0 otherwise; :math:`x_{i,j-1}` is the number of successful jumps prior to trial :math:`j`; and :math:`\pi_{i,j}` is the probability of failure.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: dogs.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
             Mean        SD        Naive SE        MCSE        ESS
        A  0.783585 0.018821132 0.00021732772 0.00034707278 2940.6978
        B  0.924234 0.010903201 0.00012589932 0.00016800825 4211.5969
    alpha -0.244165 0.024064384 0.00027787158 0.00044330447 2946.7635
     beta -0.078860 0.011812880 0.00013640339 0.00018220980 4203.0847

    Quantiles:
             2.5%       25.0%      50.0%      75.0%      97.5%
        A  0.7462511  0.7709223  0.7836328  0.7962491  0.8197714
        B  0.9018854  0.9170128  0.9245159  0.9319588  0.9445650
    alpha -0.2926932 -0.2601676 -0.2438148 -0.2278432 -0.1987298
     beta -0.1032678 -0.0866338 -0.0784850 -0.0704666 -0.0570308
