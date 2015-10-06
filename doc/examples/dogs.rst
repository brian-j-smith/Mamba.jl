.. index:: Examples; Dogs: Loglinear Model for Binary Data

.. _example-Dogs:

Dogs: Loglinear Model for Binary Data
-------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`, Lindley and Smith :cite:`lindley:1972:BEL`, and Kalbfleisch :cite:`kalbfleisch:1985:PSI` concerning the Solomon-Wynne experiment on dogs.  In the experiment, 30 dogs were subjected to 25 trials.  On each trial, a barrier was raised, and an electric shock was administered 10 seconds later if the dog did not jump the barrier.

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
             Mean         SD        Naive SE        MCSE        ESS
    alpha -0.2441654 0.024064384 0.00027787158 0.00045249868 2828.2310
     beta -0.0788601 0.011812880 0.00013640339 0.00018435162 3750.0000
        B  0.9242336 0.010903201 0.00012589932 0.00017012187 3750.0000
        A  0.7835846 0.018821132 0.00021732772 0.00035466850 2816.0882

    Quantiles:
              2.5%        25.0%       50.0%       75.0%       97.5%
    alpha -0.29269315 -0.260167645 -0.24381478 -0.22784322 -0.19872975
     beta -0.10326776 -0.086633844 -0.07848501 -0.07046663 -0.05703080
        B  0.90188545  0.917012804  0.92451592  0.93195884  0.94456498
        A  0.74625110  0.770922334  0.78363276  0.79624909  0.81977141
