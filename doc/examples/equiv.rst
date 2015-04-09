.. index:: Examples; Equiv: Bioequivalence in a Cross-Over Trial

.. _example-Equiv:

Equiv: Bioequivalence in a Cross-Over Trial
-------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning a two-treatment, cross-over trial with 10 subjects.

Model
^^^^^

Treatment responses are modelled as

.. math::

    y_{i,j} &\sim \text{Normal}(m_{i,j}, \sigma_1) \quad\quad i=1,\ldots,10; j=1,2 \\
    m_{i,j} &= \mu + (-1)^{T_{i,j} - 1} \phi / 2 + (-1)^{j-1} \pi / 2 + \delta_i \\
    \delta_i &\sim \text{Normal}(0, \sigma_2) \\
    \mu, \phi, \pi &\sim \text{Normal}(0, 1000) \\
    \sigma_1^2, \sigma_2^2 &\sim \text{InverseGamma}(0.001, 0.001) \\

where :math:`y_{i,j}` is the response for patient :math:`i` in period :math:`j`; and :math:`T_{i,j} = 1,2` is the treatment received.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: equiv.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:12500
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 5000

    Empirical Posterior Estimates:
             Mean        SD       Naive SE       MCSE        ESS
       pi -0.1788624 0.07963600 0.0007963600 0.0022537740 1248.5276
       mu  1.4432301 0.04735444 0.0004735444 0.0020359933  540.9645
    equiv  0.9835000 0.12739456 0.0012739456 0.0022978910 3073.5685
     s2_2  0.0187006 0.01466228 0.0001466228 0.0004799037  933.4580
     s2_1  0.0164033 0.01412396 0.0001412396 0.0005079176  773.2610
    theta  0.9837569 0.07936595 0.0007936595 0.0025041274 1004.5132
      phi -0.0195929 0.08005300 0.0008005300 0.0025432838  990.7533

    Quantiles:
             2.5%      25.0%     50.0%      75.0%      97.5%
       pi -0.332750 -0.2265245 -0.184088 -0.1285565 -0.0142326
       mu  1.353882  1.4111562  1.441424  1.4740774  1.5345210
    equiv  1.000000  1.0000000  1.000000  1.0000000  1.0000000
     s2_2  0.001568  0.0061461  0.016302  0.0273322  0.0530448
     s2_1  0.001035  0.0043785  0.013677  0.0242892  0.0507633
    theta  0.838904  0.9336518  0.974953  1.0320046  1.1567524
      phi -0.175659 -0.0686517 -0.025366  0.0315031  0.1456164
