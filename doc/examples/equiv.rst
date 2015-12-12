.. index:: Examples; Equiv: Bioequivalence in a Cross-Over Trial

.. _example-Equiv:

Equiv: Bioequivalence in a Cross-Over Trial
-------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex` and Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning a two-treatment, cross-over trial with 10 subjects.

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
               Mean          SD        Naive SE       MCSE         ESS
     s2_2  0.0173121833 0.014549568 0.00014549568 0.0007329722  394.02626
     s2_1  0.0184397014 0.013837972 0.00013837972 0.0005689492  591.55873
       pi -0.1874240524 0.086420302 0.00086420302 0.0032257037  717.76558
      phi -0.0035569545 0.087590520 0.00087590520 0.0035141650  621.25503
    theta  1.0002921934 0.088250458 0.00088250458 0.0036227671  593.40761
    equiv  0.9751000000 0.155828169 0.00155828169 0.0036666529 1806.14385
       mu  1.4387396416 0.042269208 0.00042269208 0.0013735876  946.96847

    Quantiles:
               2.5%         25.0%         50.0%        75.0%        97.5%
     s2_2  0.0016375061  0.0056159514  0.013968228  0.024613730  0.053154674
     s2_1  0.0017289114  0.0074958338  0.015849718  0.025963832  0.051967161
       pi -0.3579631753 -0.2432161807 -0.187946319 -0.130454165 -0.014910965
      phi -0.1723722017 -0.0623573600 -0.005681830  0.053144647  0.172913654
    theta  0.8416658455  0.9395470702  0.994334281  1.054582177  1.188763456
    equiv  1.0000000000  1.0000000000  1.000000000  1.000000000  1.000000000
       mu  1.3552569594  1.4110400018  1.438593809  1.466525521  1.519643109
