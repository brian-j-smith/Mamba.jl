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
               Mean          SD        Naive SE        MCSE         ESS
      phi -0.0065694145 0.085709855 0.00085709855 0.00237984045 1297.07723
    theta  0.9971035189 0.085511488 0.00085511488 0.00236997451 1301.85188
     s2_1  0.0210195763 0.014770588 0.00014770588 0.00064784959  519.81308
    equiv  0.9787000000 0.144389732 0.00144389732 0.00248492424 3376.34112
       pi -0.1798858210 0.088104904 0.00088104904 0.00257908779 1166.99209
       mu  1.4366641956 0.042648958 0.00042648958 0.00103397482 1701.36279
     s2_2  0.0144323595 0.013757849 0.00013757849 0.00069294275  394.19068

    Quantiles:
               2.5%         25.0%        50.0%        75.0%        97.5%
      phi -0.1766257876 -0.062578911 -0.005911288  0.050530106  0.1591732568
    theta  0.8380933508  0.939338936  0.994106149  1.051828530  1.1725410798
     s2_1  0.0025352373  0.009584715  0.018779792  0.028801243  0.0555316591
    equiv  1.0000000000  1.000000000  1.000000000  1.000000000  1.0000000000
       pi -0.3478283178 -0.238173615 -0.182191849 -0.124090215 -0.0009208467
       mu  1.3512670514  1.409186605  1.436922479  1.464276491  1.5202253288
     s2_2  0.0010675857  0.003651001  0.010235890  0.021460364  0.0490664434
