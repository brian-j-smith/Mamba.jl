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
     s2_2  0.018671728 0.013978201 0.00013978201 0.00059050068  560.35305
     s2_1  0.016001577 0.013611413 0.00013611413 0.00059059236  531.16660
       pi -0.180867259 0.085194250 0.00085194250 0.00310291679  753.84196
      phi -0.003251984 0.081750665 0.00081750665 0.00259822564  989.98556
    theta  1.000096211 0.082212440 0.00082212440 0.00260223583  998.11773
    equiv  0.978300000 0.145709413 0.00145709413 0.00312065455 2180.13720
       mu  1.437186076 0.040644361 0.00040644361 0.00095360517 1816.61716

    Quantiles:
               2.5%         25.0%         50.0%         75.0%        97.5%
     s2_2  0.0017608381  0.0074524601  0.0164784368  0.026222893  0.052165313
     s2_1  0.0012768915  0.0048312056  0.0128296755  0.023427479  0.049757645
       pi -0.3530867687 -0.2345023028 -0.1807673264 -0.126533622 -0.008848990
      phi -0.1618313902 -0.0548132403 -0.0038903765  0.046482376  0.164021207
    theta  0.8505846103  0.9466619298  0.9961171813  1.047579617  1.178239301
    equiv  1.0000000000  1.0000000000  1.0000000000  1.000000000  1.000000000
       mu  1.3586651875  1.4099902797  1.4360144213  1.463952605  1.518733806
