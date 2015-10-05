.. index:: Examples; LSAT: Item Response

.. _example-LSAT:

LSAT: Item Response
-------------------

An example from OpenBUGS :cite:`openbugs:2014:ex` and Boch and Lieberman :cite:`boch:1970:FRM` concerning a 5-item multiple choice test (32 possible response patters) given to 1000 students.

Model
^^^^^

Item responses are modelled as

.. math::

    r_{i,j} &\sim \text{Bernoulli}(p_{i,j}) \quad\quad i=1,\ldots,1000; j=1,\ldots,5 \\
    \operatorname{logit}(p_{i,j}) &= \beta \theta_i - \alpha_j \\
    \theta_i &\sim \text{Normal}(0, 1) \\
    \alpha_j &\sim \text{Normal}(0, 100) \\
    \beta &\sim \text{Flat}(0, \infty),

where :math:`r_{i,j}` is an indicator for correct response by student :math:`i` to questions :math:`j`.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: lsat.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
            Mean        SD      Naive SE        MCSE        ESS
    a[1] -1.2631137 0.1050811 0.0012133719 0.00222282916 2234.791
    a[2]  0.4781765 0.0690242 0.0007970227 0.00126502181 2977.191
    a[3]  1.2401762 0.0679403 0.0007845071 0.00155219697 1915.850
    a[4]  0.1705201 0.0716647 0.0008275129 0.00146377676 2396.962
    a[5] -0.6257592 0.0858278 0.0009910540 0.00156546941 3005.846
    beta  0.7648299 0.0567109 0.0006548411 0.00277090485  418.880

    Quantiles:
            2.5%       25.0%     50.0%      75.0%      97.5%
    a[1] -1.4797642 -1.3325895 -1.262183 -1.1898184 -1.0647215
    a[2]  0.3493393  0.4300223  0.477481  0.5253468  0.6143681
    a[3]  1.1090783  1.1936967  1.240509  1.2865947  1.3699112
    a[4]  0.0285454  0.1221579  0.170171  0.2191623  0.3100302
    a[5] -0.7970834 -0.6828463 -0.625013 -0.5685627 -0.4590098
    beta  0.6618355  0.7249738  0.761457  0.8013506  0.8839082
