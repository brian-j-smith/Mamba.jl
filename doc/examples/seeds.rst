.. index:: Examples; Seeds: Random Effect Logistic Regression

.. _example-Seeds:

Seeds: Random Effect Logistic Regression
----------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex`, Crowder :cite:`crowder:1978:BBA`, and Breslow and Clayton :cite:`breslow:1993:AIG` concerning the proportion of seeds that germinated on each of 21 plates arranged according to a 2 by 2 factorial layout by seed and type of root extract.

Model
^^^^^

Germinations are modelled as

.. math::

    r_i &\sim \text{Binomial}(n_i, p_i) \quad\quad i=1,\ldots,21 \\
    \operatorname{logit}(p_i) &= \alpha_0 + \alpha_1 x_{1i} + \alpha_2 x_{2i} + \alpha_{12} x_{1i} x_{2i} + b_i \\
    b_i &\sim \text{Normal}(0, \sigma) \\
    \alpha_0, \alpha_1, \alpha_2, \alpha_{12} &\sim \text{Normal}(0, 1000) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),

where :math:`r_i` are the number of seeds, out of :math:`n_i`, that germinate on plate :math:`i`; and :math:`x_{1i}` and :math:`x_{2i}` are the seed type and root extract.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: seeds.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:12500
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 5000

    Empirical Posterior Estimates:
                Mean         SD       Naive SE       MCSE        ESS
     alpha2  1.310728093 0.26053104 0.0026053104 0.0153996801 286.21707
     alpha1  0.088700176 0.26872879 0.0026872879 0.0128300598 438.70341
     alpha0 -0.556154341 0.17595432 0.0017595432 0.0101730837 299.15388
    alpha12 -0.746440855 0.43006756 0.0043006756 0.0251658152 292.04607
         s2  0.085705306 0.09738014 0.0009738014 0.0080848189 145.07755

    Quantiles:
                 2.5%         25.0%        50.0%       75.0%       97.5%
     alpha2  0.8040593795  1.148881827  1.309947687  1.48076318  1.82815608
     alpha1 -0.4250164771 -0.093637900  0.094390643  0.26007581  0.62353933
     alpha0 -0.9149197759 -0.666632319 -0.551292851 -0.44262420 -0.22244775
    alpha12 -1.5457041398 -1.027576522 -0.757250262 -0.49149187  0.17029702
         s2  0.0011739822  0.021117624  0.059376315  0.11140082  0.35234645
