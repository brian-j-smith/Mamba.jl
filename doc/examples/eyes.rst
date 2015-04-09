.. index:: Examples; Eyes: Normal Mixture Model

.. _example-Eyes:

Eyes: Normal Mixture Model
--------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Bowmaker :cite:`bowmaker:1985:TTT`, and Robert :cite:`robert:1994:MDI` concerning 48 peak sensitivity wavelength measurements taken on a set of monkey's eyes.

Model
^^^^^

Measurements are modelled as the mixture distribution

.. math::

    y_i &\sim \text{Normal}(\lambda_{T_i}, \sigma) \quad\quad i=1,\ldots,48 \\
    T_i &\sim \text{Categorical}(p, 1 - p) \\
    \lambda_1 &\sim \text{Normal}(0, 1000) \\
    \lambda_2 &= \lambda_1 + \theta \\
    \theta &\sim \text{Uniform}(0, 1000) \\
    \sigma^2 &\sim \text{InverseGamma}(0.001, 0.001) \\
    p &= \text{Uniform}(0, 1)

where :math:`y_i` is the measurement on monkey :math:`i`.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: eyes.jl
    :language: julia

Results
^^^^^^^

.. code-block:: julia

    Iterations = 2502:10000
    Thinning interval = 2
    Chains = 1,2
    Samples per chain = 3750

    Empirical Posterior Estimates:
    5x6 Array{Any,2}:
     ""              "Mean"    "SD"       "Naive SE"   "MCSE"        "ESS"
     "lambda[1]"  536.813     0.988348   0.0114125    0.0523235   356.801
     "lambda[2]"  548.76      1.46609    0.0169289    0.0616532   565.469
     "p"            0.596457  0.0909798  0.00105054   0.00300531  916.458
     "s2"          15.3789    7.19153    0.0830407    0.406764    312.578

    Quantiles:
    5x6 Array{Any,2}:
     ""              "2.5%"      "25.0%"     "50.0%"     "75.0%"     "97.5%"
     "lambda[1]"  535.086     536.175     536.752     537.363     538.94
     "lambda[2]"  545.199     548.079     548.87      549.63      551.174
     "p"            0.413162    0.542063    0.600803    0.656188    0.760195
     "s2"           8.63741    11.4069     13.681      16.7839     39.0676
