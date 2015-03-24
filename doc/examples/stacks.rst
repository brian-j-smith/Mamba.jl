.. index:: Examples; Stacks: Robust Regression

.. _example-Stacks:

Stacks: Robust Regression
-------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Brownlee :cite:`brownlee:1965:STM`, and Birkes and Dodge :cite:`birkes:1993:AMR` concerning 21 daily responses of stack loss, the amount of ammonia escaping, as a function of air flow, temperature, and acid concentration.

Model
^^^^^

Losses are modelled as

.. math::

	y_i &\sim \text{Laplace}(\mu_i, \sigma^2) \quad\quad i=1,\ldots,21 \\
	\mu_i &= \beta_0 + \beta_1 z_{1i} + \beta_2 z_{2i} + \beta_3 z_{3i} \\
	\beta_0, \beta_1, \beta_2, \beta_3 &\sim \text{Normal}(0, 1000) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`y_i` is the stack loss on day :math:`i`; and :math:`z_{1i}, z_{2i}, z_{3i}` are standardized predictors.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: stacks.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia

	Iterations = 2502:10000
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 3750

	Empirical Posterior Estimates:
	10x6 Array{Any,2}:
	 ""                "Mean"     "SD"      "Naive SE"   "MCSE"         "ESS"
	 "b[1]"           0.834448   0.129766  0.0014984    0.00227082  3265.53
	 "b[2]"           0.751506   0.33291   0.00384411   0.00562761  3499.49
	 "b[3]"          -0.117143   0.119602  0.00138105   0.00160727  5537.35
	 "outlier[1]"     0.0389333  0.193449  0.00223376   0.0028511   4603.71
	 "outlier[3]"     0.0558667  0.229679  0.00265211   0.00357717  4122.54
	 "outlier[4]"     0.3        0.458288  0.00529186   0.00861175  2832.01
	 "outlier[21]"    0.609867   0.487813  0.00563277   0.0110649   1943.61
	 "b0"           -38.7375     8.67973   0.100225     0.103293    7061.04
	 "sigma"          3.46651    0.85488   0.0098713    0.0275886    960.173

	Quantiles:
	10x6 Array{Any,2}:
	 ""                "2.5%"      "25.0%"     "50.0%"     "75.0%"      "97.5%"
	 "b[1]"           0.568641    0.754152    0.835908    0.91699      1.08951
	 "b[2]"           0.178021    0.525812    0.723047    0.94778      1.4807
	 "b[3]"          -0.36604    -0.189466   -0.112007   -0.0410014    0.111335
	 "outlier[1]"     0.0         0.0         0.0         0.0          1.0
	 "outlier[3]"     0.0         0.0         0.0         0.0          1.0
	 "outlier[4]"     0.0         0.0         0.0         1.0          1.0
	 "outlier[21]"    0.0         0.0         1.0         1.0          1.0
	 "b0"           -56.5036    -44.0468    -38.7457    -33.4175     -21.7374
	 "sigma"          2.17743     2.86055     3.33732     3.94188      5.50963
