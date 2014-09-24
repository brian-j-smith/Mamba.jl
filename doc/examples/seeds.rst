.. index:: Examples; Seeds: Random Effect Logistic Regression

.. _example-Seeds:

Seeds: Random Effect Logistic Regression
----------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Crowder :cite:`crowder:1978:BBA`, and Breslow and Clayton :cite:`breslow:1993:AIG` concerning the proportion of seeds that germinated on each of 21 plates arranged according to a 2 by 2 factorial layout by seed and type of root extract.

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
	6x6 Array{Any,2}:
	 ""           "Mean"     "SD"       "Naive SE"   "MCSE"         "ESS"
	 "alpha2"    1.31073    0.260531   0.00260531   0.0163913   1589.45
	 "s2"        0.0857053  0.0973801  0.000973801  0.00816715  1192.34
	 "alpha0"   -0.556154   0.175954   0.00175954   0.011366    1548.07
	 "alpha12"  -0.746441   0.430068   0.00430068   0.0267964   1604.94
	 "alpha1"    0.0887002  0.268729   0.00268729   0.0137713   1951.36

	Quantiles:
	6x6 Array{Any,2}:
	 ""           "2.5%"       "25.0%"     "50.0%"     "75.0%"    "97.5%"
	 "alpha2"    0.804059     1.14888     1.30995     1.48076    1.82816
	 "s2"        0.00117398   0.0211176   0.0593763   0.111401   0.352346
	 "alpha0"   -0.91492     -0.666632   -0.551293   -0.442624  -0.222448
	 "alpha12"  -1.5457      -1.02758    -0.75725    -0.491492   0.170297
	 "alpha1"   -0.425016    -0.0936379   0.0943906   0.260076   0.623539
