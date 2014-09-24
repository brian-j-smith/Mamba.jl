.. index:: Examples; Oxford: Smooth Fit to Log-Odds Ratios

.. _example-Oxford:

Oxford: Smooth Fit to Log-Odds Ratios
-------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Breslow and Clayton :cite:`breslow:1993:AIG` concerning the association between death from childhood cancer and maternal exposure to X-rays, for subjects partitioned into 120 age and birth-year strata. 

Model
^^^^^

Deaths are modelled as

.. math::

	r^0_i &\sim \text{Binomial}(n^0_i, p^0_i) \quad\quad i=1,\ldots,120 \\
	r^1_i &\sim \text{Binomial}(n^1_i, p^1_i) \\
	\operatorname{logit}(p^0_i) &= \mu_i \\
	\operatorname{logit}(p^1_i) &= \mu_i + \log(\psi_i) \\
	\log(\psi) &= \alpha + \beta_1 \text{year}_i + \beta_2 (\text{year}^2_i - 22) + b_i \\
	\mu_i &\sim \text{Normal}(0, 1000) \\
	b_i &\sim \text{Normal}(0, \sigma) \\
	\alpha, \beta_1, \beta_2 &\sim \text{Normal}(0, 1000) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`r^0_i` is the number of deaths among unexposed subjects in stratum :math:`i`, :math:`r^1_i` is the number among exposed subjects, and :math:`\text{year}_i` is the stratum-specific birth year (relative to 1954).
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: oxford.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia

	Iterations = 2502:12500
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 5000

	Empirical Posterior Estimates:
	5x6 Array{Any,2}:
	 ""         "Mean"      "SD"        "Naive SE"   "MCSE"          "ESS"
	 "alpha"   0.565785    0.0630051   0.000630051  0.00340073   1852.69
	 "s2"      0.026239    0.0307989   0.000307989  0.00260769   1181.08
	 "beta1"  -0.0433363   0.0161754   0.000161754  0.0011078    1460.14
	 "beta2"   0.00547712  0.00356757  3.56757e-5   0.000235842  1512.69

	Quantiles:
	5x6 Array{Any,2}:
	 ""         "2.5%"        "25.0%"      "50.0%"      "75.0%"      "97.5%"
	 "alpha"   0.443826      0.52388      0.567504     0.605143     0.695968
	 "s2"      0.000713442   0.00333527   0.0146737    0.0397133    0.118202
	 "beta1"  -0.0745152    -0.054318    -0.0434426   -0.0321216   -0.00992079
	 "beta2"  -0.0010499     0.00284892   0.00565004   0.00774736   0.0136309
