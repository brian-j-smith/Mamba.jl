.. index:: Examples; Rats: A Normal Hierarchical Model

.. _example-Rats:

Rats: A Normal Hierarchical Model
---------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and section 6 of Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning 30 rats whose weights were measured at each of five consecutive weeks.

Model
^^^^^

Weights are modeled as

.. math::

	y_{i,j} &\sim \text{Normal}\left(\alpha_i + \beta_i (x_j - \bar{x}), \sigma_c\right) \quad\quad i=1,\ldots,30; j=1,\ldots,5 \\
	\alpha_i &\sim \text{Normal}(\mu_\alpha, \sigma_\alpha) \\
	\beta_i &\sim \text{Normal}(\mu_\beta, \sigma_\beta) \\
	\mu_\alpha, \mu_\beta &\sim \text{Normal}(0, 1000) \\
	\sigma^2_\alpha, \sigma^2_\beta, \sigma^2_c &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`y_{i,j}` is repeated weight measurement :math:`j` on rat :math:`i`, and :math:`x_j` is the day on which the measurement was taken.
	
Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: rats.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia

	Iterations = 2502:10000
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 3750

	Empirical Posterior Estimates:
	4x6 Array{Any,2}:
	 ""            "Mean"   "SD"      "Naive SE"   "MCSE"         "ESS"
	 "mu_beta"    6.18623  0.107224  0.00123811   0.00215008  2486.99
	 "alpha0"   106.575    3.65326   0.0421842    0.057381    4053.45
	 "s2_c"      37.02     5.51732   0.0637085    0.189731     845.626

	Quantiles:
	4x6 Array{Any,2}:
	 ""           "2.5%"     "25.0%"     "50.0%"     "75.0%"     "97.5%"
	 "mu_beta"   5.97845    6.11535     6.18663     6.25654     6.40124
	 "alpha0"   99.3656   104.135     106.549     109.013     113.783
	 "s2_c"     27.6933    33.0952     36.5856     40.2953     49.3723
