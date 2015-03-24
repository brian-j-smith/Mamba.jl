.. index:: Examples; Blocker: Random Effects Meta-Analysis of Clinical Trials

.. _example-Blocker:

Blocker: Random Effects Meta-Analysis of Clinical Trials
--------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Carlin :cite:`carlin:1992:ME2` concerning a meta-analysis of 22 clinical trials to prevent mortality after myocardial infarction.

Model
^^^^^

Events are modelled as

.. math::

	r^c_i &\sim \text{Binomial}(n^c_i, p^c_i) \quad\quad i=1,\ldots,22 \\
	r^t_i &\sim \text{Binomial}(n^t_i, p^t_i) \\
	\operatorname{logit}(p^c_i) &= \mu_i \\
	\operatorname{logit}(p^t_i) &= \mu_i + \delta_i \\
	\mu_i &\sim \text{Normal}(0, 1000) \\
	\delta_i &\sim \text{Normal}(d, \sigma) \\
	d &\sim \text{Normal}(0, 1000) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`r^c_i` is the number of control group events, out of :math:`n^c_i`, in study :math:`i`; and :math:`r^t_i` is the number of treatment group events.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: blocker.jl
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
	 ""             "Mean"     "SD"       "Naive SE"   "MCSE"         "ESS"
	 "delta_new"  -0.250058   0.150325   0.00173581   0.00429703  1223.85
	 "s2"          0.0182219  0.0211213  0.000243887  0.00124973   285.634
	 "d"          -0.255636   0.0618419  0.000714089  0.00304573   412.272

	Quantiles:
	4x6 Array{Any,2}:
	 ""             "2.5%"        "25.0%"      "50.0%"     "75.0%"     "97.5%"
	 "delta_new"  -0.538541     -0.327996    -0.255785   -0.177588    0.0798606
	 "s2"          0.000685545   0.00416488   0.0107616   0.0244421   0.0773571
	 "d"          -0.373412     -0.295917    -0.258185   -0.218341   -0.128426
