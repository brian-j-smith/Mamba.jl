.. index:: Examples; Equiv: Bioequivalence in a Cross-Over Trial

.. _example-Equiv:

Equiv: Bioequivalence in a Cross-Over Trial
-------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning a two-treatment, cross-over trial with 10 subjects.

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
	8x6 Array{Any,2}:
	 ""         "Mean"     "SD"       "Naive SE"   "MCSE"          "ESS"
	 "pi"     -0.178862   0.079636   0.00079636   0.00225377   1248.53
	 "mu"      1.44323    0.0473544  0.000473544  0.00203599    540.964
	 "equiv"   0.9835     0.127395   0.00127395   0.00229789   3073.57
	 "s2_2"    0.0187006  0.0146623  0.000146623  0.000479904   933.458
	 "s2_1"    0.0164033  0.014124   0.00014124   0.000507918   773.261
	 "theta"   0.983757   0.079366   0.00079366   0.00250413   1004.51
	 "phi"    -0.0195929  0.080053   0.00080053   0.00254328    990.753

	Quantiles:
	8x6 Array{Any,2}:
	 ""         "2.5%"       "25.0%"      "50.0%"     "75.0%"     "97.5%"
	 "pi"     -0.33275     -0.226525    -0.184088   -0.128557   -0.0142326
	 "mu"      1.35388      1.41116      1.44142     1.47408     1.53452
	 "equiv"   1.0          1.0          1.0         1.0         1.0
	 "s2_2"    0.00156773   0.00614605   0.0163017   0.0273322   0.0530448
	 "s2_1"    0.00103501   0.00437854   0.0136769   0.0242892   0.0507633
	 "theta"   0.838904     0.933652     0.974953    1.032       1.15675
	 "phi"    -0.175659    -0.0686517   -0.0253661   0.0315031   0.145616
