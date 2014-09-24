.. index:: Examples; Dyes: Variance Components Model

.. _example-Dyes:

Dyes: Variance Components Model
-------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Davies :cite:`davies:1967:SMR`, and Box and Tiao :cite:`box:1973:BIS` concerning batch-to-batch variation in yields from six batches and five samples of dyestuff.

Model
^^^^^

Yields are modelled as

.. math::

	y_{i,j} &\sim \text{Normal}(\mu_i, \sigma_\text{within}) \quad\quad i=1,\ldots,6; j=1,\ldots,5 \\
	\mu_i &\sim \text{Normal}(\theta, \sigma_\text{between}) \\
	\theta &\sim \text{Normal}(0, 1000) \\
	\sigma^2_\text{within}, \sigma^2_\text{between} &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`y_{i,j}` is the response for batch :math:`i` and sample :math:`j`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: dyes.jl
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
	 ""                "Mean"      "SD"    "Naive SE"     "MCSE"       "ESS"
	 "s2_between"  2442.3      2281.59   26.3455       150.701     1311.15
	 "s2_within"   2778.03      917.079  10.5895        43.5067    1825.5
	 "theta"       1526.43       22.385   0.25848        0.357869  5417.07

	Quantiles:
	4x6 Array{Any,2}:
	 ""                "2.5%"      "25.0%"  .      "75.0%"       "97.5%"
	 "s2_between"   170.563     902.938        3219.21      10329.7
	 "s2_within"   1534.53     2137.37         3218.14       5062.93
	 "theta"       1480.84     1513.62         1539.09       1572.34
