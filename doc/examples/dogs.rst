.. index:: Examples; Dogs: Loglinear Model for Binary Data

.. _example-Dogs:

Dogs: Loglinear Model for Binary Data
-------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Lindley and Smith :cite:`lindley:1972:BEL`, and Kalbfleisch :cite:`kalbfleisch:1985:PSI` concerning the Solomon-Wynne experiment on dogs.  In the experiment, 30 dogs were subjected to 25 trials.  On each trial, a barrier was raised, and an electric shock was administered 10 seconds later if the dog did not jump the barrier.

Model
^^^^^

Failures to jump the barriers in time are modelled as

.. math::

	y_{i,j} &= \text{Bernoulli}(\pi_{i,j}) \quad\quad i=1,\ldots,30; j=2,\ldots,25 \\
	\log(\pi_{i,j}) &= \alpha x_{i,j-1} + \beta (j - 1 - x_{i,j-1}) \\
	\alpha, \beta &\sim \text{Flat}(-\infty, -0.00001) \\
	
where :math:`y_{i,j} = 1` if dog :math:`i` fails to jump the barrier before the shock on trial :math:`j`, and 0 otherwise; :math:`x_{i,j-1}` is the number of successful jumps prior to trial :math:`j`; and :math:`\pi_{i,j}` is the probability of failure.

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: dogs.jl
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
	 ""         "Mean"     "SD"       "Naive SE"   "MCSE"          "ESS"
	 "A"       0.783585   0.0188211  0.000217328  0.000347073  4696.3
	 "B"       0.924234   0.0109032  0.000125899  0.000168008  5620.23
	 "alpha"  -0.244165   0.0240644  0.000277872  0.000443304  4701.14
	 "beta"   -0.0788601  0.0118129  0.000136403  0.00018221   5614.55

	Quantiles:
	5x6 Array{Any,2}:
	 ""         "2.5%"     "25.0%"     "50.0%"    "75.0%"     "97.5%"
	 "A"       0.746251   0.770922    0.783633   0.796249    0.819771
	 "B"       0.901885   0.917013    0.924516   0.931959    0.944565
	 "alpha"  -0.292693  -0.260168   -0.243815  -0.227843   -0.19873
	 "beta"   -0.103268  -0.0866338  -0.078485  -0.0704666  -0.0570308
