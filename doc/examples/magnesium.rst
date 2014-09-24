.. index:: Examples; Magnesium: Meta-Analysis Prior Sensitivity

.. _example-Magnesium:

Magnesium: Meta-Analysis Prior Sensitivity
------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`.

Model
^^^^^

Number of events reported for treatment and control subjects in 8 studies is modelled as

.. math::

	r^c_j &\sim \text{Binomial}(n^c_j, p^c_{i,j}) \quad\quad i=1,\ldots,6; j=1,\ldots,8 \\
	p^c_{i,j} &\sim \text{Uniform}(0, 1) \\
	r^t_j &\sim \text{Binomial}(n^t_j, p^t_{i,j}) \\
	\operatorname{logit}(p^t_{i,j}) &= \theta_{i,j} + \operatorname{logit}(p^c_{i,j}) \\
	\theta_{i,j} &\sim \text{Normal}(\mu_i, \tau_i) \\
	\mu_i &\sim \text{Uniform}(-10, 10) \\
	\tau_i &\sim \text{Different Priors},

where :math:`r^c_j` is the number of control group events, out of :math:`n^c_j`, in study :math:`j`; :math:`r^t_j` is the number of treatment group events; and :math:`i` indexes differ prior specifications.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: magnesium.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia

	Iterations = 2502:12500
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 5000

	Empirical Posterior Estimates:
	13x6 Array{Any,2}:
	 ""         "Mean"    "SD"      "Naive SE"   "MCSE"         "ESS"
	 "OR[1]"   0.493881  0.167242  0.00167242   0.00678096  2466.35
	 "OR[2]"   0.469325  1.10491   0.0110491    0.0301488   3664.86
	 "OR[3]"   0.444941  0.214777  0.00214777   0.00542145  3961.61
	 "OR[4]"   0.47496   0.138735  0.00138735   0.00596257  2326.77
	 "OR[5]"   0.493619  0.150459  0.00150459   0.00634346  2371.87
	 "OR[6]"   0.445284  0.144975  0.00144975   0.00555468  2609.96
	 "tau[1]"  0.500508  0.380023  0.00380023   0.0194311   1955.75
	 "tau[2]"  1.1646    0.75346   0.0075346    0.0379042   1987.8
	 "tau[3]"  0.808261  0.492263  0.00492263   0.0186635   2637.57
	 "tau[4]"  0.467692  0.270444  0.00270444   0.0114467   2362.63
	 "tau[5]"  0.45235   0.347114  0.00347114   0.0172585   2011.27
	 "tau[6]"  0.57345   0.194067  0.00194067   0.00629904  3080.9

	Quantiles:
	13x6 Array{Any,2}:
	 ""         "2.5%"     "25.0%"   "50.0%"   "75.0%"   "97.5%"
	 "OR[1]"   0.209056   0.382552  0.487851  0.593326  0.813823
	 "OR[2]"   0.112234   0.270276  0.382542  0.514492  1.03727
	 "OR[3]"   0.159154   0.319832  0.420812  0.534651  0.880643
	 "OR[4]"   0.215638   0.382852  0.472237  0.560446  0.764725
	 "OR[5]"   0.208413   0.387238  0.493571  0.608529  0.765637
	 "OR[6]"   0.211006   0.344643  0.429437  0.527401  0.775947
	 "tau[1]"  0.0367253  0.206871  0.432488  0.698297  1.46017
	 "tau[2]"  0.295543   0.684773  0.983598  1.40115   3.30371
	 "tau[3]"  0.131057   0.482964  0.714225  1.02093   2.01561
	 "tau[4]"  0.0783779  0.284607  0.417649  0.592994  1.13504
	 "tau[5]"  0.020413   0.183672  0.397168  0.635079  1.28576
	 "tau[6]"  0.203636   0.438994  0.571246  0.703771  0.960932
