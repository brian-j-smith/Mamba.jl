.. index:: Examples; Pumps: Gamma-Poisson Hierarchical Model

.. _example-Pumps:

Pumps: Gamma-Poisson Hierarchical Model
---------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and George *et al.* :cite:`george:1993:CLD` concerning the number of failures of 10 power plant pumps.

Model
^^^^^

Pump failure are modelled as

.. math::

	y_i &\sim \text{Poisson}(\theta_i t_i) \quad\quad i=1,\ldots,10 \\
	\theta_i &\sim \text{Gamma}(\alpha, 1 / \beta) \\
	\alpha &\sim \text{Gamma}(1, 1) \\
	\beta &\sim \text{Gamma}(0.1, 1),
	
where :math:`y_i` is the number of times that pump :math:`i` failed, and :math:`t_i` is the operation time of the pump (in 1000s of hours).
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: pumps.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia

	## MCMC Simulations
	
	Iterations = 2502:10000
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 3750

	Empirical Posterior Estimates:
	13x6 Array{Any,2}:
	 ""            "Mean"     "SD"       "Naive SE"   "MCSE"          "ESS"
	 "theta[1]"   0.0597276  0.0251118  0.000289966  0.000357021  6091.37
	 "theta[2]"   0.0997635  0.0773339  0.000892976  0.00113173   5917.76
	 "theta[3]"   0.0889841  0.0379549  0.000438265  0.000510879  6433.99
	 "theta[4]"   0.116179   0.030564   0.000352922  0.000382914  6912.57
	 "theta[5]"   0.603654   0.32236    0.00372229   0.00612325   4559.21
	 "theta[6]"   0.608916   0.139405   0.00160972   0.00155659   7500.0
	 "theta[7]"   0.894307   0.691573   0.0079856    0.0301097    1989.13
	 "theta[8]"   0.884304   0.726039   0.00838357   0.0252293    2492.22
	 "theta[9]"   1.56235    0.755878   0.00872813   0.0238079    2749.55
	 "theta[10]"  1.98785    0.427868   0.0049406    0.0102204    3625.53
	 "alpha"      0.679247   0.267001   0.00308307   0.00708272   3264.71
	 "beta"       0.893892   0.528231   0.00609949   0.0182549    2505.97

	Quantiles:
	13x6 Array{Any,2}:
	 ""            "2.5%"      "25.0%"    "50.0%"    "75.0%"    "97.5%"
	 "theta[1]"   0.0210364   0.0413043  0.0561389  0.0742908  0.117973
	 "theta[2]"   0.00770316  0.0433944  0.0812108  0.135976   0.296096
	 "theta[3]"   0.0309151   0.0616544  0.0841478  0.110698   0.178913
	 "theta[4]"   0.0635652   0.0944482  0.113746   0.135176   0.182949
	 "theta[5]"   0.14817     0.369621   0.549603   0.775779   1.38092
	 "theta[6]"   0.366678    0.510037   0.599759   0.695918   0.915158
	 "theta[7]"   0.0781523   0.374379   0.722913   1.22739    2.6892
	 "theta[8]"   0.0734788   0.373673   0.690015   1.19953    2.72981
	 "theta[9]"   0.462534    1.00451    1.43995    1.99099    3.35956
	 "theta[10]"  1.23458     1.68077    1.95688    2.25502    2.90776
	 "alpha"      0.278598    0.485673   0.640735   0.825628   1.32444
	 "beta"       0.177017    0.502637   0.783992   1.17955    2.19996
	 
	 ## Posterior Predictive Distribution
	 
	Iterations = 2502:10000
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 3750

	Empirical Posterior Estimates:
	11x6 Array{Any,2}:
	 ""         "Mean"   "SD"     "Naive SE"   "MCSE"        "ESS"
	 "y[1]"    5.668    3.4152   0.0394353    0.0440254  6718.06
	 "y[2]"    1.53347  1.70394  0.0196754    0.0240401  6138.29
	 "y[3]"    5.56853  3.38075  0.0390376    0.0406907  7195.29
	 "y[4]"   14.6473   5.4517   0.0629508    0.0659796  7155.71
	 "y[5]"    3.1764   2.46365  0.0284478    0.0374077  5703.58
	 "y[6]"   19.1655   6.20492  0.0716483    0.0664772  7500.0
	 "y[7]"    0.9384   1.20349  0.0138968    0.0347485  2999.43
	 "y[8]"    0.9372   1.22294  0.0141213    0.0279075  3795.03
	 "y[9]"    3.28933  2.39193  0.0276196    0.0516479  4010.75
	 "y[10]"  20.8635   6.39303  0.0738203    0.111248   4976.72

	Quantiles:
	11x6 Array{Any,2}:
	 ""         "2.5%"    "25.0%"    "50.0%"    "75.0%"    "97.5%"
	 "y[1]"    1.0       3.0        5.0        8.0       14.0
	 "y[2]"    0.0       0.0        1.0        2.0        6.0
	 "y[3]"    1.0       3.0        5.0        7.0       14.0
	 "y[4]"    5.0      11.0       14.0       18.0       26.0
	 "y[5]"    0.0       1.0        3.0        4.0        9.0
	 "y[6]"    9.0      15.0       19.0       23.0       33.0
	 "y[7]"    0.0       0.0        1.0        1.0        4.0
	 "y[8]"    0.0       0.0        1.0        1.0        4.0
	 "y[9]"    0.0       2.0        3.0        5.0        9.0
	 "y[10]"  10.0      16.0       20.0       25.0       35.0
