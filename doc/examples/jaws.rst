.. index:: Examples; Jaws: Repeated Measures Analysis of Variance

.. _example-Jaws:

Jaws: Repeated Measures Analysis of Variance
--------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Elston and Grizzle :cite:`elston:1962:ETR` concerning jaw bone heights measured repeatedly in a cohort of 20 boys at ages 8, 8.5, 9, and 9.5 years. 

Model
^^^^^

Bone heights are modelled as

.. math::

	\bm{y}_i &\sim \text{Normal}(\bm{X} \bm{\beta}, \bm{\Sigma}) \quad\quad i=1,\ldots,20 \\
	\bm{X} &= \begin{bmatrix}
	  1 & 8 \\
	  1 & 8.5 \\
	  1 & 9 \\
	  1 & 9.5 \\
	\end{bmatrix} \quad
	\bm{\beta} = \begin{bmatrix}
	  \beta_0 \\
	  \beta_1 \\
	\end{bmatrix} \quad
	\bm{\Sigma} = \begin{bmatrix}
	  \sigma_{1,1} & \sigma_{1,2} & \sigma_{1,3} & \sigma_{1,4} \\
	  \sigma_{2,1} & \sigma_{2,2} & \sigma_{2,3} & \sigma_{2,4} \\
	  \sigma_{3,1} & \sigma_{3,2} & \sigma_{3,3} & \sigma_{3,4} \\
	  \sigma_{4,1} & \sigma_{4,2} & \sigma_{4,3} & \sigma_{4,4} \\
	\end{bmatrix} \\
	\beta_0, \beta_1 &\sim \text{Normal}(0, \sqrt{1000}) \\
	\bm{\Sigma} &\sim \text{InverseWishart}(4, \bm{I})
	
where :math:`\bm{y}_i` is a vector of the four repeated measurements for boy :math:`i`.  In the model specification below, the bone heights are arranged into a 1-dimensional vector on which a :ref:`section-Distribution-BDiagNormal` is specified.  Furthermore, since :math:`\bm{\Sigma}` is a covariance matrix, it is symmetric with ``M * (M + 1) / 2`` unique (upper or lower triangular) parameters, where ``M`` is the matrix dimension.  Consequently, that is the number of parameters to account for when defining samplers for :math:`\bm{\Sigma}`; e.g., ``AMWG([:Sigma], fill(0.1, int(M * (M + 1) / 2)))``.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: jaws.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia

	Iterations = 2502:10000
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 3750

	Empirical Posterior Estimates:
	13x6 Array{Any,2}:
	 ""              "Mean"   "SD"      "Naive SE"   "MCSE"         "ESS"
	 "beta0"       33.6421   1.91664   0.0221314    0.0673701   2463.79
	 "beta1"        1.87509  0.216427  0.00249909   0.00776451  2413.95
	 "Sigma[1,1]"   7.38108  2.70562   0.0312418    0.185902    1260.42
	 "Sigma[1,2]"   7.19741  2.69221   0.031087     0.190452    1224.21
	 "Sigma[1,3]"   6.78381  2.67475   0.0308853    0.192253    1204.87
	 "Sigma[1,4]"   6.53015  2.68514   0.0310053    0.192724    1206.6
	 "Sigma[2,2]"   7.53852  2.77055   0.0319916    0.19577     1225.61
	 "Sigma[2,3]"   7.20718  2.76476   0.0319247    0.198141    1208.41
	 "Sigma[2,4]"   6.96772  2.77756   0.0320725    0.198367    1212.62
	 "Sigma[3,3]"   8.03805  2.96494   0.0342362    0.205394    1250.14
	 "Sigma[3,4]"   8.00528  3.01322   0.0347936    0.206604    1263.06
	 "Sigma[4,4]"   8.58703  3.18705   0.0368009    0.210459    1311.45

	Quantiles:
	13x6 Array{Any,2}:
	 ""              "2.5%"    "25.0%"    "50.0%"    "75.0%"    "97.5%"
	 "beta0"       29.8247   32.4084    33.6549    34.9082    37.4228
	 "beta1"        1.4535    1.72672    1.87357    2.01553    2.30744
	 "Sigma[1,1]"   3.87079   5.40311    6.72138    8.7796    14.3201
	 "Sigma[1,2]"   3.69391   5.27486    6.5153     8.56786   14.0538
	 "Sigma[1,3]"   3.31306   4.88846    6.13707    8.15682   13.2664
	 "Sigma[1,4]"   3.02151   4.63627    5.91389    7.8605    13.1591
	 "Sigma[2,2]"   3.92078   5.5768     6.85708    8.95095   14.5351
	 "Sigma[2,3]"   3.61053   5.23388    6.52116    8.64111   14.0402
	 "Sigma[2,4]"   3.29096   4.99759    6.297      8.38234   13.6798
	 "Sigma[3,3]"   4.16453   5.90491    7.32127    9.61586   15.0376
	 "Sigma[3,4]"   4.04242   5.85098    7.29527    9.6252    15.276
	 "Sigma[4,4]"   4.37143   6.34044    7.84308   10.1665    16.4657
