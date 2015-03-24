.. index:: Examples; Epilepsy: Repeated Measures on Poisson Counts

.. _example-Epilepsy:

Epilepsy: Repeated Measures on Poisson Counts
---------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Thall and Vail :cite:`thall:1990:SCM` Breslow and Clayton :cite:`breslow:1993:AIG` concerning the effects of treatment, baseline seizure counts, and age on follow-up seizure counts at four visits in 59 patients.

Model
^^^^^

Counts are modelled as

.. math::

	y_{i,j} &\sim \text{Poisson}(\mu_{i,j}) \quad\quad i=1,\ldots,59; j=1,\ldots,4 \\
	\log(\mu_{i,j}) &= \alpha_0 + \alpha_\text{Base} \log(\text{Base}_i / 4) +
	  \alpha_\text{Trt} \text{Trt}_i + \alpha_\text{BT} \text{Trt}_i \log(\text{Base}_i / 4) + \\
	  & \quad\quad \alpha_\text{Age} \log(\text{Age}_i) + \alpha_\text{V4} \text{V}_4 + \text{b1}_i +
	    \text{b}_{i,j} \\
	\text{b1}_i &\sim \text{Normal}(0, \sigma_\text{b1}) \\
	\text{b}_{i,j} &\sim \text{Normal}(0, \sigma_\text{b}) \\
	\alpha_* &\sim \text{Normal}(0, 100) \\
	\sigma^2_\text{b1}, \sigma^2_\text{b} &\sim \text{InverseGamma}(0.001, 0.001),
	  
	
where :math:`y_{ij}` are the counts on patient :math:`i` at visit :math:`j`, :math:`\text{Trt}` is a treatment indicator, :math:`\text{Base}` is baseline seizure counts, :math:`\text{Age}` is age in years, and :math:`\text{V}_4` is an indicator for the fourth visit.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: epil.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia

	Iterations = 2502:15000
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 6250

	Empirical Posterior Estimates:
	9x6 Array{Any,2}:
	 ""              "Mean"     "SD"       "Naive SE"   "MCSE"        "ESS"
	 "alpha_Age"    0.458309   0.394536   0.00352884   0.0203364   376.381
	 "alpha0"      -1.35617    1.31324    0.011746     0.0721      331.755
	 "alpha_BT"     0.24217    0.190566   0.00170448   0.0107593   313.706
	 "alpha_Base"   0.91105    0.135354   0.00121065   0.00720844  352.584
	 "s2_b"         0.135238   0.0318193  0.0002846    0.00155135  420.688
	 "alpha_Trt"   -0.759314   0.397734   0.00355744   0.0234788   286.968
	 "s2_b1"        0.249119   0.0731667  0.000654423  0.00290063  636.271
	 "alpha_V4"    -0.0928793  0.0836669  0.000748339  0.00360419  538.878

	Quantiles:
	9x6 Array{Any,2}:
	 ""              "2.5%"      "25.0%"    "50.0%"    "75.0%"     "97.5%"
	 "alpha_Age"   -0.19667     0.176356   0.416087   0.696648    1.30508
	 "alpha0"      -4.16889    -2.15793   -1.26343   -0.436226    0.866196
	 "alpha_BT"    -0.0902501   0.108103   0.226562   0.358354    0.657805
	 "alpha_Base"   0.663188    0.817701   0.902682   0.997417    1.20062
	 "s2_b"         0.0715581   0.112591   0.136265   0.158033    0.193716
	 "alpha_Trt"   -1.63682    -1.01139   -0.75654   -0.480871   -0.0161134
	 "s2_b1"        0.138175    0.197135   0.237611   0.289655    0.422805
	 "alpha_V4"    -0.255045   -0.148157  -0.093136  -0.0366814   0.072099
