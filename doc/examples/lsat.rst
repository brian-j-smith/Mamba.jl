.. index:: Examples; LSAT: Item Response

.. _example-LSAT:

LSAT: Item Response
-------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Boch and Lieberman :cite:`boch:1970:FRM` concerning a 5-item multiple choice test (32 possible response patters) given to 1000 students. 

Model
^^^^^

Item responses are modelled as

.. math::

	r_{i,j} &\sim \text{Bernoulli}(p_{i,j}) \quad\quad i=1,\ldots,1000; j=1,\ldots,5 \\
	\operatorname{logit}(p_{i,j}) &= \beta \theta_i - \alpha_j \\
	\theta_i &\sim \text{Normal}(0, 1) \\
	\alpha_j &\sim \text{Normal}(0, 100) \\
	\beta &\sim \text{Flat}(0, \infty),
	
where :math:`r_{i,j}` is an indicator for correct response by student :math:`i` to questions :math:`j`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: lsat.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia

	Iterations = 2502:10000
	Thinning interval = 2
	Chains = 1,2
	Samples per chain = 3750

	Empirical Posterior Estimates:
	7x6 Array{Any,2}:
	 ""        "Mean"    "SD"       "Naive SE"   "MCSE"         "ESS"
	 "a[1]"  -1.26311   0.105081   0.00121337   0.00222283  4094.01
	 "a[2]"   0.478177  0.0690242  0.000797023  0.00126502  4725.35
	 "a[3]"   1.24018   0.0679403  0.000784507  0.0015522   3790.63
	 "a[4]"   0.17052   0.0716647  0.000827513  0.00146378  4239.95
	 "a[5]"  -0.625759  0.0858278  0.000991054  0.00156547  4748.04
	 "beta"   0.76483   0.0567109  0.000654841  0.0027709   1772.46

	Quantiles:
	7x6 Array{Any,2}:
	 ""        "2.5%"      "25.0%"    "50.0%"    "75.0%"    "97.5%"
	 "a[1]"  -1.47976    -1.33259   -1.26218   -1.18982   -1.06472
	 "a[2]"   0.349339    0.430022   0.477481   0.525347   0.614368
	 "a[3]"   1.10908     1.1937     1.24051    1.28659    1.36991
	 "a[4]"   0.0285454   0.122158   0.170171   0.219162   0.31003
	 "a[5]"  -0.797083   -0.682846  -0.625013  -0.568563  -0.45901
	 "beta"   0.661835    0.724974   0.761457   0.801351   0.883908
