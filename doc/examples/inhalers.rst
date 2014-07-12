.. index:: Examples; Inhalers: Ordered Categorical Data

.. _example-Inhalers:

Inhalers: Ordered Categorical Data
----------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Ezzet and Whitehead :cite:`ezzet:1993:REM` concerning a two-treatment, two-period crossover trial comparing salbutamol inhalation devices in 286 asthma patients. 

Model
^^^^^

Treatment responses are modelled as

.. math::

	R_{i,t} &= j \quad \text{if}\; Y_{i,t} \in [a_{j-1}, a_j) \quad\quad i=1,\ldots,4; t=1,2; j=1,\ldots,3 \\
	\operatorname{logit}(Q_{i,t,j}) &= -(a_j + \mu_{s_i,t} + b_i) \\
	\mu_{1,1} &= \beta / 2 + \pi / 2 \\
	\mu_{1,2} &= -\beta / 2 - \pi / 2 - \kappa \\
	\mu_{2,1} &= -\beta / 2 + \pi / 2 \\
	\mu_{2,2} &= \beta / 2 - \pi / 2 + \kappa \\
	b_i &\sim \text{Normal}(0, \sigma) \\
	a[1] &\sim \text{Flat}(-1000, a[2]) \\
	a[2] &\sim \text{Flat}(-1000, a[3]) \\
	a[3] &\sim \text{Flat}(-1000, 1000) \\
	\beta &\sim \text{Normal}(0, 1000) \\
	\pi &\sim \text{Normal}(0, 1000) \\
	\kappa &\sim \text{Normal}(0, 1000) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`R_{i,t}` is a 4-point ordinal rating of the device used by patient :math:`i`, and :math:`Q_{i,t,j}` is the cumulative probability of the rating in treatment period :math:`t` being worse than category :math:`j`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: inhalers.jl
	:language: julia

Results
^^^^^^^

