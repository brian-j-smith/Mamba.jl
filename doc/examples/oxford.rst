.. index:: Examples; Oxford: Smooth Fit to Log-Odds Ratios

.. _example-Oxford:

Oxford: Smooth Fit to Log-Odds Ratios
-------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Breslow and Clayton :cite:`breslow:1993:AIG` concerning the association between death from childhood cancer and maternal exposure to X-rays, for subjects partitioned into 120 age and birth-year strata. 

Model
^^^^^

Deaths are modelled as

.. math::

	r^0_i &\sim \text{Binomial}(n^0_i, p^0_i) \quad\quad i=1,\ldots,120 \\
	r^1_i &\sim \text{Binomial}(n^1_i, p^1_i) \\
	\operatorname{logit}(p^0_i) &= \mu_i \\
	\operatorname{logit}(p^1_i) &= \mu_i + \log(\psi_i) \\
	\log(\psi) &= \alpha + \beta_1 \text{year}_i + \beta_2 (\text{year}^2_i - 22) + b_i \\
	\mu_i &\sim \text{Normal}(0, 1000) \\
	b_i &\sim \text{Normal}(0, \sigma) \\
	\alpha, \beta_1, \beta_2 &\sim \text{Normal}(0, 1000) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`r^0_i` is the number of deaths among unexposed subjects in stratum :math:`i`, :math:`r^1_i` is the number among exposed subjects, and :math:`\text{year}_i` is the stratum-specific birth year (relative to 1954).
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: oxford.jl
	:language: julia

Results
^^^^^^^

