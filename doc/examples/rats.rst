.. index:: Examples; Rats: A Normal Hierarchical Model

Rats: A Normal Hierarchical Model
---------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and section 6 of Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning 30 rats whose weights were measured at each of five consecutive weeks.

Model
^^^^^

Weights are modeled as

.. math::

	y_{i,j} &\sim \text{Normal}\left(\alpha_i + \beta_i (x_j - \bar{x}), \sigma_c\right) \quad\quad i=1,\ldots,30; j=1,\ldots,5 \\
	\alpha_i &\sim \text{Normal}(\mu_\alpha, \sigma_\alpha) \\
	\beta_i &\sim \text{Normal}(\mu_\beta, \sigma_\beta) \\
	\mu_\alpha, \mu_\beta &\sim \text{Normal}(0, 1000) \\
	\sigma^2_\alpha, \sigma^2_\beta &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`y_{i,j}` is repeated weight measurement :math:`j` on rat :math:`i`, and :math:`x_j` is the day on which the measurement was taken.
	
Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: rats.jl
	:language: julia

Results
^^^^^^^
