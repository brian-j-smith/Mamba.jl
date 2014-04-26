.. index:: Examples; Rats: A Normal Hierarchical Model

Rats: A Normal Hierarchical Model
---------------------------------

An example from from section 6 of Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning 30 young rats whose weights were measured weekly for five weeks.

Model
^^^^^

.. math::

	y_{i,j} &\sim \text{Normal}\left(\alpha_i + \beta_i (x_j - \bar{x}), \sigma^2_c\right) \quad\quad i=1,\ldots,30; j=1,\ldots,5 \\
	\alpha_i &\sim \text{Normal}(\mu_\alpha, \sigma^2_\alpha) \\
	\beta_i &\sim \text{Normal}(\mu_\beta, \sigma^2_\beta) \\
	\mu_\alpha, \mu_\beta &\sim \text{Normal}(0, 1e12) \\
	\sigma^2_\alpha, \sigma^2_\beta &\sim \text{InverseGamma}(0.001, 0.001)
	
.. literalinclude:: rats.jl
	:language: julia
