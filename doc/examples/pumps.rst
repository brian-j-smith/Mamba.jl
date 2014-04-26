.. index:: Examples; Pumps: Gamma-Poisson Hierarchical Model

Pumps: Gamma-Poisson Hierarchical Model
---------------------------------------

An example from from George *et al.* :cite:`george:1993:CLD` concerning the number of failures of 10 power plant pumps.

Model
^^^^^

.. math::

	y_i | \theta_i &\sim \text{Poisson}(\theta_i t_i) \quad\quad i=1,\ldots,10 \\
	\theta_i | \alpha, \beta &\sim \text{Gamma}(\alpha, \beta) \\
	\alpha &\sim \text{Gamma}(1, 1) \\
	\beta &\sim \text{Gamma}(0.1, 1)
	
.. literalinclude:: pumps.jl
	:language: julia
