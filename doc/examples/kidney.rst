.. index:: Examples; Kidney: Weibull Regression with Random Effects

Kidney: Weibull Regression with Random Effects
----------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and McGilchrist and Aisbett :cite:`mcgilchrist:1991:RFS` concerning time to first and second infection recurrence in 38 kidney patients on dialysis. 

Model
^^^^^

Time to recurrences are modelled as

.. math::

	t_{i,j} &\sim \text{Weibull}(\tau, 1 / \mu_{i,j}) \quad\quad i=1,\ldots,38; j=1,2 \\
	\log(\mu_{i,j}) &= \bm{z}_{i,j}^\top \bm{\beta} + b_i \\
	b_i &\sim \text{Normal}(0, \sigma) \\
	\sigma^2 &\sim \text{InverseGamma}(0.0001, 0.0001) \\
	\beta_k &\sim \text{Normal}(0, 100) \\
	\tau &\sim \text{Gamma}(1, 0.001),
	
where :math:`t_{i,j}` is the time of infection :math:`j` in patient :math:`i`, :math:`\bm{z}_{i,j}` is a vector of covariates, and :math:`b_i` are patient-specific random effects.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: kidney.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia