.. index:: Examples; Mice: Weibull Regression

Mice: Weibull Regression
------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Grieve :cite:`grieve:1987:ABS`, and Dellaportas and Smith :cite:`dellaportas:1993:BIG` concerning time to death or censoring among four groups of 20 mice each. 

Model
^^^^^

Time to events are modelled as

.. math::

	t_i &\sim \text{Weibull}(\tau, \mu_i) \quad\quad i=1,\ldots,20 \\
	\mu_i &= \exp(-\bm{z}_i^\top \bm{\beta}) \\
	\beta_k &\sim \text{Normal}(0, 10) \\
	\tau &\sim \text{Gamma}(1, 0.001),
	
where :math:`t_i` is the time of death for mouse :math:`i`, and :math:`\bm{z}_i` is a vector of covariates.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: mice.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
