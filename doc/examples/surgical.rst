.. index:: Examples; Surgical: Institutional Ranking

Surgical: Institutional Ranking
-------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` concerning mortality rates in 12 hospitals performing cardiac surgery in infants.

Model
^^^^^

Number of deaths are modelled as

.. math::

	r_i &\sim \text{Binomial}(n_i, p_i) \quad\quad i=1,\ldots,12 \\
	\operatorname{logit}(p_i) &= b_i \\
	b_i &\sim \text{Normal}(\mu, \sigma) \\
	\mu &\sim \text{Normal}(0, 1e6) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0001),
	
where :math:`r_i` are the number of deaths, out of :math:`n_i` operations, at hospital :math:`i`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: surgical.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
