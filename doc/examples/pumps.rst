.. index:: Examples; Pumps: Gamma-Poisson Hierarchical Model

Pumps: Gamma-Poisson Hierarchical Model
---------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and George *et al.* :cite:`george:1993:CLD` concerning the number of failures of 10 power plant pumps.

Model
^^^^^

Pump failure are modelled as

.. math::

	y_i | \theta_i &\sim \text{Poisson}(\theta_i t_i) \quad\quad i=1,\ldots,10 \\
	\theta_i | \alpha, \beta &\sim \text{Gamma}(\alpha, \beta) \\
	\alpha &\sim \text{Gamma}(1, 1) \\
	\beta &\sim \text{Gamma}(0.1, 1),
	
where :math:`y_i` is the number of times that pump :math:`i` failed, and :math:`t_i` is the operation time of the pump (in 1000s of hours).
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: pumps.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
