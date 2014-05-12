.. index:: Examples; Leuk: Cox Regression

Leuk: Cox Regression
--------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Ezzet and Whitehead :cite:`frierich:1963:E6M` concerning survival in 42 leukemia patients treated with 6-mercaptopurine or placebo. 

Model
^^^^^

Times to death are using the Bayesian Cox proportional hazards model formulated by Clayton :cite:`clayton:1994:BAF` as

.. math::

	dN_i(t) &\sim \text{Poisson}(I_i(t) dt) \quad\quad i=1,\ldots,42 \\
	I_i(t)dt &= Y_i(t) \exp(\beta Z_i) d\Lambda_0(t) \\
	\beta &\sim \text{Normal}(0, 1000) \\
	d\Lambda_0(t) &\sim \text{Gamma}(c d\Lambda_0^*(t), c) \\
	d\Lambda_0^*(t) &= r dt \\
	c &= 0.001 \\
	r &= 0.1,
	
where :math:`dN_i(t)` is a counting process increment in time interval :math:`[t, t + dt)` for patient :math:`i`; :math:`Y_i(t)` is an indicator for whether the patient is observed at time :math:`t`; :math:`\bm{z}_i` is a vector of covariates;  and :math:`d\Lambda_0(t)` is the increment in the integrated baseline hazard function during :math:`[t, t + dt)`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: leuk.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
