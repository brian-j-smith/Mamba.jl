.. index:: Examples; Blocker: Random Effects Meta-Analysis of Clinical Trials

Blocker: Random Effects Meta-Analysis of Clinical Trials
--------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Carlin :cite:`carlin:1992:ME2` concerning a meta-analysis of 22 clinical trials to prevent mortality after myocardial infarction.

Model
^^^^^

Events are modelled as

.. math::

	r^c_i &\sim \text{Binomial}(n^c_i, p^c_i) \quad\quad i=1,\ldots,22 \\
	r^t_i &\sim \text{Binomial}(n^t_i, p^t_i) \\
	\operatorname{logit}(p^c_i) &= \mu_i \\
	\operatorname{logit}(p^t_i) &= \mu_i + \delta_i \\
	\mu_i &\sim \text{Normal}(0, 1000) \\
	\delta_i &\sim \text{Normal}(d, \sigma) \\
	d &\sim \text{Normal}(0, 1000) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`r^c_i` is the number of control group events, out of :math:`n^c_i`, in study :math:`i`; and :math:`r^t_i` is the number of treatment group events.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: blocker.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
