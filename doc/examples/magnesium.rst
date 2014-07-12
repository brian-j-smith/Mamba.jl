.. index:: Examples; Magnesium: Meta-Analysis Prior Sensitivity

.. _example-Magnesium:

Magnesium: Meta-Analysis Prior Sensitivity
------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`.

Model
^^^^^

Number of events reported for treatment and control subjects in 8 studies is modelled as

.. math::

	r^c_j &\sim \text{Binomial}(n^c_j, p^c_{i,j}) \quad\quad i=1,\ldots,6; j=1,\ldots,8 \\
	p^c_{i,j} &\sim \text{Uniform}(0, 1) \\
	r^t_j &\sim \text{Binomial}(n^t_j, p^t_{i,j}) \\
	\operatorname{logit}(p^t_{i,j}) &= \theta_{i,j} + \operatorname{logit}(p^c_{i,j}) \\
	\theta_{i,j} &\sim \text{Normal}(\mu_i, \tau_i) \\
	\mu_i &\sim \text{Uniform}(-10, 10) \\
	\tau_i &\sim \text{Different Priors},

where :math:`r^c_j` is the number of control group events, out of :math:`n^c_j`, in study :math:`j`; :math:`r^t_j` is the number of treatment group events; and :math:`i` indexes differ prior specifications.


Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: magnesium.jl
	:language: julia

Results
^^^^^^^

