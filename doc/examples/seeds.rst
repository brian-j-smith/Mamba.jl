.. index:: Examples; Seeds: Random Effect Logistic Regression

Seeds: Random Effect Logistic Regression
----------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Crowder :cite:`crowder:1978:BBA` concerning the proportion of seeds that germinated on each of 21 plates arranged according to a 2 by 2 factorial layout by seed and type of root extract.

Model
^^^^^

Germinations are modelled as

.. math::

	r_i &\sim \text{Binomial}(n_i, p_i) \quad\quad i=1,\ldots,21 \\
	\operatorname{logit}(p_i) &= \alpha_0 + \alpha_1 x_{1i} + \alpha_2 x_{2i} + \alpha_{12} x_{1i} x_{2i} + b_i \\
	b_i &\sim \text{Normal}(0, \sigma^2) \\
	\alpha_0, \alpha_1, \alpha_2, \alpha_{12} &\sim \text{Normal}(0, 1e12) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`r_i` are the number of seeds, out of :math:`n_i` that germinate on plate :math:`i`; and :math:`x_{1i}` and :math:`x_{2i}` are the seed type and root extract.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: seeds.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
