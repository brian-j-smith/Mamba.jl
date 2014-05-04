.. index:: Examples; Stacks: Robust Regression

Stacks: Robust Regression
-------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Birkes and Dodge :cite:`birkes:1993:AMR` concerning 21 daily responses of stack loss, the amount of ammonia escaping, as a function of air flow, temperature, and acid concentration.

Model
^^^^^

Losses are modelled as

.. math::

	y_i &\sim \text{Laplace}(\mu_i, \sigma) \quad\quad i=1,\ldots,21 \\
	\mu_i &= \beta_0 + \beta_1 z_{1i} + \beta_2 z_{2i} + \beta_3 z_{3i} \\
	\beta_0, \beta_1, \beta_2, \beta_3 &\sim \text{Normal}(0, 1e5) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`y_i` is the stack loss on day :math:`i`; and :math:`z_{1i}, z_{2i}, z_{3i}` are standardized predictors.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: stacks.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
