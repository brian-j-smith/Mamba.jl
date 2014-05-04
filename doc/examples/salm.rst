.. index:: Examples; Salm: Extra-Poisson Variation in a Dose-Response Study

Salm: Extra-Poisson Variation in a Dose-Response Study
------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Breslow :cite:`breslow:1984:EPV` concerning mutagenicity assay data on salmonella in three plates exposed to six doses of quinoline.

Model
^^^^^

Number of revertant colonies of salmonella are modelled as

.. math::

	y_{i,j} &\sim \text{Poisson}(\mu_{i,j}) \quad\quad i=1,\ldots,3; j=1,\ldots,6 \\
	\log(\mu_{i,j}) &= \alpha + \beta \log(x_j + 10) + \gamma x_j + \lambda_{i,j} \\
	\alpha, \beta, \gamma &\sim \text{Normal}(0, 1e6) \\
	\lambda_{i,j} &\sim \text{Normal}(0, \sigma) \\
	\sigma^2 &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`y_i` is the number of colonies in plate :math:`i` and dose :math:`j`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: salm.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
