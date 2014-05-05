.. index:: Examples; Dyes: Variance Components Model

Dyes: Variance Components Model
-------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Davies :cite:`davies:1967:SMR` concerning batch to batch variation in yields in six batches and five samples of dyestuff.

Model
^^^^^

Yields are modelled as

.. math::

	y_{i,j} &\sim \text{Normal}(\mu_i, \sigma_\text{within}) \quad\quad i=1,\ldots,6; j=1,\ldots,5 \\
	\mu_i &\sim \text{Normal}(\theta, \sigma_\text{between}) \\
	\theta &\sim \text{Normal}(0, 1000) \\
	\sigma^2_\text{within}, \sigma^2_\text{between} &\sim \text{InverseGamma}(0.001, 0.001),
	
where :math:`y_{i,j}` is the response for batch :math:`i` and sample :math:`j`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: dyes.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
