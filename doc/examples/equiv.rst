.. index:: Examples; Equiv: Bioequivalence in a Cross-Over Trial

Equiv: Bioequivalence in a Cross-Over Trial
-------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Gelfand *et al.* :cite:`gelfand:1990:IBI` concerning a two-treatment, cross-over trial with 10 subjects.

Model
^^^^^

Treatment responses are modelled as

.. math::

	y_{i,j} &\sim \text{Normal}(m_{i,j}, \sigma_1) \quad\quad i=1,\ldots,10; j=1,2 \\
	m_{i,j} &= \mu + (-1)^{T_{i,j} - 1} \phi / 2 + (-1)^{j-1} \pi / 2 + \delta_i \\
	\delta_i &\sim \text{Normal}(0, \sigma_2) \\
	\mu, \phi, \pi &\sim \text{Normal}(0, 1e6) \\
	\sigma_1^2, \sigma_2^2 &\sim \text{InverseGamma}(0.001, 0.001) \\
	
where :math:`y_{i,j}` is the response for patient :math:`i` in period :math:`j`; and :math:`T_{i,j} = 1,2` is the treatment received.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: equiv.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
