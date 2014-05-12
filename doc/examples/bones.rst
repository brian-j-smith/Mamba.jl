.. index:: Examples; Bones: Latent Trait Model for Multiple Ordered Categorical Responses

Bones: Latent Trait Model for Multiple Ordered Categorical Responses
--------------------------------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Roche *et al.* :cite:`roche:1975:SM`, and Thissen :cite:`thissen:1986:MUG` concerning skeletal age in 13 boys predicted from 34 radiograph indicators of skeletal maturity. 

Model
^^^^^

Skeletal ages are modelled as

.. math::

	\operatorname{logit}(Q_{i,j,k}) &= \delta_j (\theta_i - \gamma_{j,k}) \quad\quad i=1,\ldots,13; j=1,\ldots,34; k=1,\ldots,4 \\
	\theta_i &\sim \text{Normal}(0, 100),
	
	
where :math:`\delta_j` is a discriminability parameter for indicator :math:`j`, :math:`\gamma_{j,k}` is a threshold parameter, and :math:`Q_{i,j,k}` is the cumulative probability that boy :math:`i` with skeletal age :math:`\theta_i` is assigned a more mature grade than :math:`k`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: bones.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
