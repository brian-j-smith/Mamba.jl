.. index:: Examples; Epilepsy: Repeated Measures on Poisson Counts

Epilepsy: Repeated Measures on Poisson Counts
---------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1`, Thall and Vail :cite:`thall:1990:SCM` Breslow and Clayton :cite:`breslow:1993:AIG` concerning the effects of treatment, baseline seizure counts, and age on follow-up seizure counts at four visits in 59 patients.

Model
^^^^^

Counts are modelled as

.. math::

	y_{i,j} &\sim \text{Poisson}(\mu_{i,j}) \quad\quad i=1,\ldots,59; j=1,\ldots,4 \\
	\log(\mu_{i,j}) &= \alpha_0 + \alpha_\text{Base} \log(\text{Base}_i / 4) +
	  \alpha_\text{Trt} \text{Trt}_i + \alpha_\text{BT} \text{Trt}_i \log(\text{Base}_i / 4) + \\
	  & \quad\quad \alpha_\text{Age} \log(\text{Age}_i) + \alpha_\text{V4} \text{V}_4 + \text{b1}_i +
	    \text{b}_{i,j} \\
	\text{b1}_i &\sim \text{Normal}(0, \sigma_\text{b1}) \\
	\text{b}_{i,j} &\sim \text{Normal}(0, \sigma_\text{b}) \\
	\alpha_* &\sim \text{Normal}(0, 100) \\
	\sigma^2_\text{b1}, \sigma^2_\text{b} &\sim \text{InverseGamma}(0.001, 0.001),
	  
	
where :math:`y_{ij}` are the counts on patient :math:`i` at visit :math:`j`, :math:`\text{Trt}` is a treatment indicator, :math:`\text{Base}` is baseline seizure counts, :math:`\text{Age}` is age in years, and :math:`\text{V}_4` is an indicator for the fourth visit.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: epil.jl
	:language: julia

Results
^^^^^^^

.. code-block:: julia
