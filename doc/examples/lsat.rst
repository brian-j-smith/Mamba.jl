.. index:: Examples; LSAT: Item Response

.. _example-LSAT:

LSAT: Item Response
-------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Boch and Lieberman :cite:`boch:1970:FRM` concerning a 5-item multiple choice test (32 possible response patters) given to 1000 students. 

Model
^^^^^

Item responses are modelled as

.. math::

	r_{i,j} &\sim \text{Bernoulli}(p_{i,j}) \quad\quad i=1,\ldots,1000; j=1,\ldots,5 \\
	\operatorname{logit}(p_{i,j}) &= \beta \theta_i - \alpha_j \\
	\theta_i &\sim \text{Normal}(0, 1) \\
	\alpha_j &\sim \text{Normal}(0, 100) \\
	\beta &\sim \text{Flat}(0, \infty),
	
where :math:`r_{i,j}` is an indicator for correct response by student :math:`i` to questions :math:`j`.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: lsat.jl
	:language: julia

Results
^^^^^^^

