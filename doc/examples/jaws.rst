.. index:: Examples; Jaws: Repeated Measures Analysis of Variance

.. _example-Jaws:

Jaws: Repeated Measures Analysis of Variance
--------------------------------------------

An example from OpenBUGS :cite:`openbugs:2014:ex1` and Elston and Grizzle :cite:`elston:1962:ETR` concerning jaw bone heights measured repeatedly in a cohort of 20 boys at ages 8, 8.5, 9, and 9.5 years. 

Model
^^^^^

Bone heights are modelled as

.. math::

	\bm{y}_i &\sim \text{Normal}(\bm{X} \bm{\beta}, \bm{\Sigma}) \quad\quad i=1,\ldots,20 \\
	\bm{X} &= \begin{bmatrix}
	  1 & 8 \\
	  1 & 8.5 \\
	  1 & 9 \\
	  1 & 9.5 \\
	\end{bmatrix} \quad
	\bm{\beta} = \begin{bmatrix}
	  \beta_0 \\
	  \beta_1 \\
	\end{bmatrix} \quad
	\bm{\Sigma} = \begin{bmatrix}
	  \sigma_{1,1} & \sigma_{1,2} & \sigma_{1,3} & \sigma_{1,4} \\
	  \sigma_{2,1} & \sigma_{2,2} & \sigma_{2,3} & \sigma_{2,4} \\
	  \sigma_{3,1} & \sigma_{3,2} & \sigma_{3,3} & \sigma_{3,4} \\
	  \sigma_{4,1} & \sigma_{4,2} & \sigma_{4,3} & \sigma_{4,4} \\
	\end{bmatrix} \\
	\beta_0, \beta_1 &\sim \text{Normal}(0, \sqrt{1000}) \\
	\bm{\Sigma} &\sim \text{InverseWishart}(4, \bm{I})
	
where :math:`\bm{y}_i` is a vector of the four repeated measurements for boy :math:`i`.  In the model specification below, the bone heights are arranged into a 1-dimensional vector on which a :ref:`section-Distribution-BDiagNormal` is specified.  Furthermore, since :math:`\bm{\Sigma}` is a covariance matrix, it is symmetric with ``M * (M + 1) / 2`` unique (upper or lower triangular) parameters, where ``M`` is the matrix dimension.  Consequently, that is the number of parameters to account for when defining samplers for :math:`\bm{\Sigma}`; e.g., ``AMWG([:Sigma], fill(0.1, int(M * (M + 1) / 2)))``.
	

Analysis Program
^^^^^^^^^^^^^^^^

.. literalinclude:: jaws.jl
	:language: julia

Results
^^^^^^^

