.. _section-Supplement:

Supplement
==========

Bayesian Linear Regression Model
--------------------------------

The unnormalized posterior distribution was given for a :ref:`section-Line-Model` in the tutorial.  Additional forms of that posterior are given in the following section.

Log-Transformed Distribution and Gradient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`\mathcal{L}` denote the logarithm of a density of interest up to a normalizing constant, and :math:`\nabla \mathcal{L}` its gradient.  Then, the following are obtained for the regression example parameters :math:`\bm{\beta}` and :math:`\theta = \log(\sigma^2)`, for samplers, like NUTS, that can utilize both.

.. math::

	\mathcal{L}(\bm{\beta}, \theta | \bm{y}) &= \log(p(\bm{y} | \bm{\beta}, \theta) p(\bm{\beta}) p(\theta)) \\
	  &= (-n/2 -\alpha_\pi) \theta - \frac{1}{\exp\{\theta\}} \left(\frac{1}{2} (\bm{y} - \bm{X} \bm{\beta})^\top (\bm{y} - \bm{X} \bm{\beta}) + \beta_\pi \right) \\
	  &\quad - \frac{1}{2} (\bm{\beta} - \bm{\mu}_\pi)^\top \bm{\Sigma}_\pi^{-1} (\bm{\beta} - \bm{\mu}_\pi) \\
	\nabla \mathcal{L}(\bm{\beta}, \theta | \bm{y}) &= \begin{bmatrix}
	  \frac{1}{\exp\{\theta\}} \bm{X}^\top (\bm{y} - \bm{X} \bm{\beta}) - \bm{\Sigma}_\pi^{-1} (\bm{\beta} - \bm{\mu}_\pi) \\
	  -n/2 -\alpha_\pi + \frac{1}{\exp\{\theta\}} \left(\frac{1}{2} (\bm{y} - \bm{X} \bm{\beta})^\top (\bm{y} - \bm{X} \bm{\beta}) + \beta_\pi \right)
	\end{bmatrix}
