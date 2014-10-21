.. index:: Distributions

.. _section-Distributions:

Distributions
-------------

Listed in this section are distributions, as provided by the *Distributions* :cite:`bates:2014:DP` and *Mamba* packages, supported for the specification of :ref:`section-Stochastic` nodes.  Truncated versions of the continuous univariate distributions are also supported.


Univariate
^^^^^^^^^^

`Distributions Package (Univariate) <http://distributionsjl.readthedocs.org/en/latest/univariate.html>`_
````````````````````````````````````````````````````````````````````````````````````````````````````````

.. code-block:: julia

	Arcsine       Cosine            Hypergeometric    NegativeBinomial   Rayleigh
	Bernoulli     DiscreteUniform   InverseGamma      NoncentralBeta     Skellam
	Beta          Edgeworth         InverseGaussian   NoncentralChisq    TDist
	BetaPrime     Erlang            KSDist            NoncentralF        TriangularDist
	Binomial      Exponential       KSOneSided        NoncentralT        Uniform
	Categorical   FDist             Laplace           Normal             Weibull
	Cauchy        Gamma             Levy              NormalCanon     
	Chi           Geometric         Logistic          Pareto
	Chisq         Gumbel            LogNormal         Poisson

.. index:: Distributions; Flat

.. _section-Distribution-Flat:

Flat Distribution
`````````````````

The Flat distribution has the degenerate probability density function:

.. math::

	f(x) \propto 1, \quad -\infty < x < \infty.
	
.. code-block:: julia

	Flat()    # Flat distribution


Multivariate
^^^^^^^^^^^^

`Distributions Package (Multivariate) <http://distributionsjl.readthedocs.org/en/latest/multivariate.html>`_
````````````````````````````````````````````````````````````````````````````````````````````````````````````

.. code-block:: julia

	MvNormal     Dirichlet       DiagNormalCanon   DiagTDist
	DiagNormal   Multinomial     IsoNormalCanon    IsoTDist
	IsoNormal    MvNormalCanon   MvTDist           VonMisesFisher


.. index:: Distributions; Block-Diagonal Normal

.. _section-Distribution-BDiagNormal:

Block-Diagonal Multivariate Normal Distribution
```````````````````````````````````````````````

The	Block-Diagonal Multivariate Normal distribution has probability density function:

.. math::

	f(\bm{x}; \bm{\mu}, \bm{\Sigma}) = \frac{1}{(2 \pi)^{d/2} |\bm{\Sigma}|^{1/2}} \exp\left(-\frac{1}{2} (\bm{x} - \bm{\mu})^\top \bm{\Sigma}^{-1} (\bm{x} - \bm{\mu})\right), \quad -\infty < \bm{x} < \infty,
	
where

.. math::

	\bm{\Sigma} = \begin{bmatrix}
		\bm{\Sigma_1} & \bm{0} & \hdots & \bm{0} \\
		\bm{0} & \bm{\Sigma_2} & \hdots & \bm{0} \\
		\vdots & \vdots & \ddots & \vdots \\
		\bm{0} & \bm{0} & \hdots & \bm{\Sigma}_m \\
	\end{bmatrix}.

.. code-block:: julia

	BDiagNormal(mu, C)    # multivariate normal with mean vector mu and block-
	                      # diagonal covariance matrix Sigma such that
	                      # length(mu) = dim(Sigma), and Sigma_1 = ... = Sigma_m = C
	                      # for a matrix C or Sigma_1 = C[1], ..., Sigma_m = C[m]
	                      # for a vector of matrices C.


Matrix
^^^^^^

`Distributions Package (Matrix) <http://distributionsjl.readthedocs.org/en/latest/matrix.html>`_
````````````````````````````````````````````````````````````````````````````````````````````````

.. code-block:: julia

	InverseWishart   Wishart