.. index:: Distributions

.. _section-Distributions:

Distributions
-------------

Given in this section are distributions, as provided by the *Distributions* :cite:`bates:2014:DP` and *Mamba* packages, supported for the specification of :ref:`section-Stochastic` nodes.  Truncated versions of continuous univariate distributions are also supported.


.. index:: Distributions; Univariate

Univariate Distributions
^^^^^^^^^^^^^^^^^^^^^^^^

Distributions Package Univariate Types
``````````````````````````````````````

The following `univariate types <http://distributionsjl.readthedocs.org/en/latest/univariate.html>`_ from the *Distributions* package are supported.

.. code-block:: julia

    Arcsine       DiscreteUniform   InverseGamma       NoncentralChisq            SymTriangularDist
    Bernoulli     Edgeworth         InverseGaussian    NoncentralF                TDist
    Beta          Epanechnikov      Kolmogorov         NoncentralHypergeometric   TriangularDist
    BetaPrime     Erlang            KSDist             NoncentralT                Triweight
    Binomial      Exponential       KSOneSided         Normal                     Uniform
    Biweight      FDist             Laplace            NormalCanon                VonMises
    Categorical   Frechet           Levy               Pareto                     Weibull
    Cauchy        Gamma             Logistic           PoissonBinomial
    Chi           Geometric         LogNormal          Poisson
    Chisq         Gumbel            NegativeBinomial   Rayleigh
    Cosine        Hypergeometric    NoncentralBeta     Skellam

.. index:: Distributions; Flat

.. _section-Distribution-Flat:

Flat Distribution
`````````````````

A Flat distribution is supplied with the degenerate probability density function:

.. math::

    f(x) \propto 1, \quad -\infty < x < \infty.

.. code-block:: julia

    Flat()    # Flat distribution

.. index:: Distributions; User-Defined Univariate

User-Defined Univariate Distributions
`````````````````````````````````````

New known, unknown, or unnormalized univariate distributions can be created and added to *Mamba* as subtypes of the *Distributions* package ``ContinuousUnivariateDistribution`` or ``DiscreteUnivariateDistribution`` types.  *Mamba* requires only a partial implementation of the method functions described in the `full instructions for creating univariate distributions <http://distributionsjl.readthedocs.org/en/latest/extends.html#create-a-univariate-distribution>`_.  The specific workflow is given below.

    #. Create a ``quote`` block for the new distribution.  Assign the block a variable name, say ``extensions``, preceded by the ``@everywhere`` macro to ensure compatibility when **julia** is run in multi-processor mode.

    #. The *Distributions* package contains types and method definitions for new distributions.  Load the package and import any of its methods (indicated below) that are extended.

    #. Declare the new distribution subtype, say ``D``, within the block.  Any constructors explicitly defined for the subtype should accept un-typed or abstract-type (``Real``, ``AbstractArray``, or ``DenseArray``) arguments.  Implementing constructors in this way ensures that they will be callable with the *Mamba* ``Stochastic`` and ``Logical`` types.

    #. Extend/define the following *Distributions* package methods for the new distribution ``D``.

        .. function:: minimum(d::D)

            Return the lower bound of the support of ``d``.

        .. function:: maximum(d::D)

            Return the upper bound of the support of ``d``.

        .. function:: logpdf(d::D, x::Real)

            Return the normalized or unnomalized log-density evaluated at ``x``.

    #. Test the subtype.

    #. Add the ``quote`` block (new distribution) to *Mamba* with the following calls.

        .. code-block:: julia

            using Mamba
            @everywhere eval(extensions)

Below is a univariate example based on the linear regression model in the :ref:`section-Line`.

.. literalinclude:: newunivardist.jl
    :language: julia


.. index:: Distributions; Multivariate

Multivariate Distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

Distributions Package Multivariate Types
````````````````````````````````````````

The following `multivariate types <http://distributionsjl.readthedocs.org/en/latest/multivariate.html>`_ from the *Distributions*  package are supported.

.. code-block:: julia

    Dirichlet   MvNormal   MvTDist   Multinomial   MvNormalCanon   VonMisesFisher


.. index:: Distributions; Block-Diagonal Normal

.. _section-Distribution-BDiagNormal:

Block-Diagonal Multivariate Normal Distribution
```````````````````````````````````````````````

A Block-Diagonal Multivariate Normal distribution is supplied with the probability density function:

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

.. index:: Distributions; User-Defined Multivariate

User-Defined Multivariate Distributions
```````````````````````````````````````

New known, unknown, or unnormalized multivariate distributions can be created and added to *Mamba* as subtypes of the *Distributions* package ``ContinuousMultivariateDistribution`` or ``DiscreteMultivariateDistribution`` types.  *Mamba* requires only a partial implementation of the method functions described in the `full instructions for creating multivariate distributions <http://distributionsjl.readthedocs.org/en/latest/extends.html#create-a-multivariate-distribution>`_.  The specific workflow is given below.

    #. Create a ``quote`` block for the new distribution.  Assign the block a variable name, say ``extensions``, preceded by the ``@everywhere`` macro to ensure compatibility when **julia** is run in multi-processor mode.

    #. The *Distributions* package contains types and method definitions for new distributions.  Load the package and import any of its methods (indicated below) that are extended.

    #. Declare the new distribution subtype, say ``D``, within the block.  Any constructors explicitly defined for the subtype should accept un-typed or abstract-type (``Real``, ``AbstractArray``, or ``DenseArray``) arguments.  Implementing constructors in this way ensures that they will be callable with the *Mamba* ``Stochastic`` and ``Logical`` types.

    #. Extend/define the following *Distributions* package methods for the new distribution ``D``.

        .. function:: length(d::D)

            Return the sample space size (dimension) of ``d``.

        .. function:: insupport{T<:Real}(d::D, x::AbstractVector{T})

            Return a logical indicating whether ``x`` is in the support of ``d``.

        .. function:: _logpdf{T<:Real}(d::D, x::AbstractVector{T})

            Return the normalized or unnomalized log-density evaluated at ``x``.

    #. Test the subtype.

    #. Add the ``quote`` block (new distribution) to *Mamba* with the following calls.

        .. code-block:: julia

            using Mamba
            @everywhere eval(extensions)

Below is a multivariate example based on the linear regression model in the :ref:`section-Line`.

.. literalinclude:: newmultivardist.jl
    :language: julia


.. index:: Distributions; Matrix-Variate

Matrix-Variate Distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Distributions Package Matrix-Variate Types
``````````````````````````````````````````

The following `matrix-variate <http://distributionsjl.readthedocs.org/en/latest/matrix.html>`_ types from the *Distributions*  package are supported.

.. code-block:: julia

    InverseWishart   Wishart
