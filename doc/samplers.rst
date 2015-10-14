.. _section-Sampling-Functions:

Sampling Functions
==================

Listed below are the sampling methods for which functions are provided to simulating draws from distributions that can be specified up to constants of proportionalities.  Stand-alone functions are available as well as model-based :ref:`section-Sampler` constructors for use with the :func:`mcmc` engine.

.. toctree::
    :maxdepth: 1

    samplers/amm.rst
    samplers/amwg.rst
    samplers/bmmg.rst
    samplers/dgs.rst
    samplers/mala.rst
    samplers/miss.rst
    samplers/nuts.rst
    samplers/slice.rst
    samplers/slicesimplex.rst

The following table summarizes the (*d*-dimensional) sample spaces over which each method simulates draws, whether draws are generated univariately or multivariately, and whether transformations are applied to map parameters to the sample spaces.

.. table:: Summary of sampling methods and their characteristics.

    +--------------------------------------------+---------------------------------------+---------------------------+---------------------------------------------+
    |                                            |                                       | Stand-Alone Functions     | Model-Based Functions                       |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | Method                                     | Sample Space                          | Univariate | Multivariate | Univariate | Multivariate | Transformations |
    +============================================+=======================================+============+==============+============+==============+=================+
    | :ref:`AMM <section-AMM>`                   | :math:`\mathbb{R}^d`                  | No         | Yes          | No         | Yes          | Yes             |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | :ref:`AMWG <section-AMWG>`                 | :math:`\mathbb{R}^d`                  | Yes        | No           | Yes        | No           | Yes             |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | :ref:`BMMG <section-BMMG>`                 | :math:`\{0, 1\}^d`                    | Yes        | Yes          | Yes        | Yes          | No              |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | :ref:`DGS <section-DGS>`                   | Finite :math:`S \subset \mathbb{Z}^d` | No         | Yes          | Yes        | No           | No              |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | :ref:`MALA <section-MALA>`                 | :math:`\mathbb{R}^n`                  | No         | Yes          | No         | Yes          | Yes             |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | :ref:`MISS <section-MISS>`                 | Parameter-defined                     | No         | No           | Yes        | Yes          | No              |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | :ref:`NUTS <section-NUTS>`                 | :math:`\mathbb{R}^d`                  | No         | Yes          | No         | Yes          | Yes             |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | :ref:`Slice <section-Slice>`               | :math:`S \subseteq \mathbb{R}^d`      | Yes        | Yes          | Yes        | Yes          | Optional        |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
    | :ref:`SliceSimplex <section-SliceSimplex>` | *d*-simplex                           | No         | Yes          | No         | Yes          | No              |
    +--------------------------------------------+---------------------------------------+------------+--------------+------------+--------------+-----------------+
