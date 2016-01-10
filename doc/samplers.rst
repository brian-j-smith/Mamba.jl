.. _section-Sampling-Functions:

Sampling Functions
==================

Listed below are the sampling methods for which functions are provided to simulating draws from distributions that can be specified up to constants of proportionalities.  Model-based :ref:`section-Sampler` constructors are available for use with the :func:`mcmc` engine as well as stand-alone functions that can be used independently.

.. toctree::
    :maxdepth: 1

    samplers/abc.rst
    samplers/amm.rst
    samplers/amwg.rst
    samplers/bhmc.rst
    samplers/bmc3.rst
    samplers/bmg.rst
    samplers/dgs.rst
    samplers/hmc.rst
    samplers/mala.rst
    samplers/miss.rst
    samplers/nuts.rst
    samplers/rwm.rst
    samplers/slice.rst
    samplers/slicesimplex.rst

The following table summarizes the (*d*-dimensional) sample spaces over which each method simulates draws, whether draws are generated univariately or multivariately, and whether transformations are applied to map parameters to the sample spaces.

.. table:: Summary of sampling methods and their characteristics.

    +--------------------------------------------+---------------------------------------+---------------------------------------------+---------------------------+
    |                                            |                                       | Model-Based Constructors                    | Stand-Alone Functions     |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | Method                                     | Sample Space                          | Univariate | Multivariate | Transformations | Univariate | Multivariate |
    +============================================+=======================================+============+==============+=================+============+==============+
    | :ref:`ABC <section-ABC>`                   | :math:`\mathbb{R}^d`                  | No         | Yes          | Yes             | No         | No           |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`AMM <section-AMM>`                   | :math:`\mathbb{R}^d`                  | No         | Yes          | Yes             | No         | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`AMWG <section-AMWG>`                 | :math:`\mathbb{R}^d`                  | Yes        | No           | Yes             | Yes        | No           |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`BHMC <section-BHMC>`                 | :math:`\{0, 1\}^d`                    | No         | Yes          | No              | No         | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`BMC3 <section-BMC3>`                 | :math:`\{0, 1\}^d`                    | Yes        | Yes          | No              | Yes        | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`BMG <section-BMG>`                   | :math:`\{0, 1\}^d`                    | Yes        | Yes          | No              | Yes        | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`DGS <section-DGS>`                   | Finite :math:`S \subset \mathbb{Z}^d` | Yes        | No           | No              | No         | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`HMC <section-HMC>`                   | :math:`\mathbb{R}^d`                  | No         | Yes          | Yes             | No         | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`MALA <section-MALA>`                 | :math:`\mathbb{R}^d`                  | No         | Yes          | Yes             | No         | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`MISS <section-MISS>`                 | Parameter-defined                     | Yes        | Yes          | No              | No         | No           |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`NUTS <section-NUTS>`                 | :math:`\mathbb{R}^d`                  | No         | Yes          | Yes             | No         | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`RWM <section-RWM>`                   | :math:`\mathbb{R}^d`                  | No         | Yes          | Yes             | No         | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`Slice <section-Slice>`               | :math:`S \subseteq \mathbb{R}^d`      | Yes        | Yes          | Optional        | Yes        | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
    | :ref:`SliceSimplex <section-SliceSimplex>` | *d*-simplex                           | No         | Yes          | No              | No         | Yes          |
    +--------------------------------------------+---------------------------------------+------------+--------------+-----------------+------------+--------------+
