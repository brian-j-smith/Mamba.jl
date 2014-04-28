Variate Types
=============

.. index:: Variate

Variate
-------

``Variate`` is an abstract type that serves as the basis for several concrete types in the *MCMCsim* package.  Conceptually, it represents a data structure that stores numeric values sampled from a target distribution.  As an abstract type, ``Variate`` cannot be instantiated and cannot have fields.  It can, however, have method functions, which descendant subtypes will inherit.  Such inheritance allows one to endow a core set of functionality to all subtypes by simply defining the functionality once on the abstract type.  Accordingly, a core set of functionality is defined for the ``Variate`` type through the method functions summarized below.  Although the (abstract) type does not have fields, its method functions assume that all subtypes will be declared with the ``value`` field defined below.

Declaration
^^^^^^^^^^^

``abstract Variate{T<:Union(VariateType, VecOrMat{VariateType})}``

Aliases
^^^^^^^

.. index:: VariateType
.. index:: VariateScalar
.. index:: VariateVector
.. index:: VariateMatrix
.. index:: VariateVecOrMat

.. code-block:: julia

	typealias VariateType Float64

	typealias VariateScalar Variate{VariateType}
	typealias VariateVector Variate{Vector{VariateType}}
	typealias VariateMatrix Variate{Matrix{VariateType}}
	typealias VariateVecOrMat Union(VariateVector, VariateMatrix)

Field
^^^^^

* ``value::T`` : a scalar, vector, or matrix of ``VariateType`` values that represent samples from a target distribution.

Methods
^^^^^^^

Method functions supported on all ``Variate`` types are summarized in the following sections; and, unless otherwise specified, are detailed in `The Julia Standard Library <http://docs.julialang.org/en/release-0.2/stdlib/base>`_ documentation.

Array Functions
```````````````

.. code-block:: julia

	cummin      cumsum         maximum     prod
	cummax      cumsum_kbn     minimum     sum
	cumprod     diff           norm        sum_kbn

Collections
```````````

.. code-block:: julia

	endof      size          show
	length     getindex      showcompact
	ndims      setindex!

Distributions
`````````````

The `univariate <http://distributionsjl.readthedocs.org/en/latest/univariate.html>`_, `multivariate <http://distributionsjl.readthedocs.org/en/latest/multivariate.html>`_, and `matrix <http://distributionsjl.readthedocs.org/en/latest/matrix.html>`_ distributions found in the *Distributions* package are supported.

Linear Algebra
``````````````

.. code-block:: julia

	dot

Mathematical Operators and Elementary Functions
```````````````````````````````````````````````

The basic numerical `Mathematical Operators and Elementary Functions <http://julia.readthedocs.org/en/release-0.2/manual/mathematical-operations/>`_ of the **julia** language are supported.

Statistics
``````````

.. code-block:: julia

	cor      median     var
	cov      std        varm
	mean     stdm


Subtypes
----------

Subtypes of ``Variate`` include the :ref:`section-MCMCDependent`, :ref:`section-MCMCLogical`, and :ref:`section-MCMCStochastic` types, as well as the those defined for supplied :ref:`section-Sampling-Functions`.

.. figure:: images/variateUML.png
	:align: center

	UML relational diagram of ``Variate`` types and their fields.
