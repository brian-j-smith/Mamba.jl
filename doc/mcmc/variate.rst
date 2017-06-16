.. index:: Variate Types
.. index:: Variate Types; AbstractVariate
.. index:: Variate Types; ScalarVariate
.. index:: Variate Types; ArrayVariate
.. index:: Variate Types; VectorVariate
.. index:: Variate Types; MatrixVariate

.. _section-Variate:

Variate
-------

``ScalarVariate`` and ``ArrayVariate{N}`` are abstract types that serve as the basis for several concrete types in the *Mamba* package.  Conceptually, they represent data structures that store numeric values simulated from target distributions.  Being abstract, these variate types cannot be instantiated and cannot have fields.  They can, however, have method functions, which descendant subtypes will inherit.  Such inheritance allows one to endow a core set of functionality to all subtypes by simply defining method functions once on the abstract types (see `julia Types <http://docs.julialang.org/en/latest/manual/types/>`_).  Accordingly, a core set of functionality is defined for the variate types through the field and method functions discussed below.  Although the (abstract) types do not have fields, their method functions assume that all subtypes will be declared with a ``value`` field.

Declarations
^^^^^^^^^^^^^
.. code-block:: julia

    abstract type ScalarVariate <: Real
    abstract type ArrayVariate{N} <: DenseArray{Float64, N}

    const AbstractVariate = Union{ScalarVariate, ArrayVariate}
    const VectorVariate = ArrayVariate{1}
    const MatrixVariate = ArrayVariate{2}

Type Hierarchy
^^^^^^^^^^^^^^

Subtypes of the variate types include the :ref:`section-Dependent`, :ref:`section-Logical`, :ref:`section-Stochastic`, and :ref:`section-SamplerVariate` types.

.. figure:: ../images/variateUML.png
    :align: center

    UML relational diagram of ``Variate`` types and their fields.

Field
^^^^^

* ``value::T`` : scalar or array of ``Float64`` values that represent simulated values from a target distribution.

Methods
^^^^^^^
Methods for ``ScalarVariate`` and ``ArrayVariate`` include `mathematical operators <http://docs.julialang.org/en/latest/stdlib/math/#mathematical-operators>`_, `mathematical functions <http://docs.julialang.org/en/latest/stdlib/math/#mathematical-functions>`_, and `statistics <http://docs.julialang.org/en/latest/stdlib/math/#statistics>`_ defined in the base **julia** language for parent types ``Real`` and ``DenseArray``.  In addition, the following functions are provided.

=============== ================
Function        Description
=============== ================
``logit(x)``    log-odds
``invlogit(x)`` inverse log-odds
=============== ================
