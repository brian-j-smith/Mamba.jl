Variate Types
=============

.. index:: Variate

Variate
-------

Declaration
^^^^^^^^^^^

``abstract Variate{T<:Union(VariateType, VecOrMat{VariateType})}``

Aliases
^^^^^^^

.. code-block:: julia

	typealias VariateType Float64

	typealias VariateScalar Variate{VariateType}
	typealias VariateVector Variate{Vector{VariateType}}
	typealias VariateMatrix Variate{Matrix{VariateType}}
	typealias VariateVecOrMat Union(VariateVector, VariateMatrix)

Fields
^^^^^^

* ``data::T`` :


Supertypes
----------

.. image:: images/variateUML.png
	:align: center

