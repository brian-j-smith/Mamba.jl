.. index:: Sampler Types
.. index:: Sampler Types; Sampler

.. _section-Sampler:

Sampler
-------

The ``Sampler`` type stores model-based :ref:`section-Sampling-Functions` for use in the :ref:`figure-Gibbs`.  Developers can use it as a wrapper for calling stand-alone samplers or as a structure for implementing self-contained samplers.

Declaration
^^^^^^^^^^^

``type Sampler{T}``

Fields
^^^^^^

* ``params::Vector{Symbol}`` : stochastic nodes in the block being updated by the sampler.
* ``eval::Function`` : sampling function that updates values of the ``params`` nodes.
* ``tune::T`` : tuning parameters needed by the sampling function.
* ``targets::Vector{Symbol}`` : ``Dependent`` nodes that depend on and whose states must be updated after ``params``.  Elements of ``targets`` are topologically sorted so that a given node in the vector is conditionally independent of subsequent nodes, given the previous ones.

Constructors
^^^^^^^^^^^^

.. function:: Sampler(param::Symbol, f::Function, tune::Any=Dict())
              Sampler(params::Vector{Symbol}, f::Function, tune::Any=Dict())

    Construct a ``Sampler`` object that defines a sampling function for a block of stochastic nodes.

    **Arguments**

        * ``param/params`` : node(s) being block-updated by the sampler.
        * ``f`` : function for the ``eval`` field of the constructed sampler and whose arguments are the other model nodes upon which the sampler depends, typed argument ``model::Model`` that contains all model nodes, and/or typed argument ``block::Integer`` that is an index identifying the corresponding sampling function in a vector of all samplers for the associated model.  Through the arguments, all model nodes and fields can be accessed in the body of the function.  The function may return an updated sample for the nodes identified in its ``params`` field.  Such a return value can be a structure of the same type as the node if the block consists of only one node, or a dictionary of node structures with keys equal to the block node symbols if one or more.  Alternatively, a value of ``nothing`` may be returned.  Return values that are not ``nothing`` will be used to automatically update the node values and propagate them to dependent nodes.  No automatic updating will be done if ``nothing`` is returned.
        * ``tune`` : tuning parameters needed by the sampling function.

    **Value**

        Returns a ``Sampler{typeof(tune)}`` type object.

    **Example**

        See the :ref:`section-Line-Schemes` section of the tutorial.

Display
^^^^^^^

.. function:: show(s::Sampler)

    Write a text representation of the defined sampling function to the current output stream.

.. function:: showall(s::Sampler)

    Write a verbose text representation of the defined sampling function to the current output stream.


.. index:: Sampler Types; SamplerTune
.. index:: Sampler Types; SamplerVariate

.. _section-SamplerVariate:

SamplerVariate
--------------

The ``SamplerVariate`` type is designed to store simulated values from and tuning parameters for stand-alone :ref:`section-Sampling-Functions`.  It is a parametric type whose parameter can be any subtype of the abstract ``SamplerTune`` type and serves to identify the family of sampling functions to which the variate belongs.

Declaration
^^^^^^^^^^^

.. code-block:: julia

    abstract type SamplerTune
    type SamplerVariate{T<:SamplerTune} <: VectorVariate

Fields
^^^^^^

* ``value::Vector{Float64}`` : simulated values.
* ``tune::T`` : tuning parameters.

Constructors
^^^^^^^^^^^^

.. function:: SamplerVariate(x::AbstractVector{U<:Real}, tune::SamplerTune)
              SamplerVariate{T<:SamplerTune}(x::AbstractVector{U<:Real}, tune::T)
              SamplerVariate{T<:SamplerTune}(x::AbstractVector{U<:Real}, pargs...; kargs...)

    Construct a ``SamplerVariate`` object for storing simulated values and tuning parameters.

    **Arguments**

        * ``x`` : simulated values.
        * ``tune`` : tuning parameters.  If not specified, the tuning parameter constructor is called with the ``value`` field of the variate to instantiate the parameters.
        * ``T`` : explicit tuning parameter type for the variate.  If not specified, the type is inferred from the ``tune`` argument.
        * ``pargs...``, ``kargs...`` : variable positional and keyword arguments that are passed, along with the ``value`` field of the variate, to the tuning parameter constructor as ``T(value, pargs...; kargs...)``.  Accordingly, the arguments that this version of the ``SamplerVariate{T}`` constructor accepts are defined by the ``T`` constructors implemented for it.

    **Value**

    Returns a ``SamplerVariate{T}`` type object with fields set to the supplied ``x`` and tuning parameter values.
