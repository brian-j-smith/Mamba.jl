.. index:: MCMCDependent

.. _section-MCMCDependent:

MCMCDependent
-----------

``MCMCDependent`` is an abstract type designed to store values and attributes of model nodes, including parameters :math:`\theta_1, \ldots, \theta_p` to be simulated via MCMC, functions of the parameters, and likelihood specifications on observed data.  It extends the base ``Variate`` type with method functions defined for the fields summarized below.  Like the type it extends, values are stored in a ``value`` field and can be used with method functions that accept ``Variate`` type objects.

Since parameter values in the ``MCMCNode`` structure are stored as a scalar, vector, or matrix, objects of this type can be created for model parameters of corresponding dimensions, with the choice between the three being user and application-specific.  At one end of the spectrum, a model might be formulated in terms of parameters that are all scalars, with a separate instances of  ``MCMCDependent`` for each one.  At the other end, a formulation might be made in terms of a single parameter vector, with one corresponding instance of ``MCMCDependent``.  Whether to formulate parameters as scalar, vector, or matrix nodes depends on the application at hand.  Vector and matrix formulations should be considered for parameters and data that have multivariate distributions, or are to be used as such in numeric operations and functions.  In other cases, scalar parameterizations may be preferable.  Situations in which parameter vectors are often used include the specification of regression coefficients and random effects.

Declaration
^^^^^^^^^^^

``abstract MCMCDependent{T} <: Variate{T}``

Fields
^^^^^^

* ``value::T`` : a scalar, vector, or matrix of ``VariateType`` values that represent samples from a target distribution.
* ``names::Vector{String}`` : element-specific names corresponding to values in the ``value`` field.
* ``monitor::Vector{Bool}`` : element-specific boolean values indicating whether to include corresponding ``value`` field values in monitored MCMC sampler output.
* ``eval::Function`` : a function for updating the state of the node.
* ``sources::Vector{String}`` : names of other nodes upon whom the values of this one depends.

Methods
^^^^^^^

.. function:: invlink(d::MCMCDependent, x)

	Apply a node-specific inverse-link transformation.  In this method, the link function is defined to be the identity function.  The method function may be redefined for subtypes of ``MCMCDependent`` to implement other link functions. 
	
	**Arguments**
	
		* ``d`` : a node on which a ``link`` transformation method is defined.
		* ``x`` : an object to which to apply the inverse-link transformation.
	
	**Value**
	
		Returns the inverse-link-transformed version of ``x``.

.. function:: link(d::MCMCDependent, x)

	Apply a node-specific link transformation.  In this method, the link function is defined to be the identity function.  The method function may be redefined for subtypes of ``MCMCDependent`` to implement other link functions. 
	
	**Arguments**
	
		* ``d`` : a node on which a ``link`` transformation method is defined.
		* ``x`` : an object to which to apply the link transformation.
	
	**Value**
	
		Returns the link-transformed version of ``x``.

.. function:: logpdf(d::MCMCDependent, transform::Bool=false)

	Evaluate the log-density function for a node.  In this method, no density function is assumed for the node, and a value of 0 is thus returned.  The method function may be redefined for subtypes of ``MCMCDependent`` that have distributional specifications.
	
	**Arguments**
	
		* ``d`` : a node containing values at which to compute the log-density.
		* ``transform`` : whether to evaluate the log-density on the link-transformed scale.
		
	**Value**
	
		The resulting numeric value of the log-density.

.. function:: setmonitor!(d::MCMCDependent, monitor::Union(Bool,Vector{Bool}))

	Specify node elements to be included in monitored MCMC sampler output.
	
	**Arguments**
	
		* ``d`` : a node whose elements contain sampled MCMC values.
		* ``momitor`` : a scalar indicating whether all elements are monitored, or a vector of element-wise indicators.
		
	**Value**
	
		Returns ``d`` with its ``monitor`` field updated to reflect the specified monitoring.

.. function:: show(d::MCMCDependent)

	Write a text representation of nodal values and attributes to the current output stream.  

.. function:: showall(d::MCMCDependent)

	Write a verbose text representation of nodal values and attributes to the current output stream.  


.. index:: MCMCLogical

.. _section-MCMCLogical:

MCMCLogical
-----------

Type ``MCMCLogical`` inherits the fields and method functions from the ``MCMCDependent`` type, and adds the constructors and methods listed below.  It is designed for nodes that are deterministic functions of model parameters and data.  Stored in the field ``eval`` is an anonymous function defined as

.. code-block:: julia

	function(model::MCMCModel)

where ``model`` contains all model nodes.  The function can contain any valid **julia** expression or code block written in terms of other nodes and data structures.  It should return values with which to update the node in the same type as the ``value`` field of the node.

Declaration
^^^^^^^^^^^

``type MCMCLogical{T} <: MCMCDependent{T}``

Fields
^^^^^^

* ``value::T`` : a scalar, vector, or matrix of ``VariateType`` values that represent samples from a target distribution.
* ``names::Vector{String}`` : element-specific names corresponding to values in the ``value`` field.
* ``monitor::Vector{Bool}`` : element-specific boolean values indicating whether to include corresponding ``value`` field values in monitored MCMC sampler output.
* ``eval::Function`` : a function for updating values stored in ``value``.
* ``sources::Vector{String}`` : names of other nodes upon whom the values of this one depends.

Constructors
^^^^^^^^^^^^

.. function:: MCMCLogical(expr::Expr, monitor::Union(Bool,Vector{Bool})=true)
              MCMCLogical(length::Integer, expr::Expr, monitor::Union(Bool,Vector{Bool})=true)
              MCMCLogical(m::Integer, n::Integer, expr::Expr, monitor::Union(Bool,Vector{Bool})=true)

	Construct an ``MCMCLogical`` object that defines a logical model node.
	
	**Arguments**
	
		* ``length`` : number of vector elements in the node.
		* ``m``, ``n`` : number of matrix rows and columns in the node.
		* ``expr`` : a quoted expression or code-block defining the body of the function stored in the ``eval`` field.
		* ``monitor`` : a scalar indicating whether all elements are monitored, or a vector of element-wise indicators.
		
	**Value**
	
		Returns an ``MCMCLogical{VariateType}`` type object if no dimensional arguments are specified, a ``MCMCLogical{Vector{VariateType}}`` if ``length`` is specified, and a ``MCMCLogcial{Matrix{VariateType}}`` if ``m`` and ``n`` are specified.

Methods
^^^^^^^

.. function:: setinits!(l::MCMCLogical, m::MCMCModel, x=nothing)

	Set initial values for a logical node.
	
	**Arguments**
	
		* ``l`` : a logical node to assign initial values.
		* ``m`` : a model that contains the node.
		* ``x`` : unused.
		
	**Value**
	
		Returns the result of a call to `update!(l, m)``.

.. function:: update!(l::MCMCLogical, m::MCMCModel)

	Update the values of a logical node according to its relationship with others in a model.
	
	**Arguments**
	
		* ``l`` : a logical node to update.
		* ``m`` : a model that contains the node.
		
	**Value**
	
		Returns the node with its values updated.


.. index:: MCMCStochastic

.. _section-MCMCStochastic:

MCMCStochastic
--------------

Type ``MCMCStochastic`` inherits the fields and method functions from the ``MCMCDependent`` type, and adds the additional ones listed below.  It is designed for model parameters or data that have distributional or likelihood specifications, respectively.  Its stochastic relationship to other nodes and data structures is represented by the ``Distributions`` structure stored in field ``distr``.  Stored in the field ``eval`` is an anonymous function defined as

.. code-block:: julia

	function(model::MCMCModel)

where ``model`` contains all model nodes.  The function can contain any valid **julia** expression or code-block.  It should return a single `Distributions <http://distributionsjl.readthedocs.org/en/latest/index.html>`_ object for all node elements or a structure of the same type as the node with element-specific `Distributions` objects.

Declaration
^^^^^^^^^^^

``type MCMCStochastic{T} <: MCMCDependent{T}``

Fields
^^^^^^

* ``value::T`` : a scalar, vector, or matrix of ``VariateType`` values that represent samples from a target distribution.
* ``names::Vector{String}`` : element-specific names corresponding to values in the ``value`` field.
* ``monitor::Vector{Bool}`` : element-specific boolean values indicating whether to include corresponding ``value`` field values in monitored MCMC sampler output.
* ``eval::Function`` : a function for updating the ``distr`` field for the node.
* ``sources::Vector{String}`` : names of other nodes upon whom the distributional specification for this one depends.
* ``distr::DistributionStruct`` : the distributional specification for the node.

Aliases
^^^^^^^

.. code-block:: julia

	typealias DistributionStruct Union(Distribution, VecOrMat{Distribution})

Constructors
^^^^^^^^^^^^

.. function:: MCMCStochastic(expr::Expr, monitor::Union(Bool,Vector{Bool})=true)
              MCMCStochastic(length::Integer, expr::Expr, monitor::Union(Bool,Vector{Bool})=true)
              MCMCStochastic(m::Integer, n::Integer, expr::Expr, monitor::Union(Bool,Vector{Bool})=true)

	Construct an ``MCMCStochastic`` object that defines a stochastic model node.
	
	**Arguments**
	
		* ``length`` : number of vector elements in the node.
		* ``m``, ``n`` : number of matrix rows and columns in the node.
		* ``expr`` : a quoted expression or code-block defining the body of the function stored in the ``eval`` field.
		* ``monitor`` : a scalar indicating whether all elements are monitored, or a vector of element-wise indicators.
		
	**Value**
	
		Returns an ``MCMCStochastic{VariateType}`` type object if no dimensional arguments are specified, a ``MCMCStochastic{Vector{VariateType}}`` if ``length`` is specified, and a ``MCMCStochastic{Matrix{VariateType}}`` if ``m`` and ``n`` are specified.

Methods
^^^^^^^

.. function:: insupport(s::MCMCStochastic)

	Check whether stochastic node values are within the support of its distribution.
	
	**Arguments**
	
		* ``s`` : a stochastic node on which to perform the check.
		
	**Value**
	
		Returns ``true`` if all values are within the support, and ``false`` otherwise.

.. function:: invlink(s::MCMCStochastic, x)

	Apply an inverse-link transformation to map transformed values back to the original distributional scale of a stochastic node.
	
	**Arguments**
	
		* ``s`` : a stochastic node on which a ``link`` transformation method is defined.
		* ``x`` : an object to which to apply the inverse-link transformation.
	
	**Value**
	
		Returns the inverse-link-transformed version of ``x``.

.. function:: link(s::MCMCStochastic, x)

	Apply a link transformation to map values in a constrained distributional support to an unconstrained space. 
	
	**Arguments**
	
		* ``s`` : a stochastic node on which a ``link`` transformation method is defined.
		* ``x`` : an object to which to apply the link transformation.
	
	**Value**
	
		Returns the link-transformed version of ``x``.

.. function:: logpdf(s::MCMStochastic, transform::Bool=false)

	Evaluate the log-density function for a stochastic node.
	
	**Arguments**
	
		* ``s`` : a stochastic node containing values at which to compute the log-density.
		* ``transform`` : whether to evaluate the log-density on the link-transformed scale.
		
	**Value**
	
		The resulting numeric value of the log-density.

.. function:: setinits!(s::MCMCStochastic, m::MCMCModel, x=nothing)

	Set initial values for a stochastic node.
	
	**Arguments**
	
		* ``s`` : a stochastic node to assign initial values.
		* ``m`` : a model that contains the node.
		* ``x`` : values to assign to the node.
		
	**Value**
	
		Returns the node with its assigned initial values.

.. function:: update!(s::MCMCStochastic, m::MCMCModel)

	Update the values of a stochastic node according to its relationship with others in a model.
	
	**Arguments**
	
		* ``s`` : a stochastic node to update.
		* ``m`` : a model that contains the node.
		
	**Value**
	
		Returns the node with its values updated.
