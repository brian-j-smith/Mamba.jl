.. index:: MCMCSampler

.. _section-MCMCSampler:

MCMCSampler
-----------

Each of the :math:`\{f_j\}_{j=1}^{B}` sampling functions of the :ref:`figure-Gibbs` is implemented as an ``MCMCSampler`` type object, whose fields are summarized herein.  The ``eval`` field is an anonymous function defined as

.. code-block:: julia

	function(model::MCMCModel, block::Integer)

where ``model`` contains all model nodes and ``block`` is an index identifying the corresponding sampling function in a vector of all samplers for the associated model.  Through the arguments, all model nodes and fields can be accessed in the body of the function.  The function should return an updated sample for the nodes named in its ``params`` field.  The return value can be a structure of the same type of the node if the block consists of only one node, or a dictionary of node structures with names equal to the block nodes if one or more.

Declaration
^^^^^^^^^^^

``type MCMCSampler``

Fields
^^^^^^

* ``params::Vector{String}`` : names of the stochastic nodes in the block being updated by the sampler.
* ``eval::Function`` : a sampling function that updates values of the ``params`` nodes.
* ``tune::Dict`` : any tuning parameters needed by the sampling function.
* ``targets::Vector{String}`` : names of ``MCMCDependent`` nodes that depend on and whose states must be updated after the ``params``.  Elements of ``targets`` are topologically sorted so that a given node in the vector is conditionally independent of subsequent nodes, given the previous ones.

Constructor
^^^^^^^^^^^

.. function:: MCMCSampler(params::Vector{T<:String}, expr::Expr, tune::Dict=Dict())

	Construct a ``MCMCSampler`` object that defines a sampling function for a block of stochastic nodes.
	
	**Arguments**
	
		* ``params`` : names of nodes that are being block-updated by the sampler.
		* ``expr`` : a quoted expression that makes up the body of the sampling function whose definition is described above.
		* ``tune`` : tuning parameters needed by the sampling function.
		
	**Value**
	
		Returns an ``MCMCSampler`` type object.

Methods
^^^^^^^

.. function:: show(s::MCMCSampler)

    Write a text representation of the defined sampling function to the current output stream.

.. function:: showall(s::MCMCSampler)

    Write a verbose text representation of the defined sampling function to the current output stream.
