.. index:: MCMCModel

.. _section-MCMCModel:

MCMCModel
---------

The ``MCMCModel`` type is designed to store the set of all model nodes, including parameter set :math:`\Theta` as denoted in  the :ref:`figure-Gibbs`.  In particular, it stores ``MCMCDepNode`` type objects in its ``nodes`` dictionary field.  Valid models are ones whose nodes form directed acyclic graphs (DAGs).  Sampling functions :math:`\{f_j\}_{j=1}^{B}` are saved as ``MCMCSampler`` objects in the vector of field ``samplers``.  Vector elements :math:`j=1,\ldots,B` correspond to sampling blocks :math:`\{\Theta_j\}_{j=1}^{B}`.

Declaration
^^^^^^^^^^^

``type MCMCModel``

Fields
^^^^^^

* ``nodes::Dict{String,Any}`` : a dictionary containing all input, logical, and stochastic model nodes.
* ``links::Vector{String}`` : names of all ``MCMCDepNode`` nodes in topologically sorted order so that a given node in the vector is conditionally independent of subsequent nodes, given the previous ones.
* ``samplers::Vector{MCMCSampler}`` : sampling functions for updating blocks of stochastic nodes.
* ``iter::Integer`` : current MCMC draw from the target distribution.
* ``burnin::Integer`` : number of initial draws to discard as a burn-in sequence to allow for convergence.
* ``chain::Integer`` : current run of the MCMC simulator in a possible sequence of runs.
* ``hasinputs::Bool`` : whether values have been assigned to the input nodes.
* ``hasinits::Bool`` : whether initial values have been assigned to stochastic nodes.

Constructor
^^^^^^^^^^^

.. function:: MCMCModel(; iter::Integer=0, burnin::Integer=0, chain::Integer=1, \
				samplers::Vector{MCMCSampler}=Array(MCMCSampler, 0), nodes...)
				
	Construct an ``MCMCModel`` object that defines a model for MCMC simulation.
	
	**Arguments:**
	
		* ``iter`` : current iteration of the MCMC simulation.
		* ``burnin`` : number of initial draws to be discarded as a burn-in sequence to allow for convergence.
		* ``chain`` : current run of the MCMC simulator in a possible sequence of runs.
		* ``samplers`` : a vector of block-specific sampling functions.
		* ``nodes...`` : a variable number of user specified arguments defining logical and stochastic nodes in the model.  Argument values must be ``MCMCLogical`` or ``MCMCStochastic`` type objects.  Node names in the model will be taken from the argument names.
		
	**Value:**
	
		Returns an ``MCMCModel`` type object.

Methods
^^^^^^^

.. function:: getindex(m::MCMCModel, key::String)

	Returns a named model node.  The syntax ``m[key]`` is converted to ``getindex(m, key)``.
	
	**Arguments:**
	
		* ``m`` : a model contining the node to get.
		* ``key`` : name of the node to get.
		
	**Value:**
	
		The specified node.
	
.. function:: gradient(m::MCMCModel, block::Integer=0, transform::Bool=false, \
				dtype::Symbol=:central)
              gradient(m::MCMCModel, x::Vector{T<:Real}, block::Integer=0, \
				transform::Bool=false, dtype::Symbol=:central)
			
	Numerically approximate the gradient for stochastic nodes.
	
	**Arguments:**
	
		* ``m`` : a model containing the stochastic nodes for which to compute the gradient.
		* ``x`` : a value (other than the current one) at which to compute the gradient.
		* ``block`` : the sampling block of stochastic nodes for which to compute the gradient, if specified; otherwise, all sampling blocks are included.
		* ``transform`` : whether to compute the gradient on the link–transformed scale.
		* ``dtype`` : type of numerical approximation to use.  Options are
			* ``:central`` : central differencing.
			* ``:forward`` : forward differencing.
		
	**Value:**
	
		The resulting gradient vector.

.. function:: graph(m::MCMCModel)

	Construct a graph representation of model nodes and their relationships.
	
	**Arguments:**
	
		* ``m`` : a model for which to construct a graph.
	
	**Value:**
	
		Returns a ``GenericGraph`` type object as defined in the `Graphs <http://graphsjl-docs.readthedocs.org/en/latest/index.html>`_ package.

.. function:: graph2dot(m::MCMCModel)
              graph2dot(m::MCMCModel, filename::String)

	Construct a `GraphViz <http://www.graphviz.org/>`_ DOT-formatted graph representation of model nodes and their relationships.
	
	**Arguments:**
	
		* ``m`` : a model for which to construct a graph.
		* ``filename`` : an external file to which to save the resulting graph.
	
	**Value:**
	
		A character string represenation of the graph in DOT format.
	
.. function:: keys(m::MCMCModel, ntype::Symbol=:assigned, block::Integer=0)

	Return names of node of a specified type.
	
	**Arguments:**
	
		* ``m`` : a model containing the nodes of interest.
		* ``ntype`` : the type of node to return.  Options are
			* ``:all`` : all input, logical, and stochastic model nodes.
			* ``:assigned`` : nodes that have been assigned values.
			* ``:block`` : stochastic nodes being block-sampled.
			* ``:dep`` : logical or stochastic (dependent) nodes.
			* ``:indep`` or ``:input`` : input (independent) nodes.
			* ``:logical`` : logical nodes.
			* ``:monitor`` : stochastic nodes being monitored in MCMC sampler output.
			* ``:terminal`` : stochastic nodes upon which no other stochastic nodes depend.
			* ``:stochastic`` : stochastic nodes.
		* ``block`` : the block for which to return nodes if ``ntype = :block``, or all blocks if ``block = 0`` (default).
		
	**Value:**
	
		Node names returned as a vector of character strings.

.. function:: logpdf(m::MCMCModel, block::Integer=0, transform::Bool=false)
              logpdf(m::MCMCModel, x::Vector{T<:Real}, block::Integer=0, \
				transform::Bool=false)

	Compute the sum of log-densities for stochastic nodes.
	
	**Arguments:**
	
		* ``m`` : a model containing the stochastic nodes for which to evaluate log-densities.
		* ``x`` : a value (other than the current one) at which to evaluate densities.
		* ``block`` : the sampling block of stochastic nodes over which to sum densities, if specified; otherwise, all sampling blocks are included.
		* ``transform`` : whether to evaluate evaluate log-densities on the link–transformed scale.
		
	**Value:**
	
		The resulting numeric value of summed log-densities.
				
.. function:: mcmc(model::MCMCModel, inputs::Dict{T<:String}, \
				inits::Vector{Dict{U<:String,Any}}, iter::Integer; \
				burnin::Integer=0, thin::Integer=1, chains::Integer=1)

	Simulate MCMC draws for a specified model.
	
	**Arguments:**
	
		* ``model`` : a specified mode.
		* ``inputs`` : a dictionary of values for input model nodes.  Dictionary keys and values should be given for each input node.
		* ``inits`` : a vector of dictionaries that contain initial values for stochastic model nodes.  Dictionary keys and values should be given for each stochastic node.  Consecutive runs of the simulator will iterate through the vector's dictionary elements.
		* ``iter`` : number of draws to generate for each simulation run.
		* ``burnin`` : numer of initial draws to discard as a burn-in sequence to allow for convergence.
		* ``thin`` : step-size between draws to output.
		* ``chains`` : number of simulation runs to perform.
		
	**Value:**
	
		An ``MCMCChains`` type object of simulated draws.
		
.. function:: relist(m::MCMCModel, values::Vector{T<:Real}, block::Integer=0, \
				transform::Bool=false)
              relist(m::MCMCModel, values::Vector{T<:Real}, nkeys::Vector{U<:String}, \
				transform::Bool=false)
				
	Convert a vector of values to a set of logical and/or stochastic node values.

	**Arguments:**
	
		* ``m`` : a model with nodes to serve as the template for conversion.
		* ``values`` : values to convert.
		* ``block`` : the sampling block of nodes to which to convert ``values``.  Defaults to all blocks.
		* ``nkeys`` : a vector of names specifying the nodes to which to convert ``values``.
		* ``transform`` : whether to apply an inverse-link transformation in the conversion.
		
	**Value:**
	
		A dictionary of node names and converted values.

.. function:: relist!(m::MCMCModel, values::Vector{T<:Real}, block::Integer=0, \
				transform::Bool=false)
              relist!(m::MCMCModel, values::Vector{T<:Real}, nkeys::Vector{U<:String}, \
				transform::Bool=false)
				
	Copy a vector of values to a set of logical and/or stochastic nodes.
	
	**Arguments:**
	
		* ``m`` : a model with nodes to which values will be copied.
		* ``values`` : values to copy.
		* ``block`` : the sampling block of nodes to which to copy ``values``.  Defaults to all blocks.
		* ``nkeys`` : a vector of names specifying the nodes to which to copy ``values``.
		* ``transform`` : whether to apply an inverse-link transformation in the copy.
		
	**Value:**
	
		Returns the model with copied node values.
							
.. function:: setinits!(m::MCMCModel, inits::Dict{T<:String,Any})

	Set the initial values of stochastic model nodes.
	
	**Arguments:**
	
		* ``m`` : a model with nodes to be initialized.
		* ``inits`` : a dictionary of initial values for stochastic model nodes.  Dictionary keys and values should be given for each stochastic node.
		
	**Value:**
	
		Returns the model with initialized stochastic nodes.

.. function:: setinputs!(m::MCMCModel, inputs::Dict{T<:String,Any})

	Set the values of input model nodes.
	
	**Arguments:**
	
		* ``m`` : a model with input nodes to be assigned.
		* ``inputs`` : a dictionary of values for input model nodes.  Dictionary keys and values should be given for each input node.
		
	**Value:**
	
		Returns the model with values assigned to input nodes.

.. function:: setsamplers!(m::MCMCModel, samplers::Vector{MCMCSampler})

	Set the block-samplers for stochastic model nodes.
	
	**Arguments:**
	
		* ``m`` : a model with stochastic nodes to be sampled.
		* ``samplers`` : block-specific samplers.
		
	**Values:**
	
		Returns the model updated with the block-samplers.

.. function:: show(m::MCMCModel)

	Write a text representation of the model, nodes, and attributes to the current output stream.

.. function:: showall(m::MCMCModel)

	Write a verbose text representation of the model, nodes, and attributes to the current output stream.

.. function:: simulate!(m::MCMCModel, block::Integer=0)

	Simulate one MCMC draw from a specified model.
	
	**Argument:**
	
		* ``m`` : a model specification.
		* ``block`` : the block for which to simulate an MCMC draw, if specified; otherwise, simulate draws for all blocks (default).
		
	**Value:**
	
		Returns the model updated with the MCMC draw.

.. function:: tune(m::MCMCModel, block::Integer=0)

	Get block-sampler tuning parameters.
	
	**Arguments:**
	
		* ``m`` : a model with block-samplers.
		* ``block`` : the block for which to return the tuning parameters, if specified; otherwise, the tuning parameters for all blocks.
		
	**Value:**
	
		If ``block = 0``, a vector of dictionaries containing block-specific tuning parameters; otherwise, one block-specific dictionary.

.. function:: unlist(m::MCMCModel, block::Integer=0, transform::Bool=false)
              unlist(m::MCMCModel, nkeys::Vector{T<:String}, transform::Bool=false)
			  
	Convert a set of logical and/or stochastic node values to a vector.
	
	**Arguments:**
	
		* ``m`` : a model with nodes to be converted.
		* ``block`` : the sampling block of nodes to be converted.  Defaults to all blocks.
		* ``nkeys`` : a vector of names specifying the nodes to be converted.
		* ``transform`` : whether to apply a link transformation in the conversion.
		
	**Value:**
	
		A vector of concatenated node values.

.. function:: update!(m::MCMCModel, block::Integer=0)

	Update values of logical and stochastic model node according to their relationship with others in a model.
	
	**Arguments:**
	
		* ``m`` : a mode with nodes to be updated.
		* ``block`` : the sampling block of nodes to be updated.  Defaults to all blocks.
		
	**Value:**
	
		Returns the model with updated nodes.
