.. index:: Sampling Functions; Direct Grid Sampler

Direct Grid Sampler (DGS)
-------------------------

Implementation of a sampler for the simulation of draws from discrete univariate distributions with finite supports.  Draws are simulated directly from the individual full conditional probability mass functions for each of the specified model parameters.

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: DGS(params::Vector{Symbol})

	Construct a ``Sampler`` object for direct grid sampling.  Parameters are assumed to have discrete uniform distributions with finite supports.
	
	**Arguments**
	
		*  ``params`` : stochastic nodes to be updated with the sampler.

	**Value**
	
		Returns a ``Sampler`` type object.

	**Example**
	
		See the :ref:`Eyes <example-Eyes>` example.
