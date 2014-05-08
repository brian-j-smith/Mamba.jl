.. index:: Sampling Functions; Missing Values Sampler

Missing Values Sampler (MISS)
-----------------------------

A sampler to simulate missing output values.

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: MISS(params::Vector{T<:String})

	Construct an ``MCMCSampler`` object to sampling missing output values.
	
	**Arguments**
	
		* ``params`` : a single, named stochastic node that contains at least one missing value (``NaN``) to be updated with the sampler.

	**Value**
	
		Returns an ``MCMCSampler`` type object.
