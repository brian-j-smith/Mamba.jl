.. index:: Sampling Functions; Missing Values Sampler

Missing Values Sampler (MISS)
-----------------------------

A sampler to simulate missing output values from their likelihood distributions.

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: MISS(params::Vector{Symbol})

	Construct an ``MCMCSampler`` object to sampling missing output values.  The constructor should only be used to sample stochastic nodes upon which no other stochastic node depends.  So-called 'output nodes' can be identified with the :func:`keys` function.  Moreover, when the ``MISS`` constructor is included in a vector of ``MCMCSamplers`` to define a sampling scheme, it should be positioned at the beginning of the vector.  This ensures that missing output values are updated before any other samplers are executed.
	
	**Arguments**
	
		* ``params`` : stochastic nodes that contain missing values (``NaN``) to be updated with the sampler.

	**Value**
	
		Returns an ``MCMCSampler`` type object.
		
	**Example**
	
		See the :ref:`Bones <example-Bones>` and :ref:`Mice <example-Mice>` examples.
