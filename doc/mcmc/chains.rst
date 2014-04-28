.. index:: MCMCChains

.. _section-MCMCChains:

MCMCChains
----------

The ``MCMCChains`` type stores MCMC sampler output.

Declaration
^^^^^^^^^^^

``immutable MCMCChains``

Fields
^^^^^^

* ``value::Array{VariateType,3}`` : a 3-dimensional array of sampled values whose first, second, and third dimensions index the iterations, parameter elements, and runs of an MCMC sampler, respectively.
* ``names::Vector{String}`` : names assigned to the parameter elements.
* ``start::Integer`` : number of the iteration stored in the first row of the ``value`` array.
* ``thin::Integer`` : number of steps between consecutive iterations stored in the ``value`` array.
* ``model::MCMCModel`` : the model from which the sampled values were generated.

Constructors
^^^^^^^^^^^^

.. function:: MCMCChains(value::Array{T<:Real,2}, names::Vector{U<:String}; \
                start::Integer=1, thin::Integer=1, model::MCMCModel=MCMCModel())
		      MCMCChains(value::Array{T<:Real,3}, names::Vector{U<:String}; \
                start::Integer=1, thin::Integer=1, model::MCMCModel=MCMCModel())
              MCMCChains(iter::Integer, names::Vector{T<:String}; start::Integer=1, \
                thin::Integer=1, chains::Integer=1, model::MCMCModel=MCMCModel())
              
		
	Construct an ``MCMCChains`` object that stores MCMC sampler output.
	
	**Arguments**
	
		* ``value`` : simulated values whose first, second, and third (optional) dimensions index the iterations, parameter elements, and runs of an MCMC sampler, respectively.
		* ``iter`` : number of simulation-specific iterations to store.
		* ``names`` : names to assign to the parameter elements.
		* ``start`` : number of the first iteration to be stored.
		* ``thin`` : number of steps between consecutive iterations to be stored.
		* ``chains`` : number of simulations runs for which to store output.
		* ``model`` : the model for which the simulation was run.
		
	**Value**
	
		Returns an ``MCMCChains`` type object.

Methods
^^^^^^^

.. function:: autocor(c::MCMCChains; lags::Vector=[1,5,10,50], relative::Bool=true)

	Compute lag-k autocorrelation for MCMC sampler output.
	
	**Arguments**
	
		* ``c`` : sampler output on which to perform calculations.
		* ``lags`` : lags at which to compute autocorrelations.
		* ``relative`` : whether the lags are relative to the thinning interval of the output (``true``) or relative to the absolute iteration numbers (``false``).
		
	**Value**
	
		A ``ChainSummary`` type object of the form:
		
		.. index:: ChainSummary
		
		.. code-block:: julia
		
			immutable ChainSummary
			  value::Array{Float64,3}
			  rownames::Vector{String}
			  colnames::Vector{String}
			  header::String
			end

		with model parameters indexed by the first dimension of ``value``, lag-autocorrelations by the second, and chains by the third.

.. function:: cor(c::MCMCChains)

	Compute cross-correlations for MCMC sampler output.
	
	**Arguments**
	
		* ``c`` : sampler output on which to perform calculations.
		
	**Value**
	
		A ``ChainSummary`` type object with the first and second dimensions of the ``value`` field indexing the model parameters between which correlations.  Results are for all chains combined.

.. function:: describe(c::MCMCChains; batchsize::Integer=100, \
				q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])
				
	Compute summary statistics for MCMC sampler output.
	
	**Arguments**
	
		* ``c`` : sampler output on which to perform calculations.
		* ``batchsize`` : number of iterations to include in a partitioning of the output for calculation of batch standard errors.
		* ``q`` : probabilities at which to calculate quantiles.
		
	**Value**
	
		A tuple of results from calls to ``summarystats(c, batchsize)`` and ``quantile(c, q)``, respectively.    Results are for all chains combined.

.. function:: dic(c::MCMCChains)

	Compute the Deviance Information Criterion (DIC) of Spiegelhalter et al. :cite:`spiegelhalter:2002:BMM` and Gelman et al. :cite:`gelman:2013:bda` from MCMC sampler output.
	
	**Arguments**
	
		* ``c`` : sampler output on which to perform calculations.
		
	**Value**
	
		A ``ChainSummary`` type object with DIC results from the methods of Spiegelhalter and Gelman in the first and second rows of the ``value`` field, and the DIC value and effective numbers of parameters in the first and second columns.  Results are for all chains combined.

.. function:: gelmandiag(c::MCMCChains; alpha::Real=0.05, mpsrf::Bool=false, \
				transform::Bool=false)
	
	Compute the convergence diagnostic of Brooks, Gelman, and Rubin :cite:`brooks:1998:GMM`, :cite:`gelman:1992:IIS` for MCMC sampler output.
	
	**Arguments**
	
		* ``c`` : sampler output on which to perform calculations.
		* ``alpha`` : quantile (``1 - alpha / 2``) at which to estimate the upper limits of shrink factors.
		* ``mpsrf`` : whether to compute the multivariate potential scale reduction factor.
		
	**Value**
	
		A ``ChainSummary`` type object with parameters contained in the rows of the ``value`` field, and scale reduction factors and upper-limit quantiles in the first and second columns.

.. function:: getindex(c::MCMCChains, inds...)

	Subset MCMC sampler output.  The syntax ``c[i, j, k]`` is converted to ``getindex(c, i, j, k)``.
	
	**Arguments**
	
		* ``c`` : sampler output to subset.
		* ``inds...`` : a tuple of ``i, j, k`` indices to the iterations, parameters, and chains to be subsetted.  Indices of the form ``start:stop`` or ``start:thin:stop`` can be used to subset iterations, where ``start`` and ``stop`` define a range for the subset and ``thin`` will apply additional thinning to existing sampler output.  Indices for subsetting of parameters can be specified as a vector of strings, integers, or booleans identifying parameters to be kept.  Indices for chains can be a vector of integers or booleans.  A value of ``:`` can be specified for any of the dimensions to indicate no subsetting.
		
	**Value**
	
		Returns an ``MCMCChains`` object with the subsetted sampler output.
		
.. function:: hpd(c::MCMCChains; alpha::Real=0.05)

	Compute highest posterior density (HPD) intervals of Chen and Shao :cite:`chen:1999:MCE` for MCMC sampler output.
	
	**Arguments**
	
		* ``c`` : sampler output on which to perform calculations.
		* ``alpha`` : the ``100 * (1 - alpha)``\% interval to compute.
		
	**Value**
	
		A ``ChainSummary`` type object with parameters contained in the rows of the ``value`` field, and lower and upper intervals in the first and second columns.  Results are for all chains combined.

.. function:: quantile(c::MCMCChains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])

	Compute posterior quantiles for MCMC sampler output.
	
	**Arguments**
	
		* ``c`` : sampler output on which to perform calculations.
		* ``q`` : probabilities at which to calculate quantiles.
		
	**Value**
	
		A ``ChainSummary`` type object with parameters contained in the rows of the ``value`` field, and quantiles in the columns.  Results are for all chains combined.

.. function:: summarystats(c::MCMCChains; batchsize::Integer=100)

	Compute posterior summary statistics for MCMC sampler output.
	
	**Arguments**
	
		* ``c`` : sampler output on which to perform calculations.
		* ``batchsize`` : number of iterations to include in a partitioning of the output for calculation of batch standard errors.
		
	**Value**
	
		A ``ChainSummary`` type object with parameters in the rows of the ``value`` field; and the sample mean, standard deviation, standard error, batch standard error (estimate of Monte Carlo error), and effective sample size in the columns.  Results are for all chains combined.
