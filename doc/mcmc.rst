MCMC Types
==========

.. index:: MCMCDepNode

MCMCDepNode
-----------

Declaration
^^^^^^^^^^^

``abstract MCMCDepNode{T} <: Variate{T}``

Fields
^^^^^^

* ``data::T`` :
* ``names::Vector{String}`` : 
* ``monitor::Bool`` :
* ``eval::Function`` :
* ``deps::Vector{String}`` :

Methods
^^^^^^^

.. function:: invlink(n::MCMCDepNode, x)

.. function:: link(n::MCMCDepNode, x)

.. function:: logpdf(n::MCMCDepNode, transform::Bool=false)

.. function:: show(n::MCMCDepNode)

.. function:: showall(n::MCMCDepNode)


.. index:: MCMCLogical

MCMCLogical
-----------

Declaration
^^^^^^^^^^^

``type MCMCLogical{T} <: MCMCDepNode{T}``

Fields
^^^^^^

* ``data::T`` :
* ``names::Vector{String}`` : 
* ``monitor::Bool`` :
* ``eval::Function`` :
* ``deps::Vector{String}`` :

Constructors
^^^^^^^^^^^^

.. function:: MCMCLogical(expr::Expr, monitor::Bool=true)

.. function:: MCMCLogical(length::Integer, expr::Expr, monitor::Bool=true)

.. function:: MCMCLogical(m::Integer, n::Integer, expr::Expr, monitor::Bool=true)


Methods
^^^^^^^

.. function:: setinits!(l::MCMCLogical, m::MCMCModel, x=nothing)

.. function:: update!(l::MCMCLogical, m::MCMCModel)


.. index:: MCMCStochastic

MCMCStochastic
--------------

Declaration
^^^^^^^^^^^

``type MCMCStochastic{T} <: MCMCDepNode{T}``

Fields
^^^^^^

* ``data:T`` :
* ``names::Vector{String}`` : 
* ``monitor::Bool`` :
* ``eval::Function`` :
* ``deps::Vector{String}`` :
* ``distr::DistributionStruct`` :

Aliases
^^^^^^^

.. code-block:: julia

	typealias DistributionStruct Union(Distribution, VecOrMat{Distribution})

Constructors
^^^^^^^^^^^^

.. function:: MCMCStochastic(expr::Expr, monitor::Bool=true)

.. function:: MCMCStochastic(length::Integer, expr::Expr, monitor::Bool=true)

.. function:: MCMCStochastic(m::Integer, n::Integer, expr::Expr, monitor::Bool=true)

Methods
^^^^^^^

.. function:: insupport(s::MCMCStochastic)

.. function:: invlink(s::MCMCStochastic, x)

.. function:: link(s::MCMCStochastic, x)

.. function:: logpdf(s::MCMCStochastic, transform::Bool=false)

.. function:: setinits!(l::MCMCStochastic, m::MCMCModel, x)

.. function:: update!(l::MCMCStochastic, m::MCMCModel)


.. index:: MCMCSampler

MCMCSampler
-----------

Declaration
^^^^^^^^^^^

``type MCMCSampler``

Fields
^^^^^^

* ``params::Vector{String}`` :
* ``links::Vector{String}`` :
* ``eval::Function`` :
* ``tune::Dict`` :

Constructor
^^^^^^^^^^^

.. function:: MCMCSampler(params::Vector{T<:String}, expr::Expr, tune::Dict=Dict())

Methods
^^^^^^^

.. function:: show(s::MCMCSampler)

.. function:: showall(s::MCMCSampler)


.. index:: MCMCModel

MCMCModel
---------

Declaration
^^^^^^^^^^^

``type MCMCModel``

Fields
^^^^^^

* ``nodes::Dict{String,Any}`` :
* ``links::Vector{String}`` :
* ``samplers::Vector{MCMCSampler}`` :
* ``iter::Integer`` :
* ``burnin::Integer`` :
* ``chain::Integer`` :
* ``hasinputs::Bool`` :
* ``hasinits::Bool`` :

Constructor
^^^^^^^^^^^

.. function:: MCMCModel(; iter::Integer=0, burnin::Integer=0, chain::Integer=1, \
				samplers::Vector{MCMCSampler}=Array(MCMCSampler, 0), nodes...)

Methods
^^^^^^^

.. function:: gradient(m::MCMCModel, block::Integer=0, transform::Bool=false, \
				dtype::Symbol=:central)

.. function:: gradient(m::MCMCModel, x::Vector{T<:Real}, block::Integer=0, \
				transform::Bool=false, dtype::Symbol=:central)

.. function:: graph(m::MCMCModel)

.. function:: graph2dot(m::MCMCModel)

.. function:: graph2dot(m::MCMCModel, filename::String)
				
.. function:: keys(m::MCMCModel, ntype::Symbol=:assigned, block::Integer=0)

.. function:: logpdf(m::MCMCModel, block::Integer=0, transform::Bool=false)

.. function:: logpdf(m::MCMCModel, x::Vector{T<:Real}, block::Integer=0, \
				transform::Bool=false)
				
.. function:: mcmc(model::MCMCModel, inputs::Dict{T<:String}, \
				inits::Vector{Dict{U<:String,Any}}, iter::Integer; \
				burnin::Integer=0, thin::Integer=1, chains::Integer=1)

.. function:: relist(m::MCMCModel, values::Vector{T<:Real}, block::Integer=0, \
				transform::Bool=false)

.. function:: relist(m::MCMCModel, values::Vector{T<:Real}, nkeys::Vector{U<:String}, \
				transform::Bool=false)

.. function:: relist!(m::MCMCModel, values::Vector{T<:Real}, block::Integer=0, \
				transform::Bool=false)

.. function:: relist!(m::MCMCModel, values::Vector{T<:Real}, nkeys::Vector{U<:String}, \
				transform::Bool=false)
							
.. function:: setinits!(m::MCMCModel, inits::Dict{T<:String,Any})

.. function:: setinputs!(m::MCMCModel, inputs::Dict{T<:String,Any})

.. function:: setsamplers!(m::MCMCModel, samplers::Vector{MCMCSampler})

.. function:: show(m::MCMCModel)

.. function:: showall(m::MCMCModel)

.. function:: simulate!(m::MCMCModel, block::Integer=0)

.. function:: tune(m::MCMCModel, block::Integer=0)

.. function:: unlist(m::MCMCModel, block::Integer=0, transform::Bool=false)

.. function:: unlist(m::MCMCModel, nkeys::Vector{T<:String}, transform::Bool=false)

.. function:: update!(m::MCMCModel, block::Integer=0)


.. index:: MCMCChains

MCMCChains
----------

Declaration
^^^^^^^^^^^

``type MCMCChains``

Fields
^^^^^^

* ``data::Array{VariateType,3}`` :
* ``names::Vector{String}`` :
* ``start::Integer`` :
* ``thin::Integer`` :
* ``model::MCMCModel`` :

Constructor
^^^^^^^^^^^

.. function:: MCMCChains(names::Vector{T<:String}, iter::Integer; start::Integer=1, \
                   thin::Integer=1, chains::Integer=1, model::MCMCModel=MCMCModel())

Methods
^^^^^^^

.. function:: autocor(c::MCMCChains; lags::Vector=[1,5,10,50], relative::Bool=true)

.. function:: cor(c::MCMCChains)

.. function:: describe(c::MCMCChains; batchsize::Integer=100, \
				q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])

.. function:: dic(c::MCMCChains)

.. function:: gelmandiag(c::MCMCChains; alpha::Real=0.05, mpsrf::Bool=false, \
				transform::Bool=false)
				
.. function:: hpd(c::MCMCChains; alpha::Real=0.05)

.. function:: quantile(c::MCMCChains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])

.. function:: summarystats(c::MCMCChains; batchsize::Integer=100)
