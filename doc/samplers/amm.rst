.. index:: Sampling Functions; Adaptive Multivariate Metropolis

Adaptive Multivariate Metropolis (AMM)
--------------------------------------

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: amm(x::Vector{T<:Real}, Sigma::Cholesky{Float64}, logf::Function; \
				adapt::Bool=false)

.. function:: amm!(v::VariateAMM, Sigma::Cholesky{Float64}, logf::Function; \
				adapt::Bool=false)


.. index:: VariateAMM

VariateAMM Type
^^^^^^^^^^^^^^^

Declaration
```````````

``VariateAMM <: VariateVector``

Fields
``````

* ``data::Vector{VariateType}`` : 
* ``tune::TuneAMM`` : 

Constructors
````````````

.. function:: VariateAMM(x::Vector{T<:Real}, tune::TuneAMM)
              VariateAMM(x::Vector{T<:Real}, tune=nothing)


.. index:: TuneAMM

TuneAMM Type
^^^^^^^^^^^^

Declaration
```````````

``type TuneAMM``

Fields
``````
* ``adapt::Bool`` : 
* ``beta::Real`` : 
* ``m::Integer`` : 
* ``mu::Vector{Float64}`` : 
* ``scale::Real`` : 
* ``Sigma::Cholesky{Float64}`` : 
* ``SigmaLm::Matrix{Float64}`` : 

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SamplerAMM(params::Vector{T<:String}, Sigma::Matrix{U:<Real}; \
				adapt::Symbol=:none)

