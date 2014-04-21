.. index:: Sampling Functions; Adaptive Metropolis within Gibbs

Adaptive Metropolis within Gibbs (AMWG)
---------------------------------------

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: amwg(x::Vector{T<:Real}, sigma::Vector{Float64}, logf::Function; \
				adapt::Bool=false, batch::Integer=50, target::Real=0.44)

.. function:: amwg!(v::VariateAMWG, sigma::Vector{Float64}, logf::Function; \
				adapt::Bool=false, batch::Integer=50, target::Real=0.44)


.. index:: VariateAMWG

VariateAMWG Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``VariateAMWG <: VariateVector``

Fields
``````

* ``data::Vector{VariateType}`` : 
* ``tune::TuneAMWG`` : 

Constructors
````````````

.. function:: VariateAMWG(x::Vector{T<:Real}, tune::TuneAMWG)
              VariateAMWG(x::Vector{T<:Real}, tune=nothing)


.. index:: TuneAMWG

TuneAMWG Type
^^^^^^^^^^^^^

Declaration
```````````

``type TuneAMWG``

Fields
``````

* ``adapt::Bool`` : 
* ``accept::Vector{Integer}`` : 
* ``batch::Integer`` : 
* ``m::Integer`` : 
* ``sigma::Vector{Float64}`` : 
* ``target::Real`` : 

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SamplerAMWG(params::Vector{T<:String}, sigma::Vector{U<:Real}; \
				adapt::Symbol=:none, batch::Integer=50, target::Real=0.44)
