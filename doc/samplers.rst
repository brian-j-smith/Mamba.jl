Sampling Functions
==================

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

.. function:: VariateAMM(x::Vector{T<:Real}, tune=nothing)


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

.. function:: VariateAMWG(x::Vector{T<:Real}, tune=nothing)


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

.. index:: Sampling Functions; No-U-Turn Sampler

No-U-Turn Sampler (NUTS)
------------------------

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: nutseps(x::Vector{T<:Real}, fx::Function)

.. function:: nuts(x::Vector{T<:Real}, eps::Real, fx::Function; adapt::Bool=false, \
				target::Real=0.6)

.. function:: nuts!(v::VariateNUTS, eps::Real, fx::Function; adapt::Bool=false, \
           target::Real=0.6)


.. index:: VariateNUTS

VariateNUTS Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``VariateNUTS <: VariateVector``

Fields
``````

* ``data::Vector{VariateType}`` : 
* ``tune::TuneNUTS`` : 

Constructors
````````````

.. function:: VariateNUTS(x::Vector{T<:Real}, tune::TuneNUTS)

.. function:: VariateNUTS(x::Vector{T<:Real}, tune=nothing)


.. index:: TuneNUTS

TuneNUTS Type
^^^^^^^^^^^^^

Declaration
```````````

``type TuneNUTS``

Fields
``````
* ``adapt::Bool`` : 
* ``alpha::Float64`` : 
* ``eps::Float64`` : 
* ``epsbar::Float64`` : 
* ``gamma::Float64`` : 
* ``Hbar::Float64`` : 
* ``kappa::Float64`` : 
* ``m::Integer`` : 
* ``mu::Float64`` : 
* ``nalpha::Integer`` : 
* ``t0::Float64`` : 
* ``target::Float64`` : 

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SamplerNUTS(params::Vector{T<:String}; dtype::Symbol=:forward, \
				target::Real=0.6)

				
.. index:: Sampling Functions; Multivariate Slice

Multivariate Slice (Slice)
--------------------------

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: slice(x::Vector{T<:Real}, width::Vector{Float64}, logf::Function)

.. function:: slice!(v::VariateSlice, width::Vector{Float64}, logf::Function)


.. index:: VariateSlice

VariateSlice Type
^^^^^^^^^^^^^^^^^

Declaration
```````````

``VariateSlice <: VariateVector``

Fields
``````

* ``data::Vector{VariateType}`` : 
* ``tune::TuneSlice`` : 

Constructors
````````````

.. function:: VariateSlice(x::Vector{T<:Real}, tune::TuneSlice)

.. function:: VariateSlice(x::Vector{T<:Real}, tune=nothing)


.. index:: TuneSlice

TuneSlice Type
^^^^^^^^^^^^^^

Declaration
```````````

``type TuneSlice``

Fields
``````
* ``width::Vector{Float64}`` : 

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SamplerSlice(params::Vector{T<:String}, width::Vector{Float64})

.. index:: Sampling Functions; Slice within Gibbs

Slice within Gibbs (SliceWG)
----------------------------

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: slicewg(x::Vector{T<:Real}, width::Vector{Float64}, logf::Function)

.. function:: slicewg!(v::VariateSlice, width::Vector{Float64}, logf::Function)

MCMCSampler Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SamplerSliceWG(params::Vector{T<:String}, width::Vector{Float64})