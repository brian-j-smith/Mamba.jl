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
              VariateNUTS(x::Vector{T<:Real}, tune=nothing)


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
