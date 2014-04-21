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
              VariateSlice(x::Vector{T<:Real}, tune=nothing)


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
