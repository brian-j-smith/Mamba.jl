.. index:: Sampling Functions; Binary Hamiltonian Monte Carlo

Binary Hamiltonian Monte Carlo (BHMC)
-----------------------------------------

Implementation of the binary-state Hamiltonian Monte Carlo sampler of Pakman :cite:`pakman:2013:hmcb`.  The sampler simulates autocorrelated draws from a distribution that can be specified up to a constant of proportionality.


Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: bhmc!(v::BHMCVariate, T::Real, logf::Function)

    Simulate one draw from a target distribution using the BHMC sampler.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``T`` : Length of time over which particle paths are simulated.
        * ``logf`` : function to compute the log-transformed density (up to a normalizing constant) at ``v.value``.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-bhmc:

    **Example**

        .. literalinclude:: bhmc.jl
            :language: julia


.. index:: Sampler Types; BHMCVariate

BHMCVariate Type
^^^^^^^^^^^^^^^^

Declaration
```````````

``BHMCVariate <: VectorVariate``

Fields
``````

* ``value::Vector{Float64}`` : vector of sampled values.
* ``tune::BHMCTune`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: BHMCVariate(x::Vector{Float64}, tune::BHMCTune)
              BHMCVariate(x::Vector{Float64}, tune=nothing)

    Construct a ``BHMCVariate`` object that stores sampled values and tuning parameters for BHMC sampling.

    **Arguments**

        * ``x`` : vector of sampled values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If ``nothing`` is supplied, parameters are set to their defaults.

    **Value**

        Returns a ``BHMCVariate`` type object with fields pointing to the values supplied to arguments ``x`` and ``tune``.

.. index:: Sampler Types; BHMCTune

BHMCTune Type
^^^^^^^^^^^^^

Declaration
```````````

``type BHMCTune``

Fields
``````
* ``T::Float64`` : Length of time over which particle paths are simulated.
* ``X::Vector{Float64}`` : Initial particle positions.
* ``V::Vector{Float64}`` : Initial particle velocites. 
* ``wallhits::Int`` : Number of times particles are reflected off the 0 threshold.
* ``wallcrosses::Int`` :Number of times particles travel through 0 threshold. 

Sampler Constructor
^^^^^^^^^^^^^^^^^^^

.. function:: BHMC(params::Vector{Symbol}, T::Real)

    Construct a ``Sampler`` object for BHMC sampling.  Parameters are assumed to have binary numerical values (0 or 1).

    **Arguments**

        * ``params`` : stochastic nodes containing the parameters to be updated with the sampler.
        * ``T::Real`` : Length of time over which particle paths are simulated.
        
    **Value**

        Returns a ``Sampler`` type object.
