.. index:: Sampling Functions; Shrinkage Slice

.. _section-Slice:

Shrinkage Slice (Slice)
-----------------------

Implementation of the shrinkage slice sampler of Neal :cite:`neal:2003:SS` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: Slice(params::ElementOrVector{Symbol}, \
                    width::ElementOrVector{T<:Real}, \
                    ::Type{F<:SliceForm}=Multivariate; \
                    transform::Bool=false)

    Construct a ``Sampler`` object for Slice sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.
        * ``width`` : scaling value or vector of the same length as the combined elements of nodes ``params``, defining initial widths of a hyperrectangle from which to simulate values.
        * ``F`` : sampler type. Options are
            * ``Univariate`` : sequential univariate sampling of parameters .
            * ``Multivariate`` : joint multivariate sampling.
        * ``transform`` : whether to sample parameters on the link-transformed scale (unconstrained parameter space).  If ``true``, then constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function, and ``width`` is interpreted as being relative to the unconstrained parameter space.  Otherwise, sampling is relative to the untransformed space.

    **Value**

        Returns a ``Sampler{SliceTune{Univariate}}`` or ``Sampler{SliceTune{Multivariate}}`` type object if sampling univariately or multivariately, respectively.

    **Example**

        See the :ref:`Birats <example-Birats>`, :ref:`Rats <example-Rats>`, and other :ref:`section-Examples`.

Stand-Alone Functions
^^^^^^^^^^^^^^^^^^^^^

.. function:: sample!(v::SliceUnivariate)
              sample!(v::SliceMultivariate)

    Draw one sample from a target distribution using the Slice univariate or multivariate sampler.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-slice:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.  Also, see the :ref:`example-Line_AMWG_Slice` example.

        .. literalinclude:: slice.jl
            :language: julia


.. index:: Sampler Types; SliceUnivariate
.. index:: Sampler Types; SliceMultivariate

Slice Variate Types
^^^^^^^^^^^^^^^^^^^

Declaration
```````````

.. code-block:: julia

    typealias SliceUnivariate SamplerVariate{SliceTune{Univariate}}
    typealias SliceMultivariate SamplerVariate{SliceTune{Multivariate}}

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::SliceTune{F<:SliceForm}`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: SliceUnivariate(x::AbstractVector{T<:Real}, \
                              width::ElementOrVector{U<:Real}, logf::Function)
              SliceMultivariate(x::AbstractVector{T<:Real}, \
                                width::ElementOrVector{U<:Real}, logf::Function)

    Construct an object that stores simulated values and tuning parameters for Slice sampling.

    **Arguments**

        * ``x`` : initial values.
        * ``width`` : scaling value or vector of the same length as the combined elements of nodes ``params``, defining initial widths of a hyperrectangle from which to simulate values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).

    **Value**

        Returns an object of the same type as the constructor name for univariate or multivariate sampling, respectively, with fields set to the supplied ``x`` and tuning parameter values.

.. index:: Sampler Types; SliceForm
.. index:: Sampler Types; SliceTune

SliceTune Type
^^^^^^^^^^^^^^

Declaration
```````````

.. code-block:: julia

    typealias SliceForm Union{Univariate, Multivariate}
    type SliceTune{F<:SliceForm} <: SamplerTune

Fields
``````

* ``logf::Nullable{Function}`` : function supplied to the constructor to compute the log-transformed density, or null if not supplied.
* ``width::Union{Float64, Vector{Float64}}`` : initial widths of hyperrectangles from which to simulate values.
