.. index:: Sampling Functions; Shrinkage Slice

.. _section-Slice:

Shrinkage Slice (Slice)
-----------------------

Implementation of the shrinkage slice sampler of Neal :cite:`neal:2003:SS` for simulating autocorrelated draws from a distribution that can be specified up to a constant of proportionality.

Model-Based Constructor
^^^^^^^^^^^^^^^^^^^^^^^

.. function:: Slice(params::ElementOrVector{Symbol}, \
                    width::ElementOrVector{T<:Real}, stype::Symbol=:multivar; \
                    transform::Bool=false)

    Construct a ``Sampler`` object for shrinkage slice sampling.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``params`` : stochastic node(s) to be updated with the sampler.
        * ``width`` : scaling value or vector of the same length as the combined elements of nodes ``params``, defining initial widths of a hyperrectangle from which to simulate values.
        * ``stype`` : sampler type. Options are
            * ``:multivar`` : Joint multivariate sampling of parameters.
            * ``:univar`` : Sequential univariate sampling.
        * ``transform`` : whether to sample parameters on the link-transformed scale (unconstrained parameter space).  If ``true``, then constrained parameters are mapped to unconstrained space according to transformations defined by the :ref:`section-Stochastic` ``unlist()`` function, and ``width`` is interpreted as being relative to the unconstrained parameter space.  Otherwise, sampling is relative to the untransformed space.

    **Value**

        Returns a ``Sampler{SliceTune{Univariate}}`` or ``Sampler{SliceTune{Multivariate}}`` type object if sampling univariately or multivariately, respectively.

    **Example**

        See the :ref:`Birats <example-Birats>`, :ref:`Rats <example-Rats>`, and other :ref:`section-Examples`.

Stand-Alone Function
^^^^^^^^^^^^^^^^^^^^

.. function:: slice!(v::SamplerVariate{SliceTune{F<:SliceForm}}, \
                     width::ElementOrVector{T<:Real}, \
                     logf::Function, stype::Symbol=:multivar)

    Simulate one draw from a target distribution using a shrinkage slice sampler.  Parameters are assumed to be continuous, but may be constrained or unconstrained.

    **Arguments**

        * ``v`` : current state of parameters to be simulated.
        * ``width`` : scalar or vector of the same length as ``v``, defining initial widths of a hyperrectangle from which to simulate values.
        * ``logf`` : function that takes a single ``DenseVector`` argument of parameter values at which to compute the log-transformed density (up to a normalizing constant).
        * ``stype`` : sampler type. Options are
            * ``:multivar`` : Joint multivariate sampling of parameters.
            * ``:univar`` : Sequential univariate sampling.

    **Value**

        Returns ``v`` updated with simulated values and associated tuning parameters.

    .. _example-slice:

    **Example**

        The following example samples parameters in a simple linear regression model.  Details of the model specification and posterior distribution can be found in the :ref:`section-Supplement`.  Also, see the :ref:`example-Line_AMWG_Slice` example.

        .. literalinclude:: slice.jl
            :language: julia


.. index:: Sampler Types; SliceVariate

SliceVariate Type
^^^^^^^^^^^^^^^^^

Declaration
```````````

``type SamplerVariate{SliceTune{F<:SliceForm}}``

Fields
``````

* ``value::Vector{Float64}`` : simulated values.
* ``tune::SliceTune{F}`` : tuning parameters for the sampling algorithm.

Constructors
````````````

.. function:: SliceVariate(x::AbstractVector{T<:Real})

    Construct a ``SliceVariate`` object that stores simulated values and tuning parameters for slice sampling.

    **Arguments**

        * ``x`` : simulated values.
        * ``tune`` : tuning parameters for the sampling algorithm.  If not supplied, parameters are set to their defaults.

    **Value**

        Returns a ``SamplerVariate{SliceTune{SliceForm}}`` type object with fields set to the values supplied to argument ``x``.

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

* ``width::Union{Float64, Vector{Float64}}`` : initial widths defining hyperrectangles from which to simulate values.
