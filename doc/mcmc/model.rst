.. index:: Model

.. _section-Model:

Model
-----

The ``Model`` type is designed to store the set of all model nodes, including parameter set :math:`\Theta` as denoted in  the :ref:`figure-Gibbs`.  In particular, it stores ``Dependent`` type objects in its ``nodes`` dictionary field.  Valid models are ones whose nodes form directed acyclic graphs (DAGs).  Sampling functions :math:`\{f_j\}_{j=1}^{B}` are saved as ``Sampler`` objects in the vector of field ``samplers``.  Vector elements :math:`j=1,\ldots,B` correspond to sampling blocks :math:`\{\Theta_j\}_{j=1}^{B}.`

Declaration
^^^^^^^^^^^

``type Model``

Fields
^^^^^^

* ``nodes::Dict{Symbol,Any}`` : a dictionary containing all input, logical, and stochastic model nodes.
* ``samplers::Vector{Sampler}`` : sampling functions for updating blocks of stochastic nodes.
* ``states::Vector{Vector{Float64}}`` : states of chains at the end of a possible series of MCMC runs.
* ``iter::Int`` : current MCMC draw from the target distribution.
* ``burnin::Int`` : number of initial draws to discard as a burn-in sequence to allow for convergence.
* ``hasinputs::Bool`` : whether values have been assigned to the input nodes.
* ``hasinits::Bool`` : whether initial values have been assigned to stochastic nodes.

Constructor
^^^^^^^^^^^

.. function:: Model(; iter::Integer=0, burnin::Integer=0, \
                    samplers::Vector{Sampler}=Sampler[], nodes...)

    Construct a ``Model`` object that defines a model for MCMC simulation.

    **Arguments**

        * ``iter`` : current iteration of the MCMC simulation.
        * ``burnin`` : number of initial draws to be discarded as a burn-in sequence to allow for convergence.
        * ``samplers`` : a vector of block-specific sampling functions.
        * ``nodes...`` : an arbitrary number of user-specified arguments defining logical and stochastic nodes in the model.  Argument values must be ``Logical`` or ``Stochastic`` type objects.  Their names in the model will be taken from the argument names.

    **Value**

        Returns a ``Model`` type object.

    **Example**

        See the :ref:`section-Line-Specification` section of the tutorial.

MCMC Engine
^^^^^^^^^^^

.. function:: mcmc(m::Model, inputs::Dict{Symbol}, \
                   inits::Vector{Dict{Symbol,Any}}, iters::Integer; \
                   burnin::Integer=0, thin::Integer=1, chains::Integer=1, \
                   verbose::Bool=true)
              mcmc(mc::ModelChains, iters::Integer; verbose::Bool=true)

    Simulate MCMC draws for a specified model.

    **Arguments**

        * ``m`` : a specified model.
        * ``mc`` : chains from a previous call to ``mcmc`` for which to simulate additional draws.
        * ``inputs`` : a dictionary of values for input model nodes.  Dictionary keys and values should be given for each input node.
        * ``inits`` : a vector of dictionaries that contain initial values for stochastic model nodes.  Dictionary keys and values should be given for each stochastic node.  Consecutive runs of the simulator will iterate through the vector's dictionary elements.
        * ``iters`` : number of draws to generate for each simulation run.
        * ``burnin`` : numer of initial draws to discard as a burn-in sequence to allow for convergence.
        * ``thin`` : step-size between draws to output.
        * ``chains`` : number of simulation runs to perform.
        * ``verbose`` : whether to print sampler progress at the console.

    **Value**

        A ``ModelChains`` type object of simulated draws.

    **Example**

        See the :ref:`section-Line-Simulation` section of the tutorial.

Indexing
^^^^^^^^

.. function:: getindex(m::Model, key::Symbol)

    Returns a model node identified by its symbol.  The syntax ``m[key]`` is converted to ``getindex(m, key)``.

    **Arguments**

        * ``m`` : a model contining the node to get.
        * ``key`` : symbol of the node to get.

    **Value**

        The specified node.

.. function:: keys(m::Model)
              keys(m::Model, ntype::Symbol, at...)

    Extract the symbols (keys) for all existing nodes or for nodes of a specified type.

    **Arguments**

        * ``m`` : a model containing the nodes of interest.
        * ``ntype`` : the type of nodes to return.  Options are
            * ``:all`` : all input, logical, and stochastic model nodes.
            * ``:assigned`` : nodes that have been assigned values.
            * ``:block`` : stochastic nodes being updated by the sampling block(s) ``at::Integer=0`` (default: all blocks).
            * ``:dependent`` : logical and stochastic (dependent) nodes in topologically sorted order.
            * ``:independent`` or ``:input`` : input (independent) nodes.
            * ``:logical`` : logical nodes.
            * ``:monitor`` : stochastic nodes being monitored in MCMC sampler output.
            * ``:output`` : stochastic nodes upon which no other stochastic nodes depend.
            * ``:source`` : nodes upon which the node ``at::Symbol`` or vector of nodes ``at::Vector{Symbol}`` depends.
            * ``:stochastic`` : stochastic nodes.
            * ``:target`` : topollogically sorted nodes that depend on the node ``at::Symbol`` or vector of nodes ``at::Vector{Symbol}``.
        * ``at...`` : additional positional arguments to be passed to the ``ntype`` options, as described above.

    **Value**

        A vector of node symbols.

Display
^^^^^^^

.. function:: draw(m::Model; filename::AbstractString="")

    Draw a `GraphViz <http://www.graphviz.org/>`_ DOT-formatted graph representation of model nodes and their relationships.

    **Arguments**

        * ``m`` : a model for which to construct a graph.
        * ``filename`` : an external file to which to save the resulting graph, or an empty string to draw to standard output (default).  If a supplied external file name does not include a dot (``.``), the file extension ``.dot`` will be appended automatically.

    **Value**

        The model drawn to an external file or standard output.  Stochastic, logical, and input nodes will be represented by ellipses, diamonds, and rectangles, respectively.  Nodes that are unmonitored in MCMC simulations will be gray-colored.

    **Example**

        See the :ref:`section-Line-DAG` section of the tutorial.

.. function:: graph(m::Model)

    Construct a graph representation of model nodes and their relationships.

    **Arguments**

        * ``m`` : a model for which to construct a graph.

    **Value**

        Returns a ``GenericGraph`` type object as defined in the `Graphs <http://graphsjl-docs.readthedocs.org/en/latest/index.html>`_ package.

.. function:: graph2dot(m::Model)

    Draw a `GraphViz <http://www.graphviz.org/>`_ DOT-formatted graph representation of model nodes and their relationships.

    **Arguments**

        * ``m`` : a model for which to construct a graph.

    **Value**

        A character string representation of the graph suitable for in-line processing.  Stochastic, logical, and input nodes will be represented by ellipses, diamonds, and rectangles, respectively.  Nodes that are unmonitored in MCMC simulations will be gray-colored.

    **Example**

        See the :ref:`section-Line-DAG` section of the tutorial.

.. function:: show(m::Model)

    Write a text representation of the model, nodes, and attributes to the current output stream.

.. function:: showall(m::Model)

    Write a verbose text representation of the model, nodes, and attributes to the current output stream.

Initialization
^^^^^^^^^^^^^^

.. function:: setinits!(m::Model, inits::Dict{Symbol,Any})

    Set the initial values of stochastic model nodes.

    **Arguments**

        * ``m`` : a model with nodes to be initialized.
        * ``inits`` : a dictionary of initial values for stochastic model nodes.  Dictionary keys and values should be given for each stochastic node.

    **Value**

        Returns the model with stochastic nodes initialized and the ``iter`` field set equal to 0.

    **Example**

        See the :ref:`section-Line-Development` section of the tutorial.

.. function:: setinputs!(m::Model, inputs::Dict{Symbol,Any})

    Set the values of input model nodes.

    **Arguments**

        * ``m`` : a model with input nodes to be assigned.
        * ``inputs`` : a dictionary of values for input model nodes.  Dictionary keys and values should be given for each input node.

    **Value**

        Returns the model with values assigned to input nodes.

    **Example**

        See the :ref:`section-Line-Development` section of the tutorial.

.. function:: setsamplers!(m::Model, samplers::Vector{Sampler})

    Set the block-samplers for stochastic model nodes.

    **Arguments**

        * ``m`` : a model with stochastic nodes to be sampled.
        * ``samplers`` : block-specific samplers.

    **Values:**

        Returns the model updated with the block-samplers.

    **Example**

        See the :ref:`section-Line-Specification` and :ref:`section-Line-Simulation` sections of the tutorial.

Parameter Block Operations
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: gradlogpdf(m::Model, block::Integer=0, transform::Bool=false; \
                         dtype::Symbol=:forward)
              gradlogpdf(m::Model, x::AbstractArray{T<:Real}, block::Integer=0, \
                         transform::Bool=false; dtype::Symbol=:forward)
              gradlogpdf!(m::Model, x::AbstractArray{T<:Real}, block::Integer=0, \
                          transform::Bool=false; dtype::Symbol=:forward)

    Compute the gradient of log-densities for stochastic nodes.

    **Arguments**

        * ``m`` : a model containing the stochastic nodes for which to compute the gradient.
        * ``x`` : a value (possibly different than the current one) at which to compute the gradient.
        * ``block`` : the sampling block of stochastic nodes for which to compute the gradient (default: all stochastic nodes).
        * ``transform`` : whether to compute the gradient of block parameters on the link–transformed scale.
        * ``dtype`` : type of differentiation for gradient calculations.  Options are
            * ``:central`` : central differencing.
            * ``:forward`` : forward differencing.

    **Value**

        The resulting gradient vector.  Method ``gradlogpdf!()`` additionally updates model ``m`` with supplied values ``x``.

    **Note**

        Numerical approximation of derivatives by central and forward differencing is performed with the `Calculus` package :cite:`white:2014:CP`.

.. function:: logpdf(m::Model, block::Integer=0, transform::Bool=false)
              logpdf(m::Model, x::AbstractArray{T<:Real}, block::Integer=0, \
                     transform::Bool=false)
              logpdf!(m::Model, x::AbstractArray{T<:Real}, block::Integer=0, \
                      transform::Bool=false)

    Compute the sum of log-densities for stochastic nodes.

    **Arguments**

        * ``m`` : a model containing the stochastic nodes for which to evaluate log-densities.
        * ``x`` : a value (possibly different than the current one) at which to evaluate densities.
        * ``block`` : the sampling block of stochastic nodes over which to sum densities (default: all stochastic nodes).
        * ``transform`` : whether to evaluate evaluate log-densities of block parameters on the link–transformed scale.

    **Value**

        The resulting numeric value of summed log-densities.  Method ``logpdf!()`` additionally updates model ``m`` with supplied values ``x``.

.. function:: simulate!(m::Model, block::Integer=0)

    Simulate one MCMC draw from a specified model.

    **Argument:**

        * ``m`` : a model specification.
        * ``block`` : the block for which to simulate an MCMC draw (default: all blocks).

    **Value**

        Returns the model updated with the MCMC draw and, in the case of ``block=0``, the ``iter`` field incremented by 1.

    **Example**

        See the :ref:`section-Line-Development` section of the tutorial.

.. function:: tune(m::Model, block::Integer=0)

    Get block-sampler tuning parameters.

    **Arguments**

        * ``m`` : a model with block-samplers.
        * ``block`` : the block for which to return the tuning parameters (default: all blocks).

    **Value**

        If ``block = 0``, a vector of dictionaries containing block-specific tuning parameters; otherwise, one block-specific dictionary.

.. function:: unlist(m::Model, block::Integer=0, transform::Bool=false)
              unlist(m::Model, nodekeys::Vector{Symbol}, transform::Bool=false)
              relist(m::Model, values::AbstractArray{T<:Real}, \
                     block::Integer=0, transform::Bool=false)
              relist(m::Model, values::AbstractArray{T<:Real}, \
                     nodekeys::Vector{Symbol},transform::Bool=false)
              relist!(m::Model, values::AbstractArray{T<:Real}, \
                      block::Integer=0, transform::Bool=false)
              relist!(m::Model, values::AbstractArray{T<:Real}, \
                      nodekeys::Vector{Symbol}, transform::Bool=false)

    Convert (unlist) sets of logical and/or stochastic node values to vectors, or reverse (relist) the process.

    **Arguments**

        * ``m`` : a model containing nodes to be unlisted or relisted.
        * ``values`` : values to re-list.
        * ``block`` : the sampling block of nodes to be listed (default: all blocks).
        * ``nodekeys`` : a vector of symbols identifying the nodes to be listed.
        * ``transform`` : whether to apply a link transformation in the conversion.

    **Value**

        The ``unlist`` methods return vectors of concatenated node values, ``relist`` return dictionaries of symbol keys and values for the specified nodes, and ``relist!`` return their model argument with values copied to the nodes.

.. function:: update!(m::Model, block::Integer=0)

    Update values of logical and stochastic model node according to their relationship with others in a model.

    **Arguments**

        * ``m`` : a mode with nodes to be updated.
        * ``block`` : the sampling block of nodes to be updated (default: all blocks).

    **Value**

        Returns the model with updated nodes.
