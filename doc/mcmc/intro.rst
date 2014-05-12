Type Relationships
------------------

A Unified Modeling Language (UML) relational diagram of *MCMC* types is shown in the figure below.  In the diagram, types are represented with boxes that display their respective names in the top-most panels, and fields in the second panels.  By convention, plus signs denote fields that are publicly accessible, which is always the case for these structures in **julia**.  Hollow triangle arrows point to types that the originator extends.  Solid diamond arrows indicate that a number of instances of the type being pointed to are contained in the originator.  The undirected line between ``MCMCSampler`` and ``MCMCStochastic`` represents a bi-directional association.  Numbers on the graph indicate that there is one (1), zero or more (0..*), or one or more (1..*) instances of a type at the corresponding end of a relationship.

.. figure:: ../images/mcmcUML.png
	:align: center

	UML relational diagram of *MCMC* types and their fields.

The *MCMC* types are related as follows.  Type ``MCMCModel`` contains a dictionary field of model nodes (``Dict{String,Any}``) and a field of one or more sampling functions (``Vector{MCMCSampler}``). Nodes can be one of the three types:

	* **Stochastic nodes** (``MCMCStochastic``) are any model terms that have likelihood or prior distributional specifications.
	* **Logical nodes** (``MCMCLogical``) are terms that are deterministic functions of other nodes.
	* **Input nodes** (not shown) are any other model terms and data types that are considered to be fixed quantities in the analysis.

``MCMCStochastic`` and ``MCMCLogical`` are inherited from the base :ref:`section-Variate` type and can be used with operators and in functions defined for that type.  The sampling functions in ``MCMCModel`` each correspond to a block of one or more model parameters (stochastic nodes) to be sampled from a target distribution (e.g. full conditional) during the simulation.  Finally, ``MCMCChain`` stores simulation output for a given model.  General information about each type is provided in the subsequent sections.

