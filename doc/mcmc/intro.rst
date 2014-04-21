Type Relationships
------------------

A Unified Modeling Language (UML) relational diagram of *MCMC* types is shown in the figure below.  In the diagram, types are represented with boxes that display their respective names in the top-most panels, and fields in the second panels.  By convention, plus signs denote fields that are publicly accessible, which is always the case for these structures in **julia**.  Hollow triangle arrows point to types that the originator extends.  Solid diamond arrows indicate that a number of instances of the type being pointed to are contained in the originator.  The undirected line between ``MCMCSampler`` and ``MCMCStochastic`` represents a bi-directional association.  Numbers on the graph indicate that there is one (1) or one or more (1..*) instances of a type at the corresponding end of a relationship.

The types are related as follows.  Type ``MCMCModel`` contains a dictionary field of model nodes (``Dict{String,Any}``) and a field of one or more sampling functions (``Vector{MCMCSampler}``).  Nodes can be either stochastic (``MCMCStochastic``) or logical (``MCMCLogical``).  They are inherited from the base ``Variate`` type and can be used with operators and in functions defined for that type.  A stochastic node is a model parameter or observed data vector on which a prior or likelihood distribution is specified, respectively.  A logical node is a deterministic function of other nodes.  Finally, the sampling functions in ``MCMCModel`` each correspond to a block of one or more model parameters (stochastic nodes) to be sampled from their full conditionals during the simulation.  General information about each type is provided in subsequent sections.

.. figure:: ../images/mcmcUML.png
	:align: center

	UML relational diagram of *MCMC* types and their fields.
