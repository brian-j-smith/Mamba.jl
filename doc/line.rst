.. index::
	single: Examples; Linear Regression

.. _section-Line:

Linear Regression Example
=========================

Bayesian Model Specification
----------------------------

In the sections that follow, the Bayesian simple linear regression example from the `BUGS 0.5` manual :cite:`spiegelhalter:1996:BIU` is used to illustrate features of the package.  The example describes a regression relationship between observations :math:`\bm{x} = (1, 2, 3, 4, 5)^\top` and :math:`\bm{y} = (1, 3, 3, 3, 5)^\top` that can be expressed as

.. math::

	\bm{y} &\sim N\left(\bm{\mu}, \sigma^2 \bm{I}\right) \\
	\bm{\mu} &= \bm{X} \bm{\beta}

with prior distribution specifications

.. math::

    \bm{\beta} &\sim N\left(
      \bm{\mu}_\pi =
      \begin{bmatrix}
        0 \\
        0 \\
      \end{bmatrix},
      \bm{\Sigma}_\pi =
      \begin{bmatrix}
        1000 & 0 \\
        0 & 1000 \\
      \end{bmatrix}
    \right) \\
    \sigma^2 &\sim InvGamma(\alpha_\pi = 0.001, \beta_\pi = 0.001),

where :math:`\bm{\beta} = (\beta_0, \beta_1)^\top`, and :math:`\bm{X}` is a design matrix with an intercept vector of ones in the first column and :math:`\bm{x}` in the second.  Primary interest lies in making inference about the :math:`\beta_0`, :math:`\beta_1`, and :math:`\sigma^2` parameters, based on their posterior distribution.  A computational consideration in this example is that inference cannot be obtain from the joint posterior directly because of its nonstandard form, derived below up to a constant of proportionality.

.. math::

  p(\bm{\beta}, \sigma^2 | \bm{y}) &\propto p(\bm{y} | \bm{\beta}, \sigma^2) p(\bm{\beta}) p(\sigma^2) \\
    &\propto \left(\sigma^2\right)^{-n/2} \exp\left\{-\frac{1}{2 \sigma^2} (\bm{y} - \bm{X} \bm{\beta})^\top (\bm{y} - \bm{X} \bm{\beta}) \right\} \\
    &\quad \times \exp\left\{-\frac{1}{2} (\bm{\beta} - \bm{\mu}_\pi)^\top \bm{\Sigma}_\pi^{-1} (\bm{\beta} - \bm{\mu}_\pi) \right\}
    \left(\sigma^2\right)^{-\alpha_\pi - 1} \exp\left\{-\beta_\pi / \sigma^2\right\}

A common alternative is to make approximate inference based on parameter values simulated from the posterior with MCMC methods.


Model Implementation
--------------------

Nodes
^^^^^

In the `MCMCsim` package, terms that appear in the Bayesian model specification are referred to as *nodes*.  Nodes are classified as one of three types: stochastic, logical, or input.

	* **Stochastic nodes** are any model terms that have likelihood or prior distributional specifications.  In the regression example, :math:`\bm{y}`, :math:`\bm{\beta}`, and :math:`\sigma^2` are stochastic nodes.
	* **Logical nodes** are terms, like :math:`\bm{\mu}`, that are deterministic functions of other nodes.
	* **Input nodes** are any remaining model terms (:math:`\bm{X}`) and are considered to be fixed quantities in the analysis.

Note that the :math:`\bm{y}` node has both a distributional specification and is a fixed quantity.  It is designated a stochastic node in `MCMCsim` because of the distributional specification.  This allows for the possibility of model terms with distributional specifications that are a mix of observed and unobserved elements, as in the case of missing values in response vectors.

For model implementation, all nodes are stored in and accessible from a **julia** dictionary structure called ``model`` with the names (keys) of nodes being character strings.  The regression nodes will be named ``"y"``, ``"beta"``, ``"s2"``, ``"mu"``, and ``"xmat"`` to correspond to the stochastic, logical, and input nodes identified in the previous paragraph.  Implementation begins by instantiating the stochastic and logical nodes using the `MCMCsim`--supplied constructors ``MCMCStochastic`` and ``MCMCLogical``.  Logical and stochastic nodes for a model are defined with a call to the ``MCMCModel`` constructor.  The model constructor formally defines and assigns names to the nodes.  User-specified names are given on the left-hand sides of the arguments to which ``MCMCLogical`` and ``MCMCStochastic`` nodes are passed.

.. code-block:: julia

	using MCMCsim
	using Distributions

	## Model Implementation

	line = MCMCModel(

	  y = MCMCStochastic(5,
	    quote
	      mu = model["mu"]
	      s2 = model["s2"]
	      IsoNormal(mu, s2)
	    end,
	    false
	  ),

	  mu = MCMCLogical(5,
	    :(model["xmat"] * model["beta"]),
	    false
	  ),

	  beta = MCMCStochastic(2,
	    :(IsoNormal(2, 1000))
	  ),

	  s2 = MCMCStochastic(
	    :(InverseGamma(0.001, 0.001))
	  )

	)
	
A single integer value for the first ``MCMCStochastic`` constructor argument indicates that the node is a vector of the specified length.  Absence of an integer value implies a scalar node.  The next argument is a quoted expression that can contain any valid **julia** code.  Expressions for stochastic nodes must return a distribution object from or compatible with the `Distributions <http://distributionsjl.readthedocs.org/en/latest/>`_ package.  Such objects represent the nodes' distributional specifications.  The dimensions of a stochastic node and its distribution object must match.  An optional boolean argument after the expression can be specified to indicate whether values of the node should be monitored (saved) during MCMC simulations (default: ``true``).

In the example, nodes ``y``, ``mu``, and ``beta`` are vectors, ``s2`` is a scalar, and the first two are not being monitored.  Further, note that the model could be implemented without the logical node ``mu``.  It is created here primarily for illustrative purposes.


Sampling Schemes
^^^^^^^^^^^^^^^^

The package provides a flexible system for the specification of schemes to sample stochastic nodes.  Arbitrary blocking of nodes and designation of block-specific samplers is supported.  Furthermore, block-updating of nodes can be performed with samplers provided, defined by the user, or available from other packages.  Schemes are specified as vectors of ``MCMCSampler`` objects.  Constructors are provided for several popular sampling algorithms, including adaptive Metropolis, No-U-Turn (NUTS), and slice sampling.  Example schemes are shown below.  In the first one, NUTS is used for the sampling of ``beta`` and slice for ``s2``.  The two nodes are block together in the second scheme and sampled jointly with NUTS.

.. code-block:: julia

	## Hybrid No-U-Turn and Slice Sampling Scheme
	scheme1 = [SamplerNUTS(["beta"]),
	           SamplerSlice(["s2"], [1.0])]

	## No-U-Turn Sampling Scheme
	scheme2 = [SamplerNUTS(["beta", "s2"])]

Additionally, users are free to create their own samplers with the generic ``MCMCSampler`` constructor.  This is particularly useful in settings were full conditional distributions are of standard forms for some nodes and can be sampled from directly.  Such is the case for the full conditional of :math:`\bm{\beta}` which can be written as

.. math::
  p(\bm{\beta} | \sigma^2, \bm{y}) &\propto p(\bm{y} | \bm{\beta}, \sigma^2) p(\bm{\beta}) \\
  &\propto \exp\left\{-\frac{1}{2} (\bm{\beta} - \bm{\mu})^\top \bm{\Sigma}^{-1} (\bm{\beta} - \bm{\mu})\right\},

where :math:`\bm{\Sigma} = \left(\frac{1}{\sigma^2} \bm{X}^\top \bm{X} + \bm{\Sigma}_\pi^{-1}\right)^{-1}` and :math:`\bm{\mu} = \bm{\Sigma} \left(\frac{1}{\sigma^2} \bm{X}^\top \bm{y} + \bm{\Sigma}_\pi^{-1} \bm{\mu}_\pi\right)`, which is recognizable as multivariate normal.  Likewise, 

.. math::

	p(\sigma^2 | \bm{\beta}, \mathbf{y}) &\propto p(\bm{y} | \bm{\beta}, \sigma^2) p(\sigma^2) \\
    &\propto \left(\sigma^2\right)^{-(n/2 + \alpha_\pi) - 1} \exp\left\{-\frac{1}{\sigma^2} \left(\frac{1}{2} (\bm{y} - \bm{X} \bm{\beta})^\top (\bm{y} - \bm{X} \bm{\beta}) + \beta_\pi \right) \right\},

whose form is inverse gamma with :math:`n / 2 + \alpha_\pi` shape and :math:`(\bm{y} - \bm{X} \bm{\beta})^\top (\bm{y} - \bm{X} \bm{\beta}) / 2 + \beta_\pi` scale parameters.  A user-defined sampling scheme to generate draws from these full conditions is constructed below.

.. code-block:: julia

	## User-Defined Samplers

	Gibbs_beta = MCMCSampler(["beta"],
	  quote
	    beta = model["beta"]
	    s2 = model["s2"]
	    xmat = model["xmat"]
	    y = model["y"]
	    beta_mean = mean(beta.distr)
	    beta_invcov = invcov(beta.distr)
	    Sigma = inv(xmat' * xmat / s2 + beta_invcov)
	    mu = Sigma * (xmat' * y / s2 + beta_invcov * beta_mean)
	    rand(MvNormal(mu, Sigma))
	  end
	)

	Gibbs_s2 = MCMCSampler(["s2"],
	  quote
	    beta = model["beta"]
	    s2 = model["s2"]
	    xmat = model["xmat"]
	    y = model["y"]
	    a = length(y) / 2.0 + s2.distr.shape
	    b = sum((y - xmat * beta).^2) / 2.0 + s2.distr.scale
	    rand(InverseGamma(a, b))
	  end
	)
	
	## User-Defined Sampling Scheme
	scheme3 = [Gibbs_beta, Gibbs_s2]

When it is possible to do so, direct sampling from full conditions is often preferred in practice because it tends to be more efficient than general-purpose algorithms.  Schemes that mix the two approaches can be used if full conditionals are available for some of the parameters but not for others.

A sampling scheme can be assigned to an existing model with a call to the ``setsamplers!`` function.

.. code-block:: julia

	## Sampling Scheme Assignment
	setsamplers!(line, scheme1)

Alternative, a predefined scheme can be passed in to the ``MCMCModel`` constructor at the time of model implementation as the value to its ``samplers`` argument.


Directed Acyclic Graphs
-----------------------

One of the internal structures created by ``MCMCModel`` is a graph representation of the model nodes and their associations.  Like `OpenBUGS`, `JAGS`, and other `BUGS` clones, `MCMCsim` fits models whose nodes form a directed acyclic graph (DAG).  A *DAG* is a graph in which nodes are connected by directed edges and no node has a path that loops back to itself.  With respect to statistical models, directed edges point from parent nodes to the child nodes that depend on them.  Thus, a child node is independent of all others, given its parents.

The DAG representation of an ``MCMCModel`` can be printed out at the command line or saved to an external file in a format that can be displayed in `Graphviz <http://www.graphviz.org/>`_.

.. code-block:: julia

	## Graph Representation of Nodes

	>>> print(graph2dot(line))
	
	digraph MCMCModel {
	  "beta" [shape="ellipse"];
	  	"beta" -> "mu";
	  "mu" [fillcolor="gray85", shape="diamond", style="filled"];
	  	"mu" -> "y";
	  "xmat" [fillcolor="gray85", shape="box", style="filled"];
	  	"xmat" -> "mu";
	  "s2" [shape="ellipse"];
	  	"s2" -> "y";
	  "y" [fillcolor="gray85", shape="ellipse", style="filled"];
	}
	
	>>> graph2dot(line, "lineDAG.dot")

Either the printed or saved output can be passed to the Graphviz software to plot a visual representation of the model.  A generated plot of the regression model graph is show in the figure below.

.. figure:: images/LineDAG.png
	:align: center
	
	Directed acyclic graph representation of the example regression model nodes.

Stochastic, logical, and input nodes are represented by ellipses, diamonds, and rectangles, respectively.  Gray-colored nodes are ones designated as unmonitored in MCMC simulations.  The DAG not only allows the user to visually check that the model specification is the intended one, but is also used internally to check that nodal relationships are acyclic.


MCMC Simulation
---------------

Data
^^^^

For the example, observations :math:`(\bm{x}, \bm{y})` are stored in a **julia** dictionary defined in the code block below.  Included are predictor and response vectors ``"x"`` and ``"y"`` as well as a design matrix ``"xmat"`` corresponding to the model matrix :math:`\bm{X}`.

.. code-block:: julia

	## Data
	data = (String => Any)[
	  "x" => [1, 2, 3, 4, 5],
	  "y" => [1, 3, 3, 3, 5]
	]
	data["xmat"] = [ones(5) data["x"]]

Initial Values
^^^^^^^^^^^^^^

A **julia** vector of dictionaries containing initial values for all stochastic nodes must be created.  The dictionary keys should match the node names, and their values should be vectors whose elements are the same type of structures as the nodes.  Vector elements are cycled through to initialize nodes when multiple runs of the MCMC simulator are performed.  Initial values for the regression example are as given below.

.. code-block:: julia

	## Initial Values
	inits = [["y" => data["y"],
	          "beta" => rand(Normal(0, 1), 2),
	          "s2" => rand(Gamma(1, 1))]
	         for i in 1:3]

Initial values for ``y`` are those in the observed response vector.  Likewise, the node is not updated in the sampling schemes defined earlier and thus retains its initial values throughout MCMC simulations.  Three different sets of initial values are generated for ``beta`` from a normal distribution and for ``s2`` from a gamma distribution.


MCMC Engine
^^^^^^^^^^^

MCMC simulation of draws from the posterior distribution of a declared set of model nodes and sampling scheme is performed with the ``mcmc`` function.  As shown below, the first three arguments are an ``MCMCModel`` object, a dictionary of values for input nodes, and a dictionary of initial values.  The number of draws to generate in each simulation run is given as the fourth argument.  The remaining arguments are named such that ``burnin`` is the number of initial values to discard to allow for convergence; ``thin`` defines the interval between draws to be retained in the output; and ``chains`` specifies the number of times to run the simulator.

.. code-block:: julia

	## MCMC Simulation Engine
	sim1 = mcmc(line, data, inits, 10000, burnin=250, thin=2, chains=3)

Results are retuned as an ``MCMCChains`` object on which methods for posterior inference are defined.


Posterior Inference
-------------------

Convergence Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^

Checks of MCMC output should be performed to assess convergence of simulated draws to the posterior distribution.  One popular check is the diagnostic of Brooks, Gelman, and Rubin :cite:`brooks:1998:GMM`, :cite:`gelman:1992:IIS`.  It is available through the ``gelmandiag`` function.

.. code-block:: julia

	## Brooks, Gelman and Rubin Convergence Diagnostic
	>>> gelmandiag(sim1, mpsrf=true, transform=true)

	5x3 Array{Any,2}:
	 ""               "PSRF"      "97.5%"
	 "beta[1]"       0.558658    0.559873
	 "beta[2]"       0.58107     0.582333
	 "s2"            1.15596     1.15848 
	 "Multivariate"  1.00221   NaN       

Values of the diagnostic that are greater than 1.2 are evidence of non-convergence.  The smaller diagnostic values for the regression example suggest that its draws have converged.
 

Posterior Summaries
^^^^^^^^^^^^^^^^^^^

Once convergence has been assessed, sample statistics can be computed on the MCMC output to estimate features of the posterior distribution.  Some of the available posterior summaries are illustrated in the code block below.

.. code-block:: julia

	## Summary Statistics
	>>> describe(sim1)

	Iterations = 252:10000
	Thinning interval = 2
	Number of chains = 3
	Samples per chain = 4875

	Empirical Posterior Estimates:
	4x6 Array{Any,2}:
	 ""          "Mean"    "SD"      "Naive SE"   "Batch SE"      "ESS"
	 "beta[1]"  0.646951  1.33333   0.0110253    0.0240354    6708.61  
	 "beta[2]"  0.788261  0.400042  0.00330794   0.00638979   7571.25  
	 "s2"       1.60026   6.616     0.0547076    0.227472     3517.35  

	Quantiles:
	4x6 Array{Any,2}:
	 ""           "2.5%"     "25.0%"    "50.0%"   "75.0%"   "97.5%"
	 "beta[1]"  -1.8012     0.0516419  0.594027  1.21848   3.24066 
	 "beta[2]"   0.0369576  0.61944    0.79773   0.965215  1.53923 
	 "s2"        0.175576   0.387191   0.674083  1.32558   7.47037 

	## Highest Posterior Density Intervals
	>>> hpd(sim1)

	4x3 Array{Any,2}:
	 ""           "2.5%"     "97.5%"
	 "beta[1]"  -1.82575    3.19734 
	 "beta[2]"   0.0734925  1.57316 
	 "s2"        0.0861001  4.69749 
	
	## Cross-Correlations
	>>> cor(sim1)

	4x4 Array{Any,2}:
	 ""           "beta[1]"    "beta[2]"    "s2"    
	 "beta[1]"   1.0         -0.899921    -0.0145811
	 "beta[2]"  -0.899921     1.0          0.0227525
	 "s2"       -0.0145811    0.0227525    1.0      

	## Lag-Autocorrelations
	>>> autocor(sim1)

	4x5x3 Array{Any,3}:
	[:, :, 1] =
	 ""          "Lag 2"   "Lag 10"    "Lag 20"    "Lag 100" 
	 "beta[1]"  0.449877  0.0867905  -0.0360931  -0.000137467
	 "beta[2]"  0.371405  0.067229   -0.0305984   0.00658498 
	 "s2"       0.725702  0.245108    0.150537   -0.0479392  

	[:, :, 2] =
	 ""          "Lag 2"   "Lag 10"     "Lag 20"    "Lag 100"
	 "beta[1]"  0.395009  0.0224026   -0.0272157  -0.00879161
	 "beta[2]"  0.303902  0.00634753  -0.0344257  -0.0140071 
	 "s2"       0.829777  0.428135     0.105249   -0.0236275 

	[:, :, 3] =
	 ""          "Lag 2"   "Lag 10"    "Lag 20"     "Lag 100"
	 "beta[1]"  0.31747   0.0143544   0.000619652  0.0142481 
	 "beta[2]"  0.295017  0.0192439  -0.0149453    0.0236515 
	 "s2"       0.811755  0.273569    0.15009      0.0561878 

	## Deviance Information Criterion
	>>> dic(sim1)

	3x3 Array{Any,2}:
	 ""      "DIC"   "Effective Parameters"
	 "pD"  13.9011  0.979245               
	 "pV"  24.4464  6.2519                 

Output Subsetting
^^^^^^^^^^^^^^^^^

Additionally, sampler output can be subsetted to perform posterior inference on select iterations, parameters, and chains.

.. code-block:: julia

	## Subset Sampler Output
	>>> describe(sim1[1000:5000, ["beta[1]", "beta[2]"], :])
	
	Iterations = 1000:5000
	Thinning interval = 2
	Number of chains = 3
	Samples per chain = 2001

	Empirical Posterior Estimates:
	3x6 Array{Any,2}
	 ""          "Mean"    "SD"      "Naive SE"   "Batch SE"      "ESS"
	 "beta[1]"  0.598793  1.19409   0.0154118    0.0268244    3448.98  
	 "beta[2]"  0.799864  0.359956  0.00464585   0.00768688   3628.13  

	Quantiles:
	3x6 Array{Any,2}:
	 ""           "2.5%"     "25.0%"    "50.0%"   "75.0%"   "97.5%"
	 "beta[1]"  -1.67395    0.0202488  0.622143  1.15556   2.95496 
	 "beta[2]"   0.0709892  0.624956   0.799999  0.976372  1.49103 
	


Computational Performance
-------------------------

Computing runtimes were recorded for different sampling algorithms applied to the regression example.  Runs wer performed on a desktop computer with an Intel i5-2500 CPU @ 3.30GHz.  Results are summarized in the table below.  Note that these are only intended to measure the raw computing performance of the package, and do not account for different efficiencies in the sampling algorithms.

.. table:: Number of draws per second for select sampling algorithms in `MCMCsim`.

	+---------------------+-----------------------+--------+--------+--------+
	| Adaptive Metropolis                         |        |        |        |
	+---------------------+-----------------------+        |        |        |
	| Within Gibbs        | Multivariate          | Gibbs  | NUTS   | Slice  |
	+=====================+=======================+========+========+========+
	| 5,000               | 7,500                 | 15,000 | 1,200  | 4,300  |
	+---------------------+-----------------------+--------+--------+--------+

	
Development and Testing
-----------------------

Command-line access is provided for all package functionality to aid in the development and testing of models.  Examples of available functions are shown in the code block below.  Documentation for these and other related functions can be found in the :ref:`section-MCMC-Types` section. 

.. code-block:: julia

	## Development and Testing

	setinputs!(line, data)             # Set input node values
	setinits!(line, inits[1])          # Set initial values
	setsamplers!(line, scheme1)        # Set sampling scheme

	showall(line)                      # Show detailed node information

	logpdf(line, 1)                    # Log-density sum for block 1
	logpdf(line, 2)                    # Block 2
	logpdf(line)                       # All blocks

	simulate!(line, 1)                 # Simulate draws for block 1
	simulate!(line, 2)                 # Block 2
	simulate!(line)                    # All blocks

In this example, functions ``setinputs!``, ``setinits!``, and ``setsampler!`` allow the user to manually set the input node values, the initial values, and the sampling scheme form the ``line`` object, and would need to be called prior to ``logpdf`` and ``simulate!``.  Updated model objects should be returned when called; otherwise, a problem with the supplied values may exist.  Method ``showall`` prints a detailed summary of all model nodes, their values, and attributes; ``logpdf`` sums the log-densities over nodes associated with a specified sampling block (second argument); and ``simulate!`` generates an MCMC draw for the nodes.  Non-numeric results may indicate problems with distributional specifications in the second case or with sampling functions in the last case.  The block arguments are optional; and, if left unspecified, will cause the corresponding functions to be applied over all sampling blocks.  This allows testing of some or all of the samplers.