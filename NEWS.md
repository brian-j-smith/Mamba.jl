# News

## Version Updates

### 0.6.1
* Compatibility updates for julia 0.4 prerelease.

### 0.6.0
* Stable release for julia 0.4.
* Added support for the specification of ``MultivariateDistribution`` arrays in stochastic nodes.
* Arrays in stochastic nodes must now be declared as a ``UnivariateDistribution[]`` or as a ``MultivariateDistribution[]``.  In previous package versions, arrays could be declared as a generic ``Distribution[]`` array.  This is no longer allowable due to the need to distinguish between arrays of univariate and multivariate distributions.
* The following changes were made to internal data structures and generally do not affect the user interface (i.e., model specification, sampling function calls, and sampler output diagnostics and summaries):
    * The ``VariateType``, which aliases ``Float64``, is deprecated and will be removed in a future version.
    * The abstract ``Variate`` type was separated into ``ScalarVariate`` and ``ArrayVariate`` types which are subtypes of ``Real`` and ``DenseArray``, respectively.
    * ``AbstractVariate`` was defined as the ``Union(ScalarVariate,ArrayVariate)``.
    * The ``Logical`` type was separated into ``ScalarLogical`` and ``ArrayLogical`` types which are subtypes of ``ScalarVariate`` and ``ArrayVariate``, respectively.
    * ``AbstractLogical`` was defined as the ``Union(ScalarLogical,ArrayLogical)``.
    * The ``Stochastic`` type was separated into ``ScalarStochastic`` and ``ArrayStochastic`` types which are subtypes of ``ScalarVariate`` and ``ArrayVariate``, respectively.
    * ``AbstractStochastic`` was defined as the ``Union(ScalarStochastic,ArrayStochastic)``.
    * ``AbstractDependent`` was defined as the ``Union(AbstractLogical,AbstractStochastic)``.
    * The ``nlink`` field of logical and stochastic types was renamed to ``linklength``.
    * An abstract ``AbstractChains`` type was implemented, the ``model`` field removed from the ``Chains`` type, and a new ``ModelChains`` type created to provide the ``model`` field.
* Added an example of sampling different parameter blocks with different stand-alone samplers (``amwg!`` and ``slice!``).
* Removed the ``insupport`` method for stochastic types.

### 0.5.2
* Applied *Mamba* changes through 0.4.12.

### 0.5.1
* Applied *Mamba* changes through 0.4.7 and updated compatibility with julia 0.4.0-dev.

### 0.5.0
* Branched off of *Mamba* 0.4.4.
* Initial release for the under-development version of julia 0.4.  *Mamba* 0.5.x releases exist to incorporate changes being made in the nightly builds of julia, and should be considered unstable.  They may contain compatibility issues or serious bugs.  Most users are advised to use *Mamba* 0.4.x releases for julia 0.3 or to wait for stable *Mamba* 0.6.x releases for julia 0.4.

### 0.4.12
* Updated documentation for user-defined distributions.

### 0.4.11
* Implemented `Chains` method function ``changerate`` to calculate parameter state space change rates (per iteration).
* Updated `Truncated` constructor for `Flat` distributions for compatibility with latest *Distributions* package.
* Simplified documentation instructions for user-defined univariate distributions.

### 0.4.10
* Optimized code and improved handling of sampler output in the Gelman convergence diagnostic.
* Added ``ask`` argument to the ``draw`` plot method.

### 0.4.9
* Added Heidelberger and Welch, and Raftery and Lewis convergence diagnostics.
* Added documentation and illustrations for all included diagnostics.
* Added documentation for creating user-defined distributions.

### 0.4.8
* Fixed `BoundsError()` occurring with autocorrelation plots.

### 0.4.7
* Improved formatting and outputting of posterior summary statistics
* Improved efficiency of DOT formatting of `Model` graphs.
* Exported `graph2dot` function to allow in-line processing of `Model` graphs with GraphViz package.

### 0.4.6
* Added `verbose` argument to `mcmc` methods to suppress printing of sampler progress to the console.
* Fixed calculation of effective sample size.

### 0.4.5
* Replaced the `Scale.discrete_color` function deprecated in the *Gadfly* package with `Scale.color_discrete`.

### 0.4.4
* Require julia 0.3.4 to remove version restriction on the *Colors* package.
* Call new *Distributions* package methods to get `InverseGamma` shape and scale parameters in the tutorial example.

### 0.4.3
* Fixed `ERROR: too many parameters for type Truncated`.

### 0.4.2
* Added support for optional arguments in `Chains` plot method.
* Implemented direct grid sampling for discrete univariate stochastic nodes with finite support.

### 0.4.1
* Updated and documented `predict` function as an official part of the package.
* Reorganized `Chains` methods documentation.

### 0.4.0

* Added support for user add-on packages and functions to allow for their inclusion in `Model` specifications.
* Added experimental `predict` (posterior prediction) function.
* Required the *Cairo* package.
* Removed deprecated `MCMC*` types and `slicewg` and `SliceWG` functions.
* Fixed `ERROR: GenericMvNormal not defined`.
* Distributions `DiagNormal` and `IsoNormal` removed and replaced by `MvNormal`.
* Distributions `DiagNormalCanon` and `IsoNormalCanon` removed and replaced by `MvNormalCanon`.
* Distributions `DiagTDist` and `IsoTDist` removed and replaced by `MvTDist`.

### 0.3.8

* Updated to fix warning and work with the latest versions of the *PDMat* and *Distributions* packages.

### 0.3.7

* Extend `Chains` draw method to allow automatic outputting of multiple plot grids to different files.
* Add `Chains` plot method to accommodate vectors of plot types.
* Fix variance calculation in `gewekediag()`.

### 0.3.6

* Fix for convert errors triggered by the *Color* package beginning with its version 0.3.9.

### 0.3.5

* Documentation updates only - primarily the addition of results to examples.
* No changes made to the source code.

### 0.3.4

* Added distributions documentation.
* Added jaws repeated measures analysis of variance example.

### 0.3.3

* Fixed the `rand` method definition error (`type DataType has no field body`) that began occurring with late release candidates and final release of julia 0.3.

### 0.3.2

* Fixed tuning parameter overwrites occurring with `pmap()` in single-processor mode.

### 0.3.1

* Added `chains` field to `Chains` type for tracking purposes.
* Fixed `mcmc` to accommodate restarting of chains subsetted by parameters and/or chains.
* Fixed plot legends to properly reference the chains being displayed.
* Added support for sampling of positive-definite matrices specified with Wishart or InverseWishart distributions.
* Added a block-diagonal normal (`BDiagNormal`) distribution.

### 0.3.0

* Implemented restarting of MCMC chains.
* Deprecated `slicewg` and `SliceWG`.  Replaced with `:univar` option to `slice` and `Slice`.

### 0.2.1

* Updated documentation.
* Simplified parallel code.

### 0.2.0

* Automatically load *Distributions* package.
* Implemented parallel execution of parallel chains on multi-processor systems.

### 0.1.0

* Removed `MCMC` prefix from type names, and deprecated `MCMC*` types.

### 0.0.2

* Renamed package from *MCMCsim* to *Mamba*.

### 0.0.1

* Initial public release.
