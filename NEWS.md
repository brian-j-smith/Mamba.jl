# News

## Version Updates


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