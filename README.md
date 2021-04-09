# microsim
A simulator package in R containing a collection of basic community models that outputs count data of microbial communities, integrating various statistical and computational methodologies for minimal runtime.

### Installation

Most recent development version can be pulled from this repository and installed directly from R command line.

First load the `devtools` package.

```{R}
library(devtools)
```

Then install from Github with following command.

```{R}
devtools::install_github("gheysenemma/microsimR")
library(microsim)
```

### Content Summary

The most recent version of this package contains the following simulation models:

1) The Generalised Lotka-Voltera Model: `microsim::glv()`
2) The Hubbell Neutral Model: `microsim::hubbell()`
3) The Self Organised Instability Model: `microsim::soi()`

Other functions included in the package:

1) RMSE-based sampling function given a time series: `microsim::rmse_sample()`
2) `phyloseq` wrapper function (samples / full time series): `microsim::asPhyloseq()`
3) Power-law interaction matrix decomposition `microsim::powerlawA()`
4) Random interaction matrix by uniform distribution `microsim::randomA()`

### Acknowledgements

This package is constructed as part of a master's thesis at KU Leuven 2020-2021.
- Author: Emma Gheysen

Under supervision of:
- Promotor: Prof. Karoline Faust (KU Leuven)
- Co-promotor: Prof. Leo Lahti (Turku University)
- Mentor: Daniel Rios Garza (KU Leuven)
