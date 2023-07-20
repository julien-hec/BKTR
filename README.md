# BKTR

## Intro
This project is a R implementation of the BKTR algorithm presented by Mengying Lei, Aurélie Labbe & Lijun Sun (2023).
The article presenting the algorithm can be found [here](https://arxiv.org/abs/2109.00046).

BKTR stands for Scalable Spatiotemporally Varying Coefficient Modelling with Bayesian Kernelized Tensor Regression.
It allows to model spatiotemporally varying coefficients using a Bayesian framework.
We implemented the algorithm and more in a R package that uses [torch](https://torch.mlverse.org/) as a tensor operation backend.

For information, an alternative Python implementation of the algorithm can be found [here](https://github.com/julien-hec/pyBKTR). The Python implementation is synchronized with this repository and development is done in parallel.

An article presenting the R package in details is currently in preparation and should be available soon.

## Installation
The package is not yet on CRAN, but can be installed from GitHub.
The user can install the library simply using devtools as follows:
```r
library(devtools)
devtools::install_github("julien-hec/BKTR", ref = "main")
```
---
**NOTE**

Installing torch can be straightforward on some systems and more complicated on others. Please refer to the [torch GitHub repo](https://torch.mlverse.org/docs/installation.html) if you encounter any problems regarding its installation.

---

## Simple Example
To verify that everything is running smooth you can try to run a BKTR regression on the BIXI data presented in the package. (The data is already preloaded in the package in the `BixiData` [R6](https://r6.r-lib.org/articles/Introduction.html) class).

The following code will run a BKTR regression using sensible defaults on the BIXI data and print a summary of the results.
```r
library(BKTR)
bixi_data <- BixiData$new()
bktr_regressor <- BKTRRegressor$new(
    data_df=bixi_data$data_df,
    spatial_positions_df=bixi_data$spatial_positions_df,
    temporal_positions_df=bixi_data$temporal_positions_df
)
bktr_regressor$mcmc_sampling()
summary(bktr_regressor)
```

## Contributing
Contributions are welcome. Do not hesitate to open an issue or a pull request if you encounter any problem or have any suggestion.
