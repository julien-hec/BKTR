
<!-- README.md is a file that contains the information that will be displayed on the GitHub repository page. -->

# BKTR

## Intro
This project is a *R* implementation of the BKTR algorithm presented by Mengying Lei, Aurélie Labbe & Lijun Sun (2023).
The article presenting the algorithm can be found [here](https://arxiv.org/abs/2109.00046).

BKTR stands for Scalable Spatiotemporally Varying Coefficient Modelling with Bayesian Kernelized Tensor Regression.
It allows to model spatiotemporally varying coefficients using a Bayesian framework.
We implemented the algorithm and more in a *R* package that uses [torch](https://torch.mlverse.org/) as a tensor operation backend.

For information, an alternative *Python* implementation of the algorithm can be found [here](https://github.com/julien-hec/pyBKTR).
The *Python* implementation is synchronized with this repository and development is done in parallel.
The synchronization of features will be done at a subrevision level (x.y.0).

An article presenting the *R* package in details is currently in preparation and should be available soon.

## Installation

### CRAN Installation
```r
install.packages('BKTR')
```

### Latest Development Version Installation
The latest development version on GitHub can be installed using the `devtools` package:
```r
library(devtools)
devtools::install_github('julien-hec/BKTR', ref = 'main')
```

### Notes
If you obtain an error message when installing the package, it may be due to the installation of the `torch` package.
A common error message that can appear during `BKTR` installation is `installation of package 'BKTR' had non-zero exit status`.
The `torch` package is a dependency of the `BKTR` package and there is a good chance that the error comes from the installation of `torch`.
Because of its ability to perform tensor operations on the GPU, it can sometimes be more complicated to install than other R packages.
We provide some guidance for the installation of `torch` below.

#### Installing torch alone
A simple way to see if BKTR installation problems come from the torch installation is to try to install torch alone first:
```r
install.packages('torch')
```
If you obtain an error message, we encourage you to continue reading the following subsections.

#### Installing torch with a non-interactive R session
If you use a non-interactive R session (e.g. in a Docker container), you need to install LibTorch and LibLantern afterwards with the following command:
```r
library(torch)
torch::torch_install()
```

#### Installation of torch for CPU only
If you have a CUDA version that causes issues during the torch installation and you just want to use the CPU version of BKTR, you can install torch with the `CPU` option:
```r
Sys.setenv(CUDA='cpu')
install.packages('torch')
```

#### Installing torch with a specific CUDA version
If your `CUDA` version does not seem to be supported correctly and you obtain the following error message:
```
Error in `check_supported_version()`:
x Unsupported CUDA version "12.2"
i Currently supported versions are: "11.7" and "11.8".
```
As specified in the [prebuilt section of torch's installation guide](https://torch.mlverse.org/docs/articles/installation.html#pre-built), you can try to install from specific precompiled binaries for another `CUDA` version:
```r
options(timeout = 600) # increasing timeout since we download a 2GB file.
# For Windows and Linux: "cpu", "cu117", "cu118" are the only currently supported
# For MacOS the supported are: "cpu-intel" or "cpu-m1"
kind <- "cu118"
version <- available.packages()["torch","Version"]
options(repos = c(
  torch = sprintf("https://torch-cdn.mlverse.org/packages/%s/%s/", kind, version),
  CRAN = "https://cloud.r-project.org" # or any other from which you want to install the other R dependencies.
))
install.packages("torch")
```

#### More information on torch installation
For more information on how to install torch, please refer to the [torch installation vignette](https://cran.r-project.org/package=torch/vignettes/installation.html).

#### Getting started with Colab
If you want to get started quickly with `BKTR` on Google Colab, you can use the following examples
- BKTR on Light BIXI Data with CPU ([Colab](https://colab.research.google.com/drive/1nfhXJTrjpJWNd8q6XP-TKCjuIhoIOfrP?usp=sharing) & [GitHub](https://github.com/julien-hec/bktr-examples/blob/main/BKTR-installations/R_BKTR_CPU.ipynb))
- BKTR on Light BIXI Data with GPU ([Colab](https://colab.research.google.com/drive/1tM5nmDYPmRWhGFLOq8Qf4a7m4bKGfDfD?usp=sharing) & [GitHub](https://github.com/julien-hec/bktr-examples/blob/main/BKTR-installations/R_BKTR_GPU.ipynb))

## Simple Example
To verify that everything is running smoothly you can try to run a BKTR regression on the BIXI data presented in the package. (The data is already preloaded in the package via the `BixiData` [R6](https://r6.r-lib.org/articles/Introduction.html) class). To use a subset of the BIXI dataset as a simple example, we can also use the `is_light` argument of the `BixiData$new()` method to only run our example on 25 stations and 50 days of data.

The following code will run a BKTR regression using sensible defaults on the simplified BIXI data and print a summary of the results.
```r
library(BKTR)
bixi_data <- BixiData$new(is_light=TRUE)
bktr_regressor <- BKTRRegressor$new(
    data_df=bixi_data$data_df,
    spatial_positions_df=bixi_data$spatial_positions_df,
    temporal_positions_df=bixi_data$temporal_positions_df,
    burn_in_iter=200,
    sampling_iter=200
)
bktr_regressor$mcmc_sampling()
summary(bktr_regressor)
```

## Contributing
Contributions are welcome. Do not hesitate to open an issue or a pull request if you encounter any problem or have any suggestion.
