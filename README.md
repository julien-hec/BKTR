# BKTR
This project is a R implementation of the BKTR algorithm (by Mengying Lei, Aurelie Labbe, Lijun Sun).

BKTR stands for Scalable Spatiotemporally Varying Coefficient Modelling with Bayesian Kernelized Tensor Regression.

### Quick Start
To install and use the library it is as simple as
```r
library(devtools)
devtools::install_github(
    "julien-hec/BKTR",
    ref = "main",
    auth_token = auth_token, 'YOUR_AUTH_TOKEN_EAAEA777777'
    dependencies=FALSE
)
```
To verify that everything is running smooth you can try to run the BIXI example from the vignette. (The data is already preloaded in the package)
```r
library(BKTR)
ex_bktr_config <- BKTR::BKTRConfig$new(
    rank_decomp = 10,
    burn_in_iter = 5,
    max_iter = 10,
    decay_max_hparam_val = log(2),
    period_max_hparam_val = log(2),
    torch_device = "cpu"
)
ex_bktr_regressor <- BKTR::BKTRRegressor$new(
    bktr_config = ex_bktr_config,
    temporal_covariate_mat = bixi_weather_matrix,
    spatial_covariate_mat = bixi_station_matrix,
    spatial_distance_mat = bixi_distance_matrix,
    y = bixi_y,
    omega = bixi_omega
)
bktr_estimates <- ex_bktr_regressor$mcmc_sampling()
```

### Documentation Generation
From the project root folder run the following line to regenerate the documentation
```r
library(devtools)
```
then to build the manual
```r
devtools::document()
```
and to build the vignettes
```r
devtools::build_vignettes()
```
