# BKTR (development version)

# BKTR 0.2.0

* Fix the `torch_triangular_solve` deprection warnings using the new `linalg_solve_triangular` function (mlverse/torch#1146). Also, update the minimal `torch` package version to 0.13.0 to be able to use this new function.

* Add api keys for the `plot_spatial_betas` function to align with the new requirements of `ggmap` for stadia maps (dkahle/ggmap#350).

* Call garbage collection directly during MCMC sampling when using only CPU and the covariates are large to avoid memory issues.

* Detail and improve the `torch` installation instructions in the README file.

* Use the `is_light` argument with the `BixiData$new()` method in the README example to only run the example on 25 stations and 50 days of data.

* Update the vignette of the BKTR presentation article to also reflect all the changes of this version.


# BKTR 0.1.1

* Initial release of the package on CRAN.

* Added a `NEWS.md` file to track changes to the package.
