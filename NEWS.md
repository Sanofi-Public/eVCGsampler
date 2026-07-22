# eVCGsampler 1.1.2

* `VCG_sampler()` now shortens long covariate names in the target plot to a maximum of 11 characters, built from the first 8 characters + `_` + the last 2 characters to keep them unique. A warning is printed when names are cut.

* `VCG_sampler()` now issues sanity warnings when running: too many covariates (> 21), a small POOL size relative to the number of covariates (ratio < 5, reported), and a requested VCG size that is >= 75% of the POOL. When stratification is used, these checks are reported per stratum, plus a warning is added when stratum sizes are strongly imbalanced.

* Covariate scaling now uses a robust median/MAD approach (falling back to the standard deviation when the MAD is 0). `robust_scale()` gained an optional `group` argument (defaulting to `NULL`, using all observations) and now backs the scaling in `VCG_sampler()`, `energy_distance()`, `combine_variables()` and `find_outliers()`, replacing base `scale()`.

* Unit test added via testthat framework
