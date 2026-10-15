# eVCGsampler 1.2.0

* `launchApp()` now starts a new, modernised Shiny app with a guided five-step workflow (Data, Setup, Pre-tests, Sampling, Report), example data, an HTML report and Excel/CSV exports. The app's extra packages are listed under Suggests; `launchApp()` reports any missing ones, and passes `...` on to `shiny::runApp()`.

* `VCG_sampler()`: the size sanity warnings now count the POOL as the first level of the group variable, so they are also correct for `combine_data()` output ("POOL"/"TG") instead of reporting a POOL size of 0.

* `VCG_sampler()` now stops with a clear message when `n` is larger than the POOL with complete data, or when one of the groups is empty, instead of silently selecting the whole POOL.

* `VCG_sampler()` with stratum now returns the rows in the input order and keeps rows with a missing stratum value (with `VCG = NA`), as is already done for rows with missing covariates.

* `BestVCGsize()` no longer fails for small POOLs or when no plateau is found (it returns a single `NA`), only scans sizes up to the POOL size with complete data, and no longer repeats the `VCG_sampler()` size warnings for every tried size.

* `robust_scale()` ignores `NA` when detecting binary/categorical variables and no longer returns `Inf`/`NaN` when the reference group is constant.

* `plot_var()` now works with a single stratum level and with character stratum variables.

# eVCGsampler 1.1.2

* `VCG_sampler()` now shortens long covariate names in the target plot to a maximum of 11 characters, built from the first 8 characters + `_` + the last 2 characters to keep them unique. A warning is printed when names are cut.

* `VCG_sampler()` now issues sanity warnings when running: too many covariates (> 21), a small POOL size relative to the number of covariates (ratio < 5, reported), and a requested VCG size that is >= 75% of the POOL. When stratification is used, these checks are reported per stratum, plus a warning is added when stratum sizes are strongly imbalanced.

* Covariate scaling now uses a robust median/MAD approach (falling back to the standard deviation when the MAD is 0). `robust_scale()` gained an optional `group` argument (defaulting to `NULL`, using all observations) and now backs the scaling in `VCG_sampler()`, `energy_distance()`, `combine_variables()` and `find_outliers()`, replacing base `scale()`.

* Unit test added via testthat framework
