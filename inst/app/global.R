# =====================================================================
# global.R - eVCGsampler Shiny app (started with eVCGsampler::launchApp())
#
# Sourced once at app start. Everything under R/ is auto-loaded by Shiny
# afterwards, so only packages, options and constants belong here.
# =====================================================================

library(shiny)
library(DT)
library(ggplot2)
library(shinycssloaders)

# The statistics engine. Attached rather than only namespaced so the
# ggplot objects it returns print with the package's own theme.
library(eVCGsampler)


# The app ships with the package, so it carries the package version.
app_version <- as.character(utils::packageVersion("eVCGsampler"))

# App colour palette. The first five are black & white print
# safe, the first eight are colour-blind safe. eVCGsampler draws its own
# plots with the same first three tones (POOL, VCG, TG).
palette_colors <- c(
    "#7A00E6",  # Chasing Purple - TG
    "#FE7500",  # Orange         - VCG
    "#00A091",  # Cyan           - POOL
    "#F5C142",
    "#694A3E",
    "#1F70A6",
    "#EBA0AB",
    "#6CADD1",
    "#1A7A16",
    "#5F7786",
    "#FF6D6A",
    "#23004C"
)

# Name of the POOL / TG indicator that combine_data() creates when the two
# uploads are merged. Kept in one place because the Import, Setup and
# Report steps all refer to it.
group_indicator <- "POOL_vs_TG"

# eVCGsampler refuses more than eight strata (VCG_sampler() and
# energy_test() both stop at nlevels > 8).
max_strata <- 8

# Shared selectize options for the multi-selects (covariates, strata, outcome
# variables). The remove_button plugin puts an "x" on every selected chip -
# without it a wrongly clicked variable can only be taken out with backspace,
# which is not discoverable.
multi_select_options <- list(
    plugins = list("remove_button")
)

# Example data for trying the app out without hunting for a data set - see
# the "Load the example data" link on the Data step. Regenerate with
# source("data/make_example_data.R").
example_data <- list(
    pool = "data/example_POOL.csv",
    tg   = "data/example_TG.csv"
)

options(
    spinner.type  = 4,
    spinner.color = "#7A00E6",
    # Historical control pools are wide but rarely huge; 30 MB covers the
    # xlsx exports we see in practice.
    shiny.maxRequestSize = 30 * 1024^2
)

