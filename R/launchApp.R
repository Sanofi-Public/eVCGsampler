
#' @title Launch the Shiny App
#' @description Starts the interactive eVCGsampler app: upload the POOL and TG data,
#' choose the covariates and strata, run the pre-tests, draw the VCG and download
#' an HTML report and the sampled data.
#' @param ... Further arguments passed to \code{\link[shiny]{runApp}}, e.g. \code{port} or \code{launch.browser}.
#' @details The app needs some packages that are not required by the rest of eVCGsampler
#' (listed under Suggests). Missing ones are reported together with the command to install them.
#' Rendering the HTML report additionally needs pandoc (it ships with RStudio).
#' @examples
#' \dontrun{
#' eVCGsampler::launchApp()
#' }
#' @export


launchApp <- function(...) {
  appDir <- system.file("app", package = "eVCGsampler")
  if (appDir == "") {
    stop("Could not find app directory. Try re-installing the package.", call. = FALSE)
  }
  app_pkgs <- c("bslib", "DT", "glue", "knitr", "openxlsx", "readxl", "rmarkdown",
                "shinycssloaders", "shinyhelper", "shinyjs")
  missing  <- app_pkgs[!vapply(app_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("The app needs the following packages: ", paste(missing, collapse = ", "),
         ".\nInstall them with: install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))",
         call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal", ...)
}
