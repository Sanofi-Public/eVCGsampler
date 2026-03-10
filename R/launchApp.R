
#'  @title Launch the Shiny App
#'
#' @examples
#' \dontrun{
#' eVCGsampler::launchApp()
#' }
#' @export


launchApp <- function() {
  appDir <- system.file("app", package = "eVCGsampler")
  if (appDir == "") {
    stop("Could not find app directory. Try re-installing the package.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}

