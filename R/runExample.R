#' Title
#'
#' @return
#' @export
#'
#' @examples
#'
runExample <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "GenomeTornadoPlot")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `myshinyapp`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
