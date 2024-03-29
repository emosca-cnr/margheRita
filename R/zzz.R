#' @importFrom utils packageVersion
#' 
.onAttach <- function(libname, pkgname){
    packageStartupMessage("margheRita ", packageVersion("margheRita"), "\n")
}