#' @importFrom utils packageVersion
#' 
.onAttach <- function(libname, pkgname){
  packageStartupMessage("If you use ", pkgname, " v", utils::packageVersion("margheRita"), ", please cite:\n\n Mosca, E., et al. MargheRita: streamlining MS-DIAL output analysis and metabolite identification in R. bioRxiv 2024.06.20.599545; doi: 10.1101/2024.06.20.599545\n\n")
}