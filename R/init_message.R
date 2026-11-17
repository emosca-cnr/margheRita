#' @importFrom utils packageVersion
#' 
.onLoad <- function(libname, pkgname){
  packageStartupMessage("margheRita ", packageVersion("margheRita"), "\n")
	cat("If you use margheRita plsase cite: Mosca, E., et al. MargheRita: streamlining MS-DIAL output analysis and metabolite identification in R. bioRxiv 2024.06.20.599545; doi: 10.1101/2024.06.20.599545\n\n")
	
}