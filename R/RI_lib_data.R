#' establishing new library data by calculating relative intensity
#'
#'
#' @export

RI_lib_data = function(reference_spectra , RT_mass, acceptable_RI = acceptable_RI) {

  #calculating relative intensity of ions for mass-mass
  # selected Relative intensity > 10 with correlated mass in lib
  RI_lib1 = lapply ( 1:length(reference_spectra) ,function(x) {reference_spectra[[x]][,2]/ max(reference_spectra[[x]][,2])* 100 })
  names(RI_lib1)= names(reference_spectra)
  RI_lib = reference_spectra
  for (x in 1: length(RI_lib)) {
    RI_lib[[x]][,2] = RI_lib1[[x]]
  }
  for(e in 1:length(RI_lib)) {
    RI_lib[[e]] = RI_lib[[e]][RI_lib[[e]][,2] >= acceptable_RI, , drop=FALSE]
  }


  #establishing new library with the same length of RT_mass
  RI_lib = RI_lib[match(names(RT_mass) ,names(RI_lib))]
  all(names(RI_lib) == names(RT_mass))



  return(RI_lib)
}
