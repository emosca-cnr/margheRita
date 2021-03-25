#' margherRita - main function to run the full pipeline
#'
#'
#'

margheRita <- function(){



  #preprocessing (Marynka): read input files: data file + metadata

#read_input_files()
#PCA_pre_processing()
#M/Z_filtering()
#inputation()
#CV()


  #filtering:
  #a)	Delete features that have mass defect 0.4-0.85 (decimals)
  #b)	Remove features that have Coefficient of Variation in QC > CV than study samples
  #c)	Select features found in at least 20% samples of at least one group (check with metadata what are observation groups)


  #normalization
#pqn()


  #sample aggregation - stat x i replicati tecnici


  # statistical analysis (check Uni Alberta - metaboanalyst)
#univarite(?)

  #outputs



}
