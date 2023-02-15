#' Filter metabolite-feature pairs
#' @param  mRList mRList object
#' @export
#' 

filter_metabolite_associations <- function(mRList=NULL){
	
	
	out_levels <- mRList$metabolite_identification$associations
	
	cat("Filtering...\n")
	cat("before filtering:", nrow(out_levels), "\n")
	out_levels <- unique(out_levels[, -which(colnames(out_levels) == "ID_peaks")]) #remove the redundant ID_peak
	
	for(lev in c(1, 2)){
		lev_m <- out_levels$ID[out_levels$Level==lev]
		if(length(lev_m) > 0){
			idx_rm <- which(out_levels$ID %in% lev_m & out_levels$Level > lev)
			if(length(idx_rm) > 0){
				out_levels <- out_levels[-idx_rm, ]
			}
		}
	}
	
	#level 3 mz,rt annotated also as level 3 mz
	lev_m <- out_levels$ID[out_levels$Level==3 & out_levels$Level_note ==  "mz, rt"]
	if(length(lev_m) > 0){
		idx_rm <- which(out_levels$ID %in% lev_m & out_levels$Level_note == "mz")
		if(length(idx_rm) > 0){
			out_levels <- out_levels[-idx_rm, ]
		}
	}
	cat("filtered:", nrow(out_levels), "\n")
	
	mRList$metabolite_identification$associations <- out_levels
	
	return(mRList)
		
}