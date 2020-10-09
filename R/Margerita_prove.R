#########################     MargheRita    ################################
############################################################################


## 1) Set directory, where all your files are: xls/cvs and metadata
#Session->Set Directory -> chose directory

## 2) Libraries needed:


Margerita_prove <- function(){

  library(readxl)
  library(ggplot2)
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(kableExtra)
  library(gridExtra)
  library(knitr)
  library(RColorBrewer)
  library("rjson")

  ## 3) import data from xls or cvs, object with data is farina

  #FARINA ---------------------------------------------------

  farina_pre<-read_excel(path="Margerita Trial.xlsx",sheet=1)

  str(farina_pre)
  summary(farina_pre)

  #Eliminate useless column, convert everything as numeric,
  #create column with feature identifiers

  farina<-farina_pre %>%
    dplyr::select(-1,-4) %>%
    mutate(feat=paste(round(Average_mz,4),round(Rt_min,2),sep="_")) %>%
    dplyr::select(-Average_mz,-Rt_min) %>%
    gather(sample,value,-feat) %>%
    spread(feat,value)


  glimpse(farina)
  head(farina)


  ## 4) Factors & Metadata definition select way 1 or 2 or 3----------------------

  #ACQUA    metadata call 1: from excel ##################################

  acqua<-read_excel(path="metadata.xlsx",sheet=1)



  #ACQUA    metadata call 2: from name ###################################

  sample_name<-farina$sample
  injection_order<-substr(sample_name,start=1,stop=4)
  QC_ID<-grep("QC",sample_name)
  diet_type<-grep("NOD_STAND","BAL_STAND","NOD_OMEGA",sample_name)
  technical_rep<-grep("_rep1_","_rep2_","_rep3_","_rep4_","_rep5_",sample_name)


  #diet (QC are explicitly indicated)
  diet_type<-substr(sample_name,start=6,stop=14)
  diet_type[QC_ID]<-"QC"

  #biological replicates (QC are explicitly indicated)
  biological_rep<-substr(sample_name,start=16,stop=18)
  biological_rep[QC_ID]<-"QC"

  #Technical repllicates (QC are explicitly indicated)
  technical_rep<-substr(sample_name,start=20,stop=23)
  technical_rep[QC_ID]<-"QC"


  #ACQUA    metadata call 3: from name ###################################

  farina$diet = NA
  farina$diet[grepl("_BAL_STAND_", rownames(farina))] = "BAL_STAND"
  farina$diet[grepl("_NOD_STAND_", rownames(farina))] = "NOD_STAND"
  farina$diet[grepl("_NOD_OMEGA_", rownames(farina))] = "NOD_OMEGA"
  farina$diet[QC_ID]<-"QC"

  farina$reptech_id = NA
  farina$reptech_id[grepl("_rep1_", rownames(farina))] = "rep1"
  farina$reptech_id[grepl("_rep2_", rownames(farina))] = "rep2"
  farina$reptech_id[grepl("_rep3_", rownames(farina))] = "rep3"
  farina$reptech_id[grepl("_rep4_", rownames(farina))] = "rep4"
  farina$reptech_id[grepl("_rep5_", rownames(farina))] = "rep5"
  farina$reptech_id[QC_ID]<-"QC"

  farina$repbiol_id = NA
  farina$repbiol_id[grepl("_m01_", rownames(farina))] = "m01"
  farina$repbiol_id[grepl("_m02_", rownames(farina))] = "m02"
  farina$repbiol_id[grepl("_m03_", rownames(farina))] = "m03"
  farina$repbiol_id[grepl("_m04_", rownames(farina))] = "m04"
  farina$repbiol_id[QC_ID]<-"QC"


  write.csv(farina, "x.csv")
  farina_2 = read.csv("x.csv")

  farina_2$diet_person = paste(farina_2$diet, farina_2$person, sep="_")
  farina_2$diet_reptech_id = paste(farina_2$diet, farina_2$reptech_id, sep="_")




  ## 5) Data inspection - first PCA for general check of outliers


  row.names(farina)<-injection_order
  farina$biological_rep<-biological_rep
  farina$diet_type<-diet_type
  farina$technical_rep<-technical_rep
  farina<-farina[,-grep("sample",names(farina))]
  str(farina)

  pca<-PCA(farina,graph = FALSE,
           quali.sup = c(grep("biological_rep",colnames(farina)),
                         grep("diet_type",colnames(farina)),
                         grep("technical_rep",colnames(farina))))

  fviz_pca_ind(pca,
               #addEllipses = TRUE,
               label="none", #to remove labels from plot
               col.ind="cos2", #TO COLOR by importance for projection
               habillage = "diet_type",
               repel = FALSE) #to avoid overlapping labels: check the output


  ######## Data inspection - PCA for batch effect check

  check_batch<-injection_order %>% gsub("x","",.) %>% as.numeric
  check_batch<-ifelse(check_batch <55,"<55",">55") %>% factor

  fviz_pca_ind(pca,
               addEllipses = TRUE,
               #label="none", #to remove labels from plot
               #col.ind="cos2", #TO COLOR by importance for projection
               habillage =check_batch,
               repel = FALSE) #to avoid overlapping labels: check the output






  ##  6) Missing value imputation: treat farina as numeric->check min value for this dataset
  ##                               ->change all 0 into NA -> compute NA with random number

  farina_comp <- sapply(farina[-1],as.numeric)
  min(farina_comp, na.rm=TRUE)

  farina_comp[farina_comp==0] <- NA
  min(farina_comp, na.rm=TRUE)

  sum(is.na(farina_comp))
  set.seed(1234)
  farina_comp[is.na(farina_comp)] <- sample(0.005:0.5, size=sum(is.na(farina_comp)), replace=TRUE)
  min(farina_comp)

  farina_comp <- farina[order(rownames(farina)),]

  ##ASK ETTORE why column disappear







  # Cleaning CV-QC > CV-samples -----------------------------------------------------------

  # Definition of Coeffitient of Variation CV function

  CV <- function(x)((100*sd(x))/(mean(x)))

  # Compute the CV intensities of each feature across QC samples
  QC <- subset(farina_2, farina_2$diet_type=="QC")
  table <- data.frame(apply(QC, 2, CV))
  colnames(table) <- "QC"
  rm(QC)






  # Compute the CV intensities of each feature across study samples
  samples <- subset(farina_2, farina_2$diet_type==diet_type)
  abcd <- data.frame(apply(samples, 2, CV))
  colnames(abcd) <- "samples"
  rm(samples)

  table$id <- rownames(table)
  abcd$id <- rownames(abcd)
  table <- merge(table, abcd, by="id")
  rm(abcd,CV)

  # Select those features where CV(samples) > CV(QC)
  table$noise[table$samples>table$QC] <- "ok"
  table$noise[table$samples<table$QC]<- "noise"
  t <- subset(table, noise=="ok")
  rownames(t) <- t$id
  t <- rownames(t)
  data <- data[,t]
  CV.QC <- table
  rm(t, table)

  dim$CV <- ncol(data)
  dim.p$CV <- (dim$'m/z .4-.8'-dim$CV)*100/dim$total



















  ## CLEAN
  rm(list = ls())
  dev.off()

}
