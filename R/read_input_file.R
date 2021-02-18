#Read files (input and metadata) both excel and .csv, default .csv

#import excel files input and metadata

library(readxl)
df <- read_excel(input, col_names = T)
metadata <- read_excel(metadata)
View(df)
View(metadata)

#import csv files MSdial ID as rownames(metabolite)
df <- read.delim(file,
                 sep = ",",
                 h = T,
                 row.names = 1)
metadata <- read.delim(metadata,
                       sep = ",",
                       h = T,
                       row.names = 1)

#setmetabolites name
metabolites <- dataframe[, c(1:3)] #ID MSDial+RT+average_m/z

#Create an additional column called "Group" in the metadata in which
#according to subclass (diet, genotype, etc.) a number as class is associated
#Code is specific for Metadata used for example, if it is useful it can be generalized


subclass <- metadata$subclass
print(subclass)
class <- metadata$class
print(class)

library(dplyr)

group <- function(metadata) {
  metadata <- metadata %>%
    mutate(
      Group = case_when(
        metadata$subclass == "STAND" ~ "1",
        metadata$subclass == "OMEGA" ~ "2",
        metadata$class == "QC" ~ "3",
      )
    )
  metadata$Group <- as.factor(metadata$Group)
}
group(metadata)
View(metadata)
