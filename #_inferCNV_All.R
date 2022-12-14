
##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(300000)

##### Load package #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

##### Load data ######
# ## Load all
# load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-12-07_TrajAna_PCA_PDAC_ROGUE_Monocle3.RData")
# ## Clean up the object
# rm(list=setdiff(ls(), c("scRNA.SeuObj")))

## Small dataset
## Load all
load("D:/Dropbox/##_GitHub/##_CAESAR/MagicDisc/2022-12-13_SeuratSmall_PADC/2022-12-13_SeuratSmall_PADC_SmallData.RData")
## Clean up the object
rm(list=setdiff(ls(), str_subset(objects(), pattern = "Small")))
scRNA.SeuObj <- scRNA_Small.SeuObj

scRNA.SeuObj@meta.data$Cell_type <- gsub(" ", "_", scRNA.SeuObj@meta.data$Cell_type) #Note!!# Will make factor become character!!
scRNA.SeuObj@meta.data$Cell_type <- as.factor(scRNA.SeuObj@meta.data$Cell_type)
scRNA.SeuObj <- scRNA.SeuObj[,!scRNA.SeuObj$Cell_type %in% c(NA)]


## plot UMAP
library(Seurat)
DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA")
DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")

##### Current path and new folder setting* #####
ProjectName = paste0("inferCNV")
Sampletype = "PDAC2"

Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
Save.Path = paste0(getwd(),"/",Version)
## Create new folder
if (!dir.exists(Save.Path)){
  dir.create(Save.Path)
}

RData_Save.Path = "D:/Dropbox/#_Dataset/Cancer/PDAC/"

##### Condition setting #####
Species = "Human"

##### 013 inferCNV #####
## Create new folder
PathinferCNV <- paste0(Save.Path,"/","D01_inferCNV")
if (!dir.exists(PathinferCNV)){
  dir.create(PathinferCNV)
}

source("FUN_inferCNV.R")

# infercnv_obj <- inferCNV(scRNA.SeuObj, AnnoSet = "singleR_classic_PredbyscRNA",
#                          SpeciSet = Species,
#                          Path = PathinferCNV,
#                          RefSet = c("Acinar cell"),
#                          CreateInfercnvObject.lt = list(chr_exclude = c("chrM")))

infercnv_obj <- inferCNV(scRNA.SeuObj, AnnoSet = "Cell_type",
                         SpeciSet = Species,
                         Path = PathinferCNV,
                         RefSet = c("Acinar_cell"),
                         CreateInfercnvObject.lt = list(delim="\t",max_cells_per_group = NULL,min_max_counts_per_cell = c(100, +Inf),chr_exclude = c("chrX", "chrY", "chrM")))


##### save.image #####
save.image(paste0(RData_Save.Path,"/",Version,".RData"))
