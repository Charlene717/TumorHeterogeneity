
##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(300000)

##### Load package #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

##### Load data ######
## Load all
load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-12-07_TrajAna_PCA_PDAC_ROGUE_Monocle3.RData")
## Clean up the object
rm(list=setdiff(ls(), c("scRNA.SeuObj")))

## plot UMAP
library(Seurat)
DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA")

##### Current path and new folder setting* #####
ProjectName = paste0("inferCNV")
Sampletype = "PDAC"

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

infercnv_obj <- inferCNV(scRNA.SeuObj, AnnoSet = "singleR_classic_PredbyscRNA",
                         SpeciSet = Species,
                         Path = PathinferCNV,
                         RefSet = c("Acinar cell"),
                         CreateInfercnvObject.lt = list(chr_exclude = c("chrM")))
##### save.image #####
save.image(paste0(RData_Save.Path,"/",Version,".RData"))
