
##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(300000)

##### Load package #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

##### Load data ######
## Load all
# load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-12-05_CTAnno_singleR_RefPRJCA001063_PDAC.RData")
load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-12-07_TrajAna_PCA_PDAC_ROGUE_Monocle3.RData")

## Clean up the object
rm(list=setdiff(ls(), c("scRNA.SeuObj","scRNA_Sub.SeuObj")))

## plot UMAP
library(Seurat)

scRNA_Sub.SeuObj@meta.data$Cell_type <- gsub("_", " ", scRNA_Sub.SeuObj@meta.data$Cell_type) #Note!!# Will make factor become character!!
scRNA_Sub.SeuObj <- scRNA_Sub.SeuObj[,!scRNA_Sub.SeuObj$Cell_type %in% c("T cell","Fibroblast cell")]
scRNA_Sub.SeuObj@meta.data$Cell_type <- as.factor(scRNA_Sub.SeuObj@meta.data$Cell_type)

DimPlot(scRNA_Sub.SeuObj, reduction = "umap",group.by = "Cell_type")


##### Current path and new folder setting* #####
ProjectName = paste0("inferCNV_Sub")
Sampletype = "PDAC"

Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
Save.Path = paste0(getwd(),"/",Version)
## Create new folder
if (!dir.exists(Save.Path)){dir.create(Save.Path)}

RData_Save.Path = "D:/Dropbox/#_Dataset/Cancer/PDAC/"

##### Condition setting #####
Species = "Human"

##### 013 inferCNV #####
## Create new folder
PathinferCNV <- paste0(Save.Path,"/","D01_inferCNV")
if (!dir.exists(PathinferCNV)){dir.create(PathinferCNV)}

source("FUN_inferCNV.R")

infercnv_obj <- inferCNV(scRNA_Sub.SeuObj, AnnoSet = "Cell_type",
                         SpeciSet = Species,
                         Path = PathinferCNV,
                         RefSet = c("Acinar cell"),
                         CreateInfercnvObject.lt = list(chr_exclude = c("chrM")))
##### save.image #####
save.image(paste0(RData_Save.Path,"/",Version,".RData"))
