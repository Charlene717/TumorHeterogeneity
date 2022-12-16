
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

# scRNA_Sub.SeuObj@meta.data$Cell_type <- gsub("_", " ", scRNA_Sub.SeuObj@meta.data$Cell_type) #Note!!# Will make factor become character!!
scRNA_Sub.SeuObj@meta.data$Cell_type <- gsub(" ", "_", scRNA_Sub.SeuObj@meta.data$Cell_type) #Note!!# Will make factor become character!!
#Try small data# scRNA_Sub.SeuObj <- scRNA_Sub.SeuObj[,scRNA_Sub.SeuObj$Cell_type %in% c("Acinar_cell","Ductal_cell_type_1","Ductal_cell_type_2")]
scRNA_Sub.SeuObj <- scRNA_Sub.SeuObj[,!scRNA_Sub.SeuObj$Cell_type %in% c("T_cell","Fibroblast_cell")]
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
                         RefSet = c("Acinar_cell"),
                         inferCNVRun.lt = list(cluster_by_groups = FALSE, cluster_references= TRUE, plot_steps=FALSE, no_plot=FALSE, resume_mode = FALSE, k_nn = 30),
                         CreateInfercnvObject.lt = list(delim="\t",max_cells_per_group = NULL,min_max_counts_per_cell = c(100, +Inf),chr_exclude = c("chrX", "chrY", "chrM")))

##### save.image #####
save.image(paste0(RData_Save.Path,"/",Version,".RData"))
