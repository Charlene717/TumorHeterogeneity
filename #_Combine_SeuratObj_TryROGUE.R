## Ref: https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html
## Ref: https://github.com/PaulingLiu/ROGUE
## Ref: https://www.nature.com/articles/s41467-020-16904-3

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)

##### Load Packages #####
  #### Basic and BiocManager installation ####
  source("FUN_Package_InstLoad.R")
  FUN_Basic.set <- c("tidyverse","Seurat","ggpmisc","stringr","magrittr")
  FUN_BiocManager.set <- c("ensembldb","SeuratDisk")

  FUN_Package_InstLoad(Basic.set = FUN_Basic.set, BiocManager.set = FUN_BiocManager.set)
  rm(FUN_Basic.set, FUN_BiocManager.set)

  #### GitHub installation ####
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("PaulingLiu/ROGUE")
  library(ROGUE)


##### Function setting #####
  ## Call function
  source("FUN_Cal_Mit.R")
  source("FUN_CombineSeuObj.R")
  source("FUN_Beautify_ggplot.R")
  source("FUN_Anno_SingleR.R")


#### Load data #####
  # load("D:/Dropbox/##_GitHub/##_Charlene/TumorHeterogeneity/2022-08-01_Com_PDAC/scRNA.SeuObj_CDS_PRJCA001063_Combine_Anno_ReDR.RData")
  # load("D:/Dropbox/##_GitHub/##_Charlene/TumorHeterogeneity/2022-11-08_Com_PDAC/SeuratObject_Com.RData")
  load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-11-15_CTAnno_singleR_RefPRJCA001063_PDAC.RData")

##### Current path and new folder setting* #####
  ProjectName = "Com_ROGUE"
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Extract data #####
  ## Gene GeneExp.dfession
  ## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix

  ## Meta.data
  # Cell annotation
  Meta.data <- scRNA.SeuObj@meta.data
  Meta.data <- data.frame(ID=row.names(Meta.data), Meta.data)


##### Run ROGUE #####
  ## Filtering out low-abundance genes and low-quality cells
  GeneExp.df <- matr.filter(GeneExp.df, min.cells = 10, min.genes = 10)

  ## GeneExp.dfession entropy model
  ent.res <- SE_fun(GeneExp.df)
  head(ent.res)

  ## S-E plot
  SEplot(ent.res)

  ## ROGUE calculation
  rogue.value <- CalculateRogue(ent.res, platform = "UMI")
  rogue.value


  ##### ************************************************************************ #####
  scRNA.SeuObj@meta.data$singleR_classic_PredbyscRNA2 <- scRNA.SeuObj@meta.data$singleR_classic_PredbyscRNA  %>% as.character()
  TTT1 <- scRNA.SeuObj@meta.data
  TTT_Ref <- scRNA.SeuObj_Ref@meta.data
  TTT1[TTT1$CELL %in% TTT_Ref$CELL,]$singleR_classic_PredbyscRNA2 <- TTT_Ref$celltype %>% as.character()
  scRNA.SeuObj@meta.data <- TTT1
  scRNA.SeuObj@meta.data$singleR_classic_PredbyscRNA2 <- as.factor(scRNA.SeuObj@meta.data$singleR_classic_PredbyscRNA2)


  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA2")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "celltype")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "celltype")

  ## Meta.data
  # Cell annotation
  Meta.data <- scRNA.SeuObj@meta.data
  Meta.data <- data.frame(ID=row.names(Meta.data), Meta.data)
  ##### ************************************************************************ ####


  ## Calculate the ROGUE value of each putative cluster for each sample
  rogue.res <- rogue(GeneExp.df, labels = Meta.data$celltype, samples = Meta.data$DataSetID, platform = "UMI", span = 0.6)
  rogue.res
  rogue.boxplot(rogue.res)

  rogue.res2 <- rogue(GeneExp.df, labels = Meta.data$singleR_classic_PredbyscRNA2, samples = Meta.data$DataSetID, platform = "UMI", span = 0.6)
  rogue.res2
  View(rogue.res2)
  rogue.res2_2 <- rogue(GeneExp.df, labels = Meta.data$singleR_classic_PredbyscRNA2, samples = Meta.data$seurat_clusters, platform = "UMI", span = 0.6)
  rogue.res2_2

  rogue.res3_1 <- rogue(GeneExp.df, labels = Meta.data$seurat_clusters, samples = Meta.data$DataSetID, platform = "UMI", span = 0.6)
  rogue.res3_1

  rogue.res3_2 <- rogue(GeneExp.df, labels = Meta.data$seurat_clusters, samples = Meta.data$singleR_classic_PredbyscRNA2, platform = "UMI", span = 0.6)
  rogue.res3_2



  ## Visualize ROGUE values on a boxplot
  rogue.boxplot(rogue.res)
  rogue.boxplot(rogue.res2)
  P2 <- rogue.boxplot(rogue.res2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P2

  rogue.boxplot(rogue.res2_2)
  P2_2 <- rogue.boxplot(rogue.res2_2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P2_2

  rogue.boxplot(rogue.res3_1)
  P3_1 <- rogue.boxplot(rogue.res3_1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P3_1

  rogue.boxplot(rogue.res3_2)
  P3_2 <- rogue.boxplot(rogue.res3_2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P3_2

  ##### Save RData #####
  save.image(paste0("D:/Dropbox/#_Dataset/Cancer/PDAC/",ProjectName,"_ROGUE.RData"))

  # save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,"_ROGUE.RData"))

# ####################################################################################################
#
#
#   #### Re-dimension reduction ####
#   # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj, selection.method = "vst", nfeatures = 2000)
#   scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)
#   # Run the standard workflow for visualization and clustering
#   scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
#   scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 160, verbose = FALSE)
#   scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:160,n.neighbors = 20,min.dist = 0.3)
#   scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:160)
#   scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)
#
#   #### Save RData #####
#   save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_Ori.RData"))
#
#
#   ##### Plot #####
#   FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
#   # DimPlot(scRNA.SeuObj, reduction = "umap")
#   DimPlot(scRNA.SeuObj, reduction = "umap",label = T)
#   DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")
#   DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "DataSetID")
#
#   # #### Re-dimension reduction ####
#   # DefaultAssay(scRNA.SeuObj) <- "integrated"
#   # # # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj, selection.method = "vst", nfeatures = 2000)
#   # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)
#   # # Run the standard workflow for visualization and clustering
#   # scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
#   # scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 10, verbose = FALSE)
#   # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:10,n.neighbors = 20,min.dist = 0.3)
#   # scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:10)
#   # scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)
#
#   rm(scRNA.SeuObj_1, scRNA.SeuObj_2, scRNA.SeuObj_3, scRNA.SeuObj_Ref,
#      scRNA.SeuObj_Ori, CTFeatures.SeuObj,SingleRResult.lt)
#
#   #### Save RData #####
#   save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_Anno_ReDR.RData"))
#
#
#   ##### Session information #####
#     sessionInfo()
#     ## Ref: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
#     writeLines(capture.output(sessionInfo()), paste0(Save.Path,"/sessionInfo.txt"))
#
#
#
#
#
