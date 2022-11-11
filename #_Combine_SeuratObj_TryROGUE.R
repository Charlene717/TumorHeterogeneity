## Ref: https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)

##### Load Packages #####
  source("FUN_Package_InstLoad.R")
  FUN_Basic.set <- c("tidyverse","Seurat","ggplot2","ggpmisc",
                     "stringr","magrittr","dplyr")
  FUN_BiocManager.set <- c("fgsea","AnnotationHub","ensembldb",
                           "SeuratDisk","monocle",
                           "SingleR","scRNAseq","celldex","scran")
  ## Set the desired organism
  # organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db")
  # c(organism,"fgsea")

  FUN_Package_InstLoad(Basic.set = FUN_Basic.set, BiocManager.set = FUN_BiocManager.set)
  rm(FUN_Basic.set, FUN_BiocManager.set)

  #### GitHub installation ####
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  # library(monocle)
  # devtools::install_github("cole-trapnell-lab/garnett")
  # devtools::install_github('cole-trapnell-lab/monocle3')
  #
  # library(monocle3)
  # library(garnett)


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
  load("D:/Dropbox/##_GitHub/##_Charlene/TumorHeterogeneity/2022-11-08_Com_PDAC/SeuratObject_Com.RData")


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

  rogue.res2 <- rogue(GeneExp.df, labels = Meta.data$singleR_classic_PredbyscRNA2, samples = Meta.data$DataSetID, platform = "UMI", span = 0.6)
  rogue.res2
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
  save.image(paste0(Save.Path,"/SeuratObject_",ProjectName,"_ROGUE.RData"))

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
# ##### Cell Type Annotation #####
#   scRNA.SeuObj@meta.data[["DataSetID"]] %>% unique()
#   scRNA.SeuObj_Ori <- scRNA.SeuObj
#   scRNA.SeuObj_Ref <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["DataSetID"]] %in% "PRJCA001063"]
#   DimPlot(scRNA.SeuObj_Ref, reduction = "umap",group.by = "Cell_type")
#
#   RefName = "Cell_type"
#   RefName2 = "DataSetID"
#   Remark1 <- "PredbyCTDB"         # "PredbyscRNA_CT_GSE131886"
#   source("CellTypeAnno_SingleR.R", encoding="UTF-8")
#
#   scRNA.SeuObj@meta.data[["Cell_type"]] <- scRNA.SeuObj@meta.data[[paste0("singleR_",Remark)]]
#   scRNA.SeuObj@meta.data[["celltype"]] <- scRNA.SeuObj@meta.data[[paste0("singleR_",Remark)]]
#   scRNA.SeuObj@meta.data[[paste0("singleR_",Remark)]] <- NULL
#
#   # DefaultAssay(scRNA.SeuObj) <- "RNA"
#   FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
#   DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")
#   DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "ReCluster2")
#   DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "DataSetID")
#   DimPlot(scRNA.SeuObj, reduction = "umap")
#
#   rm(meta_ori.df, meta_ann.df, meta.df, scRNA.SeuObj_Ori,
#      SingleRResult.lt, CTFeatures.SeuObj)
#   #### Save RData #####
#   save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_Anno.RData"))
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
# ##### Export figures #####
#   ## Export TIFF
#   for (i in seq(50,400,50)) {
#     for (j in seq(0.1,0.7,0.2)) {
#       for (k in seq(20,300,40)) {
#         try({
#           set.seed(1)
#           scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:i,n.neighbors = k, min.dist= j)
#           scRNA.SeuObj@meta.data[[paste0("UMAP_PCA",i,"_NNe",k,"_MD03",j)]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
#           # scRNA.SeuObj@reductions[["umap"]]@cell.embeddings <- scRNA.SeuObj@meta.data[[paste0("UMAP_PCA",i,"_NNe",k,"_MD03",j)]]
#
#           Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["Cell_type"]]
#           p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
#             ggtitle(paste0("CellType","  PCA:",i,"  NNe:",k,"  MD:",j)) +
#             theme(plot.title = element_text(hjust = 0.5,vjust = 0))
#           tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_CellType.tiff"),
#                width = 28, height = 20, units = "cm", res = 200)
#           print(p)
#           graphics.off()
#
#           Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["ReCluster2"]]
#           p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
#             ggtitle(paste0("ReCluster2","  PCA:",i,"  NNe:",k,"  MD:",j)) +
#             theme(plot.title = element_text(hjust = 0.5,vjust = 0))
#           tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_ReCluster.tiff"),
#                width = 35, height = 20, units = "cm", res = 200)
#           print(p)
#           graphics.off()
#
#           Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["DataSetID"]]
#           p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
#             ggtitle(paste0("ReCluster","  PCA:",i,"  NNe:",k,"  MD:",j)) +
#             theme(plot.title = element_text(hjust = 0.5,vjust = 0))
#           tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_DataSetID.tiff"),
#                width = 35, height = 20, units = "cm", res = 200)
#           print(p)
#           graphics.off()
#
#           p <-  FeaturePlot(scRNA.SeuObj, features = c("TOP2A")) %>% BeautifyggPlot(LegPos = c(1.02, 0.15)) +
#             ggtitle(paste0("TOP2A","  PCA:",i,"  NNe:",k,"  MD:",j)) +
#             theme(plot.title = element_text(hjust = 0.5,vjust = 0))
#           tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_TOP2A.tiff"),
#                width = 28, height = 20, units = "cm", res = 200)
#           print(p)
#           graphics.off()
#         })
#       }
#     }
#   }
#
#
#
#   # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:100,n.neighbors = 20, min.dist=0.05)
#   # # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:100,n.neighbors = 1000, min.dist=0.1)
#   # # scRNA.SeuObj@meta.data[["UMAP_NNe1000"]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
#   # # scRNA.SeuObj@meta.data <- scRNA.SeuObj@meta.data[,!colnames(scRNA.SeuObj@meta.data)=="UMAP_NNe1000"]
#   # scRNA.SeuObj@meta.data[["UMAP_NNe20_MD03"]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
#   #
#   # DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#
#   Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["Cell_type"]]
#   DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#
#   Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["ReCluster"]]
#   DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#   FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))
#
#   #### Save RData #####
#   save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_TryCondition.RData"))
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
