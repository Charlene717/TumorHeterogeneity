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
  load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-12-05_CTAnno_singleR_RefPRJCA001063_PDAC.RData")

##### Current path and new folder setting* #####
  ProjectName = "Com_ROGUE"
  Sampletype = "PDAC"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){dir.create(Save.Path)}

##### Extract data #####
  ## Gene GeneExp.dfession
  ## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix

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
  Meta.df <- scRNA.SeuObj@meta.data
  Meta_Ref.df <- scRNA.SeuObj_Ref@meta.data
  Meta.df[Meta.df$CELL %in% Meta_Ref.df$CELL,]$singleR_classic_PredbyscRNA2 <- Meta_Ref.df$celltype %>% as.character()
  scRNA.SeuObj@meta.data <- Meta.df
  scRNA.SeuObj@meta.data$singleR_classic_PredbyscRNA2 <- as.factor(scRNA.SeuObj@meta.data$singleR_classic_PredbyscRNA2)
  rm(Meta.df,Meta_Ref.df)

  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA2")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "seurat_clusters")

  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA")

  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "celltype")
  DimPlot(scRNA.SeuObj_Ref, reduction = "umap",group.by = "celltype")

  ## Extract Meta.data
  Meta.df <- scRNA.SeuObj@meta.data
  Meta.df <- data.frame(ID=row.names(Meta.df), Meta.df)

  ##### ************************************************************************ ####
  # ### Try Ori celltype
  # ## Calculate the ROGUE value of each putative cluster for each sample
  # rogue.res <- rogue(GeneExp.df, labels = Meta.df$celltype, samples = Meta.df$DataSetID, platform = "UMI", span = 0.6)
  # rogue.res
  # rogue.boxplot(rogue.res)
  #
  # ## Visualize ROGUE values on a boxplot
  # rogue.boxplot(rogue.res)
  # P1 <- rogue.boxplot(rogue.res) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # P1


  # #### DataSetID ####
  # ## Calculate the ROGUE value of each putative cluster for each sample
  # rogue_DataSetID.res <- rogue(GeneExp.df, labels = Meta.df$singleR_classic_PredbyscRNA2, samples = Meta.df$DataSetID, platform = "UMI", span = 0.6)
  # rogue_DataSetID.res
  # View(rogue_DataSetID.res)
  #
  # ## Visualize ROGUE values on a boxplot
  # rogue.boxplot(rogue_DataSetID.res)
  # P.DataSetID <- rogue.boxplot(rogue_DataSetID.res) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # P.DataSetID


  #### ROUGE of cell type by cluster ####
  rogue_CT_CLU.res <- rogue(GeneExp.df, labels = Meta.df$singleR_classic_PredbyscRNA2, samples = Meta.df$seurat_clusters, platform = "UMI", span = 0.6)
  rogue_CT_CLU.res
  View(rogue_CT_CLU.res)

  ## Visualize ROGUE values on a boxplot
  rogue.boxplot(rogue_CT_CLU.res)
  P.CT_CLU.res <- rogue.boxplot(rogue_CT_CLU.res) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P.CT_CLU.res

  #### ROUGE of cluster by cell type ####
  rogue_CLU_CT.res <- rogue(GeneExp.df, labels = Meta.df$seurat_clusters, samples = Meta.df$singleR_classic_PredbyscRNA2, platform = "UMI", span = 0.6)
  rogue_CLU_CT.res
  View(rogue_CLU_CT.res)

  ## Visualize ROGUE values on a boxplot
  rogue.boxplot(rogue_CLU_CT.res)
  P.CLU_CT.res <- rogue.boxplot(rogue_CLU_CT.res) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P.CLU_CT.res


  #### ROUGE of cluster by Patients(SampleID) #### ## Take 30 min -> Too long
  Meta.df$Patient <- Meta.df$SampleID
  Meta_Ref.df <- scRNA.SeuObj_Ref@meta.data
  Meta.df[Meta.df$CELL %in% Meta_Ref.df$CELL,]$Patient <- Meta_Ref.df$Patient
  rm(Meta_Ref.df)

  rogue_CLU_Samp.res <- rogue(GeneExp.df, labels = Meta.df$seurat_clusters, samples = Meta.df$Patient, platform = "UMI", span = 0.6)
  rogue_CLU_Samp.res
  View(rogue_CLU_Samp.res)

  ## Visualize ROGUE values on a boxplot
  rogue.boxplot(rogue_CLU_Samp.res)
  P.CLU_Samp.res <- rogue.boxplot(rogue_CLU_Samp.res) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P.CLU_Samp.res


  ## ************************************************************************ ##
  #### ROUGE of cluster by DataSetID #### ## Take 12 min
  rogue_CLU_DataSetID.res <- rogue(GeneExp.df, labels = Meta.df$seurat_clusters, samples = Meta.df$DataSetID, platform = "UMI", span = 0.6)
  rogue_CLU_DataSetID.res
  View(rogue_CLU_DataSetID.res)

  ## Visualize ROGUE values on a boxplot
  rogue.boxplot(rogue_CLU_DataSetID.res)
  P.CLU_DataSetID.res <- rogue.boxplot(rogue_CLU_DataSetID.res) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P.CLU_DataSetID.res

  # ## avg rouge
  # av.rogue <- c()
  # for (i in 1:ncol(rogue_CLU_DataSetID.res)) {
  #   tmp.r <- rogue_CLU_DataSetID.res[,i]
  #   tmp.r <- tmp.r[!is.na(tmp.r)]
  #   av.rogue[i] <- mean(tmp.r)
  # }

  av.rogue <- mean(rogue_CLU_DataSetID.res[!is.na(rogue_CLU_DataSetID.res)])


  # #### Try Clustering Condition ####
  # scRNA.SeuObj_Test <- FindClusters(scRNA.SeuObj, resolution = 0.5)
  # DimPlot(scRNA.SeuObj_Test, reduction = "umap",group.by = "seurat_clusters")
  #
  # scRNA.SeuObj_Test1 <- FindClusters(scRNA.SeuObj, resolution = 0.5)
  # DimPlot(scRNA.SeuObj_Test1, reduction = "umap",group.by = "seurat_clusters")
  #
  # scRNA.SeuObj_Test2 <- FindClusters(scRNA.SeuObj, resolution = 0.8)
  # DimPlot(scRNA.SeuObj_Test2, reduction = "umap",group.by = "seurat_clusters")


  ### Rogue_TryCond_DataSet.df
  Rogue_TryCond_DataSet.df <- data.frame(matrix(data = NA ,nrow = 20,ncol = 3))
  colnames(Rogue_TryCond_DataSet.df) <- c("CondSet_Res","ClusterNum","av.rogue")

  for (i in seq(1:20)) {
    CondSet_Res_Set <- i*0.01
    scRNA.SeuObj_Temp <- FindClusters(scRNA.SeuObj, resolution = CondSet_Res_Set )
    Meta.df <- scRNA.SeuObj_Temp@meta.data
    Meta.df$Patient <- Meta.df$SampleID
    Meta_Ref.df <- scRNA.SeuObj_Ref@meta.data
    Meta.df[Meta.df$CELL %in% Meta_Ref.df$CELL,]$Patient <- Meta_Ref.df$Patient
    rm(Meta_Ref.df)

    Rogue_TryCond_DataSet.df$CondSet_Res[i] <- CondSet_Res_Set
    # Rogue_TryCond_DataSet.df$ClusterNum[i] <- i*0.1
    Rogue_TryCond_DataSet.df$ClusterNum[i] <- scRNA.SeuObj_Temp$seurat_clusters %>% unique() %>% length()

    ## Rogue
    rogue_Temp.res <- rogue(GeneExp.df, labels = Meta.df$seurat_clusters, samples = Meta.df$Patient, platform = "UMI", span = 0.6)
    av.rogue <- mean(rogue_Temp.res[!is.na(rogue_Temp.res)])
    Rogue_TryCond_DataSet.df$av.rogue[i] <- av.rogue

    rm(scRNA.SeuObj_Temp, rogue_Temp.res, av.rogue, CondSet_Res_Set,scRNA.SeuObj_Ref)
  }
  rm(i)


  ### Line plot
  ## Ref: http://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization
  # Rogue_TryCond_DataSet.df$av.rogue <- c(0,0.01,0.02,0.03,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2)
  p <-ggplot(Rogue_TryCond_DataSet.df, aes(x=ClusterNum, y=av.rogue)) +
      geom_line(color="#9346b3", size=0.8)+
      geom_point(color="#9346b3", size=2)+
      theme_classic()

    # geom_line(aes(color=CondSet))+
    # geom_point(aes(color=CondSet))
  p


  Rogue_TryCond_DataSet.df$av.slope.rogue <- NA
  for (i in seq(1:20)) {
    if(i==1){
      Rogue_TryCond_DataSet.df$av.slope.rogue[i] <- 0
    }else{
      Rogue_TryCond_DataSet.df$av.slope.rogue[i] <- Rogue_TryCond_DataSet.df$av.rogue[i]-Rogue_TryCond_DataSet.df$av.rogue[i-1]
    }

  }

  FinCondSet <- Rogue_TryCond_DataSet.df[which(Rogue_TryCond_DataSet.df$av.slope.rogue == max(Rogue_TryCond_DataSet.df$av.slope.rogue)),]$CondSet
  scRNA_Fin.SeuObj <- FindClusters(scRNA.SeuObj, resolution = FinCondSet*0.1)
  DimPlot(scRNA_Fin.SeuObj, reduction = "umap",group.by = "seurat_clusters")


  ##### Redefine the cluster #####
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA2")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "seurat_clusters")
  scRNA.SeuObj$celltype_Sub <- paste0(scRNA.SeuObj$singleR_classic_PredbyscRNA2,"_",scRNA.SeuObj$seurat_clusters)
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "celltype_Sub")

  ## Sub celltype
  scRNA_Duc2.SeuObj <- scRNA.SeuObj[,grepl("Duc",scRNA.SeuObj$singleR_classic_PredbyscRNA2)]
  DimPlot(scRNA_Duc2.SeuObj, reduction = "umap",group.by = "singleR_classic_PredbyscRNA2")
  DimPlot(scRNA_Duc2.SeuObj, reduction = "umap",group.by = "seurat_clusters")
  scRNA_Duc2.SeuObj <- FindClusters(scRNA_Duc2.SeuObj, resolution = CondSet_Res_Set )

  ##### Save RData #####
  save.image(paste0("D:/Dropbox/#_Dataset/Cancer/PDAC/",Version,"_ROGUE_byPatient.RData"))



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
