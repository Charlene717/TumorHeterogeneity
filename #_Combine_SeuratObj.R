##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(300000)

##### Load basic packages #####
  if(!require("Seurat")) install.packages("Seurat")
  if(!require("tidyverse")) install.packages("tidyverse")
  library(tidyverse)
  library(Seurat)

#### Load data #####
  load("D:/Dropbox/##_GitHub/#_scRNADataset/SeuratObject_CDS_PRJCA001063.RData")
  scRNA.SeuObj_1 <- scRNA.SeuObj

  load("D:/Dropbox/##_GitHub/#_scRNADataset/SeuratObject_GSE131886_PDC_SC.RData")
  scRNA.SeuObj_2 <- scRNA.SeuObj
  DefaultAssay(scRNA.SeuObj_2) <- "RNA"

  load("D:/Dropbox/##_GitHub/#_scRNADataset/SeuratObject_GSE154778_PDAC_SC.RData")
  scRNA.SeuObj_3 <- scRNA.SeuObj
  DefaultAssay(scRNA.SeuObj_3) <- "RNA"

  scRNA_SeuObj.list <- list(PRJCA001063 = scRNA.SeuObj_1,
                            GSE131886 = scRNA.SeuObj_2,
                            GSE154778 = scRNA.SeuObj_3)

  rm(list=setdiff(ls(), "scRNA_SeuObj.list"))

##### Current path and new folder setting* #####
  ProjectName = "Com"
  Sampletype = "PDAC"
  # ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }


##### Load Packages #####
  #### Basic installation ####
  Package.set <- c("tidyverse","Seurat","ggplot2","ggpmisc",
                   "stringr","magrittr","dplyr")
  ## Check whether the installation of those packages is required from basic
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### BiocManager installation ####
  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("fgsea","AnnotationHub","ensembldb",
                   "SeuratDisk","monocle",
                   "SingleR","scRNAseq","celldex","scran")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  options(stringsAsFactors = FALSE)

  #### GitHub installation ####
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  library(monocle)
  devtools::install_github("cole-trapnell-lab/garnett")
  devtools::install_github('cole-trapnell-lab/monocle3')
  #devtools::install_github("cole-trapnell-lab/monocle3", ref="develop")
  devtools::install_github("LTLA/SingleR")

  library(monocle3)
  library(garnett)
  library(SingleR)

  remotes::install_github('satijalab/seurat-wrappers')
  library(SeuratWrappers)




##### Function setting #####
  ## Call function
  source("FUN_Cal_Mit.R")
  source("FUN_CombineSeuObj.R")
  source("FUN_Beautify_ggplot.R")
  source("FUN_Anno_SingleR.R")



##### CombineSeuObj #####
  scRNA.SeuObj <- CombineSeuObj(scRNA_SeuObj.list)
  rm(scRNA_SeuObj.list)

  #### Re-dimension reduction ####
  # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj, selection.method = "vst", nfeatures = 2000)
  scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)
  # Run the standard workflow for visualization and clustering
  scRNA.SeuObj <- ScaleData(scRNA.SeuObj)

  ## PCA: Finding the right PCA conditions
  scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 200)
  VizDimLoadings(scRNA.SeuObj, dims = 1:2, reduction = "pca")
  DimPlot(scRNA.SeuObj, reduction = "pca")
  DimHeatmap(scRNA.SeuObj, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(scRNA.SeuObj, dims = 1:15, cells = 500, balanced = TRUE)
  # Determine the ‘dimensionality’ of the dataset
    # NOTE: This process can take a long time for big datasets, comment out for expediency. More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
    pbmc <- scRNA.SeuObj
    pbmc <- JackStraw(pbmc, num.replicate = 100)
    pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
    JackStrawPlot(pbmc, dims = 1:20)
    ElbowPlot(pbmc, ndims = 200)

  ## UMAP
  scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:150,n.neighbors = 20,min.dist = 0.3)
  scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:150)
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)

  #### Save RData #####
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_Ori.RData"))


  ##### Plot #####
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))

  # DimPlot(scRNA.SeuObj, reduction = "umap")
  DimPlot(scRNA.SeuObj, reduction = "umap",label = T)
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "DataSetID")

##### Cell Type Annotation #####
  scRNA.SeuObj@meta.data[["DataSetID"]] %>% unique()
  scRNA.SeuObj_Ori <- scRNA.SeuObj
  scRNA.SeuObj_Ref <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[["DataSetID"]] %in% "PRJCA001063"]
  DimPlot(scRNA.SeuObj_Ref, reduction = "umap",group.by = "Cell_type")

  RefName = "Cell_type"
  RefName2 = "DataSetID"
  Remark1 <- "PredbyCTDB"         # "PredbyscRNA_CT_GSE131886"
  source("CellTypeAnno_SingleR.R", encoding="UTF-8")

  scRNA.SeuObj@meta.data[["Cell_type"]] <- scRNA.SeuObj@meta.data[[paste0("singleR_",Remark)]]
  scRNA.SeuObj@meta.data[["celltype"]] <- scRNA.SeuObj@meta.data[[paste0("singleR_",Remark)]]
  scRNA.SeuObj@meta.data[[paste0("singleR_",Remark)]] <- NULL

  # DefaultAssay(scRNA.SeuObj) <- "RNA"
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "Cell_type")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "ReCluster2")
  DimPlot(scRNA.SeuObj, reduction = "umap",group.by = "DataSetID")
  DimPlot(scRNA.SeuObj, reduction = "umap")

  ## Cancer stem cell marker gene
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A", "CD24", "CD44", "EPCAM"))
  # https://molecular-cancer.biomedcentral.com/articles/10.1186/s12943-022-01596-8
  FeaturePlot(scRNA.SeuObj, features = c("CD133", "CXCR4", "CD44", "ABCG2","ALDH1A1","MET"))
  # https://pubmed.ncbi.nlm.nih.gov/29948612/#:~:text=Cancer%20stem%20cells%20(CSCs)%2C,resistance%20to%20chemotherapy%20and%20radiation.


  rm(meta_ori.df, scRNA.SeuObj_Ori, SingleRResult.lt)
  #### Save RData #####
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_Anno.RData"))

  # #### Re-dimension reduction ####
  # DefaultAssay(scRNA.SeuObj) <- "integrated"
  # # # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj, selection.method = "vst", nfeatures = 2000)
  # scRNA.SeuObj <- FindVariableFeatures(scRNA.SeuObj)
  # # Run the standard workflow for visualization and clustering
  # scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  # scRNA.SeuObj <- RunPCA(scRNA.SeuObj, npcs = 10, verbose = FALSE)
  # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:10,n.neighbors = 20,min.dist = 0.3)
  # scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:10)
  # scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)

  #### Save RData #####
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_Anno_ReDR.RData"))



##### Building trajectories with Monocle3 #####
  ## https://satijalab.org/signac/articles/monocle.html
  # remotes::install_github('satijalab/seurat-wrappers')
  # library(SeuratWrappers)
  #
  library(monocle3)
  scRNA.cds <- as.cell_data_set(scRNA.SeuObj)
  scRNA.cds <- cluster_cells(cds = scRNA.cds, reduction_method = "UMAP")
  scRNA.cds <- learn_graph(scRNA.cds, use_partition = TRUE)

  scRNA.cds <- order_cells(scRNA.cds)

  plot_cells(
    cds = scRNA.cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )
  plot_cells(scRNA.cds, color_cells_by = "partition")
  plot_cells(scRNA.cds, color_cells_by = "Cell_type")

  ## plot_genes_in_pseudotime
  rowData(scRNA.cds)$gene_short_name <- scRNA.cds@assays@data@listData[["counts"]]@Dimnames[[1]]
  Int_genes <- c("TOP2A", "NSUN2", "TP53","PTK2")
  Int_lineage_cds <- scRNA.cds[rowData(scRNA.cds)$gene_short_name %in% Int_genes,]
  plot_genes_in_pseudotime(Int_lineage_cds,
                           color_cells_by="pseudotime",
                           min_expr=0.5)

  cds_sub <- choose_cells(scRNA.cds,clear_cds = F)


  cds_sub <- choose_graph_segments(scRNA.cds,clear_cds = F)
  Int_lineage_cds_sub <- cds_sub[rowData(cds_sub)$gene_short_name %in% Int_genes,]

  plot_genes_in_pseudotime(Int_lineage_cds_sub,
                           color_cells_by="pseudotime",
                           min_expr=0.5)

  #### Save RData #####
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_Traj.RData"))



##### Export figures #####
  ## Export TIFF
  for (i in seq(50,400,50)) {
    for (j in seq(0.1,0.7,0.2)) {
      for (k in seq(20,300,40)) {
        try({
          set.seed(1)
          scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:i,n.neighbors = k, min.dist= j)
          scRNA.SeuObj@meta.data[[paste0("UMAP_PCA",i,"_NNe",k,"_MD03",j)]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
          # scRNA.SeuObj@reductions[["umap"]]@cell.embeddings <- scRNA.SeuObj@meta.data[[paste0("UMAP_PCA",i,"_NNe",k,"_MD03",j)]]

          Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["Cell_type"]]
          p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
            ggtitle(paste0("CellType","  PCA:",i,"  NNe:",k,"  MD:",j)) +
            theme(plot.title = element_text(hjust = 0.5,vjust = 0))
          tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_CellType.tiff"),
               width = 28, height = 20, units = "cm", res = 200)
          print(p)
          graphics.off()

          Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["ReCluster2"]]
          p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
            ggtitle(paste0("ReCluster2","  PCA:",i,"  NNe:",k,"  MD:",j)) +
            theme(plot.title = element_text(hjust = 0.5,vjust = 0))
          tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_ReCluster.tiff"),
               width = 35, height = 20, units = "cm", res = 200)
          print(p)
          graphics.off()

          Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["DataSetID"]]
          p <-  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() %>% BeautifyggPlot(LegPos = c(1.02, 0.5)) +
            ggtitle(paste0("ReCluster","  PCA:",i,"  NNe:",k,"  MD:",j)) +
            theme(plot.title = element_text(hjust = 0.5,vjust = 0))
          tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_DataSetID.tiff"),
               width = 35, height = 20, units = "cm", res = 200)
          print(p)
          graphics.off()

          p <-  FeaturePlot(scRNA.SeuObj, features = c("TOP2A")) %>% BeautifyggPlot(LegPos = c(1.02, 0.15)) +
            ggtitle(paste0("TOP2A","  PCA:",i,"  NNe:",k,"  MD:",j)) +
            theme(plot.title = element_text(hjust = 0.5,vjust = 0))
          tiff(file = paste0(Save.Path,"/",ProjectName,"_Trajectory","_PCA",i,"_NNe",k,"_MD",j,"_TOP2A.tiff"),
               width = 28, height = 20, units = "cm", res = 200)
          print(p)
          graphics.off()
        })
      }
    }
  }



  # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:100,n.neighbors = 20, min.dist=0.05)
  # # scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, dims = 1:100,n.neighbors = 1000, min.dist=0.1)
  # # scRNA.SeuObj@meta.data[["UMAP_NNe1000"]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
  # # scRNA.SeuObj@meta.data <- scRNA.SeuObj@meta.data[,!colnames(scRNA.SeuObj@meta.data)=="UMAP_NNe1000"]
  # scRNA.SeuObj@meta.data[["UMAP_NNe20_MD03"]] <- scRNA.SeuObj@reductions[["umap"]]@cell.embeddings
  #
  # DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

  Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["Cell_type"]]
  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

  Idents(scRNA.SeuObj) <- scRNA.SeuObj@meta.data[["ReCluster"]]
  DimPlot(scRNA.SeuObj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  FeaturePlot(scRNA.SeuObj, features = c("TOP2A"))

  #### Save RData #####
  save.image(paste0(Save.Path,"/scRNA.SeuObj_CDS_PRJCA001063_Combine_TryCondition.RData"))


  ##### Session information #####
    sessionInfo()
    ## Ref: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
    writeLines(capture.output(sessionInfo()), paste0(Save.Path,"/sessionInfo.txt"))





