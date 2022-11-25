## SingleRBook Ref: http://bioconductor.org/books/release/SingleRBook/
## Example Ref: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(300000)


##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

#### Basic and BiocManager installation ####
source("FUN_Package_InstLoad.R")
FUN_Basic.set <- c("tidyverse","Seurat","ggpubr")
FUN_BiocManager.set <- c("SingleR","scRNAseq","celldex","scran","scater","scuttle")

FUN_Package_InstLoad(Basic.set = FUN_Basic.set, BiocManager.set = FUN_BiocManager.set)
rm(FUN_Basic.set, FUN_BiocManager.set)


##### Function setting #####
## Call function
source("FUN_Anno_SingleR.R", encoding="UTF-8")

#### Load data #####
## Load all
load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-11-15_scRNA_SeuObj_PDAC_SC_Combine_ReBEC.RData")

## Set Ref
scRNA.SeuObj_Ref <- scRNA.SeuObj[,scRNA.SeuObj$DataSetID %in% "PRJCA001063"]
CTFeatures.SeuObj <- scRNA.SeuObj_Ref

## Set Tar
# scRNA.SeuObj_Ref


##### Parameter setting* #####
# Remark1 <- "PredbyscRNA" # c("PredbyCTDB","PredbyscRNA")
RefType <- "BuiltIn_scRNA" # c("BuiltIn_celldex","BuiltIn_scRNA")
celldexDatabase <- "HumanPrimaryCellAtlasData"
# c("BlueprintEncodeData","DatabaseImmuneCellExpressionData","HumanPrimaryCellAtlasData","ImmGenData",
#   "MonacoImmuneData","MouseRNAseqData","NovershternHematopoieticData")
de.method <- "classic"

## Parameter of classifySingleR
quantile = 0.8
tune.thresh = 0.05
sd.thresh = 1

Remark <- paste0(Remark1,"_",de.method,"_",
                 "qua",quantile,"_tun",tune.thresh,"_sd",sd.thresh)

SmallTest = F



#####  #####
if(SmallTest == TRUE){
  ## SeuObj_Ref for small test
  # CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref$CELL %in% sample(scRNA.SeuObj_Ref$CELL,1000)] ## For small test
  CTFeatures.SeuObj <- scRNA.SeuObj_Ref[,scRNA.SeuObj_Ref@meta.data[[1]] %in% sample(scRNA.SeuObj_Ref@meta.data[[1]],1000)] ## For small test
  ## SeuObj_Tar for small test
  # scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj$CELL %in% sample(scRNA.SeuObj$CELL,1000)] ## For small test
  scRNA.SeuObj <- scRNA.SeuObj[,scRNA.SeuObj@meta.data[[1]] %in% sample(scRNA.SeuObj@meta.data[[1]],1000)] ## For small test
}else{
  ## SeuObj_Ref for full data
  CTFeatures.SeuObj <- scRNA.SeuObj_Ref
}




##### Current path and new folder setting* #####
  ProjectName_Ano = paste0("CTAnno_singleR_PRJCA001063S")
  Sampletype = "PDAC"

  # Version = paste0(Sys.Date(),"_",ProjectName_Ano,"_",Sampletype)
  # Save.Path = paste0(getwd(),"/",Version)
  # ## Create new folder
  # if (!dir.exists(Save.Path)){
  #   dir.create(Save.Path)
  # }



##### Run singleR #####
  #### Presetting ####
  SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = RefType, celldexDatabase = celldexDatabase,
                                   CTFeatures.SeuObj = CTFeatures.SeuObj,
                                   quantile = quantile, tune.thresh = tune.thresh, sd.thresh = sd.thresh,
                                   de.method = de.method,
                                   Remark = Remark,
                                   RefName = RefName,
                                   RefName2 = RefName2,
                                   Save.Path = paste0(Save.Path,"/",Remark), ProjectName = "CT")
  scRNA.SeuObj <- SingleRResult.lt[["scRNA.SeuObj"]]

#   #### Try Parameter ####
#   CC_Anno.df <- as.data.frame(matrix(nrow=0, ncol=7))
#   colnames(CC_Anno.df) <- c("TestID", "Tool", "Type","Set", "quantile", "tune_Thr","SD_Thr")
#
#   #### PredbyCTDB ####
#   for (i in seq(0.6,1,0.2)) {
#     for (j in seq(0.03,0.07,0.01)) {
#       for (k in c(1,2)) {
#         Remark1 <- "PredbyCTDB"
#         de.method <- "classic"
#         RefType <- "BuiltIn_celldex"
#         Remark <- paste0(Remark1,"_",de.method,"_",
#                          "qua",i,"_tun",j,"_sd",k)
#
#         SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = RefType, celldexDatabase = "HumanPrimaryCellAtlasData",
#                                          quantile = i, tune.thresh = j, sd.thresh = k,
#                                          Remark = Remark,Save.Path = paste0(Save.Path,"/",Remark), ProjectName = "CT")
#         scRNA.SeuObj <- SingleRResult.lt[["scRNA.SeuObj"]]
#
#         CC_Anno_Temp.df <- data.frame(TestID = "Predict", Tool = "singleR", Type = "PDAC",
#                                       Set = Remark1, quantile = i, tune_Thr = j, SD_Thr = k)
#         CC_Anno.df <- rbind(CC_Anno.df, CC_Anno_Temp.df)
#       }
#     }
#   }
#   rm(i,j,k,CC_Anno_Temp.df, Remark, Remark1, de.method, RefType)
#
#     #### Create check dataframe ####
#     CC.df <- scRNA.SeuObj@meta.data[,(ncol(scRNA.SeuObj@meta.data)-nrow(CC_Anno.df)+1):ncol(scRNA.SeuObj@meta.data)]
#     CC_Anno.df$TestID <- colnames(CC.df)
#
#     # TTT <- gsub(CC.df, pattern = " ", replacement = "_")
#     # TTT <- sub(" ", "_", CC.df)
#     # TTT <- CC.df
#     CC.df <- lapply(CC.df, gsub, pattern = "_", replacement = " ", fixed = TRUE) %>%
#              lapply(., gsub, pattern = "cells", replacement = "cell", fixed = TRUE) %>%
#              lapply(., gsub, pattern = "Macrophage", replacement = "Macrophage cell", fixed = TRUE) %>%
#              lapply(., gsub, pattern = "Fibroblasts", replacement = "Fibroblast cell", fixed = TRUE) %>%
#              lapply(., gsub, pattern = "Epithelial cell", replacement = "Ductal cell type 1", fixed = TRUE) %>%
#              as.data.frame()
#
#     CC_CT.df <- data.frame(Cell_type = scRNA.SeuObj@meta.data[,"Cell_type"])
#     CC.df <- cbind(CC_CT.df, CC.df)
#     rm(CC_CT.df)
#
#     CTReplace <- function(CC.df,ColN=1, ReferCT = "Actual",Replacement="Other") {
#       CC.df[!CC.df[,ColN] %in% c(CC.df[,ReferCT] %>% unique()),][,ColN] <- Replacement
#       return(CC.df)
#     }
#
#     for (i in 2:ncol(CC.df)) {
#       CC.df <- CTReplace(CC.df,ColN=i, ReferCT = "Cell_type",Replacement="Other")
#     }
#
#
#   #### PredbyscRNA ####
#   CC_Anno2.df <- as.data.frame(matrix(nrow=0, ncol=7))
#   colnames(CC_Anno2.df) <- c("TestID", "Tool", "Type","Set", "quantile", "tune_Thr","SD_Thr")
#
#   for (i in seq(0.6,1,0.2)) {
#     for (j in seq(0.03,0.07,0.01)) {
#       for (k in c(1,2)) {
#         Remark1 <- "PredbyscRNA"
#         de.method <- "classic"
#         RefType <- "BuiltIn_scRNA"
#         Remark <- paste0(Remark1,"_",de.method,"_",
#                          "qua",i,"_tun",j,"_sd",k)
#
#         SingleRResult.lt <- Anno_SingleR(scRNA.SeuObj, RefType = RefType, celldexDatabase = "HumanPrimaryCellAtlasData",
#                                          quantile = i, tune.thresh = j, sd.thresh = k,CTFeatures.SeuObj = CTFeatures.SeuObj,
#                                          Remark = Remark, Save.Path = paste0(Save.Path,"/",Remark), ProjectName = "CT")
#         scRNA.SeuObj <- SingleRResult.lt[["scRNA.SeuObj"]]
#
#         CC_Anno_Temp.df <- data.frame(TestID = "Predict", Tool = "singleR", Type = "PDAC",
#                                       Set = Remark1, quantile = i, tune_Thr = j, SD_Thr = k)
#         CC_Anno2.df <- rbind(CC_Anno2.df, CC_Anno_Temp.df)
#       }
#     }
#   }
#   rm(i,j,k,CC_Anno_Temp.df, Remark, Remark1, de.method, RefType)
#
#     #### Create check dataframe ####
#     CC2.df <- scRNA.SeuObj@meta.data[,(ncol(scRNA.SeuObj@meta.data)-nrow(CC_Anno2.df)+1):ncol(scRNA.SeuObj@meta.data)]
#     CC_Anno2.df$TestID <- colnames(CC2.df)
#     CC_Anno.df <- rbind(CC_Anno.df,CC_Anno2.df)
#
#     CC_CT.df <- data.frame(Cell_type = scRNA.SeuObj@meta.data[,"Cell_type"])
#     CC.df <- cbind(CC.df, CC2.df)
#     rm(CC_CT.df,CC2.df,CC_Anno2.df)
#
#     CC.df$Cell_type <- as.character(CC.df$Cell_type)
#     CC.df <- rbind(CC.df,"Other")
#
# ##### Verification (CellCheck) #####
#   #### Install ####
#   ## Check whether the installation of those packages is required
#   Package.set <- c("tidyverse","caret","cvms","DescTools","devtools","ggthemes")
#   for (i in 1:length(Package.set)) {
#     if (!requireNamespace(Package.set[i], quietly = TRUE)){
#       install.packages(Package.set[i])
#     }
#   }
#   ## Load Packages
#   # library(Seurat)
#   lapply(Package.set, library, character.only = TRUE)
#   rm(Package.set,i)
#
#   ## install CellCheck
#   # Install the CellCheck package
#   detach("package:CellCheck", unload = TRUE)
#   devtools::install_github("Charlene717/CellCheck")
#   # Load CellCheck
#   library(CellCheck)
#
#   #### Run CellCheck ####
#
#   ## For one prediction
#   DisCMSet.lt = list(Mode = "One", Actual = "Cell_type", Predict = "singleR_PredbyCTDB_classic_qua0.6_tun0.03_sd1" , FilterSet1 = "Tool", FilterSet2 = "singleR" , Remark = "") # Mode = c("One","Multiple")
#   BarChartSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "Set", Group = "Tool", Remark = "")
#   LinePlotSet.lt <- list(Mode = "One", Metrics = "Balanced.Accuracy", XValue = "tune_Thr", Group = "Set", Remark = "")
#   CCR_cm_DisMult.lt <- CellCheck_DisMult(CC.df, CC_Anno.df,
#                                          DisCMSet.lt = DisCMSet.lt,
#                                          BarChartSet.lt = BarChartSet.lt,
#                                          LinePlotSet.lt = LinePlotSet.lt,
#                                          Save.Path = Save.Path, ProjectName = paste0("CellCheck_",ProjectName_Ano))
#
#   ## For multiple prediction
#   DisCMSet.lt = list(Mode = "Multiple", Actual = "Cell_type", FilterSet1 = "Tool", FilterSet2 = "singleR" , Remark = "_All") # Mode = c("One","Multiple")
#   BarChartSet.lt <- list(Mode = "Multiple", XValue = "Set", Group = "Tool", Remark = "_All")
#   LinePlotSet.lt <- list(Mode = "Multiple", XValue = "tune_Thr", Group = "Set", Remark = "_All")
#   Sum_DisMult.df <- CellCheck_DisMult(CC.df, CC_Anno.df,
#                                       DisCMSet.lt = DisCMSet.lt,
#                                       BarChartSet.lt = BarChartSet.lt,
#                                       LinePlotSet.lt=LinePlotSet.lt,
#                                       Save.Path = Save.Path, ProjectName = paste0("CellCheck_",ProjectName_Ano))
#
#

# ##### Session information #####
#   sessionInfo()
#   ## Ref: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
#   writeLines(capture.output(sessionInfo()), paste0(Save.Path,"/sessionInfo.txt"))
#
# ##### Save RData #####
#   save.image(paste0(Save.Path,"/SeuratObject_",ProjectName_Ano,".RData"))



