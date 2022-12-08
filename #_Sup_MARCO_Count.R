

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(300000)

#### Load data #####
load("D:/Dropbox/##_GitHub/##_Charlene/TrajectoryAnalysis/2022-12-07_TrajAna_PCA_PDAC/scRNA.SeuObj_CDS_PRJCA001063_TrajAna_PCA.RData")

## Sub celltype
scRNA_Mac.SeuObj <- scRNA.SeuObj[,grepl("Macrophage cell",scRNA.SeuObj$singleR_classic_PredbyscRNA2)]

## Plot UMAP
DimPlot(scRNA_Mac.SeuObj, reduction = "umap",group.by = "Type", label = TRUE, pt.size = 0.5)

## Extract MetaData
Meta_Mac.df <- scRNA_Mac.SeuObj@meta.data %>% as.data.frame()
GeneEXP_Mac.df <- GetAssayData(scRNA_Mac.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix

GeneEXPT_Mac.df <- data.frame(Cell_ID = colnames(GeneEXP_Mac.df),GeneEXP_Mac.df %>% t() %>%  as.data.frame())
Mac_Tar.df <- left_join(Meta_Mac.df,GeneEXPT_Mac.df[,c("Cell_ID","MARCO","RAP1A")] )


## Plot BarPlot
## Ref: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)

Mac_Tar.df$Type <- factor(Mac_Tar.df$Type, levels = c("Ctrl","PT","MET"))

p <- ggboxplot(Mac_Tar.df, x = "Type", y = "MARCO",
               color = "Type", #palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test")



# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Ctrl", "PT"), c("PT", "MET") , c("Ctrl", "MET"))
ggboxplot(Mac_Tar.df, x = "Type", y = "MARCO",
          color = "Type", # palette = "jco",
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6)     # Add global p-value

