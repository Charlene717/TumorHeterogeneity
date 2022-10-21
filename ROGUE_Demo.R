## Ref: https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("PaulingLiu/ROGUE")

# suppressMessages(library(ROGUE))
# suppressMessages(library(ggplot2))
# suppressMessages(library(tidyverse))


library(ROGUE)
library(ggplot2)
library(tidyverse)

#### Load the data ####
expr <- readr::read_rds(file = "D:/Dropbox/##_GitHub/##_Charlene/TumorHeterogeneity/input_ROGUE_Demo/DC.rds.gz")
meta <- readr::read_rds(file = "D:/Dropbox/##_GitHub/##_Charlene/TumorHeterogeneity/input_ROGUE_Demo/info.rds.gz")

## Filtering out low-abundance genes and low-quality cells
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)

## Expression entropy model
ent.res <- SE_fun(expr)
head(ent.res)

## S-E plot
SEplot(ent.res)

## ROGUE calculation
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value

## Calculate the ROGUE value of each putative cluster for each sample
rogue.res <- rogue(expr, labels = meta$ct, samples = meta$Patient, platform = "UMI", span = 0.6)
rogue.res

## Visualize ROGUE values on a boxplot
rogue.boxplot(rogue.res)

