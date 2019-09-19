#load in the data
library(pacman)
p_load(PCAtools,ggplot2)
load("/Users/annaleigh/Documents/data/crick_hackathon/data/scrna_seq/PRJEB11202_SCE_norm.RData")
#first using the expression of GATA3 as the cutoff
tmp = as.data.table(melt(assays(sce)$logcounts))

gat3_expression = ggplot(tmp[Var1 == "GATA3"]) + geom_histogram(aes(x = value)) + 
  geom_vline(xintercept = 4) + 
  xlab("GATA3 Log Expression")  
  

expressing_cells = tmp[Var1 == "GATA3" & value > 4, Var2]  
  
meta_data = sce@colData
meta_data$gat3 = ifelse(rownames(meta_data) %in% expressing_cells, "GATA3+","GATA3-")
meta_data$gat3numeric = ifelse(rownames(meta_data) %in% expressing_cells, 0,1)

#quick pairplot to see how it looks
pca_norm_data <- PCAtools::pca(assays(sce)$logcounts, metadata = meta_data, removeVar = 0.1)
pairsplot(pca_norm_data, colby = "gat3")
eigencorplot(pca_norm_data,metavars = "gat3numeric")
#okay so loks like PC1, PC3, PC7 and 8 aren't bad at distinguishing
biplot(pca_norm_data, x = "PC7", y = "PC8", lab = "", colby = "gat3")
biplot(pca_norm_data, x = "PC1", y = "PC2", lab = "", colby = "gat3")

#what genes are contributing the most to this?
plotloadings(pca_norm_data,components = c("PC7","PC8"))
plotloadings(pca_norm_data,components = c("PC1","PC2"))
#uhh oohhh VCX DDX3Y, these are Y linked genes...
