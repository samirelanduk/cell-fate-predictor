
library(dplyr)
library(readr)

load("results/DGEA.RData")

rnet <- read_tsv("data/annotated_regions.txt") %>% 
  filter(`Nearest PromoterID` != "N/A")

# Some parsing of the TF column
rnet <- rnet %>% 
  mutate(TF = sapply(strsplit(TF, "[.]"), `[`, 3)) %>% 
  mutate(TF = sapply(strsplit(TF, "[(]"), `[`, 1)) %>% 
  mutate(TF = toupper(TF))

# Regulatory network construction
rnet_g3_pos <- rnet %>% 
  filter(`Nearest PromoterID` %in% gata3_pos) %>% 
  select(TF, `Nearest PromoterID`) %>% 
  rename(Gene = `Nearest PromoterID`)

rnet_g3_neg <- rnet %>% 
  filter(`Nearest PromoterID` %in% gata3_neg) %>% 
  select(TF, `Nearest PromoterID`) %>% 
  rename(Gene = `Nearest PromoterID`)

# Focus on TFs that are expressed in each cell population
library(Seurat)
library(Matrix)
load("results/Seurat_Morula.RData")

calc_avg_expr <- function(gene, expr, cl){
  gene <- strsplit(gene, "::")[[1]]
  # Get the expression of the genes in cluster cl
  idx <- gene %in% rownames(expr)
  expr_mtx <- expr@assays$RNA[gene[idx], expr@active.ident == cl]
  
  avg_expr <- rowMeans(expr_mtx)
  return(mean(avg_expr))
}

rnet_g3_pos <- rnet_g3_pos %>% 
  mutate(TF_expr = sapply(TF, calc_avg_expr, sce, 1)) %>% 
  filter(!is.na(TF_expr)) %>% 
  filter(TF_expr > median(TF_expr)) %>%
  mutate(pop = "GATA+")
rnet_g3_pos <- rnet_g3_pos[!duplicated(rnet_g3_pos),]

pos_nodes <- union(rnet_g3_pos$Gene, rnet_g3_pos$TF)

rnet_g3_neg <- rnet_g3_neg %>% 
  mutate(TF_expr = sapply(TF, calc_avg_expr, sce, 0)) %>% 
  filter(!is.na(TF_expr)) %>% 
  filter(TF_expr > median(TF_expr)) %>% 
  mutate(pop = "GATA-")
rnet_g3_neg <- rnet_g3_neg[!duplicated(rnet_g3_neg),]

neg_nodes <- union(rnet_g3_neg$Gene, rnet_g3_neg$TF)

rnet_g3 <- bind_rows(rnet_g3_pos, rnet_g3_neg)
rnet_g3 <- rnet_g3[!duplicated(rnet_g3),]

write_tsv(rnet_g3_pos, "results/rnet_g3_pos.tsv")
write_tsv(rnet_g3_neg, "results/rnet_g3_neg.tsv")
write_tsv(rnet_g3, "results/rnet_g3.tsv")

save(rnet_g3, rnet_g3_pos, rnet_g3_neg,
     pos_nodes, neg_nodes, file = "results/reg_net.RData")
