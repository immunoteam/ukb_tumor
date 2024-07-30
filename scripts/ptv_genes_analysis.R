library(dplyr)
library(magrittr)
library(ggplot2)
library(clusterProfiler)
ptvb = fread("objects/PTV_2024_05/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
compareCluster(geneClusters = strsplit(ptvb$Genes[1:10], ","), fun = "enrichGO", OrgDb = "org.Hs.eg.db")

groupGO(gene = strsplit(ptvb$Genes[ptvb$Patient.ID == "1858110"], ",")[[1]], OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

bitr(strsplit(ptvb$Genes[1:10], ",")[[1]], "ENSEMBL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)

gseGO(geneList = strsplit(ptvb$Genes[1:10], ",")[[1]],ont = "BP", OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL")
