setwd("/media/balazs/WorkL/balazs/Work/UKBiobank/")
#Load PTV burden data
library(data.table)
library("clusterProfiler")
library("org.Hs.eg.db")

ptvb_wMAF = fread("PTVvars/PTV_2024_05/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
ptvb_wMAF = fread("PTVvars/dl_20250331/MAF_10-3/PTVBurden_with_Shetscores.tsv")

genes_all = data.frame(ENSEMBL = sort(unique(unlist(strsplit(ptvb_wMAF$Genes, ","), use.names = F))))
gene_symbols = bitr(geneID = genes_all$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
genes_all$SYMBOL = gene_symbols$SYMBOL[fmatch(genes_all$ENSEMBL, gene_symbols$ENSEMBL)]
rm(gene_symbols)

load("Objects/go_desc_list")
go_desc_list[["antigen processing and presentation of exogenous antigen"]]
go_desc_list = go_desc_list[grep("antigen processing", names(go_desc_list), value = T)]
go_desc_list = go_desc_list[lengths(go_desc_list) >= 5]
go_desc_list = go_desc_list[c(1:3,5:9,11:12)]

go_desc_df = lapply(1:length(go_desc_list), function(i) {
  data.frame(GO_description = names(go_desc_list)[i], ENSEMBL = genes_all$ENSEMBL[fmatch(go_desc_list[[i]], genes_all$SYMBOL)], SYMBOL = go_desc_list[[i]])
}) %>% bind_rows()
save(go_desc_df, file = "Objects/go_desc_df_ant_proc")

#Match Hugo to ENSG
go_desc_list_ensg = lapply(go_desc_list, function(x) {
  g = genes_all$ENSEMBL[fmatch(x, genes_all$SYMBOL)]
  g[!is.na(g)]
})

plan(multisession(workers = 2))
temp = future_sapply(ptvb_wMAF$Genes, function(x) {
  g = strsplit(x, ",")[[1]]
  future_sapply(go_desc_list_ensg, function(y) sum(g %fin% y))
})
dim(temp)
temp = t(temp)
colnames(temp) = gsub(",", "", gsub(" |-", "_", names(go_desc_list_ensg)))

ptvb_antproc = cbind(ptvb_wMAF[,1], temp)

ensgs_all = unique(unlist(go_desc_list_ensg, use.names = F))
temp = future_sapply(ptvb_wMAF$Genes, function(x) {
  g = strsplit(x, ",")[[1]]
  sum(g %fin% ensgs_all)
})
ptvb_antproc$antigen_processing_all = unname(temp)

apply(ptvb_antproc[,2:ncol(ptvb_antproc)], 2, Table)

save(ptvb_antproc, file = "Objects/ptvb_antproc_103")
save(ptvb_antproc, file = "Objects/ptvb_antproc_104")

load("Objects/ptvb_antproc_104")

##MAF10-3

