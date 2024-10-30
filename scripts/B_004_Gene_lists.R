#Gene lists
setwd("/media/balazs/WorkL/balazs/Work/UKBiobank/")
library(Rfast)
library(future.apply)

go_terms = fread("Objects/go_terms.txt") #https://amigo.geneontology.org/amigo/search/bioentity
go_terms %<>% filter(grepl("UniProt", V1))
go_terms$V1 = gsub("UniProtKB:", "", go_terms$V1)
go_terms$V10 = gsub("PANTHER:", "", go_terms$V10)
go_terms = go_terms[,c(2,1,6,7,5,10,11)]
colnames(go_terms) = c("Gene", "UniProtKB", "Protein", "GO_terms", "GO_descriptions", "PANTHER_family", "PANTHER_description")

go_terms_list = strsplit(go_terms$GO_descriptions, "\\|")
names(go_terms_list) = go_terms$Gene

go_desc_freq = data.frame(go_desc = sort(unique(unlist(strsplit(go_terms$GO_descriptions, "\\|"), use.names = F))))
plan(multisession(workers = 8))
go_desc_list = future_lapply(go_desc_freq$go_desc, function(x) names(go_terms_list)[sapply(go_terms_list, function(z) x %fin% z)])
names(go_desc_list) = go_desc_freq$go_desc
go_desc_freq$freq = lengths(go_desc_list)

save(go_desc_list, file = "Objects/go_desc_list")



