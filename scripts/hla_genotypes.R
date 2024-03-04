library(data.table)
library(Rfast)
library(fastmatch)
data = fread("objects/ukbiobank_HLA.txt", sep = "\t")
cn = readLines("objects/hla_imp_header_alleles.txt")
cn = unlist(strsplit(cn, "\t"))
cn = unname(sapply(cn, function(x) {ifelse(nchar(x) %in% c(5,8), gsub("_", "_0", x), x)}))
length(unlist(strsplit(data$eid[1], ",")))
genotypes_list = lapply(data$eid, function(x) {unlist(strsplit(x, ","))})
Table(lengths(genotypes_list))

genotypes_hla = lapply(genotypes_list, function(x) {
  if(length(x) == 362) {
    x = as.numeric(x)
    x[x<0.7] = NA
    temp = cn[x>=.7 & x <= 1]
    temp = temp[!is.na(temp)]
    temp2 = cn[x>1]
    temp2 = temp2[!is.na(temp2)]
    temp2 = c(temp2, temp2)
    sort(c(temp, temp2))
  }
})
Table(lengths(genotypes_hla))
save(genotypes_hla, file = "D:/CloudStation/ukbiobank_tumor/genotypes_hla")

genotypes_eid = data$eid
saveRDS(genotypes_eid, file = "objects/genotypes_eid")

###
gergo_data = read.delim("D:/CloudStation/ukbiobank_tumor/gergo_covid/eids_covid_outcome.txt", header = F)
gergo_ids = gergo_data$V1

Table(gergo_ids %in% data$V1)

genotypes_covid = genotypes_hla[fmatch(gergo_ids, data$V1)]



genotypes_hla6 = lapply(genotypes_hla, function(x) {
  if(length(x) > 0) {
    x[substr(x,1,1) %in% c("A", "B", "C") | substr(x,1,4) %in% c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1")]
  }
})
Table(lengths(genotypes_hla6))
