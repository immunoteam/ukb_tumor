#ICD10 codes
library(data.table)
var40006 = fread("D:/CloudStation/mygit_EXT/ukb_tumor/ukb_vars/40006.txt")
icd10 = sort(unique(var40006$value))
save(icd10, file = "D:/CloudStation/mygit_EXT/ukb_tumor/icd10_unique")
