#Chunk 1
ukb_data = fread("variables/table-exporter_2025-04-04_11-58-24_data.csv")
View(head(ukb_data))
ukb_data %<>% 
  filter(`Genetic ethnic grouping` == "Caucasian") %>% 
  filter(`Skin colour | Instance 0` != "Black", `Skin colour | Instance 1` != "Black", `Skin colour | Instance 2` != "Black") #409,333
ukb_data %<>% select(-c(`Genetic ethnic grouping`, `Skin colour | Instance 0`, `Skin colour | Instance 1`, `Skin colour | Instance 2`))
colnames(ukb_data)[1:4] = c("eid", "year_of_birth", "month_of_birth", "sex")
ukb_data %<>% 
  transform(month_of_birth = match(month_of_birth, month.name)) %>% 
  mutate(date_of_birth = paste(year_of_birth, month_of_birth, "15", sep = "-"), .before = year_of_birth) %>% 
  transform(date_of_birth = as.IDate(date_of_birth)) %>% 
  select(-c(year_of_birth, month_of_birth))

ukb_data %>% filter(!is.na(`Date of death | Instance 1`)) %>% mutate(eq = `Date of death | Instance 0` == `Date of death | Instance 1`, .after = sex) %>% View()
ukb_data %<>% select(-`Date of death | Instance 1`)
colnames(ukb_data)[4] = "date_of_death"

m = max(ukb_data$date_of_death, na.rm = T)
ukb_data %<>% 
  mutate(age = round(ifelse(is.na(date_of_death), (as.numeric(m) - as.numeric(date_of_birth))/365.25, (as.numeric(date_of_death) - as.numeric(date_of_birth))/365.25),1), .after = date_of_death)

#Chunk 2
icd10_unique = ukb_data %>% select(grep("ICD10", colnames(ukb_data), value = T)) %>% unlist(, use.names = F) %>% unique()
icd10_unique = icd10_unique[icd10_unique != ""]
icd10_tbl = data.frame(icd10_long_desc = icd10_unique)

plan(multisession(workers = 4))
icd10_tbl$icd10 = future_sapply(icd10_tbl$icd10_long_desc, function(x) strsplit(x, " ")[[1]][1])

icd10_cancer_tbl = icd10_tbl %>% 
  filter(substr(icd10,1,1) == "C") %>% 
  mutate(icd10_brief = substr(icd10,1,3)) %>% 
  arrange(icd10) %>% 
  mutate(group1_desc = NA, group2_desc = NA)

#https://www.icd10data.com/ICD10CM/Codes/C00-D49/C51-C58/C54-

icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% c("C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14")] = "Malignant neoplasms of lip, oral cavity and pharynx"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 15:26)] = "Malignant neoplasms of digestive organs"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 30:39)] = "Malignant neoplasms of respiratory and intrathoracic organs"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 40:41)] = "Malignant neoplasms of bone and articular cartilage"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 43:44)] = "Melanoma and other malignant neoplasms of skin"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 45:49)] = "Malignant neoplasms of mesothelial and soft tissue"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief == "C50"] = "Malignant neoplasms of breast"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 51:58)] = "Malignant neoplasms of female genital organs"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 60:63)] = "Malignant neoplasms of male genital organs"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 64:68)] = "Malignant neoplasms of urinary tract"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 69:72)] = "Malignant neoplasms of eye, brain and other parts of central nervous system"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 73:75)] = "Malignant neoplasms of thyroid and other endocrine glands"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 76:80)] = "Malignant neoplasms of ill-defined, other secondary and unspecified sites"
icd10_cancer_tbl$group1_desc[icd10_cancer_tbl$icd10_brief %in% paste0("C", 81:96)] = "Malignant neoplasms of lymphoid, hematopoietic and related tissue"




icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C00"] = "Malignant neoplasm of lip"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C01"] = "Malignant neoplasm of base of tongue"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C02"] = "Malignant neoplasm of other and unspecified parts of tongue"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C03"] = "Malignant neoplasm of gum"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C04"] = "Malignant neoplasm of floor of mouth"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C05"] = "Malignant neoplasm of palate"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C06"] = "Malignant neoplasm of other and unspecified parts of mouth"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C07"] = "Malignant neoplasm of parotid gland"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C08"] = "Malignant neoplasm of other and unspecified major salivary glands"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C09"] = "Malignant neoplasm of tonsil"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C10"] = "Malignant neoplasm of oropharynx"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C11"] = "Malignant neoplasm of nasopharynx"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C12"] = "Malignant neoplasm of pyriform sinus"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C13"] = "Malignant neoplasm of hypopharynx"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C14"] = "Malignant neoplasm of other and ill-defined sites in the lip, oral cavity and pharynx"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C15"] = "Malignant neoplasm of esophagus"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C16"] = "Malignant neoplasm of stomach"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C17"] = "Malignant neoplasm of small intestine"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C18"] = "Malignant neoplasm of colon"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C19"] = "Malignant neoplasm of rectosigmoid junction"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C20"] = "Malignant neoplasm of rectum"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C21"] = "Malignant neoplasm of anus and anal canal"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C22"] = "Malignant neoplasm of liver and intrahepatic bile ducts"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C23"] = "Malignant neoplasm of gallbladder"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C24"] = "Malignant neoplasm of other and unspecified parts of biliary tract"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C25"] = "Malignant neoplasm of pancreas"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C26"] = "Malignant neoplasm of other and ill-defined digestive organs"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C30"] = "Malignant neoplasm of nasal cavity and middle ear"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C31"] = "Malignant neoplasm of accessory sinuses"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C32"] = "Malignant neoplasm of larynx"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C33"] = "Malignant neoplasm of trachea"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C34"] = "Malignant neoplasm of bronchus and lung"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C37"] = "Malignant neoplasm of thymus"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C38"] = "Malignant neoplasm of heart, mediastinum and pleura"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C39"] = "Malignant neoplasm of other and ill-defined sites in the respiratory system and intrathoracic organs"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C40"] = "Malignant neoplasm of bone and articular cartilage of limbs"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C41"] = "Malignant neoplasm of bone and articular cartilage of other and unspecified sites"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C43"] = "Malignant melanoma of skin"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C44"] = "Other and unspecified malignant neoplasm of skin"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C45"] = "Mesothelioma"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C46"] = "Kaposi's sarcoma"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C47"] = "Malignant neoplasm of peripheral nerves and autonomic nervous system"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C48"] = "Malignant neoplasm of retroperitoneum and peritoneum"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C49"] = "Malignant neoplasm of other connective and soft tissue"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C50"] = "Malignant neoplasm of breast"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C51"] = "Malignant neoplasm of vulva"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C52"] = "Malignant neoplasm of vagina"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C53"] = "Malignant neoplasm of cervix uteri"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C54"] = "Malignant neoplasm of corpus uteri"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C55"] = "Malignant neoplasm of uterus, part unspecified"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C56"] = "Malignant neoplasm of ovary"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C57"] = "Malignant neoplasm of other and unspecified female genital organs"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C58"] = "Malignant neoplasm of placenta"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C60"] = "Malignant neoplasm of penis"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C61"] = "Malignant neoplasm of prostate"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C62"] = "Malignant neoplasm of testis"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C63"] = "Malignant neoplasm of other and unspecified male genital organs"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C64"] = "Malignant neoplasm of kidney, except renal pelvis"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C65"] = "Malignant neoplasm of renal pelvis"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C66"] = "Malignant neoplasm of ureter"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C67"] = "Malignant neoplasm of bladder"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C68"] = "Malignant neoplasm of other and unspecified urinary organs"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C69"] = "Malignant neoplasm of eye and adnexa"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C70"] = "Malignant neoplasm of meninges"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C71"] = "Malignant neoplasm of brain"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C72"] = "Malignant neoplasm of spinal cord, cranial nerves and other parts of central"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C73"] = "Malignant neoplasm of thyroid gland"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C74"] = "Malignant neoplasm of adrenal gland"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C75"] = "Malignant neoplasm of other endocrine glands and related structures"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C76"] = "Malignant neoplasm of other and ill-defined sites"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C77"] = "Secondary and unspecified malignant neoplasm of lymph nodes"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C78"] = "Secondary malignant neoplasm of respiratory and digestive organs"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C79"] = "Secondary malignant neoplasm of other and unspecified sites"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C80"] = "Malignant neoplasm without specification of site"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C81"] = "Hodgkin lymphoma"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C82"] = "Follicular lymphoma"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C83"] = "Non-follicular lymphoma"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C84"] = "Mature T/NK-cell lymphomas"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C85"] = "Other specified and unspecified types of non-Hodgkin lymphoma"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C86"] = "Other specified types of T/NK-cell lymphoma"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C88"] = "Malignant immunoproliferative diseases and certain other B-cell lymphomas"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C90"] = "Multiple myeloma and malignant plasma cell neoplasms"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C91"] = "Lymphoid leukemia"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C92"] = "Myeloid leukemia"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C93"] = "Monocytic leukemia"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C94"] = "Other leukemias of specified cell type"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C95"] = "Leukemia of unspecified cell type"
icd10_cancer_tbl$group2_desc[icd10_cancer_tbl$icd10_brief == "C96"] = "Other and unspecified malignant neoplasms of lymphoid, hematopoietic and related tissue"


icd10_cancer_tbl$cancer_group = NA
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C00"] = "lip"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief %in% c("C01", "C02", "C03", "C04", "C05", "C06")] = "oral"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C15"] = "esophagus"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C16"] = "stomach"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C17"] = "small_intestine"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C18"] = "colon"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief %in% c("C19", "C20", "C21")] = "rectum"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C22"] = "liver"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C25"] = "pancreas"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C32"] = "larynx"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C34"] = "lung"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C37"] = "thymus"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C43"] = "melanoma"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C45"] = "mesothelioma"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C50"] = "breast"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief %in% c("C53", "C54", "C55")] = "uterus"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C56"] = "ovary"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C61"] = "prostate"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C62"] = "testis"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief %in% c("C64", "C65")] = "kidney"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C67"] = "bladder"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief %in% c("C70", "C71", "C72")] = "brain_cns"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C73"] = "thyroid"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C81"] = "hodgkin_lymphoma"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief %in% c("C82", "C83", "C84", "C85", "C86", "C88")] = "non_hodgkin_lymphoma"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief == "C90"] = "multiple_myeloma"
icd10_cancer_tbl$cancer_group[icd10_cancer_tbl$icd10_brief %in% c("C91", "C92", "C93", "C94", "C95")] = "leukemia"

#Chunk 3
plan(multisession(workers = 4))
ukb_data$havecancer = future_apply(ukb_data[,8:29], 1, function(x) any(x != ""))
ukb_data_cancer = ukb_data %>% filter(havecancer == T)

ukb_data_cancerL = ukb_data_cancer %>% 
  select(c(eid, grep("Type of cancer", colnames(ukb_data_cancer), value = T))) %>% 
  pivot_longer(cols = 2:23, names_to = "instance", values_to = "icd10_long_desc")
ukb_data_cancerL$instance = gsub("Type of cancer: ICD10 \\| Instance ", "Instance ", ukb_data_cancerL$instance)

ukb_data_cancerDate = ukb_data_cancer %>% 
  select(c(eid, grep("Date of cancer", colnames(ukb_data_cancer), value = T))) %>% 
  pivot_longer(cols = 2:23, names_to = "instance", values_to = "date_of_diagnosis")
ukb_data_cancerDate$instance = gsub("Date of cancer diagnosis \\| Instance ", "Instance ", ukb_data_cancerDate$instance)
ukb_data_cancerL %<>% left_join(ukb_data_cancerDate, by = c("eid", "instance"))
ukb_data_cancerL %<>% filter(icd10_long_desc != "")
rm(ukb_data_cancerDate)

tempdf = icd10_cancer_tbl %>% select(icd10_long_desc, cancer_group)
ukb_data_cancerL %<>% left_join(tempdf, by = "icd10_long_desc")

#ICD10 long description - no. of patients
pats_icd10_long = future_lapply(icd10_cancer_tbl$icd10_long_desc, function(x) {
  ukb_data_cancerL %>% filter(icd10_long_desc == x) %>% pull(eid) %>% unique()
})
names(pats_icd10_long) = icd10_cancer_tbl$icd10_long_desc
icd10_cancer_tbl$nPat = lengths(pats_icd10_long)

#Cancer groups - no. of patients
cancers = icd10_cancer_tbl %>% pull(cancer_group) %>% unique()
cancers = sort(cancers[!is.na(cancers)])

pats_cancer_group = future_lapply(cancers, function(x) {
  icd10ld = icd10_cancer_tbl %>% filter(cancer_group == x) %>% pull(icd10_long_desc)
  ukb_data_cancerL %>% filter(icd10_long_desc %in% icd10ld) %>% pull(eid) %>% unique()
})
names(pats_cancer_group) = cancers

tempdf = data.frame(cancer_group = cancers, nPatCG = lengths(pats_cancer_group))
icd10_cancer_tbl %<>% left_join(tempdf, by = "cancer_group")
rm(tempdf)

#saveRDS(icd10_cancer_tbl, file = "variables/icd10_cancer_tbl.rds")

#Chunk 4
ukb_data %<>% mutate(death = !is.na(date_of_death), .after = date_of_death)
m = max(ukb_data$date_of_death, na.rm = T)
#Lung, male
names(pats_cancer_group)
ids = pats_cancer_group[["lung"]]

case_df = ukb_data %>% filter(eid %in% ids) %>% select(eid, date_of_birth, sex, age, date_of_death, death)
diagnosis_tbl = ukb_data_cancerL %>% filter(eid %in% ids, cancer_group == "lung") %>% arrange(date_of_diagnosis) %>% distinct(eid, .keep_all = T) %>% select(-instance)

case_df %<>% left_join(diagnosis_tbl, by = "eid")
case_df %<>% 
  mutate(time = as.numeric(difftime(date_of_diagnosis, date_of_birth))) %>% 
  mutate(status = 1) %>% 
  select(eid, sex, time, status)


nottumorous_df = ukb_data %>% filter(!eid %in% ukb_data_cancer$eid) %>% select(eid, date_of_birth, sex, age, date_of_death, death)

nottumorous_df %<>% 
  mutate(time = case_when(
    death == T ~ as.numeric(difftime(date_of_death, date_of_birth, units = "days")),
    death == F ~ as.numeric(difftime(m, date_of_birth, units = "days"))
  )) %>% 
  mutate(status = 0) %>% 
  select(eid, sex, time, status) 

maf104 = fread("variables/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
franco = readRDS("variables/55_FRANCO_BLOOD_SANOFI_PASTEUR_SA_INACTIVATED_INFLUENZA_VACCINE_CORRELATED_WITH_ANTIBODY_RESPONSE_AGE_18_40YO_14DY_POSITIVE.rds")
maf104$ptvb = franco
maf104 %<>% select(Patient.ID, ptvb) %>% transform(Patient.ID = as.numeric(Patient.ID))

case_control = rbind(case_df, nottumorous_df)
case_control %<>% left_join(maf104, by = c("eid" = "Patient.ID"))
case_control$ptvb = ifelse(case_control$ptvb == 0, 0, 1)
case_control$ptvb = as.factor(case_control$ptvb)

case_control_male = case_control %>% filter(sex == "Male")
case_control_female = case_control %>% filter(sex == "Female")

coxph(Surv(time, status) ~ ptvb, data = case_control_male, id = case_control_male$eid)
coxph(Surv(time, status) ~ ptvb, data = case_control_female, id = case_control_female$eid)


influ = readRDS("variables/influ_gp_F.rds")
case_control$vacc = case_control$eid %in% influ$eid
coxph(Surv(time, status) ~ ptvb + vacc, data = case_control_male, id = case_control_male$eid)
coxph(Surv(time, status) ~ ptvb + vacc, data = case_control_female, id = case_control_female$eid)

coxph(Surv(time, status) ~ ptvb + vacc + factor(sex), data = case_control, id = case_control$eid)
