
#Generate 
fun_TumorIncForestPlot = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads) {
  #if(c("PTVgenes") %in% colnames(dataset))
  #plan(future::cluster, workers = threads)
  if(control == "nottumorous" & gender == "all") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
  } else if(control == "nottumorous" & gender == "female") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_female.rds")
  } else if(control == "nottumorous" & gender == "male") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_male.rds")
  }
  if(tumor == "all" & gender == "all") {
    tumor_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_all.rds")
  } else if(tumor == "all" & gender != "all") {
    tumor_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_", gender, ".rds"))
  } else if(tumor != "all" & gender == "all") {
    tumor_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tumor, ".rds"))
  } else if(tumor != "all" & gender != "all") {
    tumor_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tumor, "_", gender, ".rds"))
  }
  tempdf = rbind(tumor_data, control_data)
  rm(tumor_data, control_data)
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  library(Rfast)
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tempdf$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tempdf %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tempdf$PTVb = as.factor(ifelse(tempdf$PTVb == 0, 0, 1))}
  case_n_ptvb_1 = tempdf %>% filter(status == 2, PTVb == 1) %>% nrow()
  #Select predictors
  sel_gpcas = paste0("gpca", 1:gpca_nb)
  if(table(tempdf$PTVb[tempdf$cancer_type != "No_cancer"])["1"] != 0 & gender == "all") {
    tempdf %<>% dplyr::select(eid, time, status, PTVb, sex, all_of(sel_gpcas))
  } else if(table(tempdf$PTVb[tempdf$cancer_type != "No_cancer"])["1"] != 0 & gender != "all") {
    tempdf %<>% dplyr::select(eid, time, status, PTVb, all_of(sel_gpcas))
  } else if(table(tempdf$PTVb[tempdf$cancer_type != "No_cancer"])["1"] == 0 & gender == "all") {
    tempdf %<>% dplyr::select(eid, time, status, sex, all_of(sel_gpcas))
  } else {
    tempdf %<>% dplyr::select(eid, time, status, all_of(sel_gpcas))
  }
  myformula = as.formula(paste('Surv(time, status) ~ ', paste0(colnames(tempdf)[4:ncol(tempdf)], collapse = "+")))
  res.cox <- coxph(myformula, data = tempdf, id = tempdf$eid)
  fm = forest_model(model = res.cox, return_data = T)
  if(gender == "all") {
    tt = paste0(tumor, ", BOTH. N_PTVB_1 = ", case_n_ptvb_1)
  } else {
    tt = paste0(tumor, ", ", gender, ". N_PTVB_1 = ", case_n_ptvb_1)
  }
  fm$plot + labs(title = tt)
}




fun_TumorIncForestPlot(tumor = "Melanoma", geneset = tempensgs)
fun_TumorIncForestPlot(tumor = "Melanoma", gender = "male", geneset = tempensgs)
