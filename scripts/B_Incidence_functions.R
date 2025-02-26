library(tidycmprsk)
library(gtsummary)
#Generate Forest plot
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
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_all.rds")
    } else if(tmr == "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_", gender, ".rds"))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, ".rds"))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, "_", gender, ".rds"))
    }
  }))
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
    tt = paste0(tumor, ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(tumor, ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  fm$plot + labs(title = tt)
}

#Perform a Cox model and return the P-values and the coefficients
fun_TumorIncCoxModel = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads) {
  #if(c("PTVgenes") %in% colnames(dataset))
  #plan(future::cluster, workers = threads)
  if(control == "nottumorous" & gender == "all") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
  } else if(control == "nottumorous" & gender == "female") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_female.rds")
  } else if(control == "nottumorous" & gender == "male") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_male.rds")
  }
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_all.rds")
    } else if(tmr == "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_", gender, ".rds"))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, ".rds"))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, "_", gender, ".rds"))
    }
  }))
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
  coefs = summary(res.cox)$coefficients[,"exp(coef)"]
  names(coefs) = paste0("hr_", names(coefs))
  pvalues = summary(res.cox)$coefficients[,"Pr(>|z|)"]
  names(pvalues) = paste0("p_", names(pvalues))
  c(coefs, pvalues)
  # fm = forest_model(model = res.cox, return_data = T)
  # if(gender == "all") {
  #   tt = paste0(tumor, ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  # } else {
  #   tt = paste0(tumor, ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  # }
  # fm$plot + labs(title = tt)
}

#Perform a Fisher's exact test and return the P-value and the odds ratio
fun_TumorIncFishTest = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads) {
  #if(c("PTVgenes") %in% colnames(dataset))
  #plan(future::cluster, workers = threads)
  if(control == "nottumorous" & gender == "all") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
  } else if(control == "nottumorous" & gender == "female") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_female.rds")
  } else if(control == "nottumorous" & gender == "male") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_male.rds")
  }
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_all.rds")
    } else if(tmr == "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_", gender, ".rds"))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, ".rds"))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, "_", gender, ".rds"))
    }
  }))
  tempdf = rbind(tumor_data, control_data)
  rm(tumor_data, control_data)
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  library(Rfast)
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tempdf$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tempdf %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tempdf$PTVb = as.factor(ifelse(tempdf$PTVb == 0, 0, 1))}
  c(Fishers_exact_test_P = fisher.test(table(tempdf$status, tempdf$PTVb))$p.value, Fishers_exact_test_OR = unname(fisher.test(table(tempdf$status, tempdf$PTVb))$estimate))
}


#fun_TumorIncForestPlot(tumor = "Melanoma", geneset = genes)
# fun_TumorIncForestPlot(tumor = "Melanoma", gender = "all", geneset = tempensgs)
# fun_TumorIncCoxModel(tumor = "Melanoma", gender = "all", geneset = tempensgs)
# fun_TumorIncFishTest(tumor = "Melanoma", gender = "all", geneset = tempensgs)


fun_TumorIncPlot = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads) {
  #if(c("PTVgenes") %in% colnames(dataset))
  #plan(future::cluster, workers = threads)
  if(control == "nottumorous" & gender == "all") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
  } else if(control == "nottumorous" & gender == "female") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_female.rds")
  } else if(control == "nottumorous" & gender == "male") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_male.rds")
  }
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_all.rds")
    } else if(tmr == "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_", gender, ".rds"))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, ".rds"))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, "_", gender, ".rds"))
    }
  }))
  tempdf = rbind(tumor_data, control_data)
  rm(tumor_data, control_data)
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  library(Rfast)
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tempdf$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tempdf %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tempdf$PTVb = as.factor(ifelse(tempdf$PTVb == 0, 0, 1))}
  case_n_ptvb_1 = tempdf %>% filter(status == 2, PTVb == 1) %>% nrow()
  tempdf$inc_event = as.factor(ifelse(tempdf$status == 1, 0, 1))
  grayp = tidycmprsk::cuminc(Surv(time, inc_event) ~ PTVb, data = tempdf)$cmprsk$Tests[,"pv"] %>% round(digits = 6)
  inc_obj = cuminc(Surv(time, inc_event) ~ PTVb, data = tempdf)
  
  inc_plot = 
    ggcuminc(x = inc_obj) +
    add_confidence_interval() +
    scale_ggsurvfit() +
    annotate(geom = "text", label = paste0("\nP = ", grayp), x = 10, y = Inf ) +
    xlab("Time (days)") +
    ggtitle(paste0("Tumor: ", tumor, ". Gender: ", gender)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          legend.position = "right")
  inc_table = tbl_cuminc(x = cuminc(Surv(time, inc_event) ~ PTVb, data = tempdf), label_header = "{time} days", times = seq(5000,30000,5000)) %>% 
    add_nevent(location = c("label", "level")) %>% 
    add_n(location = c("label", "level")) %>% 
    add_p() %>%
    modify_fmt_fun(p.value ~ function(x) ifelse(is.na(x), NA, format(x, digits = 3, scientific = TRUE)))
  inc_table_df = as.data.frame(inc_table)
  inc_table_grob = tableGrob(inc_table_df)
  grid.arrange(inc_plot, inc_table_grob, layout_matrix = rbind(c(1, 1), c(1,1), c(2, 2), c(2, 2)))
}

#Overall survival plot function

#Disease specific survival plot function

