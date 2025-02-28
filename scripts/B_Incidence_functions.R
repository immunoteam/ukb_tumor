Packages <- c("tidyverse", "data.table", "Rfast", "magrittr", "fastmatch", "survival", "tidycmprsk", "ggsurvfit", "cowplot", "gtsummary", "forestmodel")
lapply(Packages, require, character.only = TRUE)
rm(Packages)
options(dplyr.summarise.inform = FALSE)


fun_TumorIncForestPlot = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
  #if(c("PTVgenes") %in% colnames(dataset))
  #plan(future::cluster, workers = threads)
  if(control == "nottumorous" & gender == "all") {
    control_data = readRDS(paste0(dataset_dir, "nottumorous_all.rds"))
  } else if(control == "nottumorous" & gender == "female") {
    control_data = readRDS(paste0(dataset_dir, "nottumorous_female.rds"))
  } else if(control == "nottumorous" & gender == "male") {
    control_data = readRDS(paste0(dataset_dir, "nottumorous_male.rds"))
  }
  m = max(control_data$date_of_death, na.rm = T)
  control_data %<>% 
    mutate(time = case_when(
      death == T ~ as.numeric(date_of_death - birthdate),
      death == F ~ as.numeric(m - birthdate)
    )) %>% 
    mutate(status = 0)
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = readRDS(paste0(dataset_dir, "tumorous_all.rds"))
    } else if(tmr == "all" & gender != "all") {
      tmr_data = readRDS(paste0(dataset_dir, "tumorous_", gender, ".rds"))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = readRDS(paste0(dataset_dir, tmr, ".rds"))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = readRDS(paste0(dataset_dir, tmr, "_", gender, ".rds"))
    }
  }))
  tumor_data %<>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(time = as.numeric((diag_date - birthdate))) %>% 
    mutate(status = 1)
  tempdf = rbind(tumor_data, control_data)
  rm(tumor_data, control_data)
  ptvb_MAF104 = fread(paste0(ptv_dir, "MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv"))
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tempdf$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tempdf %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tempdf$PTVb = as.factor(ifelse(tempdf$PTVb == 0, 0, 1))}
  case_n_ptvb_1 = tempdf %>% filter(status == 1, PTVb == 1) %>% nrow()
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
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  fm$plot + labs(title = tt)
}


#########################x


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
  m = max(control_data$date_of_death, na.rm = T)
  control_data %<>% 
    mutate(time = case_when(
      death == T ~ as.numeric(date_of_death - birthdate),
      death == F ~ as.numeric(m - birthdate)
    )) %>% 
    mutate(status = 0)
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
  tumor_data %<>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(time = as.numeric((diag_date - birthdate))) %>% 
    mutate(status = 1)
  tempdf = rbind(tumor_data, control_data)
  rm(tumor_data, control_data)
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
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
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  fm$plot + labs(title = tt)
}

#Perform a Cox model and return the P-values and the coefficients
fun_TumorIncCoxModel = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads) {
  if(control == "nottumorous" & gender == "all") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
  } else if(control == "nottumorous" & gender == "female") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_female.rds")
  } else if(control == "nottumorous" & gender == "male") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_male.rds")
  }
  m = max(control_data$date_of_death, na.rm = T)
  control_data %<>% 
    mutate(time = case_when(
      death == T ~ as.numeric(date_of_death - birthdate),
      death == F ~ as.numeric(m - birthdate)
    )) %>% 
    mutate(status = 0)
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
  tumor_data %<>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(time = as.numeric((diag_date - birthdate))) %>% 
    mutate(status = 1)
  
  tempdf = rbind(tumor_data, control_data)
  rm(tumor_data, control_data)
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
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
}

#Perform a Fisher's exact test and return the P-value and the odds ratio
fun_TumorIncFishTest = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads) {
  if(control == "nottumorous" & gender == "all") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
  } else if(control == "nottumorous" & gender == "female") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_female.rds")
  } else if(control == "nottumorous" & gender == "male") {
    control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_male.rds")
  }
  m = max(control_data$date_of_death, na.rm = T)
  control_data %<>% 
    mutate(time = case_when(
      death == T ~ as.numeric(date_of_death - birthdate),
      death == F ~ as.numeric(m - birthdate)
    )) %>% 
    mutate(status = 0)
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
  tumor_data %<>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(time = as.numeric((diag_date - birthdate))) %>% 
    mutate(status = 1)
  
  tempdf = rbind(tumor_data, control_data)
  rm(tumor_data, control_data)
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tempdf$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tempdf %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tempdf$PTVb = as.factor(ifelse(tempdf$PTVb == 0, 0, 1))}
  c(Fishers_exact_test_P = fisher.test(table(tempdf$status, tempdf$PTVb))$p.value, Fishers_exact_test_OR = unname(fisher.test(table(tempdf$status, tempdf$PTVb))$estimate))
}





fun_TumorIncPlot = function(tumor = "Melanoma", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads) {
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(gender == "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, ".rds"))
    } else {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, "_", gender, ".rds"))
    }
  }))
  tumor_data$death_type[tumor_data$death == T & is.na(tumor_data$death_type)] = "other"
  tumor_data %<>% mutate(cancerCausesDeath = death_type == cancer_type)
  tumor_data %<>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(time = case_when(
      death == T ~ as.numeric(date_of_death - birthdate),
      death == F ~ as.numeric(max(date_of_death, na.rm = T) - birthdate))) %>% 
    mutate(status = case_when(
      death == T & cancerCausesDeath == T ~ 1,
      death == F ~ 2,
      death == T & cancerCausesDeath == F ~ 3
    )) %>% 
    transform(status = factor(status, levels = c(1,2,3)))
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tumor_data$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tumor_data %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tumor_data$PTVb = as.factor(ifelse(tumor_data$PTVb == 0, 0, 1))}
  
  # cuminc(Surv(time, status) ~ 1, data = tumor_data)
  # 
  # cuminc(Surv(time, status) ~ 1, data = tumor_data) %>% 
  #   ggcuminc() + 
  #   labs(
  #     x = "Days"
  #   ) + 
  #   add_confidence_interval() +
  #   add_risktable()
  
  # cuminc(Surv(time, status) ~ 1, data = tumor_data) %>% 
  #   ggcuminc(outcome = c("3")) +
  #   ylim(c(0, 1)) + 
  #   labs(
  #     x = "Days"
  #   )
  
  cuminc(Surv(time, status) ~ PTVb, data = tumor_data) %>%
    tbl_cuminc(
      times = seq(5000,30000,5000),
      label_header = "**{time/365.25}-year cuminc**") %>%
    add_p()
  # 
  cuminc(Surv(time, status) ~ PTVb, data = tumor_data) %>%
    ggcuminc() +
    labs(
      x = "Days"
    ) +
    add_confidence_interval() +
    add_risktable()
  
  #Competing risk regression
  #I. Subdistribution hazards; HR>1 - significantly associated with increased hazard of death due to melanoma
  crr(Surv(time, status) ~ PTVb, data = tumor_data)
  crr(Surv(time, status) ~ PTVb, data = tumor_data) %>%
    tbl_regression(exp = TRUE)
  #II. Cause-specific hazards
  coxph(Surv(time, ifelse(status == 1, 1, 0)) ~ PTVb, data = tumor_data) %>% tbl_regression(exp = TRUE)
  
  # tumor_data$inc_event = as.factor(ifelse(tumor_data$status == 1, 0, 1))
  # inc_obj = cuminc(Surv(time, inc_event) ~ PTVb, data = tumor_data)
  # grayp = inc_obj$cmprsk$Tests[,"pv"] %>% round(digits = 4)
  # if(gender == "all") {
  #   tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  # } else {
  #   tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  # }
  # inc_plot = ggcuminc(x = inc_obj) +
  #   add_confidence_interval() +
  #   scale_ggsurvfit() +
  #   annotate(geom = "text", label = paste0("\nP = ", grayp), x = 5000, y = Inf ) +
  #   xlab("Time (days)") +
  #   ggtitle(tt) +
  #   theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), legend.position = "right")
  # 
  # inc_table = tbl_cuminc(x = cuminc(Surv(time, inc_event) ~ PTVb, data = tumor_data), label_header = "{time} days", times = seq(5000,30000,5000)) %>% 
  #   add_nevent(location = c("label", "level")) %>% 
  #   add_n(location = c("label", "level")) %>% 
  #   add_p() %>%
  #   modify_fmt_fun(p.value ~ function(x) ifelse(is.na(x), NA, format(x, digits = 3, scientific = TRUE)))
  # inc_table_df = as.data.frame(inc_table)
  # inc_table_grob = tableGrob(inc_table_df)
  # library()
  # out_fig = ggdraw() + 
  #   draw_plot(inc_plot, x = 0, y = .5, width = 1, height = .5) +
  #   draw_plot(inc_table_grob, x = 0.15, y = 0, width = 0.7, height = .5)
  # out_fig
}

fun_TumorIncForestPlot(tumor = "Melanoma", geneset = genes)
fun_TumorIncFishTest(tumor = "Melanoma", geneset = genes)
fun_TumorIncForestPlot(tumor = "Melanoma", gender = "female", geneset = genes)
fun_TumorIncForestPlot(tumor = "Melanoma", gender = "male", geneset = genes)

fun_TumorSurvPlotOS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads) {
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
    } else if(tmr == "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_", gender, ".rds"))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, ".rds"))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, "_", gender, ".rds"))
    }
  }))
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tumor_data$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tumor_data %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tumor_data$PTVb = as.factor(ifelse(tumor_data$PTVb == 0, 0, 1))}
  m = max(tumor_data$date_of_death, na.rm = T)
  tumorous_d = tumor_data %>%
    filter(death == T) %>%
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(surv_time = as.numeric(date_of_death - diag_date))
  tumorous_alive = tumor_data %>% 
    filter(death == F) %>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(surv_time = as.numeric(m - diag_date))
  survdf = bind_rows(tumorous_d, tumorous_alive)
  survdf$surv_event = as.numeric(ifelse(survdf$death == T, 1, 0))
  case_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  survfit2(Surv(surv_time, surv_event) ~ PTVb, data = survdf) %>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall survival probability"
    ) + 
    add_confidence_interval() +
    add_pvalue(location = c("annotation")) +
    add_risktable() +
    labs(title = tt)
}

#Disease specific
fun_TumorSurvPlotDS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads) {
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
    } else if(tmr == "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/tumorous_", gender, ".rds"))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, ".rds"))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = readRDS(paste0("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", tmr, "_", gender, ".rds"))
    }
  }))
  ptvb_MAF104 = fread("/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tumor_data$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tumor_data %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tumor_data$PTVb = as.factor(ifelse(tumor_data$PTVb == 0, 0, 1))}
  m = max(tumor_data$date_of_death, na.rm = T)
  tumorous_d = tumor_data %>%
    filter(death == T & death_type == cancer_type) %>%
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(surv_time = as.numeric(date_of_death - diag_date))
  tumorous_alive = tumor_data %>% 
    filter(death == F) %>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(surv_time = as.numeric(m - diag_date))
  survdf = bind_rows(tumorous_d, tumorous_alive)
  survdf$surv_event = as.numeric(ifelse(survdf$death == T, 1, 0))
  case_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  survfit2(Surv(surv_time, surv_event) ~ PTVb, data = survdf) %>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Disease specific survival probability"
    ) + 
    add_confidence_interval() +
    add_pvalue(location = c("annotation")) +
    add_risktable() +
    labs(title = tt)
}





res_fishtest = readRDS("Res/010/res_fishtest_104.rds")
genes = res_fishtest %>% dplyr::filter(category == "Melanoma", caseY > 5) %>% arrange(desc(OR)) %>% dplyr::slice(1:20) %>% pull(ENSEMBL)

load("Objects/go_desc_list")
hugos = go_desc_list[["antigen processing and presentation"]]
library(clusterProfiler)
library(org.Hs.eg.db)
genes = unique(bitr(geneID = hugos, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")$ENSEMBL)
genes = genes[!is.na(genes)]
rm(go_desc_list, hugos)

broca_geneset = c("ALK", "APC", "ATM", "ATR", "AXIN2", "BAP1", "BARD1", "BMPR1A", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDK4", "CDK12", "CDKN2A", "CHEK2", "CTNNA1", "DICER1", "EPCAM", "FANCM", "FH", "FLCN", "GEN1", "GREM1", "HOXB13", "MEN1", "MET", "MITF", "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "MUTYH", "NBN", "NF1", "NF2", "NTHL1", "PALB2", "PHOX2B", "PIK3CA", "PMS2", "POLD1", "POLE", "POT1", "PRKAR1A", "PTCH1", "PTEN", "RAD51B", "RAD51C", "RAD51D", "RB1", "RECQL", "RET", "RNF43", "RPS20", "SDHB", "SDHC", "SDHD", "SMAD4", "SMARCA4", "TP53", "TSC1", "TSC2", "VHL")
genes = unique(bitr(geneID = broca_geneset, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")$ENSEMBL)


fun_TumorSurvPlotDS(tumor = "Melanoma", gender = "female", geneset = genes)
fun_TumorSurvPlotOS(tumor = "Melanoma", gender = "female", geneset = genes)
#Overall survival plot function

#Disease specific survival plot function

