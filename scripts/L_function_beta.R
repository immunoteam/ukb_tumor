

fun_TumorSurvPlotOS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
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
  ptvb_MAF104 = fread(paste0(ptv_dir, "MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv"))
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
  dead_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  alive_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 0, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH. Geneset size: ", length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, "Geneset size: ",  length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
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


fun_TumorSurvCoxOS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/000_Sub_dataset/", ptv_dir = "C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/PTVvars/") {
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
  ptvb_MAF104 = fread(paste0(ptv_dir, "MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv"))
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
  sel_gpcas = paste0("gpca", 1:gpca_nb)
  if(table(survdf$PTVb)["1"] != 0 & gender == "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, sex, age, all_of(sel_gpcas))
  } else if(table(survdf$PTVb)["1"] != 0 & gender != "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, age, all_of(sel_gpcas))
  } else if(table(survdf$PTVb)["1"] == 0 & gender == "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, sex, age, all_of(sel_gpcas))
  } else {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, age, all_of(sel_gpcas))
  }
  case_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
  res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
  out = summary(res.cox)$coefficients["PTVb1",c("exp(coef)","Pr(>|z|)")]
  names(out) = c("HR", "Pvalue")
  out
}


fun_TumorSurvForestOS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/000_Sub_dataset/", ptv_dir = "C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/PTVvars/") {
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
  ptvb_MAF104 = fread(paste0(ptv_dir, "MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv"))
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
  sel_gpcas = paste0("gpca", 1:gpca_nb)
  if(table(survdf$PTVb)["1"] != 0 & gender == "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, sex, age, all_of(sel_gpcas))
  } else if(table(survdf$PTVb)["1"] != 0 & gender != "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, age, all_of(sel_gpcas))
  } else if(table(survdf$PTVb)["1"] == 0 & gender == "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, sex, age, all_of(sel_gpcas))
  } else {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, age, all_of(sel_gpcas))
  }
  dead_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  alive_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 0, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH. Geneset size: ", length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, "Geneset size: ",  length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
  }
  myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
  res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
  fm = forest_model(model = res.cox, return_data = T)
  fm$plot + labs(title = tt)
}


fun_TumorSurvPlotDS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
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
  ptvb_MAF104 = fread(paste0(ptv_dir, "MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv"))
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
  dead_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  alive_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 0, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH. Geneset size: ", length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, "Geneset size: ",  length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
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


fun_TumorSurvCoxDS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
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
  ptvb_MAF104 = fread(paste0(ptv_dir, "MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv"))
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
  sel_gpcas = paste0("gpca", 1:gpca_nb)
  survdf = survdf %>% 
    dplyr::select(eid, surv_time, surv_event, PTVb, all_of(sel_gpcas))
  case_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
  res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
  out = summary(res.cox)$coefficients["PTVb1",c("exp(coef)","Pr(>|z|)")]
  names(out) = c("HR", "Pvalue")
  out
}


fun_TumorSurvForestDS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
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
  ptvb_MAF104 = fread(paste0(ptv_dir, "MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv"))
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
  sel_gpcas = paste0("gpca", 1:gpca_nb)
  if(table(survdf$PTVb)["1"] != 0 & gender == "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, sex, age, all_of(sel_gpcas))
  } else if(table(survdf$PTVb)["1"] != 0 & gender != "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, age, all_of(sel_gpcas))
  } else if(table(survdf$PTVb)["1"] == 0 & gender == "all") {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, sex, age, all_of(sel_gpcas))
  } else {
    survdf %<>% dplyr::select(eid, surv_time, surv_event, age, all_of(sel_gpcas))
  }
  dead_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  alive_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 0, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH. Geneset size: ", length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, "Geneset size: ",  length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
  }
  myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
  res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
  fm = forest_model(model = res.cox, return_data = T)
  fm$plot + labs(title = tt)
}
