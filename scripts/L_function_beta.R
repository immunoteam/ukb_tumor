

#Overall survival

fun_TumoroaSurvPlot = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads) {
  #if(c("PTVgenes") %in% colnames(dataset))
  #plan(future::cluster, workers = threads)
  #if(control == "nottumorous" & gender == "all") {
  #  control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
  #} else if(control == "nottumorous" & gender == "female") {
  #  control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_female.rds")
  #} else if(control == "nottumorous" & gender == "male") {
  #  control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_male.rds")
  #}
  tumor_data = lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = load("C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/tumor_dfs/tumorous_all")
    } else if(tmr == "all" & gender != "all") {
      tmr_data = load(paste0("C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/tumor_dfs/tumorous_", gender))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = load(paste0("C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/tumor_dfs/", tmr))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = load(paste0("C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/tumor_dfs/", gender, "_", tmr))
    }
  })
  tempdf_t = do.call(bind_rows, tumor_data)
  #tempdf = rbind(tumor_data, control_data)
  #rm(tumor_data, control_data)
  ptvb_MAF104 = fread("C:/Users/bagil/Desktop/MyGit/ukb_tumor/rawdata/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  library(Rfast)
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tempdf_t$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tempdf_t %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tempdf_t$PTVb = as.factor(ifelse(tempdf_t$PTVb == 0, 0, 1))}
  case_n_ptvb_1 = tempdf_t %>% filter(status == 2, PTVb == 1) %>% nrow()
  tumorous_d = tempdf_t %>%
    filter(death == T) %>%
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(surv_time = as.numeric(date_of_death - diag_date))
  tumorous_alive = tempdf_t %>% 
    filter(death == F) %>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(surv_time = as.numeric(max(date_of_death) - diag_date))
  survdf = bind_rows(tumorous_d, tumorous_alive)
  survdf$surv_event = as.numeric(ifelse(survdf$death == T, 1, 0))
  survplot = survfit2(Surv(surv_time, surv_event) ~ PTVb, data = survdf) %>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Overall survival probability"
    ) + 
    add_confidence_interval() +
    add_risktable()
}

#Disease specific
fun_TumordsSurvPlot = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads) {
  #if(c("PTVgenes") %in% colnames(dataset))
  #plan(future::cluster, workers = threads)
  #if(control == "nottumorous" & gender == "all") {
  #  control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_all.rds")
  #} else if(control == "nottumorous" & gender == "female") {
  #  control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_female.rds")
  #} else if(control == "nottumorous" & gender == "male") {
  #  control_data = readRDS("/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/nottumorous_male.rds")
  #}
  tumor_data = lapply(tumor, function(tmr) {
    if(tmr == "all" & gender == "all") {
      tmr_data = load("C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/tumor_dfs/tumorous_all")
    } else if(tmr == "all" & gender != "all") {
      tmr_data = load(paste0("C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/tumor_dfs/tumorous_", gender))
    } else if(tmr != "all" & gender == "all") {
      tmr_data = load(paste0("C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/tumor_dfs/", tmr))
    } else if(tmr != "all" & gender != "all") {
      tmr_data = load(paste0("C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/tumor_dfs/", gender, "_", tmr))
    }
  })
  tempdf_t = do.call(bind_rows, tumor_data)
  #tempdf = rbind(tumor_data, control_data)
  #rm(tumor_data, control_data)
  ptvb_MAF104 = fread("C:/Users/bagil/Desktop/MyGit/ukb_tumor/rawdata/MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv")
  library(Rfast)
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tempdf_t$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tempdf_t %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tempdf_t$PTVb = as.factor(ifelse(tempdf_t$PTVb == 0, 0, 1))}
  #case_n_ptvb_1 = tempdf_t %>% filter(status == 2, PTVb == 1) %>% nrow()
  tumorous_d = tempdf_t %>%
    filter(death == T & death_type == cancer_type) %>%
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(surv_time = as.numeric(date_of_death - diag_date))
  tumorous_alive = tempdf_t %>% 
    filter(death == F) %>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(surv_time = as.numeric(max(date_of_death) - diag_date))
  survdf = bind_rows(tumorous_d, tumorous_alive)
  survdf$surv_event = as.numeric(ifelse(survdf$death == T, 1, 0))
  survplot = survfit2(Surv(surv_time, surv_event) ~ PTVb, data = survdf) %>% 
    ggsurvfit() +
    labs(
      x = "Days",
      y = "Disease specific survival probability"
    ) + 
    add_confidence_interval() +
    add_pvalue(location = c("annotation")) +
    add_risktable()
}

#cox model for Overall survival
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
  survdf = survdf %>% 
    select(eid, surv_time, surv_event, PTVb, all_of(sel_gpcas))
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

##Forest plot
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
  survdf = survdf %>% 
    select(eid, surv_time, surv_event, PTVb, all_of(sel_gpcas))
  case_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
  res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
  fm = forest_model(model = res.cox, return_data = T)
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  fm$plot + labs(title = tt)
}

#cox model for Disease specific survival
fun_TumorSurvCoxDS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/000_Sub_dataset/", ptv_dir = "C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/PTVvars/") {
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
    select(eid, surv_time, surv_event, PTVb, all_of(sel_gpcas))
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

##Forest plot
fun_TumorSurvForestDS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/000_Sub_dataset/", ptv_dir = "C:/Users/bagil/Desktop/MyGit/ukb_tumor/objects/PTVvars/") {
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
    select(eid, surv_time, surv_event, PTVb, all_of(sel_gpcas))
  case_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
  res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
  fm = forest_model(model = res.cox, return_data = T)
  if(gender == "all") {
    tt = paste0(paste0(tumor, collapse = ", "), ", BOTH.\nPatients with PTV in case group: ", case_n_ptvb_1, " Geneset size: ", length(geneset))
  } else {
    tt = paste0(paste0(tumor, collapse = ", "), ", ", gender, ".\nPatients with PTV in case group: ", case_n_ptvb_1, ". Geneset size: ", length(geneset))
  }
  fm$plot + labs(title = tt)
}
