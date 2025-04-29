Packages <- c("tidyverse", "data.table", "Rfast", "magrittr", "fastmatch", "survival", "tidycmprsk", "ggsurvfit", "cowplot", "gtsummary", "forestmodel", "circlize", "ComplexHeatmap")
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
  
  sel_gpcas = paste0("gpca", 1:gpca_nb)
  tempdf %<>% dplyr::select(eid, time, status, PTVb, sex, all_of(sel_gpcas))
  if(gender %in% c("female", "male")) {
    tempdf %<>% dplyr::select(-sex)
  }
  if(length(unique(tempdf$PTVb)) == 1) {
    tempdf %<>% dplyr::select(-PTVb)
  }
  if("PTVb" %in% colnames(tempdf)) {
    case_n_ptvb_1 = tempdf %>% filter(status == 1, PTVb == 1) %>% nrow()
  } else {
    case_n_ptvb_1 = 0
  }
  tt = paste0(paste0(tumor, collapse = ", "), ", ", paste0(gender, collapse = ", "), ". Geneset size: ", length(geneset), "\nTumorous patients (PTV=1): ", case_n_ptvb_1)
  myformula = as.formula(paste('Surv(time, status) ~ ', paste0(colnames(tempdf)[4:ncol(tempdf)], collapse = "+")))
  res.cox <- coxph(myformula, data = tempdf, id = tempdf$eid)
  fm = forest_model(model = res.cox, return_data = T)
  fm$plot + labs(title = tt)
}


#Perform a Cox model and return the P-values and the coefficients
fun_TumorIncCoxModel = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
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
  
  
  sel_gpcas = paste0("gpca", 1:gpca_nb)
  tempdf %<>% dplyr::select(eid, time, status, PTVb, sex, all_of(sel_gpcas))
  if(gender %in% c("female", "male")) {
    tempdf %<>% dplyr::select(-sex)
  }
  if(length(unique(tempdf$PTVb)) == 1) {
    tempdf %<>% dplyr::select(-PTVb)
  }
  if("PTVb" %in% colnames(tempdf)) {
    myformula = as.formula(paste('Surv(time, status) ~ ', paste0(colnames(tempdf)[4:ncol(tempdf)], collapse = "+")))
    res.cox <- coxph(myformula, data = tempdf, id = tempdf$eid)
    out = summary(res.cox)$coefficients[grep("PTV", rownames(summary(res.cox)$coefficients), value = T),c("exp(coef)","Pr(>|z|)")]
  } else {
    out = c(NA, NA)
  }
  names(out) = c("HR", "Pvalue")
  out
}

#Perform a Fisher's exact test and return the P-value and the odds ratio
fun_TumorIncFishTest = function(tumor = "all", control = "nottumorous", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
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
  c(Fishers_exact_test_OR = unname(fisher.test(table(tempdf$status, tempdf$PTVb))$estimate), Fishers_exact_test_P = fisher.test(table(tempdf$status, tempdf$PTVb))$p.value)
}





fun_TumorIncStat = function(tumor = "Melanoma", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
  tumor_data = bind_rows(lapply(tumor, function(tmr) {
    if(gender == "all") {
      tmr_data = readRDS(paste0(dataset_dir, tmr, ".rds"))
    } else {
      tmr_data = readRDS(paste0(dataset_dir, tmr, "_", gender, ".rds"))
    }
  }))
  tumor_data$death_type[tumor_data$death == T & is.na(tumor_data$death_type)] = "other"
  tumor_data %<>% mutate(cancerCausesDeath = death_type == cancer_type)
  m = max(tumor_data$date_of_death, na.rm = T)
  tumor_data %<>% 
    arrange(diag_date) %>% 
    distinct(eid, .keep_all = T) %>% 
    mutate(time = case_when(
      death == T ~ as.numeric(date_of_death - birthdate),
      death == F ~ as.numeric(m - birthdate))) %>% 
    mutate(status = case_when(
      death == T & cancerCausesDeath == T ~ 2,
      death == F ~ 1,
      death == T & cancerCausesDeath == F ~ 3
    )) %>% 
    transform(status = factor(status, levels = c(1,2,3)))
  ptvb_MAF104 = fread(paste0(ptv_dir, "MAF10_4_all_retained_variants_PTVBurden_final_Shetscores.tsv"))
  ptvb_MAF104 %<>% dplyr::select(Patient.ID, Genes) %>% set_colnames(c("eid", "ptvgenes")) %>% filter(eid %fin% tumor_data$eid)
  ptvb_MAF104$PTVb = sapply(ptvb_MAF104$ptvgenes, function(x) sum(unique(unlist(strsplit(x, ","), use.names = F)) %fin% geneset))
  tumor_data %<>% left_join(ptvb_MAF104, by = "eid") %>% dplyr::select(-ptvgenes)
  if(ptv_burden_cat == T) {tumor_data$PTVb = as.factor(ifelse(tumor_data$PTVb == 0, 0, 1))}
  cuminc(Surv(time, status) ~ PTVb, data = tumor_data) %>%
    ggcuminc() +
    labs(x = "Days") +
    add_confidence_interval() +
    add_risktable() +
    add_pvalue()
}




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
  if(length(unique(survdf$PTVb)) == 1) {
    print("No tumorous patients with geneset specific PTV")
  } else {
    dead_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
    alive_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 0, PTVb == 1) %>% nrow()
    tt = paste0(paste0(tumor, collapse = ", "), ", ", paste0(gender, collapse = ", "), ". Geneset size: ", length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1): ", alive_n_ptvb_1)
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
}


fun_TumorSurvCoxOS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
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
  survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, sex, age, all_of(sel_gpcas))
  if(gender %in% c("female", "male")) {
    survdf %<>% dplyr::select(-sex)
  }
  if(length(unique(survdf$PTVb)) == 1) {
    survdf %<>% dplyr::select(-PTVb)
  }
  if("PTVb" %in% colnames(survdf)) {
    myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
    res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
    out = summary(res.cox)$coefficients[grep("PTV", rownames(summary(res.cox)$coefficients), value = T),c("exp(coef)","Pr(>|z|)")]
  } else {
    out = c(NA, NA)
  }
  names(out) = c("HR", "Pvalue")
  out
}


fun_TumorSurvForestOS = function(tumor = "all", gender = "all", geneset = c("ENSG00000104804"), ptv_burden_cat = TRUE, gpca_nb = 10, threads, dataset_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/Objects/000_Sub_datasets/", ptv_dir = "/media/balazs/WorkL/balazs/Work/UKBiobank/PTVvars/") {
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
  survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, sex, age, all_of(sel_gpcas))
  if(gender %in% c("female", "male")) {
    survdf %<>% dplyr::select(-sex)
  }
  if(length(unique(survdf$PTVb)) == 1) {
    survdf %<>% dplyr::select(-PTVb)
  }
  if("PTVb" %in% colnames(survdf)) {
    dead_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
    alive_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 0, PTVb == 1) %>% nrow()
  } else {
    dead_n_ptvb_1 = 0
    alive_n_ptvb_1 = 0
  }
  tt = paste0(paste0(tumor, collapse = ", "), ", ", paste0(gender, collapse = ", "), ". Geneset size: ", length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
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
  
  if(length(unique(survdf$PTVb)) == 1) {
    print("No tumorous patients with geneset specific PTV")
  } else {
    dead_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
    alive_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 0, PTVb == 1) %>% nrow()
    tt = paste0(paste0(tumor, collapse = ", "), ", ", paste0(gender, collapse = ", "), ". Geneset size: ", length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1): ", alive_n_ptvb_1)
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
  survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, sex, age, all_of(sel_gpcas))
  if(gender %in% c("female", "male")) {
    survdf %<>% dplyr::select(-sex)
  }
  if(length(unique(survdf$PTVb)) == 1) {
    survdf %<>% dplyr::select(-PTVb)
  }
  if("PTVb" %in% colnames(survdf)) {
    myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
    res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
    out = summary(res.cox)$coefficients[grep("PTV", rownames(summary(res.cox)$coefficients), value = T),c("exp(coef)","Pr(>|z|)")]
  } else {
    out = c(NA, NA)
  }
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
  survdf %<>% dplyr::select(eid, surv_time, surv_event, PTVb, sex, age, all_of(sel_gpcas))
  if(gender %in% c("female", "male")) {
    survdf %<>% dplyr::select(-sex)
  }
  if(length(unique(survdf$PTVb)) == 1) {
    survdf %<>% dplyr::select(-PTVb)
  }
  if("PTVb" %in% colnames(survdf)) {
    dead_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 1, PTVb == 1) %>% nrow()
    alive_n_ptvb_1 = survdf %>% dplyr::filter(surv_event == 0, PTVb == 1) %>% nrow()
  } else {
    dead_n_ptvb_1 = 0
    alive_n_ptvb_1 = 0
  }
  tt = paste0(paste0(tumor, collapse = ", "), ", ", paste0(gender, collapse = ", "), ". Geneset size: ", length(geneset), "\nDead patients (PTV=1): ", dead_n_ptvb_1, ". Alive patients (PTV=1):", alive_n_ptvb_1)
  myformula = as.formula(paste('Surv(surv_time, surv_event) ~ ', paste0(colnames(survdf)[4:ncol(survdf)], collapse = "+")))
  res.cox = survival::coxph(myformula, data = as.data.frame(survdf))
  fm = forest_model(model = res.cox, return_data = T)
  fm$plot + labs(title = tt)
}

#Draw heatmap with HR and P-values
fun_hm = function(restbl, x_category = "TS", y_category = "ptv_type", color_var = "HR", text_var = "Pvalue", text_var_cutoff = 0.1, col_clustering = T, row_clustering = T) {
  colorDataHM = restbl %>% select(all_of(x_category), all_of(y_category), all_of(color_var)) %>% pivot_wider(id_cols = all_of(y_category), names_from = all_of(x_category), values_from = all_of(color_var))
  rn = colorDataHM %>% pull(all_of(y_category))
  colorDataHM = as.matrix(colorDataHM[,2:ncol(colorDataHM)])
  rownames(colorDataHM) = rn
  textDataHM = restbl %>% select(all_of(x_category), all_of(y_category), all_of(text_var)) %>% pivot_wider(id_cols = all_of(y_category), names_from = all_of(x_category), values_from = all_of(text_var))
  rn = textDataHM %>% pull(all_of(y_category))
  textDataHM = as.matrix(textDataHM[,2:ncol(textDataHM)])
  rownames(textDataHM) = rn
  textDataHM[textDataHM > text_var_cutoff] = NA
  Heatmap(matrix = colorDataHM, 
          row_names_max_width = max_text_width(rn),
          column_names_max_height = max_text_width(colnames(colorDataHM)),
          row_gap = unit(2, "mm"),
          name = "HR", 
          na_col = "gray",
          col = colorRamp2(c(0, 1, max(colorDataHM, na.rm = T)), c("green", "white", "red")),
          cluster_columns = col_clustering, 
          cluster_rows = row_clustering, 
          rect_gp = gpar(col = "black", lwd = 1), 
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(!is.na(textDataHM[i, j]))
              grid.text(sprintf("%.3f", textDataHM[i, j]), x, y, gp = gpar(fontsize = 8))
          })
}


# res_fishtest = readRDS("Res/010/res_fishtest_104.rds")
# genes = res_fishtest %>% dplyr::filter(category == "Melanoma", caseY > 5) %>% arrange(desc(OR)) %>% dplyr::slice(1:20) %>% pull(ENSEMBL)
# 
# 
# load("Objects/go_desc_list")
# hugos = go_desc_list[["antigen processing and presentation"]]
# library(clusterProfiler)
# library(org.Hs.eg.db)
# genes = unique(bitr(geneID = hugos, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")$ENSEMBL)
# genes = genes[!is.na(genes)]
# rm(go_desc_list, hugos)
# 
# broca_geneset = c("ALK", "APC", "ATM", "ATR", "AXIN2", "BAP1", "BARD1", "BMPR1A", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDK4", "CDK12", "CDKN2A", "CHEK2", "CTNNA1", "DICER1", "EPCAM", "FANCM", "FH", "FLCN", "GEN1", "GREM1", "HOXB13", "MEN1", "MET", "MITF", "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "MUTYH", "NBN", "NF1", "NF2", "NTHL1", "PALB2", "PHOX2B", "PIK3CA", "PMS2", "POLD1", "POLE", "POT1", "PRKAR1A", "PTCH1", "PTEN", "RAD51B", "RAD51C", "RAD51D", "RB1", "RECQL", "RET", "RNF43", "RPS20", "SDHB", "SDHC", "SDHD", "SMAD4", "SMARCA4", "TP53", "TSC1", "TSC2", "VHL")
# genes = unique(bitr(geneID = broca_geneset, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")$ENSEMBL)
# 
# fun_TumorIncForestPlot(tumor = "Melanoma", gender = "all", geneset = genes)
# fun_TumorIncCoxModel(tumor = "Melanoma", gender = "all", geneset = genes)
# fun_TumorIncFishTest(tumor = "Melanoma", gender = "all", geneset = genes)
# fun_TumorIncPlot(tumor = "Melanoma", gender = "all", geneset = genes)
# 
# 
# fun_TumorSurvPlotDS(tumor = "Melanoma", gender = "female", geneset = genes)
# fun_TumorSurvPlotOS(tumor = "Melanoma", gender = "female", geneset = genes)
# #Overall survival plot function
# 
# #Disease specific survival plot function
# 
