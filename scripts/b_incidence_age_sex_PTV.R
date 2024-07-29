library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)

load("objects/ukb_data_cancer_final")
ukb_data_cancer$cancer_type[is.na(ukb_data_cancer$cancer_type)] = "No_cancer"
tumors = unique(ukb_data_cancer$cancer_type)
tumors = tumors[tumors != "No_cancer"]

age_inc_plots = lapply(X = tumors, FUN = function(x){
  tempdf_x = ukb_data_cancer %>% filter(cancer_type == x) %>% arrange(diag_date)
  tempdf_x = tempdf_x[!duplicated(tempdf_x$eid), ]
  tempdf_x %<>% mutate(inc_time = round(as.numeric((diag_date - birthdate))/365.25),
                       final_age = case_when(
                         death == T ~ round(as.numeric(date_of_death - birthdate)/365.25),
                         death == F ~ round(as.numeric(max(date_of_death, na.rm = T) - birthdate)/365.25)))
  
  tempdf_no_x = ukb_data_cancer %>% filter(cancer_type == "No_cancer")
  tempdf_no_x = tempdf_no_x[!duplicated(tempdf_no_x$eid), ]
  tempdf_no_x %<>% mutate(inc_time = NA, 
                          final_age = case_when(
                            death == T ~ round(as.numeric(date_of_death - birthdate)/365.25),
                            death == F ~ round(as.numeric(max(date_of_death, na.rm = T) - birthdate)/365.25)))
  
  minage = min(c(tempdf_x$inc_time, tempdf_no_x$final_age))
  maxage = max(c(tempdf_x$inc_time, tempdf_no_x$final_age))
  res = expand.grid(inc_time = seq(minage, maxage, 1), sex = unique(ukb_data_cancer$sex), PTV_universal_1 = unique(tempdf_x$PTV_universal_1), stringsAsFactors = F) 
  res$PTV_universal_1 = factor(res$PTV_universal_1, levels = levels(tempdf_x$PTV_universal_1))

  res_case = tempdf_x %>% group_by(inc_time, sex, PTV_universal_1) %>% summarise(case = n())
  res %<>% left_join(res_case)
  res$case[is.na(res$case)] = 0
  
  res_nocase = tempdf_no_x %>% group_by(final_age, sex, PTV_universal_1) %>% summarise(no_case = n())
  
  res$no_case = 0
  for(i in 1:nrow(res_nocase)) {
    res$no_case[res$inc_time <= res_nocase$final_age[i] & res$sex == res_nocase$sex[i] & res$PTV_universal_1 == res_nocase$PTV_universal_1[i]] = res$no_case[res$inc_time <= res_nocase$final_age[i] & res$sex == res_nocase$sex[i] & res$PTV_universal_1 == res_nocase$PTV_universal_1[i]] + res_nocase$no_case[i]
  }
  
  res %<>% mutate(inc = case/(case+no_case))
  
  # ggplot(res) + 
  #   geom_point(mapping = aes(x = inc_time, y = inc, color = PTV_universal_1, shape = as.factor(sex))) + 
  #   geom_line(aes(x = inc_time, y = inc, color = PTV_universal_1, shape = as.factor(sex))) +
  #   facet_grid(rows = vars(sex), cols = vars(PTV_universal_1))
  
  ggplot(res, mapping = aes(x = inc_time, y = inc, color = PTV_universal_1, shape = as.factor(sex))) + 
    geom_point() + 
    geom_smooth() +
    facet_wrap(~sex) + 
    labs(title = x)
})

ggarrange(plotlist = age_inc_plots, ncol = 5, nrow = 5)
ggsave(filename = "plots/age_inc_plot_sex_PTV_universal_1.jpg", width = 100, height = 60, units = "cm", dpi = 300)

