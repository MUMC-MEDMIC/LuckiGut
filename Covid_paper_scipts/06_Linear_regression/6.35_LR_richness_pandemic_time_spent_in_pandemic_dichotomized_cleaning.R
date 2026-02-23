# Load packages ====
library(tidyverse)
library(ggplot2)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model <- read.delim("Output/Files/S6.34_LR_pandemic_rich_scaled_time_spent_in_pandemic_dichotomized_raw.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# rename variables ====
model$variables <- model$Term
plyr::count(model$Term)
model$variables[grepl("breastfeeding",model$variables)]<-"Breastfeeding (ref: no)"
model$variables[grepl("formula",model$variables)]<-"Formula feeding (ref: no)"
model$variables[grepl("inf_",model$variables)]<-"Antibiotics (ref: no)"
model$variables[grepl("probio",model$variables)]<-"Probiotics (ref: no)"
model$variables[grepl("furrypet",model$variables)]<-"Pets (ref: no)"
model$variables[grepl("daycare",model$variables)]<-"Daycare (ref: no)"
model$variables[grepl("Fur_pets_pregn_yn",model$variables)]<-"Pets pregnancy (ref: no)"
model$variables[grepl("mat_weightgain",model$variables)]<-"Maternal weight gain (kg)"
model$variables[grepl("pregn_weeks",model$variables)]<-"Gestational age (wks)"
model$variables[grepl("delivery_type",model$variables)]<-"Vaginal delivery (ref: C-section)"
model$variables[grepl("^hospitalyes$",model$variables)]<-"Hospitalisation (ref: no)"
model$variables[grepl("delivery_place",model$variables)]<-"Delivery: hospital (ref: home)"
model$variables[grepl("mat_AB_use_final_preg",model$variables)]<-"Antepartum antibiotics (ref: no)"
model$variables[grepl("mat_AB_use_final_delivery",model$variables)]<-"Intrapartum antibiotics (ref: no)"
model$variables[grepl("smoke_before_preg",model$variables)]<-"Prenatal smoking (ref: no)"
model$variables[grepl("atopy_father",model$variables)]<-"Father's atopy (ref: no)"
model$variables[grepl("atopy_mother",model$variables)]<-"Mother's atopy (ref: no)"
model$variables[grepl("days_groupbelow_or_equal_median",model$variables)]<-"Time spent in pandemic (ref: below median)"
model$variables[grepl("older_siblings",model$variables)]<-"Older siblings (ref: no)"
model$variables[grepl("twins.siblings_in_lucki",model$variables)]<-"Siblings in cohort (ref: no)"
model$variables[grepl("sex",model$variables)]<-"Male (ref. female)"
model$variables[grepl("ageintrosolids",model$variables)]<-"Complementary foods age (wks)"
model$variables[grepl("birthweight",model$variables)]<-"Birth weight (g)"

# add theme ====
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Breastfeeding (ref: no)"|.$variables == "Formula feeding (ref: no)"|.$variables == "Complementary foods age (wks)"|
      .$variables == "Probiotics (ref: no)",
    true = "Dietary", 
    false = NA))#creating new variables
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Mother's atopy (ref: no)"|.$variables == "Antepartum antibiotics (ref: no)"|.$variables == "Intrapartum antibiotics (ref: no)"|
      .$variables == "Father's atopy (ref: no)"|.$variables == "Maternal weight gain (kg)"|
      .$variables == "Prenatal smoking (ref: no)",
    true = "Parental", 
    false = Exposure_categories))
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Vaginal delivery (ref: C-section)"|.$variables == "Delivery: hospital (ref: home)"|.$variables == "Birth weight (g)"|
      .$variables == "Gestational age (wks)"|.$variables == "Male (ref. female)"|.$variables == "Hospitalisation (ref: no)",
    true = "Perinatal", 
    false = Exposure_categories))
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Antibiotics (ref: no)",
    true = "Infant's health", 
    false = Exposure_categories))
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Daycare (ref: no)"|.$variables == "Pets (ref: no)"|.$variables == "Pets pregnancy (ref: no)"|
      .$variables == "Older siblings (ref: no)"|.$variables == "Siblings in cohort (ref: no)"|.$variables == "Time spent in pandemic (ref: below median)",
    true = "Environmental", 
    false = Exposure_categories))

# finding final run for each TP
model <- model %>%
  mutate(run_num = parse_number(as.character(run))) %>%   # robust numeric extraction
  group_by(Age_individual) %>%                                      
  mutate(
    last_run_num = max(run_num, na.rm = TRUE),
    last_run = if_else(run_num == last_run_num, "yes", "no")
  ) %>%
  ungroup() %>%
  select(-run_num, -last_run_num)

# adjusting p values ====
model <- model %>%
  mutate(run_num = parse_number(as.character(run))) %>%   # robust numeric extraction
  group_by(Age_individual) %>%                                      
  mutate(
    last_run_num = max(run_num, na.rm = TRUE),
    last_run = if_else(run_num == last_run_num, "yes", "no")
  ) %>%
  ungroup() %>%
  select(-run_num, -last_run_num)

model <- model %>% mutate(p_adjusted = NA_real_)
idx <- model$last_run == "yes" #Identify the rows that correspond to the final run
model$p_adjusted[idx] <- p.adjust(model$p.value[idx], method = "fdr")

#add significance
#model <- model %>% 
#  mutate(sig_adj = if_else(
#    condition = .$p_adjusted > 0.05, 
#    true = NA, 
#    false = if_else(
#      condition = .$p_adjusted < 0.001,
#      true = "***",
#      false = if_else(
#        condition = .$p_adjusted <= 0.05&.$p_adjusted >= 0.01,
#        true = "*",
#        false = "**"))))

model <- model %>% 
  mutate(sig_padj = if_else(
    condition = .$p_adjusted > 0.1, 
    true = "", 
    false = if_else(
      condition = .$p_adjusted <= 0.1 &.$p_adjusted >= 0.05,
      true = "*",
      false = "**")))

# saving data ====
#write.table (model, "S6.35_LR_pandemic_rich_scaled_time_spent_in_pandemic_dichotomized_cleaned.txt", row.names=FALSE,sep = "\t")
