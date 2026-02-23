# Load packages ====
library(tidyverse)
library(ggplot2)
library(patchwork)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# Import data ====
model <- read.delim("Output/Files/S5.1_PERMANOVA_results_raw.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# Cleaning variables ====
#rename variables
model$wrok_var <- model$variables
model$variables[grepl("breastfeeding",model$variables)]<-"Breastfeeding"
model$variables[grepl("formula",model$variables)]<-"Formula feeding"
model$variables[grepl("inf_",model$variables)]<-"Antibiotics"
model$variables[grepl("probio",model$variables)]<-"Probiotics"
model$variables[grepl("furrypet",model$variables)]<-"Pets"
model$variables[grepl("daycare",model$variables)]<-"Daycare"
model$variables[grepl("Fur_pets_pregn_yn",model$variables)]<-"Pets pregnancy"
model$variables[grepl("mat_weightgain",model$variables)]<-"Maternal weight gain"
model$variables[grepl("pregn_weeks",model$variables)]<-"Gestational age"
model$variables[grepl("delivery_type",model$variables)]<-"Delivery type"
model$variables[grepl("delivery_place",model$variables)]<-"Delivery place"
model$variables[grepl("mat_AB_use_final_preg",model$variables)]<-"Antepartum antibiotics"
model$variables[grepl("mat_AB_use_final_delivery",model$variables)]<-"Intrapartum antibiotics"
model$variables[grepl("smoke_before_preg",model$variables)]<-"Prenatal smoking"
model$variables[grepl("atopy_father",model$variables)]<-"Father's atopy"
model$variables[grepl("atopy_mother",model$variables)]<-"Mother's atopy"
model$variables[grepl("covid_270220",model$variables)]<-"Collected during pandemic"
model$variables[grepl("older_siblings",model$variables)]<-"Older siblings"
model$variables[grepl("twins.siblings_in_lucki",model$variables)]<-"Siblings in cohort"
model$variables[grepl("sex",model$variables)]<-"Sex"
model$variables[grepl("hospital",model$variables)]<-"Hospitalisation"
model$variables[grepl("ageintrosolids",model$variables)]<-"Complementary foods age"
model$variables[grepl("birthweight",model$variables)]<-"Birth weight"

model$age<- factor(model$age, levels = c("1 week",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"))

# add exposure_categories ====
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Breastfeeding"|.$variables == "Formula feeding"|.$variables == "Complementary foods age"|
      .$variables == "Probiotics",
    true = "Dietary", 
    false = NA))#creating new variables
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Mother's atopy"|.$variables == "Antepartum antibiotics"|.$variables == "Intrapartum antibiotics"|
      .$variables == "Father's atopy"|.$variables == "Maternal weight gain"|
      .$variables == "Prenatal smoking",
    true = "Parental", 
    false = Exposure_categories))
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Delivery type"|.$variables == "Delivery place"|.$variables == "Birth weight"|
      .$variables == "Gestational age"|.$variables == "Sex"|.$variables == "Hospitalisation",
    true = "Perinatal", 
    false = Exposure_categories))
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Antibiotics",
    true = "Infant's health", 
    false = Exposure_categories))
model <- model %>%
  mutate(Exposure_categories = if_else(
    condition = .$variables == "Daycare"|.$variables == "Pets"|.$variables == "Pets pregnancy"|
      .$variables == "Older siblings"|.$variables == "Siblings in cohort"|.$variables == "Collected during pandemic",
    true = "Environmental", 
    false = Exposure_categories))

model$variables <- factor(model$variables, 
                          levels = c("Complementary foods age",
                                     "Breastfeeding",
                                     "Formula feeding",
                                     "Probiotics",
                                     "Daycare",
                                     "Pets pregnancy",
                                     "Older siblings",
                                     "Pets",
                                     "Siblings in cohort",
                                     "Antibiotics",
                                     "Collected during pandemic",
                                     "Mother's atopy",
                                     "Father's atopy",
                                     "Antepartum antibiotics",
                                     "Intrapartum antibiotics",
                                     "Maternal weight gain",
                                     "Prenatal smoking",
                                     "Birth weight",
                                     "Delivery type",
                                     "Delivery place",
                                     "Hospitalisation",
                                     "Gestational age",
                                     "Sex"))
model$run <- as.numeric(model$run)

# adjusting p values ====
# finding final run for each TP
model <- model %>%
  mutate(run_num = parse_number(as.character(run))) %>%   # robust numeric extraction
  group_by(age) %>%                                      
  mutate(
    last_run_num = max(run_num, na.rm = TRUE),
    last_run = if_else(run_num == last_run_num, "yes", "no")
  ) %>%
  ungroup() %>%
  select(-run_num, -last_run_num)

model <- model %>% mutate(p_adjusted = NA_real_)
idx <- model$last_run == "yes" #Identify the rows that correspond to the final run
model$p_adjusted[idx] <- p.adjust(model$`Pr..F.`[idx], method = "fdr")

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
#write.table (model, "S5.2_PERMANOVA_results_cleaned.txt", row.names=FALSE,sep = "\t")









