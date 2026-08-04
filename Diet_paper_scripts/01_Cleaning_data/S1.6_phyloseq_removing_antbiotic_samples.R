# load in packages ====
library(tidyverse) 
library(phyloseq)
library(microViz)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# loading in MS phyloseq ====
ps_diet <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")

#cleaning data
ps_diet <- ps_diet  %>%
  ps_mutate(
    timepoint_number = factor(str_remove(ps_diet@sam_data$timepoint,"m"), levels = c("4","5","6","9","11","14")),
  )

# filtering antibiotic infants out ====
test <- as.data.frame(as.matrix(ps_diet@sam_data))

ps_diet_filtered <- ps_filter(
  ps_diet,
  !((inf_ab_ED_4m  == "yes" & Age_individual == "4 months")  |
      (inf_ab_ED_5m  == "yes" & Age_individual == "5 months")  |
      (inf_ab_ED_6m  == "yes" & Age_individual == "6 months")  |
      (inf_ab_ED_9m  == "yes" & Age_individual == "9 months")  |
      (inf_ab_ED_14m == "yes" & Age_individual == "14 months")))

# saving phyloseq ====
#saveRDS(ps_diet_filtered, 'S1.6_cleaned_phyloseq_object_MS_without_antibiotics.rds')