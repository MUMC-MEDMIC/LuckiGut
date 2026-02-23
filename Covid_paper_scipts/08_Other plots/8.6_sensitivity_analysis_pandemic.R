# load libraries ====
library(tidyverse)
library(phyloseq)
library(scales)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

# selecting TP with the most samples ====
meta <- as.data.frame(as.matrix(infant@sam_data))

meta_covid <- meta %>% filter(covid_270220 == "yes")

count <- meta_covid %>% group_by(Age_individual)
count <- plyr::count(meta_covid$Age_individual)


ages_needed <- c("1-2 weeks","4 weeks","8 weeks","4 months")

infant_selected <- infant %>% ps_filter(Age_individual %in% ages_needed)

