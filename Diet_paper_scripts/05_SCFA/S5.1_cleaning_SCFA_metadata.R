# libraries ====
library(tidyverse)
library(microViz)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
scfa_meta <- read.delim("Input/8_KEGG_codes_scfa.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning metadata ====
scfa_meta <- scfa_meta %>% mutate(
  `pathway_Pyruvate_metabolism` = ifelse(grepl("Pyruvate metabolism", pathway_name, ignore.case = TRUE), "yes", "no"),
  `pathway_Butanoate_metabolism` = ifelse(grepl("Butanoate metabolism", pathway_name, ignore.case = TRUE), "yes", "no"),
  `pathway_Propanoate_metabolism` = ifelse(grepl("Propanoate metabolism", pathway_name, ignore.case = TRUE), "yes", "no"),
  `scfa_propionate_production` = ifelse(grepl("propionate", scfa, ignore.case = TRUE), "yes", "no"),
  `scfa_butyrate_production` = ifelse(grepl("butyrate", scfa, ignore.case = TRUE), "yes", "no"),
  `scfa_acetate_production` = ifelse(grepl("acetate", scfa, ignore.case = TRUE), "yes", "no"))

# creating uniq kegg codes ====
# removing duplicates
duplicates <- scfa_meta %>% filter(duplicated(KEGG_code) | duplicated(KEGG_code, fromLast = TRUE)) 

# saving results ====
#write.table (scfa_meta, "S5.1_cleaned_scfa_meta.txt", row.names=FALSE,sep = "\t")


