# libraries ====
library(tidyverse)
library(microViz)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
baku <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
humann_output <- read.delim("Output/Files/S5.3_cleaned_reads_scfa.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
scfa_meta <- read.delim("Output/Files/S5.1_cleaned_scfa_meta.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# dividing by SCFA df ====
ace <- scfa_meta %>% filter(scfa == "acetate")
ace_kegg <- ace$KEGG_code
ace_kegg <- unique(ace_kegg)

but <- scfa_meta %>% filter(scfa == "butyrate")
but_kegg <- but$KEGG_code
but_kegg <- unique(but_kegg)

pro <- scfa_meta %>% filter(scfa == "propionate")
pro_kegg <- pro$KEGG_code
pro_kegg <- unique(pro_kegg)

# adding scfa meta, removing NAs ====
ace_reads <- full_join(ace, humann_output, by = "KEGG_code")
but_reads <- full_join(but, humann_output, by = "KEGG_code")
pro_reads <- full_join(pro, humann_output, by = "KEGG_code")
#check 
check <- ace_reads$KEGG_code
ace_reads_no_na <- ace_reads %>% na.omit(MMHP_UM_LCK_1000)
but_reads_no_na <- but_reads %>% na.omit(MMHP_UM_LCK_1000)
pro_reads_no_na <- pro_reads %>% na.omit(MMHP_UM_LCK_1000)

# turning to long format ====
ace_long <- ace_reads_no_na %>%
  pivot_longer(
    cols = starts_with("MMHP"), # Columns to pivot
    names_to = "MMHP_SampleID",      # Name for the new variable column
    values_to = "read_count"         # Name for the new value column
  )

but_long <- but_reads_no_na %>%
  pivot_longer(
    cols = starts_with("MMHP"), 
    names_to = "MMHP_SampleID",      
    values_to = "read_count" 
  )

pro_long <- pro_reads_no_na %>%
  pivot_longer(
    cols = starts_with("MMHP"), 
    names_to = "MMHP_SampleID",      
    values_to = "read_count" 
  )

rm(ace,ace_reads,ace_reads_no_na,but,but_reads,but_reads_no_na,pro,pro_reads,pro_reads_no_na)

all_scfa <- rbind(pro_long,ace_long,but_long)
rm(pro_long,ace_long,but_long)

# adding meta to scfa output ====
meta <- as.data.frame(as.matrix(baku@sam_data))
meta <- meta %>% select(MMHP_SampleID,starts_with("diet"),Age_individual,ageintrosolids,sex,Family_ID)

long_scfa_meta <- full_join(meta,all_scfa, by = "MMHP_SampleID",multiple = "all")
long_scfa_meta$diet_class <- gsub(" ", "", long_scfa_meta$diet_class)
long_scfa_meta$diet_class <- as.character(long_scfa_meta$diet_class)
long_scfa_meta$read_count <- as.numeric(long_scfa_meta$read_count)
long_scfa_meta$read_count_scaled <- scale(long_scfa_meta$read_count)

long_scfa_meta$Age_individual<- factor(long_scfa_meta$Age_individual, 
                                  levels = c("1-2 weeks",
                                             "4 weeks",
                                             "8 weeks",
                                             "4 months",
                                             "5 months",
                                             "6 months",
                                             "9 months",
                                             "11 months",
                                             "14 months"))

long_scfa_meta$diet_class <- factor(long_scfa_meta$diet_class, 
                                  levels = c("1",
                                             "2",
                                             "3"))

ages <- c("4 months","5 months","6 months","9 months","11 months","14 months")
scfa_meta_tp  <- long_scfa_meta %>% filter(Age_individual %in% ages)

# save results ====
#write.table (scfa_meta_tp, "S5.4_long_scfa_diet_paper.txt", row.names=FALSE,sep = "\t")
