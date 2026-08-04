# libraries ====
library(tidyverse)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
ps_diet <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
freezer <- read.delim("Input/7_Diepvrieslocaties_LucKi Gut_14-01_2025.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# selecting samples codes from phyloseq ====
meta <- as.data.frame(as.matrix(ps_diet@sam_data))
meta <- meta %>% filter(Age_individual == "6 months")
families_diet_6mon <- unique(sort(meta$Family_ID))

rm(meta,ps_diet)

# freezer data cleaning ====
#filtering for 6 months
freezer_6mon <- freezer %>% filter(Time.point == "6 months")
#selecting needed columns
freezer_6mon <- freezer_6mon %>% select(all_of(c("SampleID","Sample.nr", "Sample.weight..mg.", "Box", "Place","Comments.and.purpose","CB.child.code")))

#renaming columns
names(freezer_6mon)[names(freezer_6mon) == "Sample.nr"] <- "aliquot"
names(freezer_6mon)[names(freezer_6mon) == "Sample.weight..mg."] <- "original_weight"
names(freezer_6mon)[names(freezer_6mon) == "Comments.and.purpose"] <- "comments_upon_1st_freezing"
names(freezer_6mon)[names(freezer_6mon) == "CB.child.code"] <- "comments_upon_1st_freezing2"
freezer_6mon$aliquotID <- paste(freezer_6mon$SampleID, freezer_6mon$aliquot, sep = "_")

#taking family codes out 
freezer_6mon$backup <- freezer_6mon$SampleID
freezer_6mon <- freezer_6mon %>% separate("backup",into = c("Family_ID", "musor1", "musor2"),sep = "\\.")
freezer_6mon <- freezer_6mon %>% select(!(all_of(c("musor1", "musor2"))))

#taking all family codes
families_all <- unique(sort(freezer_6mon$Family_ID))

#filtering into diet and non diet groups
freezer_6mon <- freezer_6mon %>% mutate(diet = ifelse(Family_ID %in% families_diet_6mon, "yes", "no"))

rm(freezer)

# taking samples with no comments ====
freezer_6mon_present <- freezer_6mon %>% filter(is.na(comments_upon_1st_freezing))

# counting aliquotes ====
freezer_6mon_present <- freezer_6mon_present %>%
  group_by(Family_ID) %>%
  mutate(aliquot_n_theory = n()) %>%
  ungroup()

# saving results ====
#write.table (freezer_6mon_present, "S1.4_6mon_scfa_freezer_preselection.txt", row.names=FALSE,sep = "\t")
