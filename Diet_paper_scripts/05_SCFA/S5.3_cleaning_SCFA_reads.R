# libraries ====
library(tidyverse)
library(microViz)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
baku <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
humann_output <- read.delim("Output/Files//S5.2_raw_scfa_reads.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# checking that reads per kegg code are sum of reads per baku ====
test <- humann_output %>%  filter(str_detect(X..Gene.Family, "K00074"))
rownames(test) <- test$X..Gene.Family
test <- as.data.frame(as.matrix(t(test)))
test <- test %>% filter(K00074 != "K00074")
test <- as.data.frame(lapply(test, as.numeric))
test$sum <- test %>%
  select(-K00074) %>%
  rowSums()

test <- test %>% select(K00074, sum)

# creating keys to correct sample names from metadata from phyloseq ====
keys <- as.data.frame(as.matrix(baku@sam_data))
keys <- keys %>% select(MMHP_SampleID)
keys$todestroy <- keys$MMHP_SampleID
keys <- keys %>% separate(todestroy, into = c("col1", "col2", "col3", "original"), sep = "_")
keys <- keys %>% select(MMHP_SampleID,original)
keys$original <- as.numeric(keys$original)

# cleaning names in humann output ====
wrong_sample_names <- as.data.frame(colnames(humann_output))
colnames(wrong_sample_names)[colnames(wrong_sample_names) == "colnames(humann_output)"] <- "decoy"
wrong_sample_names$todestroy <- wrong_sample_names$decoy
wrong_sample_names <- wrong_sample_names %>% separate(todestroy, into = c("col1", "col2", "original", "col3"), sep = "_")
wrong_sample_names <- wrong_sample_names %>% select(decoy,original)
wrong_sample_names$original <- as.numeric(wrong_sample_names$original)
wrong_sample_names <- wrong_sample_names[!is.na(wrong_sample_names$original), ]
str(wrong_sample_names)

# merging keys
keys <- full_join(keys,wrong_sample_names, by = "original")
rm(wrong_sample_names)

# changing names ====
rownames(humann_output) <- humann_output$ID
rev_filtered_humann <- as.data.frame(as.matrix(t(humann_output)))
rev_filtered_humann$decoy <- rownames(rev_filtered_humann)
rev_filtered_humann <- full_join(keys,rev_filtered_humann,by = "decoy")

rev_filtered_humann <- rev_filtered_humann %>% filter(!(decoy %in% c("ID", "MMHP_Lucki_0780_Abundance.RPKs"))) # removing duplicate sample
rev_filtered_humann$MMHP_SampleID <- ifelse(rev_filtered_humann$decoy == "X..Gene.Family", rev_filtered_humann$decoy, rev_filtered_humann$MMHP_SampleID)
rev_filtered_humann$original <- NULL
rev_filtered_humann$decoy <- NULL

rev_filtered_humann <- rev_filtered_humann[!is.na(rev_filtered_humann$MMHP_SampleID), ]
rownames(rev_filtered_humann) <- rev_filtered_humann$MMHP_SampleID
rev_filtered_humann$MMHP_SampleID <- NULL
renamed_scfa <- as.data.frame(as.matrix(t(rev_filtered_humann)))
length(unique(renamed_scfa$X..Gene.Family)) == nrow(renamed_scfa) #check that values are unique
rm(rev_filtered_humann,humann_output,keys)

# filtering on pathway level ====
renamed_scfa <- renamed_scfa %>% separate('X..Gene.Family', into = c("KEGG_code", "musor"), sep = "\\|")
renamed_scfa <- renamed_scfa %>% separate(musor, into = c("genus", "species"), sep = "\\.")
renamed_scfa$genus <- gsub("^g__", "", renamed_scfa$genus)
renamed_scfa$species <- gsub("^s__", "", renamed_scfa$species)
renamed_scfa$species[is.na(renamed_scfa$species)] <- renamed_scfa$genus[is.na(renamed_scfa$species)]

bakut <- renamed_scfa %>% filter(!is.na(species))
scfa <- renamed_scfa %>% filter(is.na(species))
scfa <- scfa %>% select(-genus,-species)

rm(renamed_scfa)

# save results ====
#write.table (scfa, "S5.3_cleaned_reads_scfa.txt", row.names=FALSE,sep = "\t")
#write.table (bakut, "S5.3_cleaned_reads_bakut.txt", row.names=FALSE,sep = "\t")
