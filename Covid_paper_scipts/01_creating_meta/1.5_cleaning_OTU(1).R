# load libraries ====
#cleaning read count table for phyloseq
library(tidyverse)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
# read in the file with metagenomic reads after this script: combine_taxa_reads_into_one_after_host_removal.bash
# replace | with ; in original taxonomy file in taxa column
# remove empty column and rows in excel
# replace taxa with MMHP_SampleID in A1
# Replace all NA with 0 in excel df[df == 0] <- NA
# replace "-" in sample names with "_"

all_otu = read.delim("Input/16_OTU_metagenomic_reads.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# cleaning data ====
# removing duplicate sample MMHP_UM_LCK_780 -> N = 930
all_otu <- all_otu %>%  select(!(MMHP_UM_LCK_780))

#rename
names(all_otu)[names(all_otu) == "SampleID"] <- "MMHP_SampleID"

# sub-setting to a microbe group
baku <- all_otu %>% filter(!grepl('k__Viruses', MMHP_SampleID)) %>% filter(!grepl('k__Eukaryota', MMHP_SampleID))
virus <- all_otu %>% filter(grepl('k__Viruses', MMHP_SampleID))
eukaryot <- all_otu %>% filter(grepl('k__Eukaryota', MMHP_SampleID))

# checking if any sample has 0 total reads, and if removing them from otu table
rownames(all_otu) <- all_otu$MMHP_SampleID
all_otu <- all_otu %>% select(-c(MMHP_SampleID))
trans_tot_0 = t(all_otu)
trans_tot_0 <- as.data.frame(cbind(trans_tot_0, total_read_count_total = rowSums(trans_tot_0)))
trans_tot <- trans_tot_0[trans_tot_0$total_read_count_total != '0',]
all_otu = as.data.frame(t(trans_tot))
#collecting counts
counts_total <- trans_tot %>% select(total_read_count_total)
#removing total count of reads = last row from the otu table
all_otu <- all_otu[row.names(all_otu) != "total_read_count_total", , drop = FALSE]
rm(trans_tot_0,trans_tot)

# checking if any sample has 0 total eukaryotic reads, and if removing them from otu table
rownames(eukaryot) <- eukaryot$MMHP_SampleID
eukaryot <- eukaryot %>% select(-c(MMHP_SampleID))
trans_euk_0 = t(eukaryot)
trans_euk_0 <- as.data.frame(cbind(trans_euk_0, total_read_count_eukar = rowSums(trans_euk_0)))
trans_euk <- trans_euk_0[trans_euk_0$total_read_count_eukar != '0',]
eukaryot = as.data.frame(t(trans_euk))
#collecting counts
counts_euk <- trans_euk %>% select(total_read_count_eukar)
#removing total count of reads = last row from the otu table
eukaryot <- eukaryot[row.names(eukaryot) != "total_read_count_eukar", , drop = FALSE]
rm(trans_euk_0,trans_euk)

# checking if any sample has 0 total viral reads, and if removing them from otu table
rownames(virus) <- virus$MMHP_SampleID
virus <- virus %>% select(-c(MMHP_SampleID))
trans_vir_0 = t(virus)
trans_vir_0 <- as.data.frame(cbind(trans_vir_0, total_read_count_virus = rowSums(trans_vir_0)))
trans_vir <- trans_vir_0[trans_vir_0$total_read_count_virus != '0',]
virus = as.data.frame(t(trans_vir))
#collecting counts
counts_virus <- trans_vir %>% select(total_read_count_virus) 
#removing total count of reads = last row from the otu table
virus <- virus[row.names(virus) != "total_read_count_virus", , drop = FALSE]
rm(trans_vir_0,trans_vir)

# checking if any sample has 0 total bacterial reads, and if removing them from otu table
rownames(baku) <- baku$MMHP_SampleID
baku <- baku %>% select(-c(MMHP_SampleID))
trans_baku_0 = t(baku)
trans_baku_0 <- as.data.frame(cbind(trans_baku_0, total_read_count_baku = rowSums(trans_baku_0)))
trans_baku <- trans_baku_0[trans_baku_0$total_read_count_baku != '0',]
baku = as.data.frame(t(trans_baku))
#collecting counts
counts_baku <- trans_baku %>% select(total_read_count_baku)
#removing total count of reads = last row from the otu table
baku <- baku[row.names(baku) != "total_read_count_baku", , drop = FALSE]
rm(trans_baku_0,trans_baku)

#merging counts into one file
counts_total <- tibble::rownames_to_column (counts_total, var = "MMHP_SampleID")
counts_baku <- tibble::rownames_to_column (counts_baku, var = "MMHP_SampleID")
counts_virus <- tibble::rownames_to_column (counts_virus, var = "MMHP_SampleID")
counts_euk <- tibble::rownames_to_column (counts_euk, var = "MMHP_SampleID")
count_1 <- full_join(counts_total, counts_baku, by = "MMHP_SampleID")
count_2 <- full_join(counts_virus, counts_euk, by = "MMHP_SampleID")
count <- full_join(count_1, count_2, by = "MMHP_SampleID")
rm(counts_total,counts_baku,counts_virus,counts_euk,count_1,count_2)

# write data ====
#write.table (count, "S1.5_read_counts.txt", row.names=FALSE,sep = "\t")
#write.table (all_otu, "S1.5_OTU_all.txt", row.names=TRUE,sep = "\t")
#write.table (baku, "S1.5_OTU_bacteria.txt", row.names=TRUE,sep = "\t")
#write.table (virus, "S1.5_OTU_virus.txt", row.names=TRUE,sep = "\t")
#write.table (eukaryot, "S1.5_OTU_eukaruota.txt", row.names=TRUE,sep = "\t")

