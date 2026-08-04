# libraries ====
library(tidyverse)
library(microViz)
library(phyloseq)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in data ====
humann_output <- read.delim("Input/9_kegg-orthologs.tsv",sep = "\t",header = TRUE, na.strings=c("","NA"))
scfa_meta <- read.delim("Output/Files/S5.1_cleaned_scfa_meta.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

# giving unique IDs to humann entries ====
humann_output$ID <- rownames(humann_output)
humann_output$ID <- paste0("ASV_", humann_output$ID)

# filtering humann output to only SCFA genes ====
scfa_meta$KEGG_code <- gsub(" ", "", scfa_meta$KEGG_code) #remove spaces in names
kegg_codes <- sort(unique(scfa_meta$KEGG_code)) # should be 65
#keeping only the rows that have one of the needed KEGG codes
filtered_humann <- humann_output[Reduce(`|`, lapply(kegg_codes, function(x) grepl(x, humann_output$X..Gene.Family))), ]

#save results ====
#write.table (filtered_humann, "S5.2_raw_scfa_reads.txt", row.names=FALSE,sep = "\t")