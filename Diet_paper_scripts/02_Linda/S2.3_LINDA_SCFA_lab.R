# load in packages ====
library(tidyverse) 
library(phyloseq)
library(microViz)
library(plyr)
library(reshape2)
library(dplyr)
library(readxl)
library(LinDA)
library(ComplexHeatmap)
library(patchwork)
library(gridExtra)
library(svglite)

setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# loading data ====
ps_diet <- readRDS("Output/Files/S1.1_cleaned_phyloseq_object_MS.rds")
scfa_raw <- read.delim("Output/Files/S5.12_scfa_results_raw.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

#cleaning data 
ps_diet <- ps_diet  %>%
  ps_mutate(
    timepoint_number = factor(str_remove(ps_diet@sam_data$timepoint,"m"), levels = c("4","5","6","9","11","14")),
  )

# filtering abundance > 0.01 ====
ps_diet <- phyloseq_validate(ps_diet, remove_undetected = TRUE)
ps_diet_comp <- ps_diet %>% tax_transform("compositional") %>% ps_get()

#filter per TP
ps_diet_4 <- ps_diet_comp %>% ps_filter(timepoint == "4m") 
ps_diet_6 <- ps_diet_comp %>% ps_filter(timepoint == "6m") 

#extract OTU table
otu_4 <- ps_diet_4 %>% otu_table() %>% as.data.frame()
otu_6 <- ps_diet_6 %>% otu_table() %>% as.data.frame()
rm(ps_diet_4,ps_diet_6,ps_diet_comp)

#calculate row mean per taxa
otu_4$mean <- rowMeans(otu_4)
otu_6$mean <- rowMeans(otu_6)

#select those bakut that are above 0.01
bakut_4 <- otu_4 %>% filter(mean >= 0.01) %>% rownames()
bakut_6 <- otu_6 %>% filter(mean >= 0.01) %>% rownames()
rm(otu_4,otu_6)

#combine into 1 vector
unik_baku <- c(bakut_4,bakut_6)
unik_baku <- unique(unik_baku)
rm(bakut_4,bakut_6)

#remove everything below 1% from phyloseq
ps_diet_sel <- ps_diet %>% tax_select(unik_baku)
rm(unik_baku)

# Create low/high SCFA groups ====
