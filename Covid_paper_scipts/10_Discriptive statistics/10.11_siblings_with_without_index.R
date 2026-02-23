# load libraries ====
library(tidyverse)
library(ggplot2)
library(phyloseq)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
index = read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

# cleaning data =====
index$Family_ID <- rownames(index)

meta <- as.data.frame(as.matrix(infant@sam_data))

meta <- meta %>% select(all_of(c("older_siblings","covid_270220","Age_individual","Family_ID")))

meta <- full_join(index, meta, by = "Family_ID",multiple = "all")

meta$Age_individual<- factor(meta$Age_individual, levels = c("1-2 weeks",
                                         "4 weeks",
                                         "8 weeks",
                                         "4 months",
                                         "5 months",
                                         "6 months",
                                         "9 months",
                                         "11 months",
                                         "14 months"))

rm(infant,index)

# selecting infants with index data ====
meta <- meta %>% filter(!is.na(index_m_h))

index_fam <- meta %>% group_by(older_siblings, Age_individual) %>% summarise(n = n(), .groups = "drop")

# barplot ====
ggplot(index_fam, aes(x = older_siblings, y = n, fill = Age_individual)) +
  geom_col(position = "dodge") +facet_wrap(~Age_individual)+
  labs(x = "Older Siblings", 
       y = "Sample Count", 
       fill = "Age Group") +
  theme_minimal() +
  theme(legend.position = "bottom")+
  ggtitle("Index families") 
