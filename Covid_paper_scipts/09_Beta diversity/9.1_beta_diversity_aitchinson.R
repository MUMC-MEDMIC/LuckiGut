# libraries ====
library(microViz)
library(phyloseq)
library(vegan)
library(viridis)
library(ggpubr) #for manual p values
library(broom)
library(pairwiseAdonis)
library(reshape2)
library(emmeans)
library(lme4)
library(lmerTest)
library(tidyverse)

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import and clean data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
sample_data <- as.data.frame(as.matrix(infant@sam_data))

keys <- sample_data %>% select("MMHP_SampleID", "Family_ID","Age_individual")

# break TP into pairs ====
tp1 <- infant %>% ps_filter(Age_individual == "1-2 weeks" | Age_individual == "4 weeks")
tp2 <- infant %>% ps_filter(Age_individual == "4 weeks" | Age_individual == "8 weeks")
tp3 <- infant %>% ps_filter(Age_individual == "8 weeks" | Age_individual == "4 months")
tp4 <- infant %>% ps_filter(Age_individual == "4 months" | Age_individual == "5 months")
tp5 <- infant %>% ps_filter(Age_individual == "5 months" | Age_individual == "6 months")
tp6 <- infant %>% ps_filter(Age_individual == "6 months" | Age_individual == "9 months")
tp7 <- infant %>% ps_filter(Age_individual == "9 months" | Age_individual == "11 months")
tp8 <- infant %>% ps_filter(Age_individual == "11 months" | Age_individual == "14 months")

rm(infant)

# keeping only duplicates (so infant has samples for both TP) ====
tp1 <- tp1 %>%ps_filter(Family_ID %in% Family_ID[duplicated(Family_ID)])
tp2 <- tp2 %>%ps_filter(Family_ID %in% Family_ID[duplicated(Family_ID)])
tp3 <- tp3 %>%ps_filter(Family_ID %in% Family_ID[duplicated(Family_ID)])
tp4 <- tp4 %>%ps_filter(Family_ID %in% Family_ID[duplicated(Family_ID)])
tp5 <- tp5 %>%ps_filter(Family_ID %in% Family_ID[duplicated(Family_ID)])
tp6 <- tp6 %>%ps_filter(Family_ID %in% Family_ID[duplicated(Family_ID)])
tp7 <- tp7 %>%ps_filter(Family_ID %in% Family_ID[duplicated(Family_ID)])
tp8 <- tp8 %>%ps_filter(Family_ID %in% Family_ID[duplicated(Family_ID)])

sample_data1 <- as.data.frame(as.matrix(tp1@sam_data))
sample_data2 <- as.data.frame(as.matrix(tp2@sam_data))
sample_data3 <- as.data.frame(as.matrix(tp3@sam_data))
sample_data4 <- as.data.frame(as.matrix(tp4@sam_data))
sample_data5 <- as.data.frame(as.matrix(tp5@sam_data))
sample_data6 <- as.data.frame(as.matrix(tp6@sam_data))
sample_data7 <- as.data.frame(as.matrix(tp7@sam_data))
sample_data8 <- as.data.frame(as.matrix(tp8@sam_data))

rm(sample_data1,sample_data2,sample_data3,sample_data4,sample_data5,sample_data6,sample_data7,sample_data8)

# calculating aitchison distance ====
tp_aitch1 <- tp1 %>%  dist_calc("aitchison") # input data should be virgin, not transformed
tp_aitch2 <- tp2 %>%  dist_calc("aitchison") 
tp_aitch3 <- tp3 %>%  dist_calc("aitchison") 
tp_aitch4 <- tp4 %>%  dist_calc("aitchison") 
tp_aitch5 <- tp5 %>%  dist_calc("aitchison") 
tp_aitch6 <- tp6 %>%  dist_calc("aitchison") 
tp_aitch7 <- tp7 %>%  dist_calc("aitchison")
tp_aitch8 <- tp8 %>%  dist_calc("aitchison") 

rm(tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8)

distance1 <- dist_get(tp_aitch1)
distance2 <- dist_get(tp_aitch2)
distance3 <- dist_get(tp_aitch3)
distance4 <- dist_get(tp_aitch4)
distance5 <- dist_get(tp_aitch5)
distance6 <- dist_get(tp_aitch6)
distance7 <- dist_get(tp_aitch7)
distance8 <- dist_get(tp_aitch8)

rm(tp_aitch1,tp_aitch2,tp_aitch3,tp_aitch4,tp_aitch5,tp_aitch6,tp_aitch7,tp_aitch8)

# taking distances out =====
# Convert the distance matrix to a data frame
distance1 <- distance1 %>% as.matrix() %>% melt()
distance2 <- distance2 %>% as.matrix() %>% melt()
distance3 <- distance3 %>% as.matrix() %>% melt()
distance4 <- distance4 %>% as.matrix() %>% melt()
distance5 <- distance5 %>% as.matrix() %>% melt()
distance6 <- distance6 %>% as.matrix() %>% melt()
distance7 <- distance7 %>% as.matrix() %>% melt()
distance8 <- distance8 %>% as.matrix() %>% melt()

distance <- rbind(distance1,distance2,distance3,distance4,distance5,distance6,distance7,distance8)
rm(distance1,distance2,distance3,distance4,distance5,distance6,distance7,distance8)

distance <- distance %>% filter(Var1 != Var2) #removing comparisons with itself

# merging with keys ====
keys <- keys %>% rename(Var1 = MMHP_SampleID)
keys <- keys %>% rename(Family_ID1 = Family_ID)
keys <- keys %>% rename(Age_individual1 = Age_individual)

distance <- full_join(keys,distance, by = "Var1",multiple = "all")

keys <- keys %>% rename(Var2 = Var1)
keys <- keys %>% rename(Family_ID2 = Family_ID1)
keys <- keys %>% rename(Age_individual2 = Age_individual1)

distance <- full_join(keys,distance, by = "Var2",multiple = "all")
distance <- distance %>% filter(!is.na(value))
distance <- distance %>% mutate(pairs = paste(Age_individual1, Age_individual2, sep = " vs "))
distance <- distance %>% filter(Family_ID1 == Family_ID2)

needed_pairs <- c("1-2 weeks vs 4 weeks","4 weeks vs 8 weeks","8 weeks vs 4 months","4 months vs 5 months","5 months vs 6 months",
                  "6 months vs 9 months","9 months vs 11 months","11 months vs 14 months")

distance <- distance %>% filter(pairs %in% needed_pairs)

distance$pairs<- factor(distance$pairs, 
                           levels = c("1-2 weeks vs 4 weeks",
                                      "4 weeks vs 8 weeks",
                                      "8 weeks vs 4 months",
                                      "4 months vs 5 months",
                                      "5 months vs 6 months",
                                      "6 months vs 9 months",
                                      "9 months vs 11 months",
                                      "11 months vs 14 months"))

#write.table (distance, "S9.1_aitchinson_distance_between_TP_beta_diversity.txt", row.names=FALSE,sep = "\t")

rm(keys,needed_pairs)

