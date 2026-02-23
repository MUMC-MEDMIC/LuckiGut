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

meta <- meta %>% select(all_of(c("Sample_date_filled","covid_270220","Age_individual","Family_ID","MMHP_SampleID")))

meta$count_down_date <- as.Date("27/02/2020", format = "%d/%m/%Y")

meta <- full_join(index, meta, by = "Family_ID",multiple = "all")

meta$Sample_date_filled <- as.Date(meta$Sample_date_filled)

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

# selecting infants of the covid and index period ====
index_family_vec <- meta %>% filter(covid_270220 == "yes") %>% filter(!is.na(index_m_h))
index_family_vec <- sort(unique(index_family_vec$Family_ID))

meta <- meta %>% mutate(index_family = if_else(Family_ID %in% index_family_vec, "yes", "no"))

covid_family_vec <- meta %>% filter(covid_270220 == "yes")
covid_family_vec <- sort(unique(covid_family_vec$Family_ID))

meta <- meta %>% mutate(covid_family = if_else(Family_ID %in% covid_family_vec, "yes", "no"))

# counting time spent in pandemic ====
meta$time_spent <- meta$Sample_date_filled - meta$count_down_date

meta <- meta %>% mutate(summary = pmax(0, as.numeric(time_spent)))

meta_in_pandemic <- meta %>%
  filter(covid_270220 == "yes") %>%  group_by(Age_individual) %>% 
  mutate(
    median_days = median(time_spent, na.rm = TRUE),
    days_group = case_when(
      is.na(time_spent) ~ NA_character_,
      time_spent <= median_days ~ "below_or_equal_median",
      time_spent >  median_days ~ "above_median"
    )
  ) %>% ungroup()

# histogram ====
meta %>% filter(index_family == "yes") %>% ggplot(aes(x = summary)) + 
  geom_histogram(binwidth = 10, fill = "steelblue", color = "white") +
  theme_bw() + ggtitle("36 families with index")+facet_wrap(~Age_individual)+
  xlab("Days between sample collection and start of covid") + ylab("Frequency") +
  theme(
    panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
    panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))

meta %>% ggplot(aes(x = summary)) + 
  geom_histogram(binwidth = 10, fill = "steelblue", color = "white") +
  theme_bw() + ggtitle("All families")+facet_wrap(~Age_individual)+
  xlab("Days between sample collection and start of covid") + ylab("Frequency") +
  theme(
    panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
    panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))

meta %>% filter(covid_family == "yes") %>% ggplot(aes(x = summary)) + 
  geom_histogram(binwidth = 10, fill = "steelblue", color = "white") +
  theme_bw() + ggtitle("43 families that collected samples during covid")+facet_wrap(~Age_individual)+
  xlab("Days between sample collection and start of covid") + ylab("Frequency") +
  theme(
    panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
    panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))

# saving data ====
meta <- meta %>% select(all_of(c("MMHP_SampleID","time_spent")))
#write.table (meta, "S10.10_time_spent_in_pandemic.txt", row.names=FALSE,sep = "\t")

meta_in_pandemic <- meta_in_pandemic %>% select(all_of(c("MMHP_SampleID","time_spent","days_group")))
#write.table (meta_in_pandemic, "S10.10_time_spent_in_pandemic_dichotomized.txt", row.names=FALSE,sep = "\t")


