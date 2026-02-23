# load libraries ====
library(tidyverse)
library(phyloseq)
library(viridis)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# import data ====
ps <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
food_score <- read.delim("Output/Files/S1.8_food_diversity_scores.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))

ps <- ps  %>% ps_filter(covid_270220 == "yes")

# calculating number of samples and sets
meta <- as.data.frame(as.matrix(ps@sam_data))
rm(ps)

# merging meta with food scores ====
meta <- meta %>% mutate(merge = paste0(Family_ID, "@", Age_individual))
food_score <- food_score %>% filter(covid_270220 == "yes") %>% dplyr::select(-dplyr::all_of(c("Family_ID", "Age_individual", "covid_270220")))
meta <- full_join(meta,food_score, by = "merge")

# selecting variables ====
milk <- c("breastfeeding_1_2w","breastfeeding_4w","breastfeeding_8w","breastfeeding_4m","breastfeeding_5m","breastfeeding_6m","breastfeeding_9m",
          "breastfeeding_11m","breastfeeding_14m")
formula <- c("formulafeeding_1_2w","formulafeeding_4w","formulafeeding_8w","formulafeeding_4m","formulafeeding_5m","formulafeeding_6m","formulafeeding_9m",
             "formulafeeding_11m","formulafeeding_14m")
anti <- c("inf_ab_1w","inf_ab_4w","inf_ab_8w","inf_ab_4m","inf_ab_5m","inf_ab_6m","inf_ab_9m","inf_ab_11m","inf_ab_14m")

food <- c("fvs","fad","dds")

all <- c(milk,anti,formula,"Age_individual")

meta_perm <- meta %>% select(all_of(all))

rm(milk,formula,anti,all)

# counting variables ====
df_long <- meta_perm %>%
  pivot_longer(
    cols = -Age_individual,
    names_to = "variable",
    values_to = "value"
  )
rm(meta_perm)

# Step 3: Group and count
df_counts <- df_long %>%
  group_by(Age_individual, variable, value) %>%
  summarise(count = n(), .groups = "drop")

rm(df_long)

stats_all <- df_counts %>%
  filter(
    (Age_individual == "1-2 weeks"   & str_detect(variable, "_1w")) |
      (Age_individual == "1-2 weeks"   & str_detect(variable, "1_2w")) |
      (Age_individual == "4 weeks"   & str_detect(variable, "_4w")) |
      (Age_individual == "8 weeks"   & str_detect(variable, "_8w")) |
      (Age_individual == "4 months" & str_detect(variable, "_4m")) |
      (Age_individual == "5 months" & str_detect(variable, "_5m")) |
      (Age_individual == "6 months" & str_detect(variable, "_6m")) |
      (Age_individual == "9 months" & str_detect(variable, "_9m")) |
      (Age_individual == "11 months" & str_detect(variable, "_11m")) |
      (Age_individual == "14 months" & str_detect(variable, "_14m")))

rm(df_counts)

# calculating % ====
stats_all <- stats_all %>%
  group_by(Age_individual, variable) %>%
  mutate(n = sum(count)) %>%
  ungroup()

stats_all$prosent <-  round(stats_all$count * 100 / stats_all$n, 1)

stats_all <- stats_all %>% dplyr::mutate(count_prosent = paste0(count, " (", prosent, "%)"))

# adding RN ====
stats_all <- stats_all %>%
  mutate(
    RN = case_when(
      Age_individual == "1-2 weeks" ~ 1,
      Age_individual == "4 weeks" ~ 2,
      Age_individual == "8 weeks" ~ 3,
      Age_individual == "4 months" ~ 4,
      Age_individual == "5 months" ~ 5,
      Age_individual == "6 months" ~ 6,
      Age_individual == "9 months" ~ 7,
      Age_individual == "11 months" ~ 8,
      Age_individual == "14 months" ~ 9,
      TRUE ~ NA_real_  # or 0, or whatever default you want
    ))

# calculating food score ====
meta <- meta %>% select(all_of(c(food,"Age_individual"))) %>% filter(!is.na(fvs))

ages <- unique(meta$Age_individual)

#shapiro
shapiro <- meta %>%
  group_split(Age_individual) %>%
  map_dfr(function(df_age) {
    age <- unique(df_age$Age_individual)
    
    sh <- df_age %>%
      select(any_of(food)) %>%
      map(~ shapiro.test(stats::na.omit(.x)))
    
    as.data.frame(t(sapply(sh, `[`, c("statistic", "p.value")))) %>%
      rownames_to_column("score") %>%
      mutate(
        statistic = as.numeric(statistic),
        p.value   = as.numeric(p.value),
        Age_individual = age
      )
  })

# summary
summary <- meta %>%
  group_split(Age_individual) %>%
  map_dfr(function(df_age) {
    age <- unique(df_age$Age_individual)
    
    df_age %>%
      select(any_of(food)) %>%
      select(where(is.numeric)) %>%
      imap_dfr(function(x, nm) {
        s <- summary(x, na.rm = TRUE)
        
        tibble(
          score   = nm,
          `Min.`    = as.numeric(s["Min."]),
          `1st Qu.` = as.numeric(s["1st Qu."]),
          Median  = as.numeric(s["Median"]),
          Mean    = as.numeric(s["Mean"]),
          `3rd Qu.` = as.numeric(s["3rd Qu."]),
          Max     = as.numeric(s["Max."]),
          Age_individual = age
        )
      })
  })

summary$IQR <- summary$`3rd Qu.` - summary$`1st Qu.`

n_result <- meta %>%
  group_split(Age_individual) %>%
  map_dfr(function(df_age) {
    age <- unique(df_age$Age_individual)
    
    df_age %>%
      select(any_of(food)) %>%
      select(where(is.numeric)) %>%
      imap_dfr(function(x, nm) {
        tibble(
          score = nm,
          n = sum(!is.na(x)),
          Age_individual = age
        )
      })
  })

# saving data ====
#write.table (summary, "S10.5_food_diversity_scores_pandemic_median_iqr.txt", row.names=FALSE,sep = "\t")