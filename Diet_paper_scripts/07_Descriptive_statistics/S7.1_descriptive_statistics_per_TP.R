# load libraries ====
library(tidyverse)
library(phyloseq)
library(viridis)
setwd("H:/Penders_lab/Theoretical_work/Paper_2")

# read in phyloseq data ====
baku_f_I <- readRDS(file = "Output/Files/S1.1_cleaned_phyloseq_object_MS.RDS")

# extracting data ====
meta <- as.data.frame(as.matrix(baku_f_I@sam_data))
rm(baku_f_I)

# selecting variables ====
milk <- c("breastfeeding_4m","breastfeeding_5m","breastfeeding_6m","breastfeeding_9m",
          "breastfeeding_11m","breastfeeding_14m")
formula <- c("formulafeeding_4m","formulafeeding_5m","formulafeeding_6m","formulafeeding_9m",
             "formulafeeding_11m","formulafeeding_14m")
anti <- c("ch4mq5_ab_freq_MS_4m","ch5mq5_ab_freq_MS_5m","ch6mq16_ab_freq_MS_6m","inf_ab_ED_9m","inf_ab_ED_11m","inf_ab_ED_14m")

other <- c("delivery_place_word","delivery_type_bi_word","furrypet_indoors_MS_6m_word","furrypet_indoors_6m_word_ED","sex_word")

numeric <- c("birth_weigth","ageintrosolids","bqq2_pregn_weeks")

all <- c(milk,anti,formula,other,numeric,"Age_individual","diet_class")

meta_numer <- meta %>% select(all_of(all))
meta_word <- meta %>% select(all_of(c(milk,anti,formula,other,"Age_individual","diet_class")))

rm(milk,formula,anti,all, meta)

# counting variables ====
df_long <- meta_word %>%
  pivot_longer(
    cols = -Age_individual,
    names_to = "variable",
    values_to = "value"
  )
rm(meta_word)

# Step 3: Group and count
df_counts <- df_long %>%
  group_by(Age_individual, variable, value) %>%
  summarise(count = n(), .groups = "drop")

rm(df_long)

suffixes <- c("_4m", "_5m", "_6m", "_9m", "_11m", "_14m")

stats_all <- df_counts %>%
  filter(
    (Age_individual == "1-2 weeks"   & str_detect(variable, "_1w")) |
      (Age_individual == "1-2 weeks"   & str_detect(variable, "1_2w")) |
      (Age_individual == "4 weeks"   & str_detect(variable, "_4w")) |
      (Age_individual == "8 weeks"   & str_detect(variable, "_8w")) |
      (Age_individual == "4 months" & str_detect(variable, "_4m")) |
      (Age_individual == "4 months" & str_detect(variable, "ch4")) |
      (Age_individual == "5 months" & str_detect(variable, "_5m")) |
      (Age_individual == "5 months" & str_detect(variable, "ch5")) |
      (Age_individual == "6 months" & str_detect(variable, "_6m")) |
      (Age_individual == "6 months" & str_detect(variable, "ch6")) |
      (Age_individual == "9 months" & str_detect(variable, "_9m")) |
      (Age_individual == "11 months" & str_detect(variable, "_11m")) |
      (Age_individual == "14 months" & str_detect(variable, "_14m")))

df_counts_filt <- df_counts %>%
  dplyr::filter(!stringr::str_detect(variable, paste0("(", paste(suffixes, collapse = "|"), ")$")))

stats <- rbind(stats_all,df_counts_filt)

rm(df_counts,df_counts_filt,stats_all)

# calculating % ====
stats_all <- stats %>%
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

# calculating numeric values ====
meta_numer <- meta_numer %>% select(all_of(c(numeric,"Age_individual")))
ages <- unique(meta_numer$Age_individual)
meta_numer[names(meta_numer) != "Age_individual"] <- lapply(meta_numer[names(meta_numer) != "Age_individual"], as.numeric)

#shapiro
shapiro <- meta_numer %>%
  group_split(Age_individual) %>%
  map_dfr(function(df_age) {
    age <- unique(df_age$Age_individual)
    
    sh <- df_age %>%
      select(any_of(numeric)) %>%
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
summary <- meta_numer %>%
  group_split(Age_individual) %>%
  map_dfr(function(df_age) {
    age <- unique(df_age$Age_individual)
    
    df_age %>%
      select(any_of(numeric)) %>%
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

n_result <- meta_numer %>%
  group_split(Age_individual) %>%
  map_dfr(function(df_age) {
    age <- unique(df_age$Age_individual)
    
    df_age %>%
      select(any_of(numeric)) %>%
      select(where(is.numeric)) %>%
      imap_dfr(function(x, nm) {
        tibble(
          score = nm,
          n = sum(!is.na(x)),
          Age_individual = age
        )
      })
  })

# looking at the numbers ====
test <- stats_all %>% filter(grepl("formulafeeding", variable))

# saving data ====
#write.table (summary, "S7.1_variables_per_TP_median_iqr.txt", row.names=FALSE,sep = "\t")
#write.table (stats_all, "S7.1_variables_per_TP_character.txt", row.names=FALSE,sep = "\t")
