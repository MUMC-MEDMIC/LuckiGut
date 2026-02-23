# LOAD LIBRARIES ====
library(phyloseq)
library(microViz)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggpubr) #for manual p values

setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

# read in data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")
index_mh <- read.delim("Output/Files/S1.7_index_mh.txt",sep = "\t",header = TRUE, na.strings=c("","NA"))
load("H:/Penders_lab/Theoretical_work/Paper_2/R_scripts/Input_files/12_timepoints_general_4m_to_14m_familycodes.RData")
load("H:/Penders_lab/Theoretical_work/Paper_2/R_scripts/Generated_files/2_other_variables_MS.RData")

# data cleaning ====
index_mh$Family_ID <- rownames(index_mh)
fam_mh <- sort(index_mh$Family_ID)

meta <- as.data.frame(as.matrix(infant@sam_data))
meta <- full_join(index_mh,meta, by = "Family_ID",multiple = "all")

rm(index_mh,infant,fam_mh)

# food data cleaning ====
#making columns nuemric
timepoints_general_4m_to_14m <- lapply(timepoints_general_4m_to_14m, function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) return(x)
    as.numeric(as.character(x))
  })
  df
})

# adding family code 
timepoints_general_4m_to_14m[[1]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[2]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[3]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[4]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[5]]$Family_ID <- other_variables$ï..kindcode
timepoints_general_4m_to_14m[[6]]$Family_ID <- other_variables$ï..kindcode

rm(other_variables)

# making lists int df
mon4 <- timepoints_general_4m_to_14m[[1]]
mon5 <- timepoints_general_4m_to_14m[[2]]
mon6 <- timepoints_general_4m_to_14m[[3]]
mon9 <- timepoints_general_4m_to_14m[[4]]
mon11 <- timepoints_general_4m_to_14m[[5]]
mon14 <- timepoints_general_4m_to_14m[[6]]

rownames(mon4) <- mon4$Family_ID
rownames(mon5) <- mon5$Family_ID
rownames(mon6) <- mon6$Family_ID
rownames(mon9) <- mon9$Family_ID
rownames(mon11) <- mon11$Family_ID
rownames(mon14) <- mon14$Family_ID

mon4$Age_individual <- "4 months"
mon5$Age_individual <- "5 months"
mon6$Age_individual <- "6 months"
mon9$Age_individual <- "9 months"
mon11$Age_individual <- "11 months"
mon14$Age_individual <- "14 months"

rm(timepoints_general_4m_to_14m)

food <- rbind(mon4,mon5,mon6,mon9,mon11,mon14)

rm(mon4,mon5,mon6,mon9,mon11,mon14)

# Food variety score (FVS) ====
#nr of singular food items (excl. breastfeeding / formulafeeding) max 28
columns_to_exclude <- c("breastfeeding","formulafeeding","Family_ID","Age_individual")
individual_food_items <- food %>% select(-any_of(columns_to_exclude)) %>% ncol()

food <- food %>% mutate(fvs = rowSums(across(-all_of(columns_to_exclude)), na.rm = TRUE))

max(food$fvs)

# Dietary diversity score (DDS) ====
#nr of unique food groups (here: vegetables, fruits, legumes/nuts, meat, fish, eggs, diary, grains (bread, pasta, rice)) max 8 
vegetables <- c("cauliflower", "carrots", "broccoli", "potato","margarine")
fruits <- c("banana", "pear", "apple", "melon", "peach", "kiwi", "orange", "strawberry", "tomato")
legumes <- c("beans", "porr","soyprod")
dairy <- c("milk","yoghurt","butter","cheese","pudding")
grains <- c("bread", "rice", "pasta")

food <- food %>% mutate(vegetables = rowSums(across(all_of(vegetables)), na.rm = TRUE))
food <- food %>% mutate(fruits = rowSums(across(all_of(fruits)), na.rm = TRUE))
food <- food %>% mutate(legumes = rowSums(across(all_of(legumes)), na.rm = TRUE))
food <- food %>% mutate(dairy = rowSums(across(all_of(dairy)), na.rm = TRUE))
food <- food %>% mutate(grains = rowSums(across(all_of(grains)), na.rm = TRUE))

rm(vegetables,fruits,legumes,dairy,grains)

unique_food_groups <- c("vegetables","fruits","legumes","meat","fish","egg","dairy","grains")

food <- food %>%
  mutate(
    dds = rowSums(
      across(all_of(unique_food_groups), ~ as.integer(.x > 0)),
      na.rm = TRUE
    )
  )

# Food allergen diversity (FAD) ====
#nr of main food allergens such as milk, egg, wheat, fish, soy, nuts max 6
allergy <- c("milk", "egg", "fish", "soyprod", "porr", "bread")

food <- food %>%
  mutate(
    fad = rowSums(
      across(all_of(allergy), ~ as.integer(.x > 0)),
      na.rm = TRUE
    )
  )

# selecting needed columns from food ====
food <- food %>% select(all_of(c("Family_ID","Age_individual", "fvs","fad","dds")))

food <- food %>% mutate(merge = paste0(Family_ID, "@", Age_individual))

rm(allergy,individual_food_items,columns_to_exclude,unique_food_groups)

meta <- meta %>% select(all_of(c("Family_ID","Age_individual", "covid_270220")))
meta <- meta %>% mutate(merge = paste0(Family_ID, "@", Age_individual))

# merging data together ====
all <- full_join(meta,food, by = "merge")

all$Age_individual<- factor(all$Age_individual, 
                                  levels = c("4 months",
                                             "5 months",
                                             "6 months",
                                             "9 months",
                                             "11 months",
                                             "14 months"))

rm(food,meta)

#removing samples that do not exist
all <- all %>% filter(!is.na(Family_ID.x))
#removing samples with now diet data
all <- all %>% filter(!is.na(Family_ID.y))

families <- sort(unique(all$Family_ID.x))

idx <- which(!(all$Family_ID.x == all$Family_ID.y) & !(is.na(all$Family_ID.x) & is.na(all$Family_ID.y)))
idx

all$Family_ID.y <- NULL
all <- all %>% rename(Family_ID = Family_ID.x)

idx <- which(!(all$Age_individual.x == all$Age_individual.y) & !(is.na(all$Age_individual.x) & is.na(all$Age_individual.y)))
idx

all$Age_individual.y <- NULL
all <- all %>% rename(Age_individual = Age_individual.x)

rm(idx)

all$Age_individual<- factor(all$Age_individual, levels = c("4 months",
                                                             "5 months",
                                                             "6 months",
                                                             "9 months",
                                                             "11 months",
                                                             "14 months"))

# histograms ====
all %>% ggplot(aes(x = fvs)) + 
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  theme_bw() + ggtitle("Food variety score (max 28)")+facet_wrap(~Age_individual)+
  xlab("Score values") + ylab("Frequency") +
  theme(
    panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
    panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))

all %>% ggplot(aes(x = dds)) + 
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  theme_bw() + ggtitle("Dietary diversity score (max 8)")+facet_wrap(~Age_individual)+
  xlab("Score values") + ylab("Frequency") +
  theme(
    panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
    panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))

all %>% ggplot(aes(x = fad)) + 
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  theme_bw() + ggtitle("Food allergen diversity (max 6)")+facet_wrap(~Age_individual)+
  xlab("Score values") + ylab("Frequency") +
  theme(
    panel.grid.major = element_line(color = adjustcolor("grey95", alpha.f = 0.3)),
    panel.grid.minor = element_line(color = adjustcolor("grey95", alpha.f = 0.2)))

# calculating stats ====
#how many samples left per TP
plyr::count(all$Age_individual)

table <- all %>%
  dplyr::group_by(covid_270220, Age_individual) %>%
  dplyr::summarise(
    N      = dplyr::n(),
    q1  = quantile(fvs, 0.25),
    q3  = quantile(fvs, 0.75),
    median = median(fvs),
    IQR    = IQR(fvs),
    score  = "fvs",
    .groups = "keep"   # keeps grouping by both variables after summarise
  )

temp <- all %>%
  dplyr::group_by(covid_270220, Age_individual) %>%
  dplyr::summarise(
    N      = dplyr::n(),
    q1  = quantile(dds, 0.25),
    q3  = quantile(dds, 0.75),
    median = median(dds),
    IQR    = IQR(dds),
    score  = "dds",
    .groups = "keep"   # keeps grouping by both variables after summarise
  )

table <- rbind(table,temp)

temp <- all %>%
  dplyr::group_by(covid_270220, Age_individual) %>%
  dplyr::summarise(
    N      = dplyr::n(),
    q1  = quantile(fad, 0.25),
    q3  = quantile(fad, 0.75),
    median = median(fad),
    IQR    = IQR(fad),
    score  = "fad",
    .groups = "keep"   # keeps grouping by both variables after summarise
  )

table <- rbind(table,temp)

rm(temp)

# saving data ====
#write.table (all, "S1.8_food_diversity_scores.txt", row.names=FALSE,sep = "\t")

# wilcoxon rank-sum test ====
mon4 <- all %>% filter(Age_individual == "4 months")
mon5 <- all %>% filter(Age_individual == "5 months")
mon6 <- all %>% filter(Age_individual == "6 months")
mon9 <- all %>% filter(Age_individual == "9 months")
mon11 <- all %>% filter(Age_individual == "11 months")
mon14 <- all %>% filter(Age_individual == "14 months")

fvs <- wilcox.test(fvs ~ covid_270220, data = mon4, exact = FALSE)
stat.test <- data.frame(
  p_value = fvs$p.value,
  Z = qnorm(fvs$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fvs$statistic,
  alternative = fvs$alternative,
  method = fvs$method,
  score = "fvs",
  Age_individual = "4 months"
)

dds <- wilcox.test(dds ~ covid_270220, data = mon4, exact = FALSE)
s <- data.frame(
  p_value = dds$p.value,
  Z = qnorm(dds$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = dds$statistic,
  alternative = dds$alternative,
  method = dds$method,
  score = "dds",
  Age_individual = "4 months"
)
stat.test <- rbind(stat.test, s)

fad <- wilcox.test(fad ~ covid_270220, data = mon4, exact = FALSE)
s <- data.frame(
  p_value = fad$p.value,
  Z = qnorm(fad$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fad$statistic,
  alternative = fad$alternative,
  method = fad$method,
  score = "fad",
  Age_individual = "4 months"
)
stat.test <- rbind(stat.test, s)
rm(s,fad,dds,fvs)

fvs <- wilcox.test(fvs ~ covid_270220, data = mon5, exact = FALSE)
s <- data.frame(
  p_value = fvs$p.value,
  Z = qnorm(fvs$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fvs$statistic,
  alternative = fvs$alternative,
  method = fvs$method,
  score = "fvs",
  Age_individual = "5 months"
)
stat.test <- rbind(stat.test, s)

dds <- wilcox.test(dds ~ covid_270220, data = mon5, exact = FALSE)
s <- data.frame(
  p_value = dds$p.value,
  Z = qnorm(dds$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = dds$statistic,
  alternative = dds$alternative,
  method = dds$method,
  score = "dds",
  Age_individual = "5 months"
)
stat.test <- rbind(stat.test, s)

fad <- wilcox.test(fad ~ covid_270220, data = mon5, exact = FALSE)
s <- data.frame(
  p_value = fad$p.value,
  Z = qnorm(fad$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fad$statistic,
  alternative = fad$alternative,
  method = fad$method,
  score = "fad",
  Age_individual = "5 months"
)
stat.test <- rbind(stat.test, s)
rm(s,fad,dds,fvs)

fvs <- wilcox.test(fvs ~ covid_270220, data = mon6, exact = FALSE)
s <- data.frame(
  p_value = fvs$p.value,
  Z = qnorm(fvs$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fvs$statistic,
  alternative = fvs$alternative,
  method = fvs$method,
  score = "fvs",
  Age_individual = "6 months"
)
stat.test <- rbind(stat.test, s)

dds <- wilcox.test(dds ~ covid_270220, data = mon6, exact = FALSE)
s <- data.frame(
  p_value = dds$p.value,
  Z = qnorm(dds$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = dds$statistic,
  alternative = dds$alternative,
  method = dds$method,
  score = "dds",
  Age_individual = "6 months"
)
stat.test <- rbind(stat.test, s)

fad <- wilcox.test(fad ~ covid_270220, data = mon6, exact = FALSE)
s <- data.frame(
  p_value = fad$p.value,
  Z = qnorm(fad$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fad$statistic,
  alternative = fad$alternative,
  method = fad$method,
  score = "fad",
  Age_individual = "6 months"
)
stat.test <- rbind(stat.test, s)
rm(s,fad,dds,fvs)

fvs <- wilcox.test(fvs ~ covid_270220, data = mon9, exact = FALSE)
s <- data.frame(
  p_value = fvs$p.value,
  Z = qnorm(fvs$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fvs$statistic,
  alternative = fvs$alternative,
  method = fvs$method,
  score = "fvs",
  Age_individual = "9 months"
)
stat.test <- rbind(stat.test, s)

dds <- wilcox.test(dds ~ covid_270220, data = mon9, exact = FALSE)
s <- data.frame(
  p_value = dds$p.value,
  Z = qnorm(dds$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = dds$statistic,
  alternative = dds$alternative,
  method = dds$method,
  score = "dds",
  Age_individual = "9 months"
)
stat.test <- rbind(stat.test, s)

fad <- wilcox.test(fad ~ covid_270220, data = mon9, exact = FALSE)
s <- data.frame(
  p_value = fad$p.value,
  Z = qnorm(fad$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fad$statistic,
  alternative = fad$alternative,
  method = fad$method,
  score = "fad",
  Age_individual = "9 months"
)
stat.test <- rbind(stat.test, s)
rm(s,fad,dds,fvs)

fvs <- wilcox.test(fvs ~ covid_270220, data = mon11, exact = FALSE)
s <- data.frame(
  p_value = fvs$p.value,
  Z = qnorm(fvs$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fvs$statistic,
  alternative = fvs$alternative,
  method = fvs$method,
  score = "fvs",
  Age_individual = "11 months"
)
stat.test <- rbind(stat.test, s)

dds <- wilcox.test(dds ~ covid_270220, data = mon11, exact = FALSE)
s <- data.frame(
  p_value = dds$p.value,
  Z = qnorm(dds$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = dds$statistic,
  alternative = dds$alternative,
  method = dds$method,
  score = "dds",
  Age_individual = "11 months"
)
stat.test <- rbind(stat.test, s)

fad <- wilcox.test(fad ~ covid_270220, data = mon11, exact = FALSE)
s <- data.frame(
  p_value = fad$p.value,
  Z = qnorm(fad$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fad$statistic,
  alternative = fad$alternative,
  method = fad$method,
  score = "fad",
  Age_individual = "11 months"
)
stat.test <- rbind(stat.test, s)
rm(s,fad,dds,fvs)

fvs <- wilcox.test(fvs ~ covid_270220, data = mon14, exact = FALSE)
s <- data.frame(
  p_value = fvs$p.value,
  Z = qnorm(fvs$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fvs$statistic,
  alternative = fvs$alternative,
  method = fvs$method,
  score = "fvs",
  Age_individual = "14 months"
)
stat.test <- rbind(stat.test, s)

dds <- wilcox.test(dds ~ covid_270220, data = mon14, exact = FALSE)
s <- data.frame(
  p_value = dds$p.value,
  Z = qnorm(dds$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = dds$statistic,
  alternative = dds$alternative,
  method = dds$method,
  score = "dds",
  Age_individual = "14 months"
)
stat.test <- rbind(stat.test, s)

fad <- wilcox.test(fad ~ covid_270220, data = mon14, exact = FALSE)
s <- data.frame(
  p_value = fad$p.value,
  Z = qnorm(fad$p.value / 2),
  group1 = 'yes',
  group2 = 'no',
  W_value = fad$statistic,
  alternative = fad$alternative,
  method = fad$method,
  score = "fad",
  Age_individual = "14 months"
)
stat.test <- rbind(stat.test, s)
rm(s,fad,dds,fvs)

rm(mon4,mon5,mon6,mon9,mon11,mon14)

stat.test$p_value_test_adj <- p.adjust(stat.test$p_value, method = "BH")

# wilcox results cleaning ====
#function that will round up all numeric values in p-value table
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

#rounding up the values in p-value table
stat.test <- round_df(stat.test, digits=3)

stat.test <- stat.test %>% 
  mutate(p_sig_test_adj = if_else(
    condition = .$p_value_test_adj > 0.05, 
    true = NA_character_, 
    false = if_else(
      condition = .$p_value_test_adj < 0.001,
      true = "***",
      false = if_else(
        condition = .$p_value_test_adj <= 0.05&.$p_value_test_adj > 0.01,
        true = "*",
        false = "**"))))

stat.test <- stat.test %>% 
  mutate(p_sig = if_else(
    condition = .$p_value > 0.05, 
    true = NA_character_, 
    false = if_else(
      condition = .$p_value < 0.001,
      true = "***",
      false = if_else(
        condition = .$p_value <= 0.05&.$p_value > 0.01,
        true = "*",
        false = "**"))))

rm(round_df)

stat.test$Age_individual<- factor(stat.test$Age_individual, 
                           levels = c("4 months",
                                      "5 months",
                                      "6 months",
                                      "9 months",
                                      "11 months",
                                      "14 months"))

# saving data ====
#write.table (stat.test, "S1.8_food_diversity_scores_wilcox_test.txt", row.names=FALSE,sep = "\t")

# barplots ====
wilcox_fvs <- stat.test %>% filter(score == "fvs")
wilcox_dds <- stat.test %>% filter(score == "dds")
wilcox_fad <- stat.test %>% filter(score == "fad")

#samples per TP 
table %>% filter(score == "fvs") %>% ggplot(aes(y = N, x = covid_270220)) +
  geom_bar(stat = "identity") + 
  theme_bw()+ggtitle("Number of samples with diet data")+facet_wrap(~Age_individual)+
  xlab("collected in pandemic")+ylab("N")+theme(legend.position = "none",
                                axis.text.x = element_text(size = 9),
                                axis.text.y = element_text(size = 11),
                                axis.title.x = element_text(size = 13), 
                                axis.title.y = element_text(size = 13))+scale_color_viridis(discrete = TRUE)

#boxplots
fvs <- all  %>% ggplot(aes(y = fvs, x = covid_270220, colour = covid_270220)) +
  geom_boxplot() + 
  theme_bw()+ggtitle("Food variety score, sig p values (not adj)")+facet_wrap(~Age_individual)+
  xlab("collected in pandemic")+ylab("score")+theme(legend.position = "none",
                                                     axis.text.x = element_text(size = 9),
                                                     axis.text.y = element_text(size = 11),
                                                     axis.title.x = element_text(size = 13), 
                                                     axis.title.y = element_text(size = 13))+scale_color_viridis(discrete = TRUE)

#adding p-values
fvs <- fvs + stat_pvalue_manual(wilcox_fvs,y.position = 28, label = "p_sig",hide.ns = TRUE,
                                                size = 5)
fvs

dds <- all  %>% ggplot(aes(y = dds, x = covid_270220, colour = covid_270220)) +
  geom_boxplot() + 
  theme_bw()+ggtitle("Dietary diversity score, sig p values (not adj)")+facet_wrap(~Age_individual)+
  xlab("collected in pandemic")+ylab("score")+theme(legend.position = "none",
                                                    axis.text.x = element_text(size = 9),
                                                    axis.text.y = element_text(size = 11),
                                                    axis.title.x = element_text(size = 13), 
                                                    axis.title.y = element_text(size = 13))+scale_color_viridis(discrete = TRUE)

#adding p-values
dds <- dds + stat_pvalue_manual(wilcox_dds,y.position = 3, label = "p_sig",hide.ns = TRUE,
                                size = 5)
dds

fad <- all  %>% ggplot(aes(y = fad, x = covid_270220, colour = covid_270220)) +
  geom_boxplot() + 
  theme_bw()+ggtitle("Food allergen diversity, sig p values (not adj)")+facet_wrap(~Age_individual)+
  xlab("collected in pandemic")+ylab("score")+theme(legend.position = "none",
                                                    axis.text.x = element_text(size = 9),
                                                    axis.text.y = element_text(size = 11),
                                                    axis.title.x = element_text(size = 13), 
                                                    axis.title.y = element_text(size = 13))+scale_color_viridis(discrete = TRUE)

#adding p-values
fad <- fad + stat_pvalue_manual(wilcox_fad,y.position = 3, label = "p_sig",hide.ns = TRUE,
                                size = 5)
fad

