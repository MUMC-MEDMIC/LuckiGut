# libraries ====
#This script uses ENS calculated on unfiltered reads! ENS was calculated in 4_creating_phyloseq_object_for_bacteria.R
#Alpha diversity, ENS boxplot
library(tidyverse)
library(microViz)
library(phyloseq)
library(vegan)
library(viridis)
library(ggpubr) #for manual p values
library(broom)
library(svglite)
setwd("H:/Penders_lab/Theoretical_work/Paper_1/")

#old used test_1w_4w <- pairwise.wilcox.test(wil_1w_4w$ENS, wil_1w_4w$Age_individual)

# import data ====
infant <- readRDS("Output/Files/S1.6_phyloseq_object_bacteria_filtered_infant.rds")

# cleaning data for analysis ====
infant@sam_data$Age_individual <- gsub("weeks", "wks", infant@sam_data$Age_individual)
infant@sam_data$Age_individual <- gsub("months", "mos", infant@sam_data$Age_individual)
infant@sam_data$Age_individual <- gsub("1-2 wks", "1 wk", infant@sam_data$Age_individual)

infant@sam_data$Age_individual<- factor(infant@sam_data$Age_individual, levels = c("1 wk", "4 wks", "8 wks",
                                                                                   "4 mos", "5 mos", "6 mos", "9 mos", "11 mos", "14 mos"))

stat <- as.matrix(infant@sam_data)
stat_test <- as.data.frame(stat)
rm(stat,infant)

stat_test$ENS<- as.numeric(stat_test$ENS)
stat_test$observed_richnes<- as.numeric(stat_test$observed_richnes)

stat_test$Age_individual<- factor(stat_test$Age_individual, levels = c("1 wk", "4 wks", "8 wks",
                                                                                   "4 mos", "5 mos", "6 mos", "9 mos", "11 mos", "14 mos"))

# dividing samples into pairs ====
#as samples are connected to one child, to do pairwise comparison we need to do manually as groups n are not equal
#sorting timepoints so they are paired
wil_1w_4w <- stat_test %>%filter(Age_individual == "1 wk" | Age_individual == "4 wks") %>% subset(duplicated(Family_ID) | duplicated(Family_ID, fromLast=TRUE))
wil_4w_8w <- stat_test %>%filter(Age_individual == "4 wks" | Age_individual == "8 wks") %>% subset(duplicated(Family_ID) | duplicated(Family_ID, fromLast=TRUE))
wil_8w_4m <- stat_test %>%filter(Age_individual == "8 wks" | Age_individual == "4 mos") %>% subset(duplicated(Family_ID) | duplicated(Family_ID, fromLast=TRUE))
wil_4m_5m <- stat_test %>%filter(Age_individual == "4 mos" | Age_individual == "5 mos") %>% subset(duplicated(Family_ID) | duplicated(Family_ID, fromLast=TRUE))
wil_5m_6m <- stat_test %>%filter(Age_individual == "5 mos" | Age_individual == "6 mos") %>% subset(duplicated(Family_ID) | duplicated(Family_ID, fromLast=TRUE))
wil_6m_9m <- stat_test %>%filter(Age_individual == "6 mos" | Age_individual == "9 mos") %>% subset(duplicated(Family_ID) | duplicated(Family_ID, fromLast=TRUE))
wil_9m_11m <- stat_test %>%filter(Age_individual == "9 mos" | Age_individual == "11 mos") %>% subset(duplicated(Family_ID) | duplicated(Family_ID, fromLast=TRUE))
wil_11m_14m <- stat_test %>%filter(Age_individual == "11 mos" | Age_individual == "14 mos") %>% subset(duplicated(Family_ID) | duplicated(Family_ID, fromLast=TRUE))

# ENS normality check ====
# ENS histogramm
ggplot(stat_test, aes(x=ENS)) + 
  geom_histogram(color="black", fill="white")+ facet_wrap(~Age_individual)+
  ggtitle("ENS distribution per TP")

# ENS shapiro
# Perform Shapiro-Wilk test for each age group
shapiro_test_results <- stat_test %>%
  group_by(Age_individual) %>%
  summarise(shapiro_test = list(shapiro.test(ENS))) %>%
  mutate(test_results = map(shapiro_test, tidy))
# Create the result dataframe
shapiro_ens <- shapiro_test_results %>%
  unnest(test_results) %>%
  select(Age_individual, statistic, p.value) %>%
  rename(W = statistic, p_value = p.value)
shapiro_ens$alpha_diversity <- "ENS"
rm(shapiro_test_results)
#checks
#test <- stat_test %>% filter(Age_individual == "1-2 weeks")
#shapiro.test(test$ENS)

# ENS qqline and qq test
ggplot(stat_test, aes(sample = ENS)) + 
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~Age_individual) +
  ggtitle("ENS Q-Q Plot per TP")

# Observed richness normality check ====
# Observed richness histogram
ggplot(stat_test, aes(x=observed_richnes)) + 
  geom_histogram(color="black", fill="white")+ facet_wrap(~Age_individual)+
  ggtitle("Observed richnes distribution per TP")

# Observed richness shapiro
# Perform Shapiro-Wilk test for each age group
shapiro_test_results <- stat_test %>%
  group_by(Age_individual) %>%
  summarise(shapiro_test = list(shapiro.test(observed_richnes))) %>%
  mutate(test_results = map(shapiro_test, tidy))
# Create the result dataframe
shapiro_rich <- shapiro_test_results %>%
  unnest(test_results) %>%
  select(Age_individual, statistic, p.value) %>%
  rename(W = statistic, p_value = p.value)
shapiro_rich$alpha_diversity <- "Observed richness"
rm(shapiro_test_results)
#checks
#test <- stat_test %>% filter(Age_individual == "1-2 weeks")
#shapiro.test(test$observed_richnes)

# Observed richness qqline and qq test
ggplot(stat_test, aes(sample = observed_richnes)) + 
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~Age_individual) +
  ggtitle("Observed richness Q-Q Plot per TP")

# Shapiro all results ====
shapiro_all <- rbind(shapiro_ens,shapiro_rich)
rm(shapiro_ens,shapiro_rich)
shapiro_all <- shapiro_all %>% 
  mutate(sig_padj = if_else(
    condition = .$p_value < 0.05, 
    true = "not normally distributed", 
    false = "normally distributed"))
# conclusions, majority of tha data is not normally distributed

#write.table (shapiro_all, "S8_300_shapiro_test_for_alpha_diversity_for_each_TP.txt", row.names=FALSE,sep = "\t")

# ENS wilcox paired test ====
#paired: a logical value specifying that we want to compute a paired Wilcoxon test
test_1w_4w <- wilcox.test(ENS ~ Age_individual, data = wil_1w_4w, alternative="two.sided",paired=T)
stat.test <- data.frame(
  p_value = test_1w_4w$p.value,
  Z = qnorm(test_1w_4w$p.value / 2),
  group1 = '1 wk',
  group2 = '4 wks',
  W_value = test_1w_4w$statistic,
  alternative = test_1w_4w$alternative,
  method = test_1w_4w$method
)

test_4w_8w <- wilcox.test(ENS ~ Age_individual, data = wil_4w_8w, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_4w_8w$p.value,
  Z = qnorm(test_4w_8w$p.value / 2),
  group1 = '4 wks',
  group2 = '8 wks',
  W_value = test_4w_8w$statistic,
  alternative = test_4w_8w$alternative,
  method = test_4w_8w$method
)
stat.test <- rbind(stat.test, s) 

test_8w_4m <- wilcox.test(ENS ~ Age_individual, data = wil_8w_4m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_8w_4m$p.value,
  Z = qnorm(test_8w_4m$p.value / 2),
  group1 = '8 wks',
  group2 = '4 mos',
  W_value = test_8w_4m$statistic,
  alternative = test_8w_4m$alternative,
  method = test_8w_4m$method
)
stat.test <- rbind(stat.test, s) 

test_4m_5m <- wilcox.test(ENS ~ Age_individual, data = wil_4m_5m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_4m_5m$p.value,
  Z = qnorm(test_4m_5m$p.value / 2),
  group1 = '4 mos',
  group2 = '5 mos',
  W_value = test_4m_5m$statistic,
  alternative = test_4m_5m$alternative,
  method = test_4m_5m$method
)
stat.test <- rbind(stat.test, s) 

test_5m_6m <- wilcox.test(ENS ~ Age_individual, data = wil_5m_6m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_5m_6m$p.value,
  Z = qnorm(test_5m_6m$p.value / 2),
  group1 = '5 mos',
  group2 = '6 mos',
  W_value = test_5m_6m$statistic,
  alternative = test_5m_6m$alternative,
  method = test_5m_6m$method
)
stat.test <- rbind(stat.test, s) 

test_6m_9m <- wilcox.test(ENS ~ Age_individual, data = wil_6m_9m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_6m_9m$p.value,
  Z = qnorm(test_6m_9m$p.value / 2),
  group1 = '6 mos',
  group2 = '9 mos',
  W_value = test_6m_9m$statistic,
  alternative = test_6m_9m$alternative,
  method = test_6m_9m$method
)
stat.test <- rbind(stat.test, s) 

test_9m_11m <- wilcox.test(ENS ~ Age_individual, data = wil_9m_11m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_9m_11m$p.value,
  Z = qnorm(test_9m_11m$p.value / 2),
  group1 = '9 mos',
  group2 = '11 mos',
  W_value = test_9m_11m$statistic,
  alternative = test_9m_11m$alternative,
  method = test_9m_11m$method
)
stat.test <- rbind(stat.test, s) 

test_11m_14m <- wilcox.test(ENS ~ Age_individual, data = wil_11m_14m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_11m_14m$p.value,
  Z = qnorm(test_11m_14m$p.value / 2),
  group1 = '11 mos',
  group2 = '14 mos',
  W_value = test_11m_14m$statistic,
  alternative = test_11m_14m$alternative,
  method = test_11m_14m$method
)
stat.test <- rbind(stat.test, s) 

rm(test_1w_4w,test_4w_8w,test_8w_4m,test_4m_5m,test_5m_6m,test_6m_9m,test_9m_11m,test_11m_14m, s)

#renaming columns in p-value df
wilcox_ens = stat.test
names(wilcox_ens)[names(wilcox_ens) == "Freq"] <- "p_value" 
wilcox_ens$alpha_diversity <- "ENS"
wilcox_ens$p_value_test_adj <- p.adjust(wilcox_ens$p_value, method = "BH")
rm(stat.test)

# observed wilcox paired test ====
#paired: a logical value specifying that we want to compute a paired Wilcoxon test
test_1w_4w <- wilcox.test(observed_richnes ~ Age_individual, data = wil_1w_4w, alternative="two.sided",paired=T)
stat.test <- data.frame(
  p_value = test_1w_4w$p.value,
  Z = qnorm(test_1w_4w$p.value / 2),
  group1 = '1 wk',
  group2 = '4 wks',
  W_value = test_1w_4w$statistic,
  alternative = test_1w_4w$alternative,
  method = test_1w_4w$method
)

test_4w_8w <- wilcox.test(observed_richnes ~ Age_individual, data = wil_4w_8w, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_4w_8w$p.value,
  Z = qnorm(test_4w_8w$p.value / 2),
  group1 = '4 wks',
  group2 = '8 wks',
  W_value = test_4w_8w$statistic,
  alternative = test_4w_8w$alternative,
  method = test_4w_8w$method
)
stat.test <- rbind(stat.test, s) 

test_8w_4m <- wilcox.test(observed_richnes ~ Age_individual, data = wil_8w_4m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_8w_4m$p.value,
  Z = qnorm(test_8w_4m$p.value / 2),
  group1 = '8 wks',
  group2 = '4 mos',
  W_value = test_8w_4m$statistic,
  alternative = test_8w_4m$alternative,
  method = test_8w_4m$method
)
stat.test <- rbind(stat.test, s) 

test_4m_5m <- wilcox.test(observed_richnes ~ Age_individual, data = wil_4m_5m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_4m_5m$p.value,
  Z = qnorm(test_4m_5m$p.value / 2),
  group1 = '4 mos',
  group2 = '5 mos',
  W_value = test_4m_5m$statistic,
  alternative = test_4m_5m$alternative,
  method = test_4m_5m$method
)
stat.test <- rbind(stat.test, s) 

test_5m_6m <- wilcox.test(observed_richnes ~ Age_individual, data = wil_5m_6m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_5m_6m$p.value,
  Z = qnorm(test_5m_6m$p.value / 2),
  group1 = '5 mos',
  group2 = '6 mos',
  W_value = test_5m_6m$statistic,
  alternative = test_5m_6m$alternative,
  method = test_5m_6m$method
)
stat.test <- rbind(stat.test, s) 

test_6m_9m <- wilcox.test(observed_richnes ~ Age_individual, data = wil_6m_9m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_6m_9m$p.value,
  Z = qnorm(test_6m_9m$p.value / 2),
  group1 = '6 mos',
  group2 = '9 mos',
  W_value = test_6m_9m$statistic,
  alternative = test_6m_9m$alternative,
  method = test_6m_9m$method
)
stat.test <- rbind(stat.test, s) 

test_9m_11m <- wilcox.test(observed_richnes ~ Age_individual, data = wil_9m_11m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_9m_11m$p.value,
  Z = qnorm(test_9m_11m$p.value / 2),
  group1 = '9 mos',
  group2 = '11 mos',
  W_value = test_9m_11m$statistic,
  alternative = test_9m_11m$alternative,
  method = test_9m_11m$method
)
stat.test <- rbind(stat.test, s) 

test_11m_14m <- wilcox.test(observed_richnes ~ Age_individual, data = wil_11m_14m, alternative="two.sided",paired=T)
s <- data.frame(
  p_value = test_11m_14m$p.value,
  Z = qnorm(test_11m_14m$p.value / 2),
  group1 = '11 mos',
  group2 = '14 mos',
  W_value = test_11m_14m$statistic,
  alternative = test_11m_14m$alternative,
  method = test_11m_14m$method
)
stat.test <- rbind(stat.test, s) 

rm(test_1w_4w,test_4w_8w,test_8w_4m,test_4m_5m,test_5m_6m,test_6m_9m,test_9m_11m,test_11m_14m, s)

#renaming columns in p-value df
wilcox_rich = stat.test
names(wilcox_rich)[names(wilcox_rich) == "Freq"] <- "p_value" 
wilcox_rich$alpha_diversity <- "Observed richness"
wilcox_rich$p_value_test_adj <- p.adjust(wilcox_rich$p_value, method = "BH")
rm(stat.test)

# wilcox results cleaning ====
wilcox_all <- rbind(wilcox_ens,wilcox_rich)
rm(wil_1w_4w,wil_4w_8w,wil_8w_4m,wil_4m_5m,wil_5m_6m,wil_6m_9m,wil_9m_11m,wil_11m_14m)

#function that will round up all numeric values in p-value table
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

#rounding up the values in p-value table
wilcox_all <- round_df(wilcox_all, digits=3)

wilcox_all <- wilcox_all %>% 
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

wilcox_all$group1<- factor(wilcox_all$group1, 
                      levels = c("1 wk",
                                 "4 wks",
                                 "8 wks",
                                 "4 mos",
                                 "5 mos",
                                 "6 mos",
                                 "9 mos",
                                 "11 mos"))
wilcox_all$group2<- factor(wilcox_all$group2, 
                      levels = c("4 wks",
                                 "8 wks",
                                 "4 mos",
                                 "5 mos",
                                 "6 mos",
                                 "9 mos",
                                 "11 mos",
                                 "14 mos"))

# saving data ====
#write.table (wilcox_all, "S2.1_comparison_between_TP_for_alpha_diversity_wilcox_paired_test.txt", row.names=FALSE,sep = "\t")

rm(wilcox_ens, wilcox_rich)

# wilcox plots ====
wilcox_ens <- wilcox_all %>% filter(alpha_diversity == "ENS")
wilcox_rich <- wilcox_all %>% filter(alpha_diversity == "Observed richness")

#ENS
plot_infant <- stat_test %>% ggplot(aes(y = ENS, x = Age_individual, color = Age_individual)) +
  geom_boxplot(width = 0.3) + 
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5,size = 0.75)+ 
  theme_bw()+ggtitle("")+
  xlab("Age")+ylab("ENS")+theme(legend.position = "none",
                                axis.text.x = element_text(size = 9),
                                axis.text.y = element_text(size = 11),
                                axis.title.x = element_text(size = 13), 
                                axis.title.y = element_text(size = 13))+scale_color_viridis(discrete = TRUE)

#adding p-values
boxplot_ENS <- plot_infant + stat_pvalue_manual(wilcox_ens,
                                                y.position = 10, step.increase = 0.1,label = "p_sig_test_adj",hide.ns = TRUE,
                                                size = 5)
boxplot_ENS

#Observed richness
plot_infant <- stat_test %>% ggplot(aes(y = observed_richnes, x = Age_individual, color = Age_individual)) +
  geom_boxplot(width = 0.3,outlier.size = 0.5) + 
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5,size = 0.15)+ 
  ggtitle("")+
  xlab("Age")+ylab("Observed richness")+theme(legend.position = "none",
                                              axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
                                              axis.text.y = element_text(size = 7),
                                              axis.title.x = element_text(size = 7), 
                                              axis.title.y = element_text(size = 7),
                                              plot.title = element_text(size = 1),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank(),
                                              panel.background = element_rect(fill = "white", colour = NA),
                                              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
                                              )+scale_color_viridis(discrete = TRUE)

#adding p-values
boxplot_rich <- plot_infant + stat_pvalue_manual(wilcox_rich,
                                                    y.position = 60, step.increase = 0.1,size = 3,
                                                    label = "p_sig_test_adj",hide.ns = TRUE)
boxplot_rich

# saving figures for the paper ====
ggsave(
  filename = "S2.1_richness_barplot_stats.svg",
  plot     = boxplot_rich,
  width    = 100,
  height   = 60,
  units    = "mm",
  device   = svglite
)

